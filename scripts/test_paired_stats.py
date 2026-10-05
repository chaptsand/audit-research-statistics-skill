#!/usr/bin/env python3
"""Synthetic unit tests; no research inputs or files are changed."""
import unittest
import numpy as np
from scipy import stats
from paired_stats import (
    exact_signed_rank, rank_biserial, hodges_lehmann, bootstrap_hl,
    bh_adjust, analyse_pairs, analyse_family,
)


class PairedStatsTests(unittest.TestCase):
    def test_exact_matches_scipy_no_ties(self):
        d = np.array([1., -2., 3., 4., -5., 6., -7.])
        for alt in ("two-sided", "greater", "less"):
            w, p = exact_signed_rank(d, alt)
            ref = stats.wilcoxon(d, alternative=alt, method="exact")
            self.assertEqual(w, ref.statistic)
            self.assertEqual(p, ref.pvalue)

    def test_ties_and_zeros(self):
        self.assertEqual(exact_signed_rank([0, 1, 1, 1], "greater"), (6.0, 0.125))
        self.assertEqual(exact_signed_rank([0, 1, 1, 1], "two-sided"), (0.0, 0.25))
        self.assertEqual(exact_signed_rank([0, 0]), (0.0, 1.0))
        self.assertEqual(rank_biserial([0, 0]), 0.0)

    def test_direction(self):
        d = [1, 2, -3, 4]
        self.assertAlmostEqual(rank_biserial(d), -rank_biserial(np.negative(d)))
        self.assertAlmostEqual(hodges_lehmann(d), -hodges_lehmann(np.negative(d)))

    def test_hl_is_walsh_not_median(self):
        self.assertEqual(hodges_lehmann([1, 2, 10]), 3.75)

    def test_bh_independent_reference(self):
        p = [0.431640625, 0.76953125, 0.845703125, 1.]
        np.testing.assert_equal(bh_adjust(p), np.ones(4))
        p = [0.04, 0.01, 0.2, 0.03]
        np.testing.assert_allclose(bh_adjust(p), stats.false_discovery_control(p, method="bh"))

    def test_bootstrap_deterministic(self):
        self.assertEqual(bootstrap_hl([0, 0, 0]), (0.0, 0.0))
        self.assertEqual(bootstrap_hl([1, -2, 3]), bootstrap_hl([1, -2, 3]))

    def test_block_bootstrap(self):
        # Homogeneous within each block: preserves blocks rather than treating folds independently.
        block = np.array([[0., 0.], [2., 2.]])
        self.assertEqual(bootstrap_hl(block), (0.0, 2.0))

    def test_reject_bad_family(self):
        with self.assertRaises(ValueError):
            analyse_family({"family_id": "test", "expected_members": ["a", "b"], "p_by_id": {"a": 0.2}})

    def test_helper_records_zero_pairs(self):
        spec = {
            "comparison_id": "synthetic", "left_name": "L", "right_name": "R",
            "left": [1., 2., 3.], "right": [1., 1., 2.], "pair_ids": ["a", "b", "c"],
            "unit": "fixed_run", "alternative": "two-sided", "p_method": "exact_signflip",
            "pairing_evidence": "synthetic test mapping",
        }
        result = analyse_pairs(spec)
        self.assertEqual(result["n_pairs"], 3)
        self.assertEqual(result["n_nonzero"], 2)
        self.assertFalse(result["pairing_verified_by_helper"])

    def test_reject_nonfinite_or_shapes(self):
        with self.assertRaises(ValueError):
            bh_adjust([np.nan])
        with self.assertRaises(ValueError):
            bootstrap_hl([1, 2], n_resamples=100)


if __name__ == "__main__":
    unittest.main(verbosity=2)
