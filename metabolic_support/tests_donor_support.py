#!/usr/bin/env python3
"""Synthetic invariants only; no COVID group outcomes are loaded."""
import unittest
import numpy as np
from compute_donor_support import SupportKernel, sample_matched, bootstrap_mean, bh


class SupportTests(unittest.TestCase):
    def setUp(self):
        self.ids = ["a", "b", "c", "d"]
        self.sets = {"a": {"requirements": ["x", "y"], "providers": ["x"]},
                     "b": {"requirements": ["y"], "providers": ["x", "y"]},
                     "c": {"requirements": [], "providers": ["y"]},
                     "d": {"requirements": [], "providers": []}}
        self.kernel = SupportKernel(self.ids, self.sets, ["a", "b"])

    def test_self_excluded_with_absent_reference_recipient(self):
        result = self.kernel.scores([1, 0, 0, 0])
        self.assertAlmostEqual(result["mss_reference"][0], 0)
        result = self.kernel.scores([0, 1, 0, 0])
        self.assertAlmostEqual(result["mss_reference"][0], .5)

    def test_union_and_recipient_denominator(self):
        result = self.kernel.scores([0, 0, 1, 0])
        self.assertAlmostEqual(result["mss_reference"][0], .75)
        self.assertAlmostEqual(result["redundancy_ge2"][0], 0)
        result = self.kernel.scores([1, 1, 1, 0])
        self.assertAlmostEqual(result["mss_reference"][0], 1)
        self.assertAlmostEqual(result["redundancy_mean"][0], 1.25)
        self.assertAlmostEqual(result["redundancy_ge2"][0], .25)

    def test_batch_matches_set_oracle(self):
        supplier_sets = np.random.default_rng(42).integers(0, 2, (100, 4)).astype(bool)
        expected = []
        for present in supplier_sets:
            score = []
            for recipient in ["a", "b"]:
                providers = set().union(*(set(self.sets[k]["providers"])
                                          for k, yes in zip(self.ids, present) if yes and k != recipient))
                needs = set(self.sets[recipient]["requirements"])
                score.append(len(needs & providers) / len(needs))
            expected.append(np.mean(score))
        np.testing.assert_allclose(self.kernel.scores(supplier_sets)["mss_reference"], expected)

    def test_addition_removal_monotonic(self):
        rng = np.random.default_rng(9)
        base = rng.integers(0, 2, (50, 4)).astype(bool)
        target = rng.integers(0, 2, (50, 4)).astype(bool)
        original = self.kernel.scores(base)["mss_reference"]
        self.assertTrue(np.all(self.kernel.scores(base | target)["mss_reference"] >= original))
        self.assertTrue(np.all(self.kernel.scores(base & ~target)["mss_reference"] <= original))

    def test_null_exact_strata_without_replacement(self):
        target = np.array([1, 0, 1, 0, 1, 0], bool)
        pool = np.ones(6, bool)
        strata = np.array(["x", "x", "y", "y", "z", "z"])
        draws, movable, audit = sample_matched(target, pool, strata, np.random.default_rng(42), 100)
        self.assertEqual(movable, 3)
        for label in ["x", "y", "z"]:
            self.assertTrue(np.all(draws[:, strata == label].sum(axis=1) == 1))
        self.assertTrue(np.any(np.all(draws == target, axis=1)))
        self.assertEqual(sum(x["target_n"] for x in audit), 3)

    def test_degenerate_null_recorded(self):
        target = np.array([1, 1, 0, 0], bool)
        draws, movable, _ = sample_matched(target, np.ones(4, bool),
                                           np.array(["a", "a", "b", "b"]), np.random.default_rng(1), 10)
        self.assertEqual(movable, 0)
        self.assertTrue(np.all(draws == target))

    def test_invalid_target_rejected(self):
        with self.assertRaises(ValueError):
            sample_matched(np.array([1, 0], bool), np.array([0, 1], bool),
                           np.array(["x", "x"]), np.random.default_rng(1), 5)

    def test_bootstrap_source_counts_preserved(self):
        boot = bootstrap_mean(np.array([0., 0., 3.]), np.array(["A", "A", "B"]), np.random.default_rng(42), 100)
        np.testing.assert_array_equal(boot, np.ones(100))

    def test_bh_missing_and_order(self):
        q = bh([.01, np.nan, .04, .03])
        np.testing.assert_allclose(q[[0, 2, 3]], [.03, .04, .04])
        self.assertTrue(np.isnan(q[1]))

    def test_family_single_support_and_resident(self):
        result = self.kernel.detail(np.array([0, 0, 1, 0], bool), np.array([0, 0, .1, 0]),
                                    np.array(["A", "B", "C", "D"]))
        self.assertAlmostEqual(result["single_family_share"], 1)
        self.assertTrue(np.isnan(result["mss_resident"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
