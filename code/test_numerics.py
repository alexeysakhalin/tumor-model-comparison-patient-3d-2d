"""Regression tests for ranks, matched controls and cluster resampling."""

import unittest

import numpy as np

from build_layer_scores import healthy_layer, site_mask
from build_summary_tables import ORGANS, point_and_ci
from level_matched_scores import control_sets, controls_are_nearest, percentile_ranks


class NumericalTests(unittest.TestCase):
    def test_ties_and_gene_axis_permutation(self):
        values = np.array([[0, 0, 4, 9], [8, 3, 3, 3]], dtype=float)
        expected = np.array([[37.5, 37.5, 75, 100], [100, 50, 50, 50]])
        np.testing.assert_array_equal(percentile_ranks(values), expected)
        permutation = [2, 0, 3, 1]
        np.testing.assert_array_equal(
            percentile_ranks(values[:, permutation]), expected[:, permutation]
        )

    def test_nonfinite_ranks_fail(self):
        with self.assertRaises(ValueError):
            percentile_ranks(np.array([[0, np.nan]]))

    def test_control_matching_and_tie_break(self):
        ranks = np.array([50, 49, 51, 40, 60, 50])
        index = dict(zip(["target", "B", "A", "C", "D", "excluded"], range(6)))
        actual = control_sets(ranks, ["target"], index, {"target", "excluded"}, 3)
        self.assertEqual(actual["target"].tolist(), [2, 1, 3])
        with self.assertRaises(ValueError):
            control_sets(ranks, ["target"], index, {"target", "excluded"}, 5)

    def test_scored_gene_is_never_its_own_control(self):
        # The gene is absent from the exclusion set, as the 200 background genes were: without
        # self-exclusion it sits at distance zero from itself and takes the first control slot.
        ranks = np.array([50, 49, 51, 40, 60, 50])
        index = dict(zip(["target", "B", "A", "C", "D", "other"], range(6)))
        selected = control_sets(ranks, ["target"], index, set(), 3)["target"].tolist()
        self.assertNotIn(index["target"], selected)
        self.assertEqual(len(selected), 3)
        self.assertEqual(len(set(selected)), 3)
        # "other" shares the target's rank exactly, so it is the nearest eligible control.
        self.assertEqual(selected[0], index["other"])
        self.assertTrue(controls_are_nearest(ranks, index, set(), "target", n_ctrl=3))

    def test_self_exclusion_inside_a_tie_block(self):
        # Every candidate has the target's rank: the set must still be full, distinct and self-free.
        ranks = np.zeros(6)
        index = {name: i for i, name in enumerate(["G0", "G1", "G2", "G3", "G4", "G5"])}
        for target in index:
            selected = control_sets(ranks, [target], index, set(), 4)[target].tolist()
            self.assertNotIn(index[target], selected)
            self.assertEqual(sorted(selected), sorted(selected))
            self.assertEqual(len(set(selected)), 4)

    def test_pool_too_small_after_self_exclusion(self):
        # Five eligible genes, one of them the target: four controls remain, five cannot be served.
        ranks = np.arange(5.0)
        index = {f"G{i}": i for i in range(5)}
        self.assertEqual(len(control_sets(ranks, ["G2"], index, set(), 4)["G2"]), 4)
        with self.assertRaises(ValueError):
            control_sets(ranks, ["G2"], index, set(), 5)

    def test_anatomical_selection(self):
        result = site_mask(
            ["STAD", "ESCA", "STAD", "LUAD"],
            ["Adenocarcinoma", "Carcinoma", "Non-Cancerous", "Adenocarcinoma"],
            ["stomach", "stomach", "stomach", "lung"],
        )
        self.assertEqual(result.tolist(), [True, False, False, True])

    def test_missing_healthy_manifest_fails(self):
        import tempfile

        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(FileNotFoundError):
                healthy_layer(directory)

    def test_weighted_bootstrap_matches_explicit_cluster_duplication(self):
        # Unequal cluster sizes and shared donors distinguish cluster resampling
        # from independent sample resampling and from equal weighting of donors.
        by = {}
        for i, organ in enumerate(ORGANS):
            by[organ, "healthy"] = (
                np.array([i - 8.0, i + 2.0, i + 13.0, i + 22.0]),
                np.array(["A", "A", "B", "C"] if i else ["A", "A", "B", "B"]),
            )
        reference_rng = np.random.default_rng(17)
        values, rejected = [], 0
        while len(values) < 120:
            draw = reference_rng.choice(["A", "B", "C"], 3, replace=True)
            medians = []
            for organ in ORGANS:
                scores, labels = by[organ, "healthy"]
                duplicated = [
                    score for cluster in draw for score in scores[labels == cluster]
                ]
                if not duplicated:
                    rejected += 1
                    break
                medians.append(np.median(duplicated))
            if len(medians) == 6:
                values.append(np.median(medians))
        log = []
        actual = point_and_ci(
            by, lay="healthy", n=120, rng=np.random.default_rng(17), reject_log=log
        )
        np.testing.assert_array_equal(actual[1:3], np.percentile(values, [2.5, 97.5]))
        self.assertEqual(log[0]["rejected_draws"], rejected)
        self.assertEqual(actual[3], 6)
        by[ORGANS[0], "healthy"] = (np.array([]), np.array([]))
        with self.assertRaises(ValueError):
            point_and_ci(by, lay="healthy", n=10)


if __name__ == "__main__":
    unittest.main(verbosity=2)
