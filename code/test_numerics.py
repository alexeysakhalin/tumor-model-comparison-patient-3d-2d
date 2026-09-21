"""Regression tests for ranks, matched controls and cluster resampling."""

import unittest

import numpy as np

from build_layer_scores import healthy_layer, site_mask
from build_summary_tables import ORGANS, point_and_ci
from level_matched_scores import control_sets, percentile_ranks


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
