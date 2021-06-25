import unittest
import pandas as pd
import numpy as np

from anndata import AnnData
from descartes_rpa.analyze.analyze import scanpy_format


class TestFormatAnnData(unittest.TestCase):
    def setUp(self):
        """Setup attributes for each test iteration object"""
        df = pd.DataFrame(
            {
                "test": ["just", "for", "test"],
                "gene_short_name": ["a", "a", "b"]
            }
        )
        obs = pd.DataFrame(
            {
                "Main_cluster_umap_1": list(np.ones((2, 3))),
                "Main_cluster_umap_2": list(np.zeros((2, 3)))
            }
        )
        self.adata = AnnData(
            np.ones((2, 3)),
            var=df,
            obs=obs
        )

    def test_scanpy_format_make_unique(self):
        """Test if scanpy formatting can create unique index"""
        scanpy_format(self.adata)
        expected_names = ["a", "a-1", "b"]
        var_names = list(self.adata.var_names)
        self.assertEqual(var_names, expected_names)

    def test_scanpy_format_x_umap(self):
        """Test if scanpy formatting can create
        right X_umap from obs to obsm
        """
        x_umap = self.adata.obs[
            ["Main_cluster_umap_1", "Main_cluster_umap_2"]
        ].to_numpy()
        scanpy_format(self.adata)
        self.assertNotIn(
            False,
            np.equal(self.adata.obsm["X_umap"], x_umap, dtype=np.object)
        )
