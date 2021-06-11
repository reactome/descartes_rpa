import unittest
import os
import loompy
import anndata
import numpy as np

from tempfile import NamedTemporaryFile
from descartes_rpa.convert.loom import loom_to_anndata


class TestConvertLoom(unittest.TestCase):
    """Class to test convert.loom functions"""
    def setUp(self):
        """Setup attributes for each test iteration object"""
        loom_file = NamedTemporaryFile(suffix=".loom",  delete=False)
        with loompy.new(loom_file.name) as ds:
            m = np.zeros((20, 100))
            ra = {"Gene": [x for x in "ABCDEFGHILJKLMNOPQRS"]}
            ca = {"Cell": np.arange(100)}
            ds.add_columns(m, ca, row_attrs=ra)
            ds.add_columns(m, ca, row_attrs=ra)

        self.loom_file = loom_file.name

    def tearDown(self):
        """Remove setup attributes after each test iteration"""
        os.unlink(self.loom_file)

    def test_convert_to_anndata(self) -> None:
        """Test if function converts loom to AnnData"""
        adata = loom_to_anndata(self.loom_file)
        self.assertIsInstance(adata, anndata._core.anndata.AnnData)
