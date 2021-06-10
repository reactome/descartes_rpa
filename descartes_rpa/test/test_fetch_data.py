import unittest


from descartes_rpa.fetch.ensebml import ensembl_to_gene


class TestFetchEnsebml(unittest.TestCase):
    """Class to test Ensebml functions"""
    def setUp(self):
        """Method to start each test iteration object"""
        self.ensembl_ids = [
            "ENSG00000223972",
            "ENSG00000237613",
            "ENSG00000243485"
        ]

    def test_correct_fetch(self):
        """Test if fetch can correctly fetch gene names"""
        gene_names = ensembl_to_gene(self.ensembl_ids)
        expected_names = ["DDX11L1", "FAM138A", "MIR1302-2HG"]
        for name in gene_names:
            self.assertIn(name, expected_names)


if __name__ == '__main__':
    unittest.main()
