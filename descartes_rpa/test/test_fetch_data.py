import unittest


from descartes_rpa.fetch.ensebml import ensembl_to_gene
from descartes_rpa.fetch.descartes import fetch_de_genes_for_cell_type


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


class TestFetchDescartes(unittest.TestCase):
    """Class to test descartes data fetching"""

    def test_fetch_de_genes(self):
        """Method to test fetch of Descartes
        Differentially Expressed genes for each
        Main Cell type
        """
        de_mapping = fetch_de_genes_for_cell_type()

        self.assertEqual(len(de_mapping), 77)
        self.assertEqual(len(de_mapping["Microglia"]), 782)


if __name__ == '__main__':
    unittest.main()
