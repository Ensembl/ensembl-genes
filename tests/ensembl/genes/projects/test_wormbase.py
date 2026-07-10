import unittest
from unittest.mock import patch, Mock

from ensembl.genes.projects.wormbase_ftp import WormBaseHTTP
from ensembl.genes.projects.wormbase_renderer import WormBaseRenderer


class TestWormBaseHTTP(unittest.TestCase):
    @patch("requests.get")
    def test_get_species(self, mock_get):
        mock_response = Mock()
        mock_response.text = '<a href="species_a/">species_a/</a><a href="species_b/">species_b/</a>'
        mock_response.raise_for_status = Mock()
        mock_get.return_value = mock_response

        client = WormBaseHTTP("WBPS19")
        species = client.get_species()
        self.assertCountEqual(species, ["species_a", "species_b"])

    @patch("requests.get")
    def test_get_bioprojects(self, mock_get):
        mock_response = Mock()
        mock_response.text = '<a href="PRJEB1/">PRJEB1/</a><a href="other/">other/</a>'
        mock_response.raise_for_status = Mock()
        mock_get.return_value = mock_response

        client = WormBaseHTTP("WBPS19")
        bps = client.get_bioprojects("species_a")
        self.assertEqual(bps, ["PRJEB1"])

    @patch("requests.get")
    def test_get_files(self, mock_get):
        mock_response = Mock()
        mock_response.text = '<a href="file.txt">file.txt</a>'
        mock_response.raise_for_status = Mock()
        mock_get.return_value = mock_response

        client = WormBaseHTTP("WBPS19")
        files = client.get_files("species_a", "PRJEB1")
        self.assertEqual(files, ["file.txt"])


class TestWormBaseRenderer(unittest.TestCase):
    def setUp(self):
        self.renderer = WormBaseRenderer("WBPS19")

    def test_derive_species_name(self):
        self.assertEqual(self.renderer._derive_species_name("acanthocheilonema_viteae"), "Acanthocheilonema viteae")
        self.assertEqual(self.renderer._derive_species_name("caenorhabditis_elegans"), "Caenorhabditis elegans")

    @patch("ensembl.genes.projects.wormbase_renderer.check_url_status")
    def test_render_row(self, mock_check_url):
        mock_check_url.return_value = True

        files = [
            "species_a.PRJEB1.WBPS19.genomic.fa.gz",
            "species_a.PRJEB1.WBPS19.annotations.gff3.gz",
        ]
        
        row, missing = self.renderer.render_row("species_a", "PRJEB1", files)
        
        self.assertEqual(row["species"], "Species a")
        self.assertEqual(row["bioproject"], "PRJEB1")
        self.assertEqual(row["bioproject_link"], "https://www.ncbi.nlm.nih.gov/bioproject/PRJEB1")
        self.assertTrue(row["genome"].endswith("species_a.PRJEB1.WBPS19.genomic.fa.gz"))
        self.assertTrue(row["annotation_gff3"].endswith("species_a.PRJEB1.WBPS19.annotations.gff3.gz"))
        self.assertEqual(row["species_page"], "https://parasite.wormbase.org/species_a_prjeb1")
        
        self.assertIn("species_a.PRJEB1.WBPS19.protein.fa.gz", missing)
        self.assertNotIn("species_a.PRJEB1.WBPS19.genomic.fa.gz", missing)

    def test_is_valid_row(self):
        self.assertTrue(self.renderer.is_valid_row({"genome": "url"}))
        self.assertTrue(self.renderer.is_valid_row({"annotation_gtf": "url"}))
        self.assertFalse(self.renderer.is_valid_row({"species": "name", "bioproject": "PRJEB1"}))
