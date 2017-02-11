from unittest import TestCase

from applications.GeneticModeling.modules.chromosome import ChromosomeRegion


class TestChromosomeRegion(TestCase):
    def test_reverse(self):
        region = ChromosomeRegion('ATTGTGT')
        reversed_content = 'ACACAAT' # reversed and each amino acid is mapped to theidr pair
        region.reverse()
        self.assertEqual(reversed_content, region.represent())

    def test_represent(self):
        region = ChromosomeRegion('ATTGTGT')
        self.assertEqual(region.content, region.represent())
