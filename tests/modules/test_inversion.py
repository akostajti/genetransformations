from unittest import TestCase

from applications.GeneticModeling.modules.chromosome import IntergenicRegion, Gene, Chromosome
from applications.GeneticModeling.modules.transformation import Inversion


class TestInversion(TestCase):
    def test_transform_with_regions(self):
        regions = self._create_regions()
        chromosome = Chromosome(regions)

        Inversion(chromosome).transform_with_regions(left_region=regions[2],
                                                     right_region=regions[-1],
                                                     left_region_breaking_point=5,
                                                     right_region_breaking_point=7)

        # original value: 'TCGTTCATCTGACGGTT.CCCCGTGTGTGCCCGTCCCGCCCTCCCGGGG.TCT' (breaking at the dots)
        # expected result 'TCGTTCATCTGACGGTT.CCCCGGGAGGGCGGGACGGGCACACACGGGG.TCT' (without the dot)

        expected = 'TCGTTCATCTGACGGTTCCCCGGGAGGGCGGGACGGGCACACACGGGGTCT'
        self.assertEqual(expected, chromosome.represent())

    def test_transform(self):
        regions = self._create_regions()
        chromosome = Chromosome(regions)

        original_representation = chromosome.represent()

        Inversion(chromosome).transform()

        self.assertEqual(len(original_representation), len(chromosome.represent()))

    def _create_regions(self):
        regions = [
            IntergenicRegion('TCGTTC'),
            Gene('ATCTGA'),
            IntergenicRegion('CGGTTC'),
            Gene('CCCGT'),
            IntergenicRegion('GTGTG'),
            Gene('CCCGT'),
            Gene('CCCGCCCT'),
            IntergenicRegion('CCCGGGGTCT')
        ]

        return regions

