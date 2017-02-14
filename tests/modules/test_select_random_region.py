from unittest import TestCase

from applications.GeneticModeling.modules.chromosome import IntergenicRegion, Gene
from applications.GeneticModeling.modules.transformation import select_random_region


class TestSelect_random_region(TestCase):
    def test_select_random_region(self):
        # create 10 random intergenic regions
        regions = [IntergenicRegion('AGTTCG') for i in range(1, 10)]

        choices = select_random_region(regions)
        self.assertEqual(len(choices), 1)

        # check if all choices are in the original list and they are not identical
        self.assertTrue(regions.index(choices[0]) >= 0)

        # selecting more items
        choices = select_random_region(regions, count=3)
        self.assertEqual(len(choices), 3)

        # check if all choices are in the original list and they are not identical
        uniq = set([choice for choice in choices if regions.index(choice) >= 0])
        self.assertEqual(len(uniq), 3)

        # create a list of non-breakable regions
        non_breakable = [Gene('AGTCCCCC') for i in range(1, 10)]
        choices = select_random_region(non_breakable)
        self.assertEqual(len(choices), 0)

        # test with random error == 1.0
        choices = select_random_region(non_breakable, random_error=1.0)
        self.assertEqual(len(choices), 1)