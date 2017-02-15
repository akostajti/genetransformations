from unittest import TestCase

from applications.GeneticModeling.modules.chromosome import IntergenicRegion, Gene
from applications.GeneticModeling.modules.transformation import select_random_region_with_constraints


class TestSelect_random_region_with_constraints(TestCase):
    def test_select_random_region_with_constraints(self):
        regions = [IntergenicRegion('ATCTA'),
                   Gene('TTTATG'),
                   IntergenicRegion('GAGGCT'),
                   Gene('TCTGTC'),
                   IntergenicRegion('TATATC'),
                   Gene('AGTCGTA'),
                   IntergenicRegion('TTTGAT')]

        # no essentialgenes, everything works as the default
        selected = select_random_region_with_constraints(regions,
                                                         essential_genes_in_window=1,
                                                         essential_genes_window_size=2)
        self.assertEqual(1, len(selected))
        self.assertTrue(selected[0] in regions)

        selected = select_random_region_with_constraints(regions,
                                                         count=2,
                                                         essential_genes_in_window=1,
                                                         essential_genes_window_size=2)

        self.assertEqual(2, len(selected))
        self.assertTrue(len([s for s in selected if s in regions]))
