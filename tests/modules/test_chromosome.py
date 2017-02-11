from unittest import TestCase

from applications.GeneticModeling.modules.chromosome import Chromosome
from applications.GeneticModeling.modules.chromosome import Gene
from applications.GeneticModeling.modules.chromosome import IntergenicRegion

class TestChromosome(TestCase):
    def test_represent(self):
        regions = [
            Gene('ADAC'),
            Gene('GGGA'),
            IntergenicRegion('TTGA'),
            IntergenicRegion('TCAG')
        ]
        chromosome = Chromosome(regions)

        self.assertEqual('ADACGGGATTGATCAG', chromosome.represent())

    def test_reverse(self):
        regions = [
            Gene('ADAC'),
            Gene('GGGA'),
            IntergenicRegion('TTGA'),
            IntergenicRegion('TCAG')
        ]
        chromosome = Chromosome(regions)

        chromosome.reverse()
        self.assertEqual(''.join(reversed('ADACGGGATTGATCAG')), chromosome.represent())

    def test_get_breakable_regions(self):
        regions = [
            Gene('ADAC'),
            Gene('GGGA'),
            IntergenicRegion('TTGA'),
            IntergenicRegion('TCAG')
        ]
        chromosome = Chromosome(regions)

        # expected: the two intergenic regions
        breakable_regions = chromosome.get_breakable_regions()
        self.assertEqual(len(breakable_regions), 2)
        self.assertEqual(len([region for region in breakable_regions if region.can_break]), 2)

        # set the first gene region breakable. expected: 3 regions
        regions[0].can_break = True
        breakable_regions = chromosome.get_breakable_regions()
        self.assertEqual(len(breakable_regions), 3)
        self.assertEqual(len([region for region in breakable_regions if region.can_break]), 3)

        # set the last intergenic region  non-breakable. expected: 2 regions
        regions[3].can_break = False
        breakable_regions = chromosome.get_breakable_regions()
        self.assertEqual(len(breakable_regions), 2)
        self.assertEqual(len([region for region in breakable_regions if region.can_break]), 2)

    def test_parse(self):
        chromosome_descripton = '''<ADCGTGGG>AAAGT(DAC)TTTGACU(UUTGAAA)AGT<CCCGTU>'''
        chromosome = Chromosome.parse(chromosome_descripton)

        # expected: 7. 2 non-breakable regions at the beginning and the end and the other regions
        self.assertEqual(len(chromosome.regions), 7)

        expectations = [
            (False, 'ADCGTGGG'),
            (True, 'AAAGT'),
            (False, 'DAC'),
            (True, 'TTTGACU'),
            (False, 'UUTGAAA'),
            (True, 'AGT'),
            (False, 'CCCGTU')
        ]

        for region, expectation in zip(chromosome.regions, expectations):
            self.assertEqual(region.can_break, expectation[0])
            self.assertEqual(region.content, expectation[1], region.represent())