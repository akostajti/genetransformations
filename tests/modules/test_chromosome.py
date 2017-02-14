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
            Gene('ATAC'),
            Gene('GGGA'),
            IntergenicRegion('TTGA'),
            IntergenicRegion('TCAG')
        ]
        chromosome = Chromosome(regions)

        chromosome.reverse()

        expected = 'CTGATCAATCCCGTAT'
        self.assertEqual(expected, chromosome.represent())

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

    def test_parse_coexpression(self):
        chromosome_descripton = '''<ADCGTGGG>{AAAGT(DAC)TTTGACU(UUTGAAA)}AGT<CCCGTU>'''

        chromosome = Chromosome.parse(chromosome_descripton)

        # expected: the intergenic regions between the curly brackets cannot break
        self.assertEqual(len(chromosome.regions), 7)

        expectations = [
            (False, 'ADCGTGGG'),
            (False, 'AAAGT'),
            (False, 'DAC'),
            (False, 'TTTGACU'),
            (False, 'UUTGAAA'),
            (True, 'AGT'),
            (False, 'CCCGTU')
        ]

        for region, expectation in zip(chromosome.regions, expectations):
            self.assertEqual(region.can_break, expectation[0])
            self.assertEqual(region.content, expectation[1], region.represent())

    def test_parse_essential_gene(self):
        chromosome_descripton = '''<ADCGTGGG>AAAGT(DAC)TTTGACU(UUTGAAA;)AGT<CCCGTU>'''

        chromosome = Chromosome.parse(chromosome_descripton)

        # expected: the gene UUTGAAA is essential
        self.assertEqual(len(chromosome.regions), 7)

        self.assertTrue(chromosome.regions[4].is_essential)

    def test_get_gene_ordinals(self):
        chromosome = Chromosome(
            [
                IntergenicRegion('TCTA'),
                Gene('AAAGTA'),
                IntergenicRegion('CCCCCGTG'),
                Gene('ATCTGA'),
                Gene('TTTATTTA')
            ]
        )

        original_ordinals = chromosome.get_gene_ordinals()[:]

        # change the order
        chromosome.regions[0] = chromosome.regions[-1]
        chromosome.regions = chromosome.regions[:-1]
        self.assertEqual([original_ordinals[2], original_ordinals[0], original_ordinals[1]], chromosome.get_gene_ordinals())