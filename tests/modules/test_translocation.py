import random
from unittest import TestCase

from applications.GeneticModeling.modules.chromosome import IntergenicRegion, Chromosome, Gene
from applications.GeneticModeling.modules.transformation import Translocation


class TestTranslocation(TestCase):
    def test_transform(self):
        # create two chromosomes and transform them. this is a random event so we can test only the lengths
        left = Chromosome([
            IntergenicRegion('ATCG'),
            Gene('TTGT'),
            IntergenicRegion('ATGCTG'),
            Gene('CCCTCG'),
            IntergenicRegion('CTCTCG')
        ])

        right = Chromosome([
            IntergenicRegion('GTCCG'),
            Gene('CCTA'),
            IntergenicRegion('GTCGGGAG'),
            Gene('TCCCTG'),
            IntergenicRegion('CCCCTGAAAA'),
            Gene('GTCCCC'),
            IntergenicRegion('CGTGCTG')
        ])

        print left.describe()
        print right.describe()

        original_left = left.represent()
        original_right = right.represent()

        Translocation(left_chromosome=left, right_chromosome=right).transform()

        self.assertNotEqual(original_left, left.represent())
        self.assertNotEqual(original_right, right.represent())

        self.assertEqual(len(original_left) + len(original_right), len(left.represent()) + len(right.represent()))

        print left.describe()
        print right.describe()

    def test_multiple_transform(self):
        left = Chromosome.parse('<GTACG>AGCG(AGTC)ATAA(AGTA)GGTGC(GTCGC)GCAAC<ACTC>')
        right = Chromosome.parse('<GACGG>TATT(CCAA)TATAC(ATTG)GTAG(AGCG)CCGAC<GTCA>')

        original_length = len(left.represent()) + len(right.represent())

        # execute several transformations and assert the length
        for i in range(40):
            print left.describe()
            print right.describe()
            print
            transformation = Translocation(left_chromosome=left, right_chromosome=right)
            transformation.transform()
            new_length = len(left.represent()) + len(right.represent())

            self.assertEqual(original_length, new_length)

    def test_transform_with_parameters(self):
        """
        source: AT.CGTTGTATGCTGCCCTCGCT.CTCG
        target GTCCGCCTAGTCGGGAGTCCCTGCCCCTG.AAAAGTCCCCCGTGCTG
        expected AT.CTCG
                GTCCGCCTAGTCGGGAGTCCCTGCCCCTG.CGTTGTATGCTGCCCTCGCT.AAAAGTCCCCCGTGCTG
        """
        source = Chromosome([
            IntergenicRegion('ATCG'),
            Gene('TTGT'),
            IntergenicRegion('ATGCTG'),
            Gene('CCCTCG'),
            IntergenicRegion('CTCTCG')
        ])

        target = Chromosome([
            IntergenicRegion('GTCCG'),
            Gene('CCTA'),
            IntergenicRegion('GTCGGGAG'),
            Gene('TCCCTG'),
            IntergenicRegion('CCCCTGAAAA'),
            Gene('GTCCCC'),
            IntergenicRegion('CGTGCTG')
        ])

        source_original_regions = source.regions[:]
        target_original_regions = target.regions[:]

        left_source_region = source.regions[0]  # ATCG
        right_source_region = source.regions[-1]  # CTCTCG
        target_insertion_region = target.regions[4]  # CCCCTGAAAA
        split_left_source_region_at = 1  # AT.CG
        split_right_source_region_at = 1  # CT.CTCG
        split_target_region_at = 5  # CCCTG.AAAA

        translocation = Translocation(source, target)
        translocation.transform_with_parameters(source, target=target,
                                                left_source_region=left_source_region,
                                                right_source_region=right_source_region,
                                                target_insertion_region=target_insertion_region,
                                                split_left_source_region_at=split_left_source_region_at,
                                                split_right_source_region_at=split_right_source_region_at,
                                                split_target_region_at=split_target_region_at,
                                                reverse=False)
        self.assertEqual('ATCTCG', source.represent())
        self.assertEqual('GTCCGCCTAGTCGGGAGTCCCTGCCCCTGCGTTGTATGCTGCCCTCGCTAAAAGTCCCCCGTGCTG', target.represent())

        # test with reverse == True
        source = Chromosome(source_original_regions)
        target = Chromosome(target_original_regions)

        left_source_region = source.regions[0]  # ATCG
        right_source_region = source.regions[-1]  # UTCTCG
        target_insertion_region = target.regions[4]  # UUUUTGAAAA

        translocation = Translocation(source, target)
        translocation.transform_with_parameters(source, target=target,
                                                left_source_region=left_source_region,
                                                right_source_region=right_source_region,
                                                target_insertion_region=target_insertion_region,
                                                split_left_source_region_at=split_left_source_region_at,
                                                split_right_source_region_at=split_right_source_region_at,
                                                split_target_region_at=split_target_region_at,
                                                reverse=True)
        self.assertEqual('ATCTCG', source.represent())
        self.assertEqual('GTCCGCCTAGTCGGGAGTCCCTGCCCCTGAGCGAGGGCAGCATACAACGAAAAGTCCCCCGTGCTG', target.represent())


    def test_transform_with_special_parameters(self):
        """
        Tests a case where the regions break in such a way that two genes are put one after th other
        """
        left = Chromosome.parse("<GTACG>AGTACGGTCGG(CGCT)CTA(AGTA)(GTCGC)GCAGC<ACTC>")
        right = Chromosome.parse("<GACGG>TATCG(AGTC)ATAT(CCAA)TATAC(ATTG)GTC<GTCA>")

        translocation = Translocation(left_chromosome=left, right_chromosome=right)

    def test__split_region_and_insert(self):
        regions = [
            IntergenicRegion('ATCT'),
            IntergenicRegion('GGGTGGG'),
            IntergenicRegion('AAAGA')
        ]

        chromosome = Chromosome(regions)

        original_regions = regions[:]

        Translocation._split_region_and_insert(regions[1], chromosome)

        self.assertEqual(4, len(chromosome.regions))
        self.assertEqual(original_regions[0].content, chromosome.regions[0].content)
        self.assertEqual(original_regions[-1].content, chromosome.regions[-1].content)

        self.assertEqual(original_regions[1].content, chromosome.regions[1].content + chromosome.regions[2].content)

    def test__split_region_random(self):
        for i in range(10):
            length = random.randint(a=3, b=11)
            random_sequence = "".join([c for i in range(length) for c in random.choice(list('ATCG'))])
            new_regions = Translocation._split_region(IntergenicRegion(random_sequence))

            self.assertEqual(random_sequence, new_regions[0].content + new_regions[1].content)

    def test__split_region(self):
        region = IntergenicRegion('ATCGTGCG')

        new_regions = Translocation._split_region(region, breaking_point=3)

        self.assertEqual('ATCG', new_regions[0].content)
        self.assertEqual('TGCG', new_regions[1].content)

