from unittest import TestCase

from applications.GeneticModeling.modules.chromosome import IntergenicRegion
from applications.GeneticModeling.modules.transformation import compute_probabilities_based_on_length


class TestCompute_probabilities_based_on_length(TestCase):
    def test_compute_probabilities_based_on_length(self):
        regions = [
            IntergenicRegion('ACGTCGT'),
            IntergenicRegion('AAACTCT')
        ]

        probabilities = compute_probabilities_based_on_length(regions)
        self.assertEqual([0.5, 0.5], probabilities)

        regions = [
            IntergenicRegion('A'),
            IntergenicRegion('AC'),
            IntergenicRegion('ACT'),
            IntergenicRegion('ACTG')
        ]

        probabilities = compute_probabilities_based_on_length(regions)
        self.assertEqual([0.1, 0.2, 0.3, 0.4], probabilities)
