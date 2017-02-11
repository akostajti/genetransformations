"""
Defines the chromosome transformations
"""
import random

from applications.GeneticModeling.modules.chromosome import DnaBaseMappings


def select_random_region(regions,
                         count=1,
                         only_breakable=True):
    """
    Randomly selects some regions from the regions list.

    :param regions: the region list from which to select
    :param count: the number of regions to select. If greater than len(regions) throws an error.
    :param only_breakable: if True then only breakable regions are returned

    :returns a list of size count containing regions from the regions list. If there are no matching regions
    than returns an empty list. The elements in the resulting list are in the same order as in the original.
    """
    if only_breakable:
        filtered = [region for region in regions if region.can_break]
    else:
        filtered = regions[:]  # copy the list because the function modifies it

    if len(filtered) == 0:
        return []

    result = []
    for i in range(0, count):
        choice = random.choice(filtered)
        filtered.remove(choice)
        result.append(choice)

    # sort the items to their original order
    tuples = sorted(zip([regions.index(region) for region in result], result))
    return map(lambda (k, v): v, tuples)


class Inversion:
    """
    The simplest transformation.

    Breakc a chromosome at two points and reverses the region between the two points. The breaking
    points are always inside breaking regions (can_break == True). In most cases this means breaking inside
    intergenic regions but when the random error parameter is not zero then the genes can also break with
    some probability.

    Example:
        the chromosome 'AGT.TCGAA.GT' breaks at the points marked with dots. After the inversion the chromosome
        looks like this: 'AGT.TTCGA.GT'. Note that during the reversal the amino acid characters are replaced
        with their pairs.

    Attributes:
        chromosome: the chromosome to transform
    """
    def __init__(self, chromosome):
        self.chromosome = chromosome

    def transform_with_regions(self,
                               regions,
                               left_region,
                               right_region,
                               left_region_breaking_point,
                               right_region_breaking_point):
        """
        Transforms the chromosome
        """
        regions = self.chromosome.regions

        left_index = regions.index(left_region)
        right_index = regions.index(right_region)

        # reverse regions
        to_reverse = regions[left_index + 1:right_index]
        [region.reverse() for region in to_reverse]
        reversed_regions = list(reversed(to_reverse))
        regions[left_index + 1:right_index] = reversed_regions

        # and also update the breaking regions
        old_left_content = left_region.content
        old_right_content = right_region.content

        new_left_content = old_left_content[:left_region_breaking_point] \
                           + DnaBaseMappings.map_string("".join(reversed(old_right_content[:right_region_breaking_point])))
        new_right_content = DnaBaseMappings.map_string("".join(reversed(old_left_content[left_region_breaking_point:]))) \
                            + old_right_content[right_region_breaking_point:]

        left_region.content = new_left_content
        right_region.content = new_right_content

    def transform(self):
        """
        Transforms the chromosome
        """
        regions = self.chromosome.regions
        breaking_points = select_random_region(self.chromosome.regions, count=2)

        left_region_breaking_point = random.randint(a=0, b=len(breaking_points[0].content))
        right_region_breaking_point = random.randint(a=0, b=len(breaking_points[1].content))

        self.transform_with_regions(regions,
                                    left_region=breaking_points[0],
                                    right_region=breaking_points[1],
                                    left_region_breaking_point=left_region_breaking_point,
                                    right_region_breaking_point=right_region_breaking_point)


class Translocation:
    pass