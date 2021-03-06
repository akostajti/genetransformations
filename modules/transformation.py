"""
Defines the chromosome transformations
"""
from numpy import random, math

from applications.GeneticModeling.modules.chromosome import DnaBaseMappings, IntergenicRegion

import logging

logger = logging.getLogger("web2py.app.GeneticModeling")
logger.setLevel(logging.DEBUG)


def select_random_region_with_constraints(regions,
                                          essential_genes_window_size,
                                          essential_genes_in_window,
                                          count=1,
                                          longer_breaks_often=True,
                                          random_error=None):
    """
    Selects a random region segment based on the essential gene constraints.

    Essential genes are special genes. The assumption is that if in a given window size (expressed in number of genes)
    there are more than x essential genes then the intergenic regions in that segment cannot break. So for example
    if E1,E2,E3 are essential genes and the window size is 3 with 2 essential genes required so that the segment is
    non-breakable then the segment in the below example is non-breakable.

     I1 G1 I2 G2 E1 [G1 I3 E2 I4 E3] I5 G3

     Window size 3 means that we check the segments of 3 consecuting genes. essential_genes_in_window == 2 means that if
     there are at lest 2 essential genes in the segment than it cannot break. This function returns a list of segments
     that can break respecting these constraints.
    """
    MAX_TRIALS = 100

    regions = regions[:]  # copy the list because this function modifies the list

    result = []
    for i in range(count):
        # randomly find a sublist that matches the constraint
        genes = [region for region in regions if not region.can_break]

        if len(genes) < essential_genes_window_size:
            """
            if there are less genes than the essential window size then just randomly select th remaining regions
            """
            selected = select_random_region(regions, count=count - i,
                                            longer_breaks_often=longer_breaks_often,
                                            random_error=random_error)
            result += selected
            return result

        found = False
        trials = 0
        while not found and trials < MAX_TRIALS:
            trials += trials
            # the last index where we can start the sublist to stay inside the boundaries
            max_start = len(genes) - essential_genes_window_size
            sublist_start = random.randint(0, high=max_start)
            gene_sublist = genes[sublist_start:sublist_start + essential_genes_window_size]

            # count the number of essential genes in the section
            essentials = sum(1 for gene in gene_sublist if gene.is_essential)
            if essentials < essential_genes_in_window:
                found = True

                # randomly select a region from the subsegment containing the genes
                subsegment = regions[regions.index(gene_sublist[0]):regions.index(gene_sublist[-1]) + 1]
                selected = select_random_region(subsegment, count=1,
                                     longer_breaks_often=longer_breaks_often,
                                     random_error=random_error)

                result += selected

                # remove the subsegment elements from the regions to prevent them from being selected again
                regions = [region for region in regions if region not in subsegment]

        if trials == MAX_TRIALS:
            raise RuntimeError("Reached maximum number of trials (" + MAX_TRIALS + ")")

    return result


def select_random_region(regions,
                         count=1,
                         longer_breaks_often=True,
                         random_error=None,
                         essential_genes_window_size=None,
                         essential_genes_in_window=None):
    """
    Randomly selects some regions from the regions list.

    :param regions: the region list from which to select
    :param count: the number of regions to select. If greater than len(regions) throws an error.
    :param longer_breaks_often: if true then the regions that are longer break with a higher probability
    :param random_error the error used when selecting the regions that can break. Normally only intergenic regions
        can break. If random_error is greater than 0 then with some chance also non-breaking regions can be selected
        for break.
    :param essential_genes_window_size the essential gene window size. see select_random_region_with_constraints
    :param essential_genes_in_window the number of essential genes in the window required for the subsegment
        for being non-breakable. see select_random_region_with_constraints

    :returns a list of size count containing regions from the regions list. If there are no matching regions
    than returns an empty list. The elements in the resulting list are in the same order as in the original.
    """

    if essential_genes_in_window and essential_genes_window_size:
        return select_random_region_with_constraints(regions,
                                                     essential_genes_window_size=essential_genes_window_size,
                                                     essential_genes_in_window=essential_genes_in_window,
                                                     count=count,
                                                     longer_breaks_often=longer_breaks_often,
                                                     random_error=random_error)

    def error():
        return random.choice([True, False], p=[random_error, 1.0 - random_error])

    if random_error is None:
        filtered = [region for region in regions if region.can_break]
    else:
        filtered = [region for region in regions if region.can_break or error()]

    if len(filtered) == 0:
        return []

    if longer_breaks_often:
        probabilities = compute_probabilities_based_on_length(filtered)

    else:
        probabilities = None

    if len(filtered) < count:
        return None

    result = random.choice(filtered, p=probabilities, size=count, replace=False)

    # sort the items to their original order
    tuples = sorted(enumerate(result))
    return map(lambda (k, v): v, tuples)


def compute_probabilities_based_on_length(regions):
    """
    Returns a weight between 0.0 and 1.0 for each region in the list. The weight is higher for longer
    regions
    """
    total_length = sum([len(region.content) for region in regions])
    probabilities = [float(len(region.content)) / total_length for region in regions]

    return probabilities


class Transformation:
    """
    Base class for all transformations.

    Arguments:
        longer_breaks_often: if true then the longer intergenic regions break with a higher probability
    """
    def __init__(self,
                 longer_breaks_often=True,
                 random_error=None,
                 essential_genes_window_size=None,
                 essential_genes_in_window=None):
        self.longer_breaks_often = longer_breaks_often
        self.random_error = random_error
        self.essential_genes_window_size=essential_genes_window_size
        self.essential_genes_in_window=essential_genes_in_window

    def select_random_region(self, regions, count=1):
        return select_random_region(regions, count=count,
                             longer_breaks_often=self.longer_breaks_often,
                             essential_genes_in_window=self.essential_genes_in_window,
                             essential_genes_window_size=self.essential_genes_window_size,
                             random_error=self.random_error)



class Inversion(Transformation):
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
    def __init__(self,
                 chromosome,
                 longer_breaks_often=True,
                 random_error=None,
                 essential_genes_window_size=None,
                 essential_genes_in_window=None):
        Transformation.__init__(self,
                                longer_breaks_often=longer_breaks_often,
                                random_error=random_error,
                                essential_genes_window_size=essential_genes_window_size,
                                essential_genes_in_window=essential_genes_in_window)
        self.chromosome = chromosome

    def transform_with_regions(self,
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
        breaking_points = self.select_random_region(self.chromosome.regions, count=2)

        left_region_breaking_point = random.randint(0, high=len(breaking_points[0].content))
        right_region_breaking_point = random.randint(0, high=len(breaking_points[1].content))

        self.transform_with_regions(left_region=breaking_points[0],
                                    right_region=breaking_points[1],
                                    left_region_breaking_point=left_region_breaking_point,
                                    right_region_breaking_point=right_region_breaking_point)


class Translocation(Transformation):
    """
    A more complex transformation that involves two chromosomes.

    In the target chromosome we break a selected region and move a part of the source chromosome
    between those points. Eventually the moved part can be reversed (this is totally random).

    Note that the target chromosome may be the same as the source # TODO: implement this

    Attributes:
        left_chromosome: the left chromosome
        right_chromosome: the right chromosome
    """
    def __init__(self,
                 left_chromosome,
                 right_chromosome,
                 longer_breaks_often=True,
                 random_error=None,
                 essential_genes_window_size=None,
                 essential_genes_in_window=None):
        Transformation.__init__(self, longer_breaks_often=longer_breaks_often,
                                random_error=random_error,
                                essential_genes_window_size=essential_genes_window_size,
                                essential_genes_in_window=essential_genes_in_window)
        self.left_chromosome = left_chromosome
        self.right_chromosome = right_chromosome

    def transform_with_parameters(self,
                                  source,
                                  target,
                                  left_source_region,
                                  right_source_region,
                                  target_insertion_region,
                                  split_left_source_region_at=None,
                                  split_right_source_region_at=None,
                                  split_target_region_at=None,
                                  reverse=False):
        """
        Translocation based on the input parameters

        Example:
            source: ATC.UTTTCG.CT (breaks at the dots)
            target UUUGTA.CTGGG (breaks at the dot)

            result
            ATC..CT
            UUUGTA.UTTTCG.CTGGG

        :param source The source chromosome. One section from this chromosome wil be moved to target. Note that target
            may eventually be the same as source
        :param target The section from source will be moved to this chromosome
        :param left_source_region The first region that will break in source
        :param right_source_region The second region that will break in source
        :param target_insertion_region This region will break in target and the section from source will be added to the
            breaking point
        :param split_left_source_region_at The index where the left source region breaks
        :param split_right_source_region_at The index where the right source region breaks
        :param split_target_region_at The target region will break at this index
        :param reverse If this parameter is True then after moving the section to target the method executes an Inversion
            on the same section (reverses it)
        """
        append_to_same = source == target

        logger.debug("transform with parameters: %s, %s, %s, %d, %d, %d, %s", left_source_region.content,
                     right_source_region.content, target_insertion_region.content, split_left_source_region_at,
                     split_right_source_region_at, split_target_region_at, str(reverse))

        print "transform with parameters: {0}, {1}, {2}, {3}, {4}, {5}, {6}".format(left_source_region.content,
                      right_source_region.content, target_insertion_region.content, split_left_source_region_at,
                      split_right_source_region_at, split_target_region_at, str(reverse))
        print

        # TODO: implement the case when the source and the target are the same

        # in the target create two new intergenic regions from the broken one
        new_regions = Translocation._split_region_and_insert(target_insertion_region, target, split_target_region_at)

        # insert the parts from the source chromosome AFTER the first new region created
        insertion_index = target.regions.index(new_regions[1])

        # move the sections from between the source splitting regions
        to_move_left_index = source.regions.index(left_source_region) + 1
        to_move_right_index = source.regions.index(right_source_region)
        to_move = source.regions[to_move_left_index:to_move_right_index]

        target.regions = target.regions[:insertion_index] + to_move + target.regions[insertion_index:]
        source.regions[to_move_left_index:to_move_right_index] = []

        # finally break the breaking intergenic regions in the source and move their parts to target
        left_source_region_split = self._split_region(left_source_region, breaking_point=split_left_source_region_at)
        right_source_region_split = self._split_region(right_source_region, breaking_point=split_right_source_region_at)

        # in the source merge the two broken intergenic region parts
        merged = IntergenicRegion(left_source_region_split[0].content + right_source_region_split[1].content)
        source.regions[source.regions.index(left_source_region)] = merged
        del source.regions[source.regions.index(right_source_region)]

        prefix_length = len(new_regions[0].content)
        postfix_length = len(right_source_region_split[0].content)

        new_regions[0].content = new_regions[0].content + left_source_region_split[1].content
        new_regions[1].content = right_source_region_split[0].content + new_regions[1].content

        # TODO if there's a random inversion then call Infersion on the moved part in the target chromosome
        if reverse:
            inversion = Inversion(target)
            inversion.transform_with_regions(left_region=new_regions[0],
                                             right_region=new_regions[1],
                                             left_region_breaking_point=prefix_length,
                                             right_region_breaking_point=postfix_length)
        logger.debug("left chromosome after transformation: " + self.left_chromosome.describe())
        logger.debug("right chromosome after transformation: " + self.right_chromosome.describe())

    def transform(self):
        """
        Random translocation
        """
        chromosomes = [self.left_chromosome, self.right_chromosome]
        source = random.choice(chromosomes)

        chromosomes.remove(source)
        target = chromosomes[0]

        # these are the tw regions that will break in the source
        selected_regions = self.select_random_region(source.regions, count=2)
        if selected_regions is None:
            logger.debug("Not enough breakable regions, quitting translocation " + source.describe())
            return

        left_source_region, right_source_region = selected_regions

        # if the left region comes affter the right then swap the two regions
        if source.regions.index(right_source_region) < source.regions.index(left_source_region):
            tmp = left_source_region
            left_source_region = right_source_region
            right_source_region = tmp

        # this is the insertion point in the target region
        target_insertion_region = self.select_random_region(target.regions)[0]

        reverse = random.choice([True, False])

        if left_source_region is None or \
                        right_source_region is None or \
                        target_insertion_region is None:
            logger.debug("left_source_region, right_source_region or target_insertion_region is empty, "
                         "returning from translocation")
            return
        self.transform_with_parameters(source, target=target,
                       left_source_region=left_source_region, right_source_region=right_source_region,
                       target_insertion_region=target_insertion_region, reverse=reverse)

    @staticmethod
    def _split_region_and_insert(region, target_chromosome, breaking_point=None):
        """
        Splits the region into two parts and inserts the resulting regions back in to the
        chromosome at the original position
        """
        index = target_chromosome.regions.index(region)

        new_regions = Translocation._split_region(region, breaking_point=breaking_point)

        target_chromosome.regions = target_chromosome.regions[:index] + new_regions + \
                                    target_chromosome.regions[index + 1:]

        return new_regions

    @staticmethod
    def _split_region(region, breaking_point=None):
        content = region.content

        logger.debug("splitting region. content: " + content)
        if breaking_point is None:
            breaking_point = random.randint(0, high=len(content) - 1)

        new_regions = [IntergenicRegion(content[:breaking_point + 1]),
                       IntergenicRegion(content[breaking_point + 1:])]

        return new_regions

