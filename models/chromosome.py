"""
All chromosome related class and function definitions
"""


class Chromosome:
    """
    Represents a chromosome

    Attributes:
        regions: the ChromosomeRegions this chromosome consists of.
    """
    RIGHT_NOBREAK_BOUNDARY = '>'
    LEFT_NOBREAK_BOUNDARY = '<'
    RIGHT_GENE_BOUNDARY = ')'
    LEFT_GENE_BOUNDARY = '('
    LEFT_COEXPRESSION_BOUNDARY = '{'
    RIGHT_COEXPRESSION_BOUNDARY = '}'

    def __init__(self):
        self.regions = []

    def represent(self):
        """
        Returns the string representation of this chromosome

        :returns the joined list of the regions
        """
        return "".join(self.regions)

    def reverse(self):
        """
        Reverses the region list and the regions themselves
        """
        self.regions = [region.reverse() for region in self.regions].reverse()

    def diff(self, chromosome, diff_method):
        """
        Computes the difference between this chromosome and the other using the diff method
        """
        pass

    def get_breakable_regions(self):
        """
        Returns the chromosome regions that can break

        :returns a list of ChromosomeRegion instances
        """
        return [region for region in self.regions if region.canBreak]

    @staticmethod
    def parse(chromosome_string):
        """
        Parses the string representation of a chromosome anr returns the chromosome

        :param chromosome_string the textual representaion of the chromosome. This string
         contains the amino acid initials that make up the chromosome and some special separator
         characters marking the gene boundaries, the coexpression boundaries and some other non breakable
         segments.

        :returns a Chromosome instance
        """
        regions = []  # the list of chromosome regions parsed
        char_buffer = []
        can_break = True

        for next_char in chromosome_string:
            if next_char == Chromosome.LEFT_GENE_BOUNDARY:
                # create intergenic region and clear the buffer
                Chromosome.create_intergenic_region(char_buffer, can_break)
                char_buffer = []
            elif next_char == Chromosome.RIGHT_GENE_BOUNDARY:
                # create a gene
                gene = Gene("".join(char_buffer))
                regions.append(gene)
                char_buffer = []
                pass
            elif next_char == Chromosome.LEFT_COEXPRESSION_BOUNDARY or next_char == Chromosome.LEFT_NOBREAK_BOUNDARY:
                # create intergenic region and clear the buffer
                Chromosome.create_intergenic_region(char_buffer, can_break)
                char_buffer = []
                can_break = False
            elif next_char == Chromosome.RIGHT_COEXPRESSION_BOUNDARY or next_char == Chromosome.RIGHT_NOBREAK_BOUNDARY:
                # create intergenic region and clear the buffer
                Chromosome.create_intergenic_region(char_buffer, can_break)
                char_buffer = []
                can_break = True
            else:
                char_buffer.append(next_char)

    @staticmethod
    def create_intergenic_region(char_buffer, can_break):
        intergenic_region = IntergenicRegion("".join(char_buffer))
        intergenic_region.can_break = can_break


class ChromosomeRegion:
    """
    Base class for chromosome regions

    Attributes:
        can_break: if this chromosome is breakable
        content: the string representation of the region (a sequence of amino acid characters)
        reversed: if the region was reversed during a transformation
        ordinal: a unique number that is assigned to the region on creation
    """
    nextOrdinal = 1

    def __init__(self, content):
        self.content = content
        self.can_break = False
        self.reversed = False
        self.ordinal = ChromosomeRegion.nextOrdinal

        ChromosomeRegion.nextOrdinal += 1
    
    def reverse(self):
        """
        Reverses the region
        """
        self.content = self.content.reverse()


class IntergenicRegion(ChromosomeRegion):
    """
    A sequence that is not a gene but is between two genes. The only chromosome region
    that is breakable.
    """
    def __init__(self, content):
        ChromosomeRegion.__init__(self, content)
        self.can_reak = True


class Gene(ChromosomeRegion):
    """
    Represents a gene (a chromosome region that is not breakable)
    """
    def __init__(self, content):
        ChromosomeRegion.__init__(self, content)


