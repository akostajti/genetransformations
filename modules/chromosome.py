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
    ESSENTIAL_GENE_MARKER = ';'  # this marker must be at the end of the gene description (just before the closing bracket)

    def __init__(self, regions):
        self.regions = regions

    def represent(self):
        """
        Returns the string representation of this chromosome

        :returns the joined list of the regions
        """
        return "".join([region.represent() for region in self.regions])

    def reverse(self):
        """
        Reverses the region list and the regions themselves
        """
        for region in self.regions:
            region.reverse()
        self.regions.reverse()

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
        return [region for region in self.regions if region.can_break]

    def get_gene_ordinals(self):
        """
        Returns the original gene ordinals in their current order.
        """
        return [gene.ordinal for gene in self.regions if isinstance(gene, Gene)]

    def describe(self):
        """
        Returns the chromosome description (parseable by Chromosome.parse)
        """
        result = '<' + self.regions[0].content + '>'
        for region in self.regions[1:-1]:
            if region.is_gene:
                result += '(' + region.content
                if region.is_essential:
                    result += ';'
                result += ')'
            else:
                result += region.content

        result += '<' + self.regions[-1].content + '>'
        return result

    @staticmethod
    def parse(chromosome_string, use_coexpression=False):
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

        # if the current gene is essential or not
        is_essential = False

        """
        inside a coexpression segment neither
        the intergenic regions can break.
        """
        for next_char in chromosome_string:
            if next_char == Chromosome.LEFT_GENE_BOUNDARY:
                # create intergenic region and clear the buffer
                Chromosome.create_intergenic_region(char_buffer, can_break, regions)
                char_buffer = []
            elif next_char == Chromosome.RIGHT_GENE_BOUNDARY:
                # create a gene
                gene = Gene("".join(char_buffer), is_essential=is_essential)
                is_essential = False
                regions.append(gene)
                char_buffer = []
            elif next_char == Chromosome.LEFT_COEXPRESSION_BOUNDARY or next_char == Chromosome.LEFT_NOBREAK_BOUNDARY:
                # create intergenic region and clear the buffer
                Chromosome.create_intergenic_region(char_buffer, can_break, regions)
                char_buffer = []

                if next_char == Chromosome.LEFT_NOBREAK_BOUNDARY or use_coexpression:
                    can_break = False
            elif next_char == Chromosome.RIGHT_COEXPRESSION_BOUNDARY or next_char == Chromosome.RIGHT_NOBREAK_BOUNDARY:
                # create intergenic region and clear the buffer
                Chromosome.create_intergenic_region(char_buffer, can_break, regions)
                char_buffer = []

                can_break = True
            elif next_char == Chromosome.ESSENTIAL_GENE_MARKER:
                is_essential = True
            else:
                char_buffer.append(next_char)
        return Chromosome(regions)

    @staticmethod
    def create_intergenic_region(char_buffer, can_break, regions):
        if len(char_buffer):
            intergenic_region = IntergenicRegion("".join(char_buffer))
            intergenic_region.can_break = can_break
            regions.append(intergenic_region)


class ChromosomeRegion:
    """
    Base class for chromosome regions

    Attributes:
        can_break: if this chromosome is breakable
        content: the string representation of the region (a sequence of amino acid characters)
        reversed: if the region was reversed during a transformation
    """
    def __init__(self, content):
        self.content = content
        self.can_break = False
        self.reversed = False
        self.is_gene = False

    def reverse(self):
        """
        Reverses the region
        """
        self.content = "".join(reversed(DnaBaseMappings.map_string(self.content)))

    def represent(self):
        return self.content


class IntergenicRegion(ChromosomeRegion):
    """
    A sequence that is not a gene but is between two genes. The only chromosome region
    that is breakable.
    """
    def __init__(self, content):
        ChromosomeRegion.__init__(self, content)
        self.can_break = True


class Gene(ChromosomeRegion):
    """
    Represents a gene (a chromosome region that is not breakable)

    Attributes:
        ordinal: a unique number that is assigned to the region on creation
        is_essential: If the gene is an essential gene
    """
    next_ordinal = 0

    def __init__(self, content, is_essential=False):
        ChromosomeRegion.__init__(self, content)
        self.is_essential = is_essential
        self.ordinal = Gene.next_ordinal
        Gene.next_ordinal += 1
        self.is_gene = True


class DnaBaseMappings:
    """
    A simple class that converts amino acid characters to their pairs (A to T, G to C etc)
    """
    MAPPINGS = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G'}

    @staticmethod
    def map_string(string):
        return "".join([DnaBaseMappings.MAPPINGS[c.upper()] for c in string])



