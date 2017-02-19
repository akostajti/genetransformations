"""
Controller for generating random chromosomes
"""

from numpy import random


def index():
    """
    Returns the form and generates the random chromosome
    """
    form = _get_generator_form().process()

    if form.accepted:
        session.flash = T('Generating the chromosome')
        filename = generate_chromosome_description(form)
        redirect(URL('generate', 'download', args=[filename]))

    return dict(form=form)


def download():
    """
    Downloads a simulation file
    """
    file_name = request.args[0]

    response.stream(os.path.join(request.folder, 'uploads/' + file_name), attachment=True, filename=file_name)


def _get_generator_form():
    form = SQLFORM.factory(
        Field('number_of_genes', 'integer',
              requires=IS_INT_IN_RANGE(minimum=1, maximum=1000),
              default=100,
              label=T('Number of genes')),
        Field('minimum_gene_length', 'integer',
              requires=IS_INT_IN_RANGE(minimum=1, maximum=1000),
              default=10,
              label=T('Minimum gene length')),
        Field('maximum_gene_length', 'integer',
              requires=IS_INT_IN_RANGE(minimum=1, maximum=1000),
              default=100,
              label=T('Maximum gene length')),
        Field('minimum_intergenic_region_length', 'integer',
              requires=IS_INT_IN_RANGE(minimum=1, maximum=1000),
              default=10,
              label=T('Minimum intergenic region length')),
        Field('maximum_intergenic_region_length', 'integer',
              requires=IS_INT_IN_RANGE(minimum=1, maximum=1000),
              default=100,
              label=T('Maximum intergenic region length')),
        Field('generate_essential_genes', 'boolean',
              default=False,
              label=T('Generate essential genes'))
    )

    return form


def generate_chromosome_description(form):
    """
    Generates a chromosome description based on the form parameters. The chromosome description is in a format parseable
    by the Chromosome class.
    """
    generate_essential_genes = form.vars['generate_essential_genes']
    number_of_genes = form.vars['number_of_genes']
    maximum_gene_length = form.vars['maximum_gene_length']
    minimum_gene_length = form.vars['minimum_gene_length']
    minimum_intergenic_region_length = form.vars['minimum_intergenic_region_length']
    maximum_intergenic_region_length = form.vars['maximum_intergenic_region_length']

    index = random.randint(1, high=10000)
    filename = 'generated-chromosome' + str(index) + '.txt'
    with open(os.path.join(request.folder, 'uploads/' + filename), 'w') as output:
        # generate the opening non-breakable region
        start = _generate_random_region(random.randint(minimum_intergenic_region_length, high=maximum_intergenic_region_length))
        output.write('<' + start + '>')

        # generate the genes and intergenics in sequence
        next_is_gene = True
        count = number_of_genes

        while count > 0:
            count -= 1
            if next_is_gene:
                gene =  _generate_random_region(random.randint(minimum_gene_length, high=maximum_gene_length))

                if generate_essential_genes:
                    essential = random.choice([True, False], p=[0.2, 0.8])
                else:
                    essential = False

                output.write('(' + gene)
                if essential:
                    output.write(';')
                output.write(')')
            else:
                intergenic = _generate_random_region(random.randint(minimum_intergenic_region_length, high=maximum_intergenic_region_length))
                output.write(intergenic)
            next_is_gene = not next_is_gene

        # generate the opening non-breakable region
        end = _generate_random_region(
            random.randint(minimum_intergenic_region_length, high=maximum_intergenic_region_length))
        output.write('<' + end + '>')

    return filename


def _generate_random_region(length):
    return "".join(random.choice(['A', 'T', 'C', 'G'], size=length))

