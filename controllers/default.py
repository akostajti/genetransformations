
import os

def index():
    """
    Shows the simulation form page
    """
    form = _get_simulation_form().process()

    if form.accepted:
        session.flash=T('The simulation starts')
        redirect(URL(c='simulation', vars=form.vars))

    return dict(form=form)


def _get_simulation_form():
    form = SQLFORM.factory(
        Field('left_chromosome_file', 'upload', uploadfolder=os.path.join(request.folder,'uploads'),
              requires=IS_NOT_EMPTY(),
              label=T('Left chromosome file')),
        Field('right_chromosome_file', 'upload', uploadfolder=os.path.join(request.folder,'uploads'),
              requires=IS_NOT_EMPTY(),
              label=T('Right chromosome file')),
        Field('number_of_transformations', 'integer',
              requires=IS_NOT_EMPTY(),
              label=T('Number of transformations'),
              default=50),
        Field('rate_of_translocations', 'integer',
              default=0,
              requires=IS_INT_IN_RANGE(minimum=0, maximum=100),
              label=T('Percentage of translocations')),
        Field('random_error', 'double',
              default=0.0,
              requires=IS_DECIMAL_IN_RANGE(minimum=0.0, maximum=1.0),
              label=T('Probability of random error')),
        Field('diff_strategy', 'string',
              requires=IS_IN_SET(_get_diff_strategies(), zero=None),
              label=T('Diff Strategy')),
        Field('sequence_patterns', 'text',
              label=T('Sequence patterns used for breaking (comma separated)')),
        Field('longer_breaks_often', 'boolean',
              label=T('Longer intergenic regions break with higher probability'),
              default=True),
        Field('use_essential_gene_pairs', 'boolean',
              label=T('Use essential gene pairs')),
        Field('essential_genes_window_size', 'integer',
              label=T('Essential genes window size')),
        Field('essential_genes_in_window', 'integer',
              label=T('Essential genes in window to prevent breaking')),
        Field('use_coexpression', 'boolean',
              label=T('Use coexpression')),
        Field('compute_diff', 'boolean',
              label=T('Compute diffs'),
              default=False)
    )

    return form


def _get_diff_strategies():
    return {
        'GeneLevelLevenshtein': T('Gene Level Levenshtein'),
        'EditDistance': T('Levenshtein distance'),
        'Hamming': T('Hamming distance'),
        'SMC': T('SMC'),
        'Entropy': T('Entropy')
    }


@cache.action()
def download():
    """
    allows downloading of uploaded files
    http://..../[app]/default/download/[filename]
    """
    return response.download(request, db)


def call():
    """
    exposes services. for example:
    http://..../[app]/default/call/jsonrpc
    decorate with @services.jsonrpc the functions to expose
    supports xml, json, xmlrpc, jsonrpc, amfrpc, rss, csv
    """
    return service()


