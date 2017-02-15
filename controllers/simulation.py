import os


def index():
    """
    Shows the simulation page
    """
    # get the simulation data from the request
    number_of_transformations = int(request.vars['number_of_transformations'])
    rate_of_translocations = float(request.vars['rate_of_translocations'])
    random_error = float(request.vars['random_error'])

    longer_breaks_often = 'True' == request.vars['longer_breaks_often']
    use_coexpression = 'True' == request.vars['use_coexpression']
    use_essential_gene_pairs = 'True' == request.vars['use_essential_gene_pairs']

    sequence_patterns = request.vars['sequence_patterns']

    essential_genes_window_size = None
    essential_genes_in_window = None
    if use_essential_gene_pairs:
        if 'essential_genes_window_size' in request.vars:
            essential_genes_window_size = int(request.vars['essential_genes_window_size'])

        if 'essential_genes_in_window' in request.vars:
            essential_genes_in_window = int(request.vars['essential_genes_in_window'])

    scheduler, task = queue_simulation(left_chromosome_file=request.vars['left_chromosome_file'],
                     right_chromosome_file=request.vars['right_chromosome_file'],
                     rate_of_translocations=rate_of_translocations,
                     number_of_transformations=number_of_transformations,
                     sequence_patterns=sequence_patterns,
                     use_essential_gene_pairs=use_essential_gene_pairs,
                     use_coexpression=use_coexpression,
                     longer_breaks_often=longer_breaks_often,
                     random_error=random_error,
                                       essential_genes_in_window=essential_genes_in_window,
                                       essential_genes_window_size=essential_genes_window_size)

    redirect(URL('simulation', 'result', vars=dict(task_id=task.id)))


def result():
    task_id = int(request.vars['task_id'])
    status = get_task_status(task_id)
    traceback = status.scheduler_run.traceback

    result = status.result

    data = dict(traceback=BEAUTIFY(traceback),
                status=status)
    if status.scheduler_run.status == 'COMPLETED':
        # write the gene orders to files
        left_ordinals_file = "left-chromosome-ordinals.txt"
        _write_ordinals_to_file(result['history'], 'left', left_ordinals_file)

        right_ordinals_file = "right-chromosome-ordinals.txt"
        _write_ordinals_to_file(result['history'], 'left', right_ordinals_file)

        left_result_file = 'left-chromosome-result.txt'
        _write_chromosome_to_file(result['final_left'], left_result_file)

        right_result_file = 'right-chromosome-result.txt'
        _write_chromosome_to_file(result['final_right'], right_result_file)

        data['left_result_download_url'] = _get_download_url(left_result_file)
        data['right_result_download_url'] = _get_download_url(right_result_file)

        data['left_ordinals_download_url'] = _get_download_url(left_ordinals_file)
        data['right_ordinals_download_url'] = _get_download_url(right_ordinals_file)

        # create the dataset for the chart
        dataset = [history_item['distance'] for history_item in result['history'] if 'distance' in history_item]
        data['dataset'] = dataset

    return data


def download():
    """
    Downloads a simulation file
    """
    file_name = request.args[0]

    response.stream(os.path.join(request.folder, 'uploads/' + file_name), attachment=True, filename=file_name)


def _write_ordinals_to_file(history, key, left_result_file):
    with open(os.path.join(request.folder, 'uploads/' + left_result_file), 'w') as output:
        for history_item in history:
            ordinals = history_item[key]['ordinals']
            output.write(', '.join([str(ordinal) for ordinal in ordinals]) + '\n')


def _write_chromosome_to_file(chromosome_description, file_name):
    with open(os.path.join(request.folder, 'uploads/' + file_name), 'w') as output:
        output.write(chromosome_description)


def _get_download_url(file_name):
    return URL(c='simulation', f='download', args=[file_name])