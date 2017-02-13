import os

def index():
    """
    Shows the simulation page
    """
    # get the simulation data from the request
    number_of_transformations = int(request.vars['number_of_transformations'])
    rate_of_translocations = float(request.vars['rate_of_translocations'])
    random_error = float(request.vars['random_error'])

    longer_breaks_often = bool(request.vars['longer_breaks_often'])
    use_coexpression = bool(request.vars['use_coexpression'])
    use_essential_gene_pairs = bool(request.vars['use_essential_gene_pairs'])

    sequence_patterns = request.vars['sequence_patterns']

    scheduler, task = queue_simulation(left_chromosome_file=request.vars['left_chromosome_file'],
                     right_chromosome_file=request.vars['right_chromosome_file'],
                     rate_of_translocations=rate_of_translocations,
                     number_of_transformations=number_of_transformations,
                     sequence_patterns=sequence_patterns,
                     use_essential_gene_pairs=use_essential_gene_pairs,
                     use_coexpression=use_coexpression,
                     longer_breaks_often=longer_breaks_often,
                     random_error=random_error)

    redirect(URL('simulation', 'result', vars=dict(task_id=task.id)))


def result():
    task_id = int(request.vars['task_id'])
    status = get_task_status(task_id)
    traceback = status.scheduler_run.traceback

    result = status.result

    # write the gene orders to files
    left_result_file = "left-chromosome-result.txt"
    write_ordinals_to_file(result['history'], 'left', left_result_file)

    right_result_file = "right-chromosome-result.txt"
    write_ordinals_to_file(result['history'], 'left', right_result_file)

    return dict(traceback=BEAUTIFY(traceback),
                status=status,
                left_result_download_url=URL(c='default', f='download', args=[left_result_file]),
                right_result_download_url=URL(c='default', f='download', args=[right_result_file]))
    """
    result['history']
    return dict(traceback=result['history'], status=status, left_result_download_url='GG', right_result_download_url='RR')
    """


def write_ordinals_to_file(history, key, left_result_file):
    with open(os.path.join(request.folder, 'uploads/' + left_result_file), 'a') as output:
        for history_item in history:
            ordinals = history_item[key]['ordinals']
            output.write(', '.join([str(ordinal) for ordinal in ordinals]) + '\n')



