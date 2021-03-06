import copy
import random

import editdistance

from applications.GeneticModeling.modules.chromosome import Chromosome
from applications.GeneticModeling.modules.transformation import Translocation, Inversion
from gluon.scheduler import Scheduler

import os

import logging

logger = logging.getLogger("web2py.app.GeneticModeling")
logger.setLevel(logging.DEBUG)


def simulation(number_of_transformations,
                   rate_of_translocations,
                   random_error,
                   longer_breaks_often,
                   use_coexpression,
                   use_essential_gene_pairs,
                   sequence_patterns,
                   left_chromosome_file,
                   right_chromosome_file,
               essential_genes_window_size,
               essential_genes_in_window,
               compute_diffs=False):

    # build the chromosomes from the files
    uploads_folder = os.path.join(request.folder, 'uploads')
    left_chromosome = _build_chromosome_from_file(os.path.join(uploads_folder, left_chromosome_file),
                                                  use_coexpression=use_coexpression)
    right_chromosome = _build_chromosome_from_file(os.path.join(uploads_folder, right_chromosome_file),
                                                   use_coexpression=use_coexpression)

    number_of_translocations = int((float(rate_of_translocations)/100) * number_of_transformations)
    number_of_inversions = number_of_transformations - number_of_translocations

    logger.debug("number of translocations: %d, number of inversions: %d", number_of_translocations, number_of_inversions)

    # this sequence contains the list if transformations to execute
    transformation_sequence = ['Translocation'] * number_of_translocations + ['Inversion'] * number_of_inversions
    random.shuffle(transformation_sequence)

    # history contains tuples: each tuple contains the left and the right chromosome in the given step
    history = []

    # combine the two chromosome representations for later diffs
    original_combination = left_chromosome.represent() + right_chromosome.represent()

    if not use_essential_gene_pairs:
        essential_genes_in_window =None
        essential_genes_window_size = None

    for step in transformation_sequence:
        if step == 'Translocation':
            transformation = Translocation(left_chromosome,
                                          right_chromosome,
                                          random_error=random_error,
                                          longer_breaks_often=longer_breaks_often,
                                            essential_genes_in_window=essential_genes_in_window)
        elif step == 'Inversion':
            selected = random.choice([left_chromosome,
                                      right_chromosome])
            transformation = Inversion(selected,
                                  random_error=random_error,
                                  longer_breaks_often=longer_breaks_often,
                                    essential_genes_window_size=essential_genes_window_size,
                                    essential_genes_in_window=essential_genes_in_window)

        transformation.transform()

        levenshtein_distance = 0
        if compute_diffs:
            combination = left_chromosome.represent() + right_chromosome.represent()
            levenshtein_distance = editdistance.eval(original_combination, combination)

        history.append(dict(left=_create_history_item(left_chromosome),
                            right=_create_history_item(right_chromosome),
                            transformation=step,
                            distance=levenshtein_distance))

    result = dict(history=history,
                  final_left=left_chromosome.represent(),
                  final_right=right_chromosome.represent())
    return result


def _create_history_item(chromosome):
    ordinals = chromosome.get_gene_ordinals()

    return dict(ordinals=ordinals)


def _build_chromosome_from_file(file_name, use_coexpression=False):
    with open(file_name, 'r') as file:
        chromosome_description = file.read()
        chromosome = Chromosome.parse(chromosome_description)

        return chromosome


scheduler = Scheduler(db, dict(simulation=simulation))


def queue_simulation(**args):
    task = scheduler.queue_task(simulation, pvars=args, immediate=True, timeout=3000)

    return scheduler, task


def get_task_status(task_id):
    return scheduler.task_status(task_id, output=True)
