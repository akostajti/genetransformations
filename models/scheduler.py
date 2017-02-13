import copy
import random

import editdistance

from applications.GeneticModeling.modules.chromosome import Chromosome
from applications.GeneticModeling.modules.transformation import Translocation, Inversion
from gluon.scheduler import Scheduler

import os

def simulation(number_of_transformations,
                   rate_of_translocations,
                   random_error,
                   longer_breaks_often,
                   use_coexpression,
                   use_essential_gene_pairs,
                   sequence_patterns,
                   left_chromosome_file,
                   right_chromosome_file):

    # build the chromosomes from the files
    uploads_folder = os.path.join(request.folder, 'uploads')
    left_chromosome = _build_chromosome_from_file(os.path.join(uploads_folder, left_chromosome_file))
    right_chromosome = _build_chromosome_from_file(os.path.join(uploads_folder, right_chromosome_file))

    number_of_translocations = int(rate_of_translocations * number_of_transformations)
    number_of_inversions = number_of_transformations - number_of_translocations

    # this sequence contains the list if transformations to execute
    transformation_sequence = ['Translocation'] * number_of_translocations + ['Inversion'] * number_of_inversions
    random.shuffle(transformation_sequence)

    # history contains tuples: each tuple contains the left and the right chromosome in the given step
    history = []
    original_left = copy.deepcopy(left_chromosome)
    original_right = copy.deepcopy(right_chromosome)

    for step in transformation_sequence:
        if step == 'Translocation':
            translocation = Translocation(left_chromosome, right_chromosome)
            translocation.transform()
        elif step == 'Inversion':
            selected = random.choice([left_chromosome, right_chromosome])
            inversion = Inversion(selected)
            inversion.transform()
        history.append(dict(left=_create_history_item(left_chromosome, original_left),
                            right=_create_history_item(right_chromosome, original_right),
                            transformation=step))

    result = dict(history=history,
                  final_left=left_chromosome.represent(),
                  final_right=right_chromosome.represent())
    return result


def _create_history_item(chromosome, previous):
    ordinals = chromosome.get_gene_ordinals()

    levenshtein_distance = editdistance.eval(previous.represent(), chromosome.represent())

    return dict(ordinals=ordinals, distance=levenshtein_distance)


def _build_chromosome_from_file(file_name):
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
