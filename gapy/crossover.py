import copy
import datetime
import numpy as np

from .individual import Individual


def uniform_crossover(parent_1: Individual, parent_2: Individual, mixing_ratio: float = 0.5):

    """Obtains an offspring individual through uniform crossover.

    Arguments:
        parent_1 (Individual): The first parent individual.
        parent_2 (Individual): The second parent individual.
        mixing_ratio: float: The ratio at which to mix the genomes.

    Returns:
        dict: The child.
    """

    child_genome = []
    for i in range(len(parent_1.genome)):
        if np.random.rand() < mixing_ratio:
            child_genome.append(parent_1.genome[i])
        else:
            child_genome.append(parent_2.genome[i])

    meta = copy.deepcopy(parent_1.meta)
    meta['creation_date'] = str(datetime.datetime.now())
    return Individual(child_genome, meta=meta)
