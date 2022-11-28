import numpy as np

from individual import Individual


def uniform_integer_mutation(individual: Individual, mutation_space: int, mutation_rate: float = 0.1):

    """Randomly mutates an individuals genome within the given mutation space.

    Arguments:
        individual (Individual): The individual to be mutated.
        mutation_space (int): .
        mutation_rate (float): The mutation rate.

    Returns:
        Individual: The mutated individual.
    """

    for i in range(len(individual.genome)):
        if np.random.rand() < mutation_rate:
            individual.genome[i] = np.random.randint(0, mutation_space)

    return individual

