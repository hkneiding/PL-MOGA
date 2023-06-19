import numpy as np

from .individual import Individual


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

def swap_mutation(individual: Individual, mutation_rate: float = 0.1):

    """Randomly picks two genes of the individuals chromosome and swaps them.

    Arguments:
        individual (Individual): The individual to be mutated.
        mutation_rate (float): The mutation rate.

    Returns:
        Individual: The mutated individual.    
    """

    if np.random.rand() < mutation_rate:
        swap_indices = np.random.choice(list(range(len(individual.genome))), size=2, replace=False)

        swap_gene_1 = individual.genome[swap_indices[0]]
        swap_gene_2 = individual.genome[swap_indices[1]]

        individual.genome[swap_indices[0]] = swap_gene_2
        individual.genome[swap_indices[1]] = swap_gene_1

    return individual