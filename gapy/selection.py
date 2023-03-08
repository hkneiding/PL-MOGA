import numpy as np

from .individual import Individual


def select_by_fitness(individuals: list, n_selected: int):

    """Selects n individuals from a list of individuals by their absolute fitness.

    Arguments:
        individuals (list[Individual]): List of individuals.
        n_selected (int): The number individuals to select.

    Returns:
        list[Individual]: The selected individuals.
    """

    return sorted(individuals, key=lambda x: x.fitness_sum, reverse=True)[0:n_selected]

def select_by_rank(individuals: list, n_selected: int, rank_function: callable):

    """Selects n individuals from a list of individuals by their rank.

        Arguments:
        individuals (list[Individual]): List of individuals.
        n_selected (int): The number individuals to select.

    Returns:
        list[Individual]: The selected individuals.
    """

    rankings = rank_function(individuals)
    return [x for _, x in sorted(zip(rankings, individuals), key=lambda pair: pair[0], reverse=False)][0:n_selected]

def roulette_wheel_rank(individuals: list, n_selected: int, rank_function: callable):

    """Selects individuals from a list of individuals by randomly drawing based on
    probabilities proportional to their relative rank.

    Arguments:
        individuals (list[Individual]): List of individuals.
        n_selected (int): The number individuals to select.

    Returns:
        list[Individual]: The selected individuals.
    """

    rankings = rank_function(individuals)
    inverse_rankings = [max(rankings) + 1 - (rank) for rank in rankings]
    selection_probabilities = [inverse_ranking/sum(inverse_rankings) for inverse_ranking in inverse_rankings]
    
    selection = np.random.choice(individuals, size=n_selected, p=selection_probabilities).tolist()
    return selection


def roulette_wheel_fitness(individuals: list, n_selected: int):

    """Selects individuals from a list of individuals by randomly drawing based on
    probabilities proportional to their absolute fitness.

    Arguments:
        individuals (list[Individual]): List of individuals.
        n_selected (int): The number individuals to select.

    Returns:
        list[Individual]: The selected individuals.
    """

    fitness_sum = sum([individual.fitness_sum for individual in individuals])
    relative_fitnesses = [individual.fitness_sum/fitness_sum for individual in individuals]
    
    selection = np.random.choice(individuals, size=n_selected, p=relative_fitnesses).tolist()
    return selection
