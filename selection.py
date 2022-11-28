import numpy as np

from individual import Individual


def select_by_fitness(individuals: list, n_selected: int):

    """Selects individuals from a list of individuals by their absolute fitness.

    Arguments:
        individuals (list[Individual]): List of individuals.
        n_selected (int): The number individuals to select.

    Returns:
        list[Individual]: The selected individuals.
    """

    return sorted(individuals, key=lambda x: x.fitness_sum, reverse=True)[0:n_selected]

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

def pareto_domination_rank(individuals: list, n_selected):

    """Selects individuals from a list of individuals by their pareto domination rank.

    Arguments:
        individuals (list[Individual]): List of individuals.
        n_selected (int): The number individuals to select.

    Returns:
        list[Individual]: The selected individuals.
    """

    rankings = [0 for _ in individuals]
    for i, individual_1 in enumerate(individuals):
        for j, individual_2 in enumerate(individuals):

            if i == j:
                continue
                
            if individual_1.dominates(individual_2):
                rankings[i] += 1
            
    return [x for _, x in sorted(zip(rankings, individuals), key=lambda pair: pair[0])][0:n_selected]
