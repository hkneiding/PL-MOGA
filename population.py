from typing import Callable

from individual import Individual


class Population:

    def __init__(self, individuals: list):

        """Constructor

        Arguments:
            individuals (list[Individual]): List of Individual objects.
        """

        self.individuals = individuals

    def select(self, selection_function: Callable, n_elitism: int = 0):
        
        """Selects individuals based on specified rule.

        Arguments:
            selection_function (Callable): The function used to select individuals.
            n_elitism (int): The number of elitism to apply.
        """

        # take out n_elitism individuals
        selected_individuals = self.individuals[:n_elitism]
        # perform selection on remaining
        selected_individuals += selection_function(self.individuals[n_elitism:])

        # update individuals
        self.individuals = selected_individuals

    def add_individual(self, individual: Individual):

        """Adds an individual to the population.

        Arguments:
            individual (Individual): The individual to be added.
        """

        self.individuals.append(individual)

    def extend(self, population):

        """Extends the population by another population object.

        Arguments:
            population (Population): The population to be added.
        """

        self.individuals = self.individuals + population.individuals

    def calculate_population_fitness(self, fitness_function: Callable):

        """Calculates the fitnesses of all individuals.

        Arguments:
            fitness_function (Callable): The fitness function to be used.
        """

        for individual in self.individuals:
            individual.calculate_fitness(fitness_function)

    def rank_population(self, ranking_function: Callable):

        """Ranks and sorts the population according to a specified rule.

        Arguments:
            ranking_function (Callable): The function used to rank the individuals.
        """

        self.individuals = ranking_function(self.individuals) 