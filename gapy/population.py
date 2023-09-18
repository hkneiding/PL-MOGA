from typing import Callable
import multiprocessing as mp

from .individual import Individual


class Population:

    def __init__(self, individuals: list):

        """Constructor

        Arguments:
            individuals (list[Individual]): List of Individual objects.
        """

        self.individuals = individuals

    def as_dict(self):

        """Returns the class fields as a dict.

        Returns:
            dict: The dict of class fields.
        """

        return {
            'individuals': [_.as_dict() for _ in self.individuals]
        }

    def get_selected_individuals(self, selection_function: Callable, n_elitism: int = 0):
        
        """Selects individuals based on specified rule.

        Arguments:
            selection_function (Callable): The function used to select individuals.
            n_elitism (int): The number of elitism to apply.
        """

        # take out n_elitism individuals
        selected_individuals = self.individuals[:n_elitism]
        # perform selection on remaining
        selected_individuals += selection_function(self.individuals[n_elitism:])

        return selected_individuals

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
  
        with mp.Pool(int(mp.cpu_count() - 2)) as p:    
            individuals = p.map(fitness_function, self.individuals)
            p.close()
            p.join()

        self.individuals = individuals

    def get_masked_population_fitness(self, masking_function: Callable):
        
        """Gets the masked fitnesses of all individuals according to some masking function.

        Arguments:
            masking_function (Callable): The masking function to be used.
        """

        for individual in self.individuals:
            individual._masked_fitness = masking_function(individual, self.individuals)

    def get_dominating_individuals(self, query_individual):

        """Gets a list of individuals that a query individual is dominated by.

        Returns:
            list[Individuals]: The list of dominating individuals.
        """

        dominating_individuals = []
        for individual in self.individuals:
            if query_individual.is_dominated_by(individual):
                dominating_individuals.append(individual)

        return dominating_individuals

    def get_dominated_individuals(self, query_individual):

        """Gets a list of individuals that a query individual dominates.

        Returns:
            list[Individuals]: The list of dominated individuals.
        """

        dominated_individuals = []
        for individual in self.individuals:
            if query_individual.dominates(individual):
                dominated_individuals.append(individual)

        return dominated_individuals
