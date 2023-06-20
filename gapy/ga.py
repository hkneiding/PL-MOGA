import copy
import numpy as np

from .population import Population


class GA:
        
    def __init__(self, fitness_function, parent_selection, survivor_selection, crossover, mutation, n_offspring, n_allowed_duplicates, solution_constraints, genome_equivalence_function):

        """Constructor for the GA class.

        Arguments:
            n_allowed_duplicates (int): The number of allowed duplicates in the population.
        """ 

        self.fitness_function = fitness_function

        self.parent_selection = parent_selection
        self.survivor_selection = survivor_selection

        self.crossover = crossover
        self.mutation = mutation

        self.n_offspring = n_offspring
        self.n_allowed_duplicates = n_allowed_duplicates

        self.solution_constraints = solution_constraints

        self.genome_equivalence_function = genome_equivalence_function

    def run(self, n_epochs: int, initial_population: Population = None):

        """Runs the genetic algorithm with the specified parameters and functions.

        Arguments:
            n_epochs (int): The max number of epochs to run before terminatings.
            initial_population (Population): The initial population.

        Returns:
            population: The final population.
            dict: Dictionary containing the history of populations and all investigated individuals.
        """

        if initial_population is not None:
            population = initial_population
        else:
            population = self._initialise_population()

        # set up history of generations
        log = {
            'history': [population.as_dict()],
            'explored_individuals': [_.as_dict() for _ in population.individuals]    
        }

        # calculate fitness of initial population
        print('Calculating fitness of initial population..')
        population.calculate_population_fitness(self.fitness_function)

        # start GA iterations
        print('Starting iterations..')
        for epoch in range(1, n_epochs + 1):
            
            # get offspring
            offspring = self._get_offspring(population, self.n_offspring)

            # calculate fitness of all offspring individuals
            offspring.calculate_population_fitness(self.fitness_function)

            # extend existing population with offspring
            population.extend(offspring)

            # select survivors
            population.select(self.survivor_selection)

            # update log
            log['history'].append(population.as_dict())
            log['explored_individuals'].extend([_.as_dict() for _ in offspring.individuals])
            
            print('Epoch ' + str(epoch) + ' | Average fitness: ' + str(np.round(np.mean([individual.fitness for individual in population.individuals], axis=0), decimals=3)))

        return population, log

    def _get_offspring(self, population: list, n_offspring: int):

        """Generates an offspring population.

        Arguments:
            population (Population): The current population.
            n_offspring (int): The number of offspring individuals to generate.

        Returns:
            Population: The generated offspring population.
        """

        # select parents
        parent_individuals = population.get_selected_individuals(self.parent_selection)

        offspring = Population([])
        for i in range(n_offspring):
            
            # perform crossover and mutation until a valid candidate is found
            while(True):

                # generate new child through recombination of two random parents
                child = self.crossover(*np.random.choice(parent_individuals, size=2))
                # mutate offspring
                child = self.mutation(child)

                if self.n_allowed_duplicates is not None:

                    # check if genome exists in current population and whether it exceeds allowed
                    # number of duplicates

                    n_duplicates = 0
                    for individual in population.individuals:

                        if self.genome_equivalence_function(child.genome, individual.genome):
                            n_duplicates += 0

                    if n_duplicates > self.n_allowed_duplicates:
                        continue                

                are_constraints_satified = True
                # iterate through additional solution constraints
                for solution_constraint in self.solution_constraints:
                    # if they are not satisfied request new trial solution
                    if not solution_constraint(child):
                        are_constraints_satified = False
                        break
                
                if are_constraints_satified:
                    break

            offspring.add_individual(child)
        
        return offspring

    def _initialise_population(self):

        pass
