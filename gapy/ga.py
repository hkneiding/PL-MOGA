import numpy as np

from .population import Population


class GA:
        
    def __init__(self, fitness_function, parent_selection, survivor_selection, crossover, mutation, n_offspring, n_allowed_duplicates, solution_constraints):

        self.fitness_function = fitness_function

        self.parent_selection = parent_selection
        self.survivor_selection = survivor_selection

        self.crossover = crossover
        self.mutation = mutation

        self.n_offspring = n_offspring
        self.n_allowed_duplicates = n_allowed_duplicates

        self.solution_constraints = solution_constraints

    def run(self, initial_population: Population, n_epochs: int):

        population = initial_population

        # calculate fitness of initial population
        print('Calculating fitness of initial population..')
        population.calculate_population_fitness(self.fitness_function)

        # start GA iterations
        print('Starting iterations..')
        for epoch in range(1, n_epochs + 1):
            
            # get offspring
            offspring = self.get_offspring(population, self.n_offspring, n_allowed_duplicates=self.n_allowed_duplicates)

            # calculate fitness of all offspring individuals
            offspring.calculate_population_fitness(self.fitness_function)

            # extend existing population with offspring
            population.extend(offspring)

            # select survivors
            population.select(self.survivor_selection)

            print(epoch)
            # print('Epoch ' + str(epoch) + ' | Average fitness: ' + str(np.round(np.mean([individual['fitness'] for individual in population]), decimals=3)))
        
            for individual in population.individuals:
                print(individual.genome)
        return population

    def get_offspring(self, population: list, n_offspring: int, n_allowed_duplicates: int = None):

        # select parents
        parent_individuals = self.parent_selection(population.individuals)

        offspring = Population([])
        for i in range(n_offspring):
            
            # perform crossover and mutation until a valid candidate is found
            while(True):

                # generate new child through recombination of two random parents
                child = self.crossover(*np.random.choice(parent_individuals, size=2))
                # mutate offspring
                child = self.mutation(child)

                if n_allowed_duplicates is not None:

                    #  check if genome exists in current population and whether it exceeds allowed
                    # number of duplicates

                    n_duplicates = [individual.genome for individual in population.individuals].count(child.genome) + [individual.genome for individual in population.individuals].count(child.genome)

                    if n_duplicates > n_allowed_duplicates:
                        continue                

                # iterate through additional solution constraints
                for solution_constraint in self.solution_constraints:
                    # if they are not satisfied request new trial solution
                    if not solution_constraint(child):
                        continue

                break

            offspring.add_individual(child)
        
        return offspring