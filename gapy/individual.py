from typing import Callable


class Individual:

    def __init__(self, genome: list, fitness=None, meta={}):

        """Constructor

        Arguments:
            genome (list): List representation of the individuals genome.
            fitness (float): The fitness value.
            meta (dict): Additional meta data for the individual.
        """

        self.genome = genome
        self._fitness = fitness
        self.meta = meta
    
    def as_dict(self):

        """Returns the class fields as a dict.

        Returns:
            dict: The dict of class fields.
        """

        return {
            'genome': self.genome,
            'fitness': self._fitness,
            'meta': self.meta
        }

    @property
    def fitness(self):
        return self._fitness

    @property
    def fitness_sum(self):
        return sum(self._fitness)

    def calculate_fitness(self, fitness_function: Callable, force: bool = False):
        
        """Calculates and updates the fitness of the individual.

        Arguments:
            fitness_function (Callable): The fitness function to be used.
            force (bool): Flag used to force fitness calculate even when fitness is defined.
        """

        # check if 
        if self._fitness is None or force is True:
            self._fitness = fitness_function(self)

    def dominates(self, individual):

        """Determines whether this instance dominates another individual.

        Arguments:
            individual (Individual): The individual to compare with.

        Returns:
            bool: The flag denoting whether it dominates or not.
        """

        if self.fitness == individual.fitness:
            return False

        for i in range(len(self.fitness)):
            if self.fitness[i] < individual.fitness[i]:
                return False

        return True

    def is_dominated_by(self, individual):

        """Determines whether this instance is dominated by another individual.

        Arguments:
            individual (Individual): The individual to compare with.

        Returns:
            bool: The flag denoting whether it is dominated or not.
        """

        if self.fitness == individual.fitness:
            return False

        for i in range(len(self.fitness)):
            if self.fitness[i] > individual.fitness[i]:
                return False

        return True