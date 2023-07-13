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
        self._masked_fitness = None
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

        if self._masked_fitness is not None:
            return self._masked_fitness

        return self._fitness

    @property
    def fitness_sum(self):

        if self._masked_fitness is not None:
            return sum(self._masked_fitness)
        
        return sum(self._fitness)

    def set_fitness(self, fitness):

        """Sets the fitness of the individual.

        Arguments:
            fitness (list): The fitness.
        """

        self._fitness = fitness
        self._masked_fitness = None

    def calculate_fitness(self, fitness_function: Callable, force: bool = False):
        
        """Calculates and updates the fitness of the individual.

        Arguments:
            fitness_function (Callable): The fitness function to be used.
            force (bool): Flag used to force fitness calculate even when fitness is defined.
        """

        # check if 
        if self._fitness is None or force is True:
            self._fitness = fitness_function(self)

    def get_dominating_features(self, individual):

        """Determines which features of this instance are dominating compared to another indiviudal.
        
        Arguments:
            individual (Individual): The individual to compare with.

        Returns:
            list[bool]: List of flags to determine which features are dominating.
        """

        dominating_features = []        
        for i in range(len(self.fitness)):
            if self.fitness[i] <= individual.fitness[i]:
                dominating_features.append(False)
            else:
                dominating_features.append(True)

        return dominating_features

    def get_dominated_features(self, individual):

        """Determines which features of this instance are dominated compared to another indiviudal.
        
        Arguments:
            individual (Individual): The individual to compare with.

        Returns:
            list[bool]: List of flags to determine which features are dominated.
        """

        dominated_features = []        
        for i in range(len(self.fitness)):
            if self.fitness[i] >= individual.fitness[i]:
                dominated_features.append(False)
            else:
                dominated_features.append(True)

        return dominated_features

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