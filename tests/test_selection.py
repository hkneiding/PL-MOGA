import unittest
from parameterized import parameterized

from .utils import Utils
from gapy.selection import *


class TestTools(unittest.TestCase):

    # @parameterized.expand([

    #     [
    #         [
    #             Individual(genome=[], fitness=[1]),
    #             Individual(genome=[], fitness=[2]),
    #             Individual(genome=[], fitness=[3]),
    #             Individual(genome=[], fitness=[4])
    #         ],
    #         1,
    #         [
    #             Individual(genome=[], fitness=[4]),
    #         ]
    #     ],

    #     [
    #         [
    #             Individual(genome=[], fitness=[4]),
    #             Individual(genome=[], fitness=[4]),
    #             Individual(genome=[], fitness=[4]),
    #             Individual(genome=[], fitness=[4])
    #         ],
    #         1,
    #         [
    #             Individual(genome=[], fitness=[4]),
    #         ]
    #     ],

    #     [
    #         [
    #             Individual(genome=[], fitness=[-4]),
    #             Individual(genome=[], fitness=[-4]),
    #             Individual(genome=[], fitness=[-4]),
    #             Individual(genome=[], fitness=[0])
    #         ],
    #         2,
    #         [
    #             Individual(genome=[], fitness=[0]),
    #             Individual(genome=[], fitness=[-4]),
    #         ]
    #     ],

    #     [
    #         [
    #             Individual(genome=[], fitness=[1, 2]),
    #             Individual(genome=[], fitness=[2, 7]),
    #             Individual(genome=[], fitness=[3, 5]),
    #             Individual(genome=[], fitness=[4, 2])
    #         ],
    #         2,
    #         [
    #             Individual(genome=[], fitness=[2, 7]),
    #             Individual(genome=[], fitness=[3, 5]),
    #         ]
    #     ],

    # ])
    # def test_select_by_fitness(self, individuals, n_selected, expected):

    #     result = select_by_fitness(individuals, n_selected)

    #     Utils.assert_are_almost_equal(result, expected)

    # @parameterized.expand([

    #     [
    #         [
    #             Individual(genome=[], fitness=[1]),
    #             Individual(genome=[], fitness=[2]),
    #             Individual(genome=[], fitness=[3]),
    #             Individual(genome=[], fitness=[4])
    #         ],
    #         1,
    #         [
    #             Individual(genome=[], fitness=[4]),
    #         ]
    #     ],

    #     [
    #         [
    #             Individual(genome=[], fitness=[4]),
    #             Individual(genome=[], fitness=[4]),
    #             Individual(genome=[], fitness=[4]),
    #             Individual(genome=[], fitness=[4])
    #         ],
    #         1,
    #         [
    #             Individual(genome=[], fitness=[4]),
    #         ]
    #     ],

    #     [
    #         [
    #             Individual(genome=[], fitness=[-4]),
    #             Individual(genome=[], fitness=[-4]),
    #             Individual(genome=[], fitness=[-4]),
    #             Individual(genome=[], fitness=[0])
    #         ],
    #         2,
    #         [
    #             Individual(genome=[], fitness=[0]),
    #             Individual(genome=[], fitness=[-4]),
    #         ]
    #     ],

    #             [
    #         [
    #             Individual(genome=[], fitness=[0, 0]),
    #             Individual(genome=[], fitness=[1, 0]),
    #             Individual(genome=[], fitness=[1, 1]),
    #             Individual(genome=[], fitness=[2, 0])
    #         ],
    #         2,
    #         [
    #             Individual(genome=[], fitness=[1, 1]),
    #             Individual(genome=[], fitness=[2, 0]),
    #         ]
    #     ],

    # ])
    # def test_pareto_domination_rank(self, individuals, n_selected, expected):

    #     result = pareto_domination_rank(individuals, n_selected)

    #     Utils.assert_are_almost_equal(result, expected)

    pass