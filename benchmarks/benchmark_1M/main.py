import os
import sys
import yaml
import datetime
import shutil
import pickle
import functools
import subprocess
import numpy as np
import pandas as pd
from contextlib import contextmanager
from scipy.sparse.csgraph import connected_components

import uxtbpy

element_identifiers = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                        'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',
                        'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',
                        'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
                        'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',
                        'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs',
                        'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
                        'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
                        'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb',
                        'Bi', 'Po', 'At', 'Rn']

transition_metal_atomic_numbers = [21, 22, 23, 24, 25, 26, 27, 28, 29, 30,                           # first block
                                    39, 40, 41, 42, 43, 44, 45, 46, 47, 48,                          # second block
                                    57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,      # lanthanides
                                    72, 73, 74, 75, 76, 77, 78, 79, 80,                              # third block
                                    89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103,  # actinides
                                    104, 105, 106, 107, 108, 109, 110, 111, 112]                     # fourth block

@contextmanager
def change_directory(destination: str):
    try:
        cwd = os.getcwd()
        os.chdir(destination)
        yield
    finally:
        os.chdir(cwd)

def compose(*functions):

    """Composes a given list of functions.

    Returns:
        callable: The composed function.
    """

    return functools.reduce(lambda f, g: lambda x: f(g(x)), functions)

def flatten_list(list: list):

    """Flattens a given list.

    Returns:
        list: The flattened list.
    """

    return [item for sublist in list for item in sublist]

def charge_range(individual, charges, allowed_charges):

    """Checks if a given individual is within the allowed charge range.

    Returns:
        bool: The flag indicating whether the individual is within the allowed charge range.
    """

    if calculate_total_charge(individual, charges) in allowed_charges:
        return True
    
    return False

def are_rotation_equivalents(l1: list, l2: list):

    """Checks if two lists are rotational equivalents of each other.

    Arguments:
        l1 (list): The first list.
        l2 (list): The second list.

    Returns:
        bool: The flag indicating whether or not the lists are rotation equivalent.
    """

    if len(l1) != len(l2):
        return False

    l1_extended = l1 * 2
    for i in range(len(l1_extended) - len(l1)):
        if l1_extended[i:i + len(l1)] == l2:
            return True
    
    return False

def zero_mask_target_by_population_median(individual, individuals, target_indices, scaling=1):

    """Masks one fitness target of the individual by zeroing below median of population.

    Arguments:
        individual (Individual): The current individual.
        individuals (list[Individual]: The list of individuals.
        target_indices (list[int]): The target indices to zero mask.
        scaling: (list[float]): The scalings to apply to each of the target indices.

    Returns:
        list: The masked fitness
    """

    for target_idx in target_indices:
      
        target_population = []
        for _ in individuals:
            target_population.append(_._fitness[target_idx])
        
        if individual._fitness[target_idx] < scaling[target_idx] * np.median(target_population):
            return [0, 0]

    return individual._fitness

def weight_mask(individual, individuals, weights):

    individual_masked_fitness = [None, None]

    if individual._fitness == [0,0]:
        return individual._fitness

    # min max scaling of properties
    # obtained from the calculated ground truth data
    individual_masked_fitness[0] = (individual._fitness[0] - 55.632277) / (493.76163 - 55.632277)
    individual_masked_fitness[1] = (individual._fitness[1] - 0.044745372489) / (4.949142873462 - 0.044745372489)

    for i in range(len(individual_masked_fitness)):
        individual_masked_fitness[i] = weights[i] * individual_masked_fitness[i]

    return individual_masked_fitness

def parse_xyz(xyz: str):

    """Parses a given xyz into two lists for atoms and positions respectively.

    Returns:
        list[str]: The list of atom identifiers.
        list[list[float]]: The list of atomic positions.
    """

    atoms = []
    positions = []

    lines = xyz.split('\n')
    for i in range(2, len(lines), 1):

        line_split = lines[i].split()

        if len(line_split) != 4:
            break

        atoms.append(line_split[0])
        positions.append([float(line_split[i]) for i in [1, 2, 3]])

    return atoms, positions

def get_metal_connecting_indices(xyz: str, radius_cutoff: float):

    """Gets the indices of atoms connecting to metal based on a cutoff value.

    Returns:
        list[int]: The list of indices connecting to metal.
    """

    atoms, positions = parse_xyz(xyz)

    metal_indices = []
    for i, atom in enumerate(atoms):
        if (element_identifiers.index(atom) + 1) in transition_metal_atomic_numbers:
            metal_indices.append(i)

    connecting_indices = []
    for metal_index in metal_indices:

        for i, position in enumerate(positions):

            if metal_index == i:
                continue

            distance = np.linalg.norm(np.array(position) - np.array(positions[metal_index]))
            if distance <= radius_cutoff:
                connecting_indices.append(i)

    return connecting_indices

def get_radius_adjacency_matrix(positions: list, radius_cutoff: float):

    """Gets an adjacency matrix based on a radius cutoff.

    Returns:
        list[list[int]]: The adjacency matrix.
    """

    adjacency_matrix = np.zeros((len(positions), len(positions)))

    for i in range(len(positions)):
        for j in range(i+1, len(positions), 1):
            
            if np.linalg.norm(np.array(positions[i]) - np.array(positions[j])) <= radius_cutoff:
                adjacency_matrix[i,j] = 1
                adjacency_matrix[j,i] = 1

    return adjacency_matrix

def radius_graph_is_connected(positions: list, radius_cutoff: float):

    """Checks if a radius graph based on a cutoff is connected or not.

    Returns:
        bool: The flag saying whether the graph is connected or not.
    """

    adj = get_radius_adjacency_matrix(positions, radius_cutoff)
    n_connected_components = connected_components(adj)[0]

    if n_connected_components == 1:
        return True
    return False

def get_electron_count_from_xyz(xyz: str, charge: int = 0):

    """Get the number of electrons from a xyz file.

    Returns:
        int: The electron count.
    """

    lines = xyz.split('\n')

    n_electrons = 0
    for line in lines:

        line_split = line.split()
        if len(line_split) == 4:
            n_electrons += element_identifiers.index(line_split[0]) + 1

    return n_electrons - charge

def calculate_total_charge(individual, charges: list):

    """Calculates the total charge of an individual based on the metal oxidation state
    and the charges of ligands.

    Returns:
        int: The total charge.
    """

    return individual.meta['oxidation_state'] + sum([int(charges[_]) for _ in individual.genome])

def fitness_function(individual, target_properties, key_mapping, charges):

    """Calculates the fitness of a organometallic compound using xTB based on quantum properties.

    Returns:
        float: The fitness.
    """

    # make unique run directory
    tmp_dir = str(id(individual)) + individual.meta['creation_date'].replace(' ', '') + '/'
    os.mkdir(tmp_dir)

    with change_directory(tmp_dir):
        # prepare molsimplify parameters
        parameters = ['-skipANN True',
                    '-core ' + individual.meta['metal_centre'],
                    '-geometry ' + individual.meta['coordination_geometry'],
                    '-lig ' + ','.join([key_mapping[idx] for idx in individual.genome]),
                    '-coord ' + str(len(individual.genome)),
                    '-ligocc ' + ','.join(['1' for _ in individual.genome]),
                    '-name run',
                    '-oxstate ' + str(individual.meta['oxidation_state'])
        ]
        # run molsimplify
        result = subprocess.run(['molsimplify', *parameters], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    try:
        # read the xyz from molsimplify output
        with open(tmp_dir + '/Runs/run/run/run.xyz') as fh:
            xyz = fh.read()
            individual.meta['initial_xyz'] = xyz

        # setup xTB runner
        xtb_runner = uxtbpy.XtbRunner(xtb_directory=tmp_dir, output_format='dict')

        charge = calculate_total_charge(individual, charges)
        xtb_parameters = ['--opt tight --uhf 0 --norestart -v -c ' + str(charge)]
        result = xtb_runner.run_xtb_from_xyz(xyz, parameters=xtb_parameters)
        individual.meta['optimised_xyz'] = result['optimised_xyz']

    except FileNotFoundError:
        print('molSimplify failure. Genome: ' + ','.join([key_mapping[i] for i in individual.genome]))
        individual.set_fitness([0,0])
        return individual

    except RuntimeError:
        print('xTB failure. Genome: ' + ','.join([key_mapping[i] for i in individual.genome]))
        individual.set_fitness([0,0])
        return individual

    except Exception:
        print('Other error')
        print(charge)
        print(xyz)
        individual.set_fitness([0,0])
        return individual

    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        
    # check that connecting atoms are the same between initial and optimised xyz
    initial_connecting_indices = get_metal_connecting_indices(individual.meta['initial_xyz'], 2.5)
    optimised_connecting_indices = get_metal_connecting_indices(individual.meta['optimised_xyz'], 2.5)
    if initial_connecting_indices != optimised_connecting_indices:
        print('Connecting indices different')
        individual.set_fitness([0,0])
        return individual

    # check for disconnected graphs
    is_connected = radius_graph_is_connected(parse_xyz(individual.meta['optimised_xyz'])[1], 3.0)
    if not is_connected:
        print('Graph not connected')
        individual.set_fitness([0,0])
        return individual

    # update fitness and return
    fitness_vector = []
    for target_property in target_properties:
        fitness_vector.append(result[target_property])
    individual.set_fitness(fitness_vector)

    return individual 

def parse_config_file(file_path: str):

    """Parses a config file and returns a dictionary with all relevant parameters.

    Returns:
        dict: The parameter dictionary.
    """

    with open(file_path, 'r') as fh:
        config = yaml.safe_load(fh)

    # specify ligand space and charges
    config['ligands_names'], config['ligands_charges'] = get_ligand_names_and_charges(config['ligand_space'])

    if config['parent_selection'] == 'roulette_wheel_rank':
        config['parent_selection'] = roulette_wheel_rank
    elif config['parent_selection'] == 'select_by_rank':
        config['parent_selection'] = select_by_rank
    elif config['parent_selection'] == 'roulette_wheel_fitness':
        config['parent_selection'] = roulette_wheel_fitness
    elif config['parent_selection'] == 'select_by_fitness':
        config['parent_selection'] = select_by_fitness
    else:
        print('Config not found:', config['parent_selection'])

    if config['parent_rank'] == 'rank_dominate':
        config['parent_rank'] = rank_dominate
    elif config['parent_rank'] == 'rank_is_dominated':
        config['parent_rank'] = rank_is_dominated
    elif config['parent_rank'] == 'rank_non_dominated_fronts':
        config['parent_rank'] = rank_non_dominated_fronts
    else:
        print('Config not found:', config['parent_rank'])

    if config['survivor_selection'] == 'roulette_wheel_rank':
        config['survivor_selection'] = roulette_wheel_rank
    elif config['survivor_selection'] == 'select_by_rank':
        config['survivor_selection'] = select_by_rank
    elif config['survivor_selection'] == 'roulette_wheel_fitness':
        config['survivor_selection'] = roulette_wheel_fitness
    elif config['survivor_selection'] == 'select_by_fitness':
        config['survivor_selection'] = select_by_fitness
    else:
        print('Config not found:', config['survivor_selection'])

    if config['survivor_rank'] == 'rank_dominate':
        config['survivor_rank'] = rank_dominate
    elif config['survivor_rank'] == 'rank_is_dominated':
        config['survivor_rank'] = rank_is_dominated
    elif config['survivor_rank'] == 'rank_non_dominated_fronts':
        config['survivor_rank'] = rank_non_dominated_fronts
    else:
        print('Config not found:', config['survivor_rank'])

    if config['crossover'] == 'uniform_crossover':
        config['crossover'] = uniform_crossover
    else:
        print('Config not found:', config['crossover'])

    for i, _ in enumerate(config['mutations']):
        if config['mutations'][i] == 'swap_mutation':
            config['mutations'][i] = functools.partial(
                swap_mutation,
                mutation_rate=config['mutation_rates'][i])
        elif config['mutations'][i] == 'uniform_integer_mutation':
            config['mutations'][i] = functools.partial(
                uniform_integer_mutation, 
                mutation_space=len(config['ligands_names']), 
                mutation_rate=config['mutation_rates'][i]
            )
        else:
            print('Config not found:', config['mutations'][i])

    # compose mutations
    config['composed_mutation'] = config['mutations'][0]
    for i in range(1, len(config['mutations']), 1):
        config['composed_mutation'] = compose(config['composed_mutation'], config['mutations'][i])

    return config
    

if __name__ == "__main__":

    from gapy.individual import Individual
    from gapy.population import Population
    from gapy.mutation import uniform_integer_mutation, swap_mutation
    from gapy.crossover import uniform_crossover
    from gapy.selection import select_by_rank, roulette_wheel_rank, roulette_wheel_fitness, select_by_fitness
    from gapy.rank import rank_dominate, rank_is_dominated, rank_non_dominated_fronts
    from gapy.ga import GA

    from ligands_info import get_ligand_names_and_charges

    config = parse_config_file(sys.argv[1])

    # fix random seed
    np.random.seed(config['seed'])

    print('Using ' + str(len(config['ligands_names'])) + ' ligands.')

    # set up GA
    ga = GA(fitness_function=functools.partial(fitness_function, target_properties=config['target_properties'], key_mapping=config['ligands_names'], charges=config['ligands_charges']),
            parent_selection=functools.partial(config['parent_selection'], n_selected=config['n_parents']),
            survivor_selection=functools.partial(config['survivor_selection'], n_selected=config['n_population']),
            crossover=functools.partial(config['crossover'], mixing_ratio=config['crossover_mixing']),
            mutation=config['composed_mutation'],
            n_offspring=config['n_offspring'],
            n_allowed_duplicates=config['n_allowed_duplicates'],
            solution_constraints=[functools.partial(charge_range, charges=config['ligands_charges'], allowed_charges=config['allowed_charges'])],
            genome_equivalence_function=are_rotation_equivalents,
            masking_function=functools.partial(weight_mask, weights=[1,1])
    )

    # random initial population
    initial_individuals = []
    for i in range(config['n_population']):

        neutral_choice = np.random.randint(
            0, 
            high=len(config['ligands_names']) // 2,
            size=config['n_ligands'] // 2
        )
        anionic_choice = np.random.randint(
            len(config['ligands_names']) // 2, 
            high=len(config['ligands_names']), 
            size=config['n_ligands'] - (config['n_ligands'] // 2)
        )

        genome = np.random.permutation(np.concatenate((neutral_choice, anionic_choice))).tolist()
        initial_individuals.append(Individual(genome=genome, meta={
            'metal_centre': 'Pd',
            'oxidation_state': 2,
            'coordination_geometry': 'sqp',
            'creation_date': str(datetime.datetime.now())}
        ))
    initial_population = Population(initial_individuals)
    
    # run ga
    final_pop, log = ga.run(n_epochs=config['n_generations'], initial_population=initial_population)

    # save log
    with open('log.pickle', 'wb') as fh:
        pickle.dump(log, fh)
