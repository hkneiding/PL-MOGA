import shutil
import pickle
import functools
import subprocess
import numpy as np
import pandas as pd
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

transition_metal_atomic_numbers = [21, 22, 23, 24, 25, 26, 27, 28, 29, 30,                          # first block
                                    39, 40, 41, 42, 43, 44, 45, 46, 47, 48,                          # second block
                                    57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,      # lanthanides
                                    72, 73, 74, 75, 76, 77, 78, 79, 80,                              # third block
                                    89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103,  # actinides
                                    104, 105, 106, 107, 108, 109, 110, 111, 112]                     # fourth block

def compose(*functions):
    return functools.reduce(lambda f, g: lambda x: f(g(x)), functions)

def load_data(src_dir: str):

    data_dict = {}

    with open(src_dir + 'ligands_xyz.xyz') as fh:
        raw_xyz = fh.read()

    xyzs = raw_xyz.split('\n\n')
    for xyz in xyzs:

        id = xyz.split('\n')[1]
        data_dict[id] = {'xyz': xyz}
    
    fingerprints = pd.read_csv(src_dir + 'ligands_fingerprints.csv', sep=';')
    for fingerprint in fingerprints.to_dict(orient='records'):

        data_dict[fingerprint['name']]['charge'] = fingerprint['charge']
        data_dict[fingerprint['name']]['n_metal_bound'] = fingerprint['n_metal_bound']
        data_dict[fingerprint['name']]['n_atoms'] = fingerprint['n_atoms']

    df = pd.DataFrame(data_dict).T

    return df

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

def zero_mask_target_by_population_average(individual, individuals, target_idx):

    """Masks one fitness target of the individual by zeroing below average of population.

    Returns:
        list: The masked fitness
    """

    target_population = []
    for _ in individuals:
        target_population.append(_._fitness[target_idx])

    target_population_average = np.median(target_population)

    print(target_population)
    print('Population average:', target_population_average)
    
    if individual._fitness[target_idx] < target_population_average:
        return [0, 0]

    return individual._fitness

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


def fitness_function(individual, key_mapping, charges):

    """Calculates the fitness of a organometallic compound using xTB based on quantum properties.

    Returns:
        float: The fitness.
    """

    # setup xTB runner
    xtb_runner = uxtbpy.XtbRunner(output_format='dict')

    # remove existing molsimplify run directory
    shutil.rmtree('./Runs/run', ignore_errors=True)
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
        with open('./Runs/run/run/run.xyz') as fh:
            xyz = fh.read()
            individual.meta['initial_xyz'] = xyz
    except FileNotFoundError:
        print('molSimplify Failed: genome: ' + ','.join([key_mapping[i] for i in individual.genome]))
        return [0,0]

    try:
        charge = calculate_total_charge(individual, charges)
        xtb_parameters = ['--opt tight --uhf 0 --norestart -v -c ' + str(charge)]
        result = xtb_runner.run_xtb_from_xyz(xyz, parameters=xtb_parameters)
        individual.meta['optimised_xyz'] = result['optimised_xyz']
    except RuntimeError:
        print('xTB Failed: genome: ' + ','.join([key_mapping[i] for i in individual.genome]))
        return [0,0]

    # check that connecting atoms are the same between initial and optimised xyz
    initial_connecting_indices = get_metal_connecting_indices(individual.meta['initial_xyz'], 2.5)
    optimised_connecting_indices = get_metal_connecting_indices(individual.meta['optimised_xyz'], 2.5)
    if initial_connecting_indices != optimised_connecting_indices:
        print('Connecting indices different')
        return [0,0]

    # check for disconnected graphs
    is_connected = radius_graph_is_connected(parse_xyz(individual.meta['optimised_xyz'])[1], 3.0)
    if not is_connected:
        print('Graph not connected')
        return [0,0]

    # stabilisation_energy = -result['energy'] / get_electron_count_from_xyz(result['optimised_xyz'], calculate_total_charge(individual, charges))

    return [result['polarisability'], result['homo_lumo_gap']]

def flatten_list(list: list):
    return [item for sublist in list for item in sublist]

def charge_range(individual, charges, allowed_charges):

    if calculate_total_charge(individual, charges) in allowed_charges:
        return True
    
    return False
                
if __name__ == "__main__":

    from gapy.individual import Individual
    from gapy.population import Population
    from gapy.mutation import uniform_integer_mutation, swap_mutation
    from gapy.crossover import uniform_crossover
    from gapy.selection import select_by_fitness, roulette_wheel_fitness, select_by_rank, roulette_wheel_rank
    from gapy.rank import rank_dominate, rank_is_dominated, rank_non_dominated_fronts, rank_dominate_by_feature, rank_is_dominated_by_feature
    from gapy.ga import GA

    # fix random seed
    np.random.seed(2023)

    # specify ligand space and charges
    ligands_names = [
        'RUCBEY-subgraph-1', 'WECJIA-subgraph-3', 'KEYRUB-subgraph-1', 'NURKEQ-subgraph-2', 'MEBXUN-subgraph-1',
        'BIFMOV-subgraph-1', 'CUJYEL-subgraph-2', 'EZEXEM-subgraph-1', 'FOMVUB-subgraph-2', 'EFIHEJ-subgraph-3',
        'LETTEL-subgraph-1', 'KAKKIR-subgraph-3', 'BICRIQ-subgraph-3', 'UPEGAZ-subgraph-2', 'CEVJAP-subgraph-2',
        'BABTUT-subgraph-3', 'ZEJJEF-subgraph-3', 'KULGAZ-subgraph-2', 'CIGDAA-subgraph-1', 'HOVMIP-subgraph-3',
        'ULUSIE-subgraph-1', 'IBEKUV-subgraph-1', 'REQSUD-subgraph-2', 'BOSJIF-subgraph-1', 'GUVMEP-subgraph-0',
        'MAZJIJ-subgraph-0', 'OBONEA-subgraph-1', 'CORTOU-subgraph-2', 'LEVGUO-subgraph-2', 'REBWEB-subgraph-2',
        'DOGPAS-subgraph-1', 'IJIMIX-subgraph-1', 'PEJGAN-subgraph-1', 'BIFZEX-subgraph-0', 'IRIXUC-subgraph-3',
        'SAYGOO-subgraph-0', 'UROGIS-subgraph-1', 'MAQKEX-subgraph-1', 'LUQWUQ-subgraph-1', 'QAYDID-subgraph-2',
        'MOYDOV-subgraph-3', 'NIZQUK-subgraph-1', 'SAYHIJ-subgraph-1', 'CIQGOY-subgraph-0', 'VUFZUT-subgraph-1',
        'ZOQFIU-subgraph-0', 'GUQBUQ-subgraph-0', 'LEZYUM-subgraph-2', 'RAJXUX-subgraph-2', 'QEWZOH-subgraph-3'
    ]
    ligands_charges = [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
    ]
    print('Using ' + str(len(ligands_names)) + ' ligands.')

    # GA parameters
    n_parents = 65
    n_offspring = 2 * n_parents
    n_population = 130

    sub_mutation_1 = functools.partial(uniform_integer_mutation, mutation_space=len(ligands_names), mutation_rate=0.5)
    sub_mutation_2 = functools.partial(swap_mutation, mutation_rate=0.5)

    ga = GA(fitness_function=functools.partial(fitness_function, key_mapping=ligands_names, charges=ligands_charges),
            parent_selection=functools.partial(roulette_wheel_rank, n_selected=n_parents, rank_function=rank_non_dominated_fronts),
            survivor_selection=functools.partial(select_by_rank, n_selected=n_population, rank_function=rank_is_dominated),
            crossover=functools.partial(uniform_crossover, mixing_ratio=0.5),
            mutation=compose(sub_mutation_2, sub_mutation_1),
            n_offspring=n_offspring,
            n_allowed_duplicates=0,
            solution_constraints=[functools.partial(charge_range, charges=ligands_charges, allowed_charges=[-1, 0, 1])],
            genome_equivalence_function=are_rotation_equivalents,
            masking_function=lambda x,y: None #functools.partial(zero_mask_target_by_population_average, target_idx=0)
    )

    # ga = GA(fitness_function=functools.partial(fitness_function, key_mapping=ligands_names, charges=ligands_charges),
    #         parent_selection=functools.partial(roulette_wheel_rank, n_selected=n_parents, rank_function=functools.partial(rank_is_dominated_by_feature, weighting=[1,0.5])),
    #         survivor_selection=functools.partial(select_by_rank, n_selected=n_population, rank_function=rank_is_dominated),
    #         crossover=functools.partial(uniform_crossover, mixing_ratio=0.5),
    #         mutation=functools.partial(uniform_integer_mutation, mutation_space=len(ligands_names), mutation_rate=0.5),
    #         n_offspring=n_offspring,
    #         n_allowed_duplicates=0,
    #         solution_constraints=[functools.partial(charge_range, charges=ligands_charges, allowed_charges=[-1, 0, 1])]
    # )

    # random initial population
    initial_individuals = []
    for i in range(n_population):

        neutral_choice = np.random.randint(0, high=25, size=2)
        anionic_choice = np.random.randint(25, high=50, size=2)

        genome = np.random.permutation(np.concatenate((neutral_choice, anionic_choice))).tolist()
        initial_individuals.append(Individual(genome=genome, meta={'metal_centre': 'Pd', 'oxidation_state': 2, 'coordination_geometry': 'sqp'}))
    initial_population = Population(initial_individuals)

    # run ga
    final_pop, log = ga.run(n_epochs=100, initial_population=initial_population)

    for i, individual in enumerate(final_pop.individuals):

        with open('.temp/mol-' + str(i) + '-msinit.xyz', 'w') as fh:
            fh.write(individual.meta['initial_xyz']) 
        with open('.temp/mol-' + str(i) + '-xtbopt.xyz', 'w') as fh:
            fh.write(individual.meta['optimised_xyz']) 

    with open('log.pickle', 'wb') as fh:
        pickle.dump(log, fh)
