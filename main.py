import shutil
import functools
import subprocess
import numpy as np
import pandas as pd
from scipy.sparse.csgraph import connected_components

import uxtbpy
from gapy.individual import Individual
from gapy.population import Population
from gapy.mutation import uniform_integer_mutation
from gapy.crossover import uniform_crossover
from gapy.selection import select_by_fitness, roulette_wheel_fitness, pareto_domination_rank
from gapy.ga import GA

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

    lines = xyz.split('\n')

    n_electrons = 0
    for line in lines:

        line_split = line.split()
        if len(line_split) == 4:
            n_electrons += element_identifiers.index(line_split[0]) + 1

    return n_electrons - charge

def calculate_total_charge(individual: Individual, charges: list):

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
    print(1)
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

    print(xyz)
    try:
        xtb_parameters = ['--opt normal --uhf 0 --norestart -v -c ' + str(calculate_total_charge(individual, charges))]
        result = xtb_runner.run_xtb_from_xyz(xyz, parameters=xtb_parameters)
        individual.meta['optimised_xyz'] = result['optimised_xyz']
    except RuntimeError:
        print('xTB Failed: genome: ' + ','.join([key_mapping[i] for i in individual.genome]))
        return [0,0]

    # check that connecting atoms are the same between initial and optimised xyz
    initial_connecting_indices = get_metal_connecting_indices(individual.meta['initial_xyz'], 3.0)
    optimised_connecting_indices = get_metal_connecting_indices(individual.meta['optimised_xyz'], 3.0)
    if initial_connecting_indices != optimised_connecting_indices:
        return [0,0]

    # check for disconnected graphs
    is_connected = radius_graph_is_connected(parse_xyz(individual.meta['optimised_xyz'])[1], 3.0)
    if not is_connected:
        return [0,0]

    # stabilisation_energy = -result['energy'] / get_electron_count_from_xyz(result['optimised_xyz'], calculate_total_charge(individual, charges))

    return [np.exp(result['polarisability']), np.exp(result['homo_lumo_gap'])]

def flatten_list(list: list):
    return [item for sublist in list for item in sublist]

def charge_range(individual, charges, allowed_charges):

    if calculate_total_charge(individual, charges) in allowed_charges:
        return True
    
    return False
                
if __name__ == "__main__":

    # load ligand library
    ligands_fingerprints = pd.read_csv('/home/hkneiding/Documents/UiO/ligands-test/ligands/ligands_fingerprints.csv', sep=';')
    # get selection of ligands to allow
    ligands_selection = ligands_fingerprints[ligands_fingerprints['n_metal_bound'] == 1]
    ligands_selection = ligands_selection[ligands_selection['n_atoms'] <= 15]
    ligands_selection = ligands_selection[ligands_selection['charge'] >= -2]

    # build list of ligand names for mapping into integers (can add 'x' for empty coordination site)
    ligands_names = ligands_selection['name'].tolist()
    # build list of ligand charges (can add 0 for empty coordination site)
    ligands_charges = ligands_selection['charge'].tolist()

    print('Using ' + str(len(ligands_names)) + ' ligands.')

    n_parents = 5
    n_population = 5

    ga = GA(fitness_function=functools.partial(fitness_function, key_mapping=ligands_names, charges=ligands_charges),
            parent_selection=functools.partial(roulette_wheel_fitness, n_selected=n_parents),
            survivor_selection=functools.partial(pareto_domination_rank, n_selected=n_population),
            crossover=functools.partial(uniform_crossover, mixing_ratio=0.5),
            mutation=functools.partial(uniform_integer_mutation, mutation_space=len(ligands_names), mutation_rate=0.5),
            n_offspring=5,
            n_allowed_duplicates=0,
            solution_constraints=[functools.partial(charge_range, charges=ligands_charges, allowed_charges=[0])]
    )

    # palladium 2 coord 4 square planar
    # in xtb: oxidation state of palladium? 2
    # final charge should be in range -1,0,+1

    initial_individuals = [
        Individual([0,0,0,0], meta={'metal_centre': 'Pd', 'oxidation_state': 2, 'coordination_geometry': 'sqp'}),
        Individual([1,1,1,1], meta={'metal_centre': 'Pd', 'oxidation_state': 2, 'coordination_geometry': 'sqp'}),
        Individual([2,2,2,2], meta={'metal_centre': 'Pd', 'oxidation_state': 2, 'coordination_geometry': 'sqp'}),
        Individual([3,3,3,3], meta={'metal_centre': 'Pd', 'oxidation_state': 2, 'coordination_geometry': 'sqp'}),
        Individual([4,4,4,4], meta={'metal_centre': 'Pd', 'oxidation_state': 2, 'coordination_geometry': 'sqp'})
    ]
    initial_population = Population(initial_individuals)

    final_pop = ga.run(n_epochs=10, initial_population=initial_population)

    for i, individual in enumerate(final_pop.individuals):

        # for a in [ligands_names[idx] for idx in individual.genome]:
        #     print(a)

        # print(individual.meta['initial_xyz'])
        # print(individual.meta['optimised_xyz'])
        # print('')

        with open('.temp/mol-' + str(i) + '-msinit.xyz', 'w') as fh:
            fh.write(individual.meta['initial_xyz']) 
        with open('.temp/mol-' + str(i) + '-xtbopt.xyz', 'w') as fh:
            fh.write(individual.meta['optimised_xyz']) 
