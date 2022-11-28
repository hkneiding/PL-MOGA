import shutil
import functools
import subprocess
import numpy as np
import pandas as pd

import uxtbpy
from individual import Individual
from population import Population
from mutation import uniform_integer_mutation
from crossover import uniform_crossover
from selection import select_by_fitness, roulette_wheel_fitness
from ga import GA


def load_csv(file_path: str):

    dtypes = {
        'xyz': str,
        'smiles': str,
        'stoichiometry': str,
        'occurrence': int,
        'is_alternative_charge': int,
        'charge': int
    }

    converters={
        'metal_bond_node_idx_groups': eval,
        'parent_metal_counts': eval,
        'csd_origin': eval,
    }

    ligands = pd.read_csv(file_path, sep=',', dtype=dtypes, converters=converters)

    return ligands

def get_electron_count_from_xyz(xyz: str, charge: int = 0):

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


def fitness_function(individual, key_mapping, xyz, charges):

    """Calculates the fitness of a organometallic compound using xTB based on quantum properties.

    Returns:
        float: The fitness.
    """

    # setup xTB runner
    xtb_runner = uxtbpy.XtbRunner(output_format='dict')

    # calculate metal center
    try:
        xtb_parameters = ['--norestart -v -c ' + str(individual.meta['oxidation_state'])]
        result = xtb_runner.run_xtb_from_xyz('1\n\n' + individual.meta['metal_centre'] + ' 0.0 0.0 0.0', parameters=xtb_parameters)
    except RuntimeError:
        print('xTB Failed: metal center')
        return [0]

    building_block_energy = result['energy']
    # calculate stability
    for allel in individual.genome:
        
        try:
            xtb_parameters = ['--opt --norestart -v -c ' + str(charges[allel])]
            result = xtb_runner.run_xtb_from_xyz(xyz[allel], parameters=xtb_parameters)
            building_block_energy += result['energy']
        except RuntimeError:
            print('xTB Failed: ligand ' + key_mapping[allel])
            return [0]
    
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
        return [0]

    try:
        xtb_parameters = ['--opt tight --uhf 0 --norestart -v -c ' + str(calculate_total_charge(individual, charges))]
        result = xtb_runner.run_xtb_from_xyz(xyz, parameters=xtb_parameters)
        individual.meta['optimised_xyz'] = result['optimised_xyz']
    except RuntimeError:
        print('xTB Failed: genome: ' + ','.join([key_mapping[i] for i in individual.genome]))
        return [0]

    stabilisation_energy = np.exp(- building_block_energy + result['energy'])
    # return fitness vector
    return [stabilisation_energy]
    # return [result['homo_lumo_gap']]

def flatten_list(list: list):
    return [item for sublist in list for item in sublist]

def charge_range(individual, charges, allowed_charges):

    if calculate_total_charge(individual, charges) in allowed_charges:
        return True
    
    return False
                
if __name__ == "__main__":

    # load ligand library
    ligands_data = load_csv('/home/hkneiding/Documents/UiO/ligands-test/ligands.csv')

    #(len(flatten_list(ligands_data['metal_bond_node_idx_groups'])) == 1)
    # get selection of ligands to allow

    ligands_selection = ligands_data[ligands_data['metal_bond_node_idx_groups'].apply(lambda x: len(flatten_list(x)) == 1)]

    ligands_selection = ligands_selection[(ligands_selection['xyz'].str.split('\n').str.len() <= 15)]

    ligands_selection = ligands_selection[(ligands_selection['charge'] >= -2)]

    # build list of ligand names for mapping into integers (can add 'x' for empty coordination site)
    ligands_names = ligands_selection['name'].tolist()
    # build list of ligand charges (can add 0 for empty coordination site)
    ligands_charges = ligands_selection['charge'].tolist()
    # build list of ligand xyzs
    ligands_xyz = ligands_selection['xyz'].tolist()


    print('Using ' + str(len(ligands_names)) + ' ligands.')

    n_parents = 5
    n_population = 5

    ga = GA(fitness_function=functools.partial(fitness_function, key_mapping=ligands_names, xyz=ligands_xyz, charges=ligands_charges),
            parent_selection=functools.partial(roulette_wheel_fitness, n_selected=n_parents),
            survivor_selection=functools.partial(select_by_fitness, n_selected=n_population),
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

    final_pop = ga.run(n_epochs=5, initial_population=initial_population)

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
