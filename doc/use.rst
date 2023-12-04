Usage
=====

We use a configuration file (`config.yml`) to set the hyperparameters of the genetic algorithm. Then running the model is simply done with

.. code-block:: console

  $ python main.py config.yml

where the `main.py` script will automatically parse the information given in the configuration file.

.. note::

   Make sure you have activated the correct conda environment before running.

==========
Parameters
==========

The following will discuss all the different parameters that can be changed in the configuration file.

 - seed: The random seed to initialize on.
 - target_properties: List of the target properties to optimize for. This is based on xTB output and available options are: homo_energy, lumo_energy, homo_lumo_gap, energy, enthalpy_energy, free_energy, enthalpy, heat_capacity, entropy, dipole_moment, and polarisability.
 - zeromask_scaling_factors: List of scaling factors used in the Pareto Lighthouse algorithm. This list should have as many entries as there are properties to optimize for.
 - ligand_space: The ligand space to investigate. Here, we provide the different spaces discussed in the publication: 1M, 1B, and 1B_rand.
 - n_generations: The number of generations to run the MOGA for.
 - n_population: The population size.
 - n_parents: The number of parents.
 - n_offspring: The number of offspring individuals to generate each generation.
 - parent_selection: The parent selection mechanism. Possible choice are roulette_wheel_rank and select_by_rank.
 - parent_rank: The ranking funtion used in the parent selection. Possible choices are rank_dominate, rank_is_dominated, and rank_non_dominated_fronts.
 - survivor_selection: The survivor selection mechanism. Possible choice are roulette_wheel_rank and select_by_rank.
 - surivor_rank: The ranking funtion used in the survivor selection. Possible choices are rank_dominate, rank_is_dominated, and rank_non_dominated_fronts.
 - crossover: The crossover operation. Only uniform_crossover is implemented.
 - crossover_mixing: The mixing parameter for the crossover operation.
 - mutations: List of mutations that are applied separately. Possible choices are swap_mutation and uniform_integer_mutation.
 - mutation_rates: List of mutation rates for the different mutations specified.
 - n_allowed_duplicates: The number of allowed duplicate individuals in a population.
 - n_ligands: The number of ligand substitution sites of the TMC scaffold.
 - metal_center: The metal center element.
 - oxidation_state: The oxidation state of the metal center.
 - coordination_geometry: The coordination geometry. Possible choices are li (linear). tpl (trigonal planar), sqp (square planar), thd (tetrahedral), spy (square pyramidal), tbp (trigonal_bipyramidal), oct (octahedral), tpr (trigonal prismatic), pbp (pentagonal bipyramidal), sqapsqp (square antiprismatic)
 - allowed_charges: List of the allowed overall TMC charges.

All of these parameters have to be set in order for the code to function properly.

==============
Making changes
==============

Experienced Python users can also directly make changes to the code. For example you might want to add additional selection mechanisms, crossover and mutation operations, or ranking schemes to the `gapy` backbone. If you have implemented custom functions and want them to be available via the configuration file you need to modify the `parse_config_file` function in the `main.py` file to correctly parse the given strings into the actual functions.
