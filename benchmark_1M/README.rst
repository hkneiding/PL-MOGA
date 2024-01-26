=================================
MOGA benchmark on the 1.37M space
=================================

This directory contains the code and data associated to the benchmarking of the Pareto-Lighthouse MOGA (PL-MOGA) and the Weighted-sum MOGA (WS-MOGA) on the 1.37M space. In particular, the files ``main.py`` and ``config.yml`` contain the modified code and configuration file to run the WS-MOGA while the file ``ground_truth_fitness_values.csv`` holds the property labels (polarizability and HOMO-LUMO gap) for each unique member in the investigated space. The folders ``pl-moga`` and ``ws-moga`` hold the population fitness values per generation for the PL-MOGA and the WS-MOGA respectively. For both alogrihtms there are three benchmarking files: one baseline run and two extreme runs for the polarizability and HOMO-LUMO gap respectively. The following tables show how the scaling factors and weights have been chosen for the different runs.

.. list-table:: PL-MOGA parameters
   :widths: 25 25 25
   :header-rows: 1

   * - Run
     - Scaling factor polarizability
     - Scaling factor HOMO-LUMO gap
   * - Normal
     - 0.0
     - 0.0
   * - Extreme polarizability
     - 1.0
     - 0.0
   * - Extreme HOMO-LUMO gap
     - 0.0
     - 1.0

|

.. list-table:: WS-MOGA parameters
   :widths: 25 25 25
   :header-rows: 1

   * - Run
     - Weight polarizability
     - Weight HOMO-LUMO gap
   * - Normal
     - 1.0
     - 1.0
   * - Extreme polarizability
     - 0.8
     - 0.2
   * - Extreme HOMO-LUMO gap
     - 0.2
     - 0.8

