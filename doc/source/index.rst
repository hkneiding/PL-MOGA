Welcome to the PL-MOGA documentation!
===================================

`PL-MOGA <https://github.com/hkneiding/PL-MOGA>`_ (Pareto Lighthouse Multiobjective Genetic Algorithm) is a multiobjective genetic algorithm for the de novo design of transition metal complexes. It makes use of multiobjective selection functions and fitness masks to push the population towards specific regions of the Pareto front. Details of the implementation can be found in the corresponding publication `Directional Multiobjective Optimization of Metal Complexes at the Billion-Scale with the tmQMg-L Dataset and PL-MOGA Algorithm <https://doi.org/10.26434/chemrxiv-2023-k3tf2-v2>`_

There are multiple parts to PL-MOGA repository:

 - The backbone of the package `gapy/` provides a mostly generic multiobjective genetic algorithm framework.
 - Its functionalities are used in the `main.py` file that contains the code for the PL-MOGA.
 - The file `ligands_info.py` contains the names and charges of ligands used for the investigations carried out in the respective paper.
 - The directory `ligands/` contains an xyz file of all the ligands used in the corresponding publication and a script to add them to the molSimplify ligand library.
 - The directory `xyz_examples` contains examples of obtained TMCs through the use of the PL-MOGA.

The section `Get started` of this documentation discusses the required packages and gives a detailed description of the installation process. Furthermore, it explains how to run the MOGA with your custom ligands. The section `Usage` gives instructions on how to setup a MOGA run using the configuration file and explains the relevance of each parameter.

.. note::

   If you encounter any problems, errors or bugs please do not hesitate to contact me via mail (hanneskn@uio.no).

.. toctree::
   :caption: Contents

   setup
   use

