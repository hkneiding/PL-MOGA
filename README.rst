===========================================================
Pareto Lighthouse Multiobjective Genetic Algorithm (PL-MOGA)
===========================================================

``PL-MOGA`` is a multiobjective genetic algorithm for the *de novo* design of transition metal complexes. It makes use of multiobjective selection functions and fitness masks to push the population towards specific regions of the Pareto front. Details of the implementation can be found in the corresponding publication `Directional multiobjective optimization of metal complexes at the billion-system scale <https://www.nature.com/articles/s43588-024-00616-5>`_.

Requirements
-----------

This package requires a Python (>3.7.x) installation with the following packages:

- `molSimplify <https://github.com/hjkgrp/molSimplify>`_
- `xTB <https://github.com/grimme-lab/xtb>`_
- `uxtbpy <https://github.com/hkneiding/uxtbpy>`_

Note that ``molSimplify`` has a lot of dependencies on its own and I recommend installing it first before ``xTB`` and ``uxtbpy``.

How to use
-----------

The code can be obtained by running::
    
    $ git clone https://github.com/hkneiding/PL-MOGA

which copies the full project into your current working directory.

Afterwards, runs can be started directly from the projects root directory using the command line::

    $ python main.py config.yml

The ``config.yml`` file contains entries for all relevant PL-MOGA parameters and is used to configure PL-MOGA runs. 

For a detailed documentation, including comprehensive installation instructions and explainations of all parameters `read the docs <https://pl-moga.readthedocs.io/en/latest/setup.html>`_.

Difficulties?
-----------

If you encounter any problems, errors or bugs please do not hesitate to open an issue or directly contact me via mail (hanneskn@uio.no).
