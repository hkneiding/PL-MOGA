===========================================================
Pareto Lighthouse Multiobjective Genetic Algorithm (PL-MOGA)
===========================================================

``PL-MOGA`` is a multiobjective genetic algorithm for the *de novo* design of transition metal complexes. It makes use of multiobjective selection functions and fitness masks to push the population towards specific regions of the Pareto front. Details of the implementation can be found in the corresponding publication `Directional Multiobjective Optimization of Metal Complexes at the Billion-Scale with the tmQMg-L Dataset and PL-MOGA Algorithm <https://chemrxiv.org/engage/chemrxiv/article-details/651051d4ed7d0eccc32252ea>`_.

Requirements
-----------

This package requires a Python (>3.7.x) installation with the following packages:

- `molSimplify <https://github.com/hjkgrp/molSimplify>`_
- `xTB <https://github.com/grimme-lab/xtb>`_
- `uxtbpy <https://github.com/hkneiding/uxtbpy>`_

Note that ``molSimplify`` has a lot of dependencies on its own and I recommend installing it first before ``xTB`` and ``uxtbpy``.

How to use
-----------

The backbone of the package can be installed directly from this repository with ``pip``::
    
    pip install git+https://github.com/hkneiding/PL-MOGA

which installs ``gapy`` as a library to your Python installation or virtual environment.

Afterwards you can import the library with:

>>> import gapy

which provides mostly generic (MO)GA functions and a flexible framework to use your own selection, fitness and masking functions.

To get a quick start into production you can use the code provided in ``main.py``. The only thing you need to supply are the names and charges of the ligands you want to use, which should be stored in ``ligands_info.py``. Note that the names of the ligands should be the same as the ones used in your molSimplify installation. Information on how to add your own custom ligands to molSimplify can be found `here <http://hjkgrp.mit.edu/tutorials/2018-05-09-molsimplify-tutorial-10-adding-ligands-molsimplify>`_

Difficulties?
-----------

If you encounter any problems, errors or bugs please do not hesitate to open an issue or directly contact me via mail (hanneskn@uio.no).
