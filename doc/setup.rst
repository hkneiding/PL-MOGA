Get started
===========

This page discusses the prerequisites and outlines the installation procedure.

============
Requirements
============

This package is developed and tested on a Python 3.7.x version but any higher version should also work without any problems. The package requires installations of the following packages:

- `molSimplify <https://github.com/hjkgrp/molSimplify>`_
- `xTB <https://github.com/grimme-lab/xtb>`_
- `uxtbpy <https://github.com/hkneiding/uxtbpy>`_

The packages should be installed in the order they are listed. Detailed installation instructions can be found in the next section.

============
Installation
============

.. note::

   This installation guide uses the Python environment manager `conda <https://anaconda.org/>`_. If you are not familiar with conda you can familiarite yourself with it `here <https://docs.conda.io/en/latest/>`_.

1. Start by creating a conda environment and activating it.

.. code-block:: console

  $ conda create -n pl-moga
  $ conda activate pl-moga

2. Next, we will install all prerequisites starting with molSimplify. For this we need to clone their repository and navigate to the projects root directory.

.. code-block:: console

  $ git clone https://github.com/hjkgrp/molSimplify.git
  $ cd molSimplify

3. Then we will install molSimplify's prerequisites. The authors of the package provide a environment file that we can use to update our environment with all necessary dependencies. Afterwards we navigate back to our starting directory.

.. code-block:: console

  $ conda env update --file devtools/conda-envs/mols.yml
  $ cd ..

4. And install molSimplify itself.

.. code-block:: console

  $ pip install -e . --no-deps

5. Then we install xTB with conda.

.. code-block:: console

  $ conda install -c conda-forge xtb

6. Then we install uxtbpy with pip.

.. code-block:: console

  $ pip install git+https://github.com/hkneiding/uxtbpy

7. Finally we clone the PL-MOGA repository and navigate to the projects root directory.

.. code-block:: console

  $ git clone https://github.com/hkneiding/PL-MOGA.git
  $ cd PL-MOGA

The cloned repository is ready to be used as is and does not need to be installed.

.. note::

   `molSimplify` has a lot of dependencies of its own and their installation can take a long time.

==============
Adding ligands
==============

In order to use molSimplify to build transition metal complexes during the execution of the genetic algorithm we need to add our own custom ligands to its internal library. For this we need their respective xyz structure and metal coordinating atom indices. All ligands utilized in the corresponding publication are stored in the `ligands.xyz` file in the `ligands` directory. The comment line of each xyz is structured as follows:

```
<name> | charge: <charge> | conn atom: [<conn atoms lists>]
```

In order to add them to molSimplify's library you can use the `add_ligands_to_molsimplify.py` script as follows

.. code-block:: console

  $ python add_ligands_to_molsimplify.py ligands.xyz

Note that adding of the ligands can take some time. Once completed we are ready to setup runs in the chemical spaces introduced in the publication.

.. note::

  If you want to add your custom ligands to the library you should generate an xyz file with the same structure as in `ligands.xyz` and use the `add_ligands_to_molsimplify.py` script passing it as first argument.
