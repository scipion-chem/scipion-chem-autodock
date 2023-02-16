================================
Autodock scipion plugin
================================

**Documentation under development, sorry for the inconvenience**

Scipion framework plugin for the use of tools provided by Autodock software suite

===================
Install this plugin
===================

You will need to use `3.0.0 <https://github.com/I2PC/scipion/releases/tag/v3.0>`_ version of Scipion
to run these protocols. To install the plugin, you have two options:

- **Stable version**  

.. code-block:: 

      scipion installp -p scipion-chem-autodock
      
OR

  - through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**
      
- **Developer's version** 

1. **Download repository**:

.. code-block::

            git clone https://github.com/scipion-chem/scipion-chem-autodock.git

2. **Switch to the desired branch** (master or devel):

Scipion-chem-autodock is constantly under development.
If you want a relatively older an more stable version, use master branch (default).
If you want the latest changes and developments, user devel branch.

.. code-block::

            cd scipion-chem-autodock
            git checkout devel

3. **Install**:

.. code-block::

            scipion installp -p path_to_scipion-chem-autodock --devel

- **Binary files** 

Atom_struct_utils plugin is a pure Python module, no binary files are required. 

- **Tests**

To check the installation, simply run the following Scipion test:

===============
Buildbot status
===============

Status devel version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/bioinformatics_dev.svg

Status production version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/bioinformatics_prod.svg
