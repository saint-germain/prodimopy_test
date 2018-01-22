Create simple |prodimo| grids
*****************************

This module provides routines to generate simple (and small) |prodimo| model 
grids.

This is a working and usable implementation but still very limited. 
More features and improvements are planned and might also require to change 
the interfaces (e.g. how the grid routines are used) significantly.

Usage example
-------------

The following example reads defines the parameters and their ranges used for 
the grid, creates the grid and runs the grid by submitting a jobscript on 
a cluster. 

.. code-block:: python

  import prodimopy.grid as pgrid

  gridname="grid1"
  indir="IN"

  params=list()
  params.append(["Lstar",100,300,3,"lin"])
  params.append(["env_Mdot",1.e-6,1.e-4,3,"log"])
  params.append(["Rout",3000,6000,4,"lin"])


  modelnames=pgrid.make_grid(gridname, params, indir)
  pgrid.run_grid(gridname,modelnames,
    "sbatch -J $MODELNAME$ -p short ../../prodimo.slurm ParameterGrid.in")

Source documentation
--------------------

.. automodule:: prodimopy.grid


.. some usefull replacements for units

.. |gcm^-3| replace:: g cm\ :sup:`-3`
.. |cm^-3| replace:: cm\ :sup:`-3`
.. |cm^-2| replace:: cm\ :sup:`-2`
.. |cm^-1| replace:: cm\ :sup:`-1`
.. |s^-1| replace:: s\ :sup:`-1`
.. |sr^-1| replace:: sr\ :sup:`-1`
.. |Hz^-1| replace:: Hz\ :sup:`-1`