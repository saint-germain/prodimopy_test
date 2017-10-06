Read routines for a molecular cloud |prodimo| Model (0D chemistry)
******************************************************************

Routines to read the output of a (time-dependent) molecular cloud (mc) |prodimo| model. 
All the data belonging to a mc |prodimo| model is put into a hierachical data structure (:class:`~prodimopy.read_mc.Data_mc`).
Those kind of |prodimo| models are 0D chemistry models which provide the abundances for a given set of parameters (e.g. density,temperature etc.)

The module provides routines to read only the final abundances (last age or steady-state model) or the
abundances for all ages (as a function of time). Also the according ages and species names are read.

Usage example
-------------
Reads the time-dependent results of a molecular cloud |prodimo| model from the current working directory..

.. code-block:: python

   import prodimopy.read_mc as pread_mc
   
   model=pread_mc.read_mc("MC_Results.out")
   
   print(model)

Source documentation
--------------------
.. automodule:: prodimopy.read_mc


.. some useful replacements for units

.. |gcm^-3| replace:: g cm\ :sup:`-3`
.. |cm^-3| replace:: cm\ :sup:`-3`
.. |cm^-2| replace:: cm\ :sup:`-2`
.. |cm^-1| replace:: cm\ :sup:`-1`
.. |s^-1| replace:: s\ :sup:`-1`
.. |sr^-1| replace:: sr\ :sup:`-1`
.. |Hz^-1| replace:: Hz\ :sup:`-1`