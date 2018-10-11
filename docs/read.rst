Read routines for a ProDiMo Model
*********************************

This module provides several routines to read the output of a |prodimo| model. 
All the data belonging to a |prodimo| model is put into an hierachical data structure (:class:`~prodimopy.read.Data_ProDiMo`).

The module provides a routine to read "all" the data of a |prodimo| model and spezialized routines  
to read only distinct model data (e.g. only the SED of a model).

Usage example
-------------

The following example reads the output files of a |prodimo| model from the current working directory.
The :func:`~prodimopy.read.read_prodimo` function tries to read (nearly)all |prodimo| output data which can found.
There are also more spezialed read routines (see :mod:`prodimopy.read`). 

.. code-block:: python

   import prodimopy.read as pread
   
   model=pread.read_prodimo()
   
   print(model)


Source documentation
--------------------

.. automodule:: prodimopy.read


.. some usefull replacements for units

.. |gcm^-3| replace:: g cm\ :sup:`-3`
.. |gcm^-2| replace:: g cm\ :sup:`-2`
.. |ergcm^-3| replace:: erg cm\ :sup:`-3`
.. |kms^-1| replace:: km s\ :sup:`-1`
.. |cm^-3| replace:: cm\ :sup:`-3`
.. |cm^-2| replace:: cm\ :sup:`-2`
.. |cm^-1| replace:: cm\ :sup:`-1`
.. |s^-1| replace:: s\ :sup:`-1`
.. |sr^-1| replace:: sr\ :sup:`-1`
.. |Hz^-1| replace:: Hz\ :sup:`-1`