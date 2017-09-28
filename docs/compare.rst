.. _sec_compare:

Compare two |prodimo| models
****************************

This module provides several routines to compare two |prodimo| models.

Currently only certain output fields are compared. In each comparison a certain tolerance 
is considered. Those values depend on which quantity is compared (e.g. temperatures, abundances). 

Currently it is possible to compare full (standard) |prodimo| models and molecular cloud (0D chemistry) 
models.

This module is mainly used for automatic tests. However, it is also useful during development to compare
two models with one simple command.

.. automodule:: prodimopy.compare


.. some usefull replacements for units

.. |gcm^-3| replace:: g cm\ :sup:`-3`
.. |cm^-3| replace:: cm\ :sup:`-3`
.. |cm^-2| replace:: cm\ :sup:`-2`
.. |cm^-1| replace:: cm\ :sup:`-1`
.. |s^-1| replace:: s\ :sup:`-1`
.. |sr^-1| replace:: sr\ :sup:`-1`
.. |Hz^-1| replace:: Hz\ :sup:`-1`