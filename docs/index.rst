.. prodimopy documentation master file, created by
   sphinx-quickstart on Mon Aug 14 16:48:46 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

*********
prodimopy
*********

The python package provides routines for reading and plotting the output of one or more ProDiMo model.


Installation
============

The package can be installed from a bitbucket repository. 
For details see https://bitbucket.org/cheesyog/prodimopy


Source Documentation
====================

.. toctree::
   :maxdepth: 2
   :caption: Modules:

   read
   read_mc
   plot
   plot_models
   read_casasim
   plot_casasim
   compare
   grid


Command-line utilities
======================

prodimopy also installs a few command line utils which can be used without lauching a python interpreter or 
writing any python code.

pplot
-----
Produces plots for a single prodimo model. 
Useful to check the prodimopy installation or to take a quick look on a |prodimo| model 

For details see :ref:`sec_plot`. 

pplot_models
------------
Produces plots for a given set of prodimo models. 
Useful to quickly compare visualy different |prodimo| models. 

For details see :ref:`sec_plot_models`. 

pcompare
--------
Compares the results of two |prodimo| models.

For details see :ref:`sec_compare`. 


Using Jupyter
=============

Most prodimopy routines can also be used within a Jupyter notebook.
In this :download:`example notebook <./jupyter_example.ipynb>` it is shown how to use 
the reading a plotting routines for a single model. 

This feature is not well tested yet. It worked with python 3.6 as part of the anaconda distribution using 
the local Jupyter server.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. TODO's
.. ======
.. .. todolist::