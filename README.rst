prodimopy
=========

Python package for reading and plotting ProDiMo results.

Any bug reports or feature requests are very welcome.
If you want to contribute some code please contact me (Christian Rab)


Notebook examples
*****************
If you want to take a look before installing you can try prodimopy
on the web in a binder environment:

.. image:: https://mybinder.org/badge.svg 
   :target: https://mybinder.org/v2/git/https%3A%2F%2Fbitbucket.org%2Fcheesyog%2Fprodimopy/bdb789b71c61d6c55f263de57a8ff95c1e7236c8?filepath=notebooks&urlpath=lab/tree/notebooks

On your left hand side you will see the notebooks (currently only one), just open one and try it!

Requirements
************
prodimopy uses several additional python packages which are commonly used in the astronomical community. 
If you use anaconda_ all this packages should be available in your python distribution. 
The following packages are required

* *matplotib* required for the plotting part only, version>=2 is recommended  
* *astropy*     version >=1.3 is recommended
* *numpy*       no known special requirements
* *scipy*       no known special requirements

If you use the setup script (see Installation) those packages will be installed automatically if 
they are not included in your python distribution. We recommend to use python3 but python2 should
also still work.

Installation
************
Currently the easiest way to use it is to clone this repository and install the package directly from the source:

* change into a directory of your choice and 
* clone the repository 

::

  git clone https://cheesyog@bitbucket.org/cheesyog/prodimopy.git
 
* change into the newly created prodimopy directory and type:

::

  python setup.py develop

This will install the package in your current python environment (should be the one you want to use for ProDiMo). 
The develop options allows to update the python code (e.g. via git) without the need to reinstall the package.

If you do not have root access to install python packages, this should work

::

  python setup.py develop --user


Code Update
***********
Simply type 

::

  git pull 

in the prodimopy directory. You can directly use the updated code (no reinstall required).

Documentation
*************
Please check out the documentation! Click on the badge!

.. image:: https://readthedocs.org/projects/prodimopy/badge/?version=latest
  :target: https://prodimopy.readthedocs.io/en/latest/?badge=latest
  :alt: Documentation Status


.. _anaconda: https://www.anaconda.com/distribution/