prodimopy
=========

Python package for reading and plotting ProDiMo results.

This is still under development, but it is usable.
Any bug reports or feature requests are very welcome.


Notebook examples
*****************
If you want to take a look before installing something you can try prodimopy
on the web in a binder environment:

.. image:: https://mybinder.org/badge.svg :target: https://mybinder.org/v2/git/https%3A%2F%2Fbitbucket.org%2Fcheesyog%2Fprodimopy/997e05a5ea66dfdf4d01be523180d8156963c576?filepath=notebooks&urlpath=lab

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

* copy the git clone command at the top of this page
* change into a directory of your choice and execute the command
* change into the prodimopy directory and type:

::

  python setup.py develop

This will install the package in your current python environment (should be the one you want to use for ProDiMo). 
The develop options allows to update the python code (e.g. via git) without the need to reinstall the package.

If you need root access to install python packages but you do not have it, you can use this

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
Some preliminary documentation can be found here http://prodimopy.readthedocs.io


.. _anaconda: https://www.anaconda.com/distribution/