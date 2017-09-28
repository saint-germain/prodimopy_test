prodimopy
=========

Python package for plotting ProDiMo results.

This is still under development, but you can use it already, if you want. 

Requirements
************
prodimopy uses several additional python packages which are commonly used in the astronomical community. 
If you use anaconda_ all this packages should be available in your python distribution. 
The following packages are required

* *matplotib* required for the plotting part only, version>=2 is recommended  
* astropy     version >=1.3 is recommended
* numpy       no known special requirements
* scipy       no known special requirements

If you use the setup script (see Installation) those packages will be installed automatically if 
they are not included in your python distribution.

Installation
************
Currently the easiest way to use it is to clone this repository and install the package directly from the source:

* click on the Clone link in the menu and copy the git command
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