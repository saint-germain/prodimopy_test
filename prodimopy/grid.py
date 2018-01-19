from __future__ import print_function
from __future__ import division 
from __future__ import unicode_literals

import numpy
import math

import os
import shutil
import collections


def genvalues(param):

  start=param[1]
  end=param[2]
  steps=param[3]
  linlog=param[4]
  
  if linlog=="lin":
    return numpy.linspace(start,end,steps)
  elif linlog=="log":
    return numpy.logspace(math.log10(start),math.log10(end),steps,base=10.0)
  else:
    print("Unknown spacing type")
    return None


def run_grid(gridname,modeldirs,runProDiMo):
  os.chdir(gridname)
  print("Grid: ",os.getcwd())
  for modeldir in modeldirs:
    if isinstance(runProDiMo, collections.Callable):    
      runProDiMo(modeldir)
    else:
      runProDiMo=runProDiMo.replace("$MODELNAME$",modeldir)
      os.chdir(modeldir)
      os.system(runProDiMo)
      os.chdir("..")


def make_grid(gridname,params,indir=None):
  """
  Produces a grid of prodimo models. 
  
  A directory with the name gridname is created. Within this directory for 
  each parameter combination a directory is created. This directory is initially
  copied from the directory given by indir and additionally the new parameters 
  are added to a ParameterGrid.in file. If ParameterGrid.in does not exist
  it will be created.
  
  The routine returns the list of modelnames.
  """

  # make the list for the meschgrid
  values=list()
  for param in params:
    values.append(genvalues(param))
  
  grid=numpy.meshgrid(*values)

  # create the directory for the grid
  os.mkdir(gridname)
  
  # create an iterator the loops over all indixes of the firt grid array (all combinations)
  it = numpy.nditer(grid[0], flags=['multi_index'])
  imodel=0
  modelnames=list()
  while not it.finished:
    # create the model directory
    modelname="model"+"{:07d}".format(imodel)
    modeldir=gridname+"/"+modelname
    if indir is not None:
      shutil.copytree(indir,modeldir)
    else: 
      os.mkdir(modeldir)
    fparam=open(modeldir+"/ParameterGrid.in","a+")
    for iparam in range(len(params)):
      fparam.write(str(grid[iparam][it.multi_index])+"  ! "+params[iparam][0] +" \n")
    fparam.close()
    imodel=imodel+1
    modelnames.append(modelname)  
    it.iternext()

  return modelnames