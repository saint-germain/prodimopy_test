"""
.. module:: read 
   :synopsis: Read routines and data structure for molecular cloud (0D chemistry) |prodimo| models.

.. moduleauthor:: Ch. Rab


"""
from __future__ import division 
from __future__ import print_function

import numpy as np

class Data_mc(object):
  """
  Data structure for molecular cloud (0D chemistry) |prodimo| models.

  Can be used for time-dependent abundances or for steady-state (or final) abundances.
  """
  def __init__(self, name):  
    """
    Parameters
    ----------
    
    name : string
      The name of the model.
      
    
    Attributes
    ----------
    
    """
    self.name = name
    """ string :
    The name of the model (can be empty)
    """
    self.species = None
    """ array_like(string,ndim=1) :
    an ordered list of species names.
    """
    self.ages = None
    """ array_like(float,ndim=1) :
    the output ages of the model.
    """
    self.abundances = None
    """ array_like(float,ndim=2) :
    the abundances for each species and each age `DIMS:` (number of ages,number of species).
    """
    
def _read_ages(filename):
  #f=open(filename,"r")
  
  #print(f.readlines())
  ages=np.loadtxt(filename)
  # insert age zero initial abundance
  ages=np.insert(ages,0, 0.0)

  
  #f.close()
  
  return ages

def _read_species(filename):
  f=open(filename,"r")
  
  species=[str.strip() for str in f.readlines()]
  species[species.index("HN2+")]="N2H+"
  f.close()
  
  return species


def read_mc_final(directory,filename="Molecular_cloud.out",name=None):
  """
  Reads the final (last timestep) molecular cloud abundances.
  
  Parameters
  ----------
  
  directory : string
    The model directory.
    
  filename : string
    The name of the file containing the abundances (default: `Molecular_cloud.out`). 
    
  name : string
    The name of the model. Will be shown in the plots (default: `None`). 
  
  
  FIXME: ist not consistent with read_mc, e.g. the species names such as N2H+ are not adapted here
  FIXME: use numpy arrays for the abundances such as for time-dependent models.
  
  """
  if name == None:
    dirfields = directory.split("/")
    name = dirfields[len(dirfields) - 1]

  mc=Data_mc(name) 
  
  f=open(directory+"/"+filename)
  lines = f.readlines()
  f.close()
  
  species=list()
  abun=list()
  
  for line in lines:
    fields=line.strip().split()
    species.append(fields[0])
    abun.append(float(fields[2]))
    
  mc.species=species
  mc.abundances=np.array(abun)
  
  return mc
  
def read_mc(directory,filename,agesfile="mc_ages.txt",speciesfile="mc_species.txt",name=None):
  """
  Read routine for molecular cloud |prodimo| models including the ages and the species list. 

  Parameters
  ----------
  
  directory : string
    The model directory.
    
  filename: string
    The name of the file containing the abundances for all ages. 
    
  agesfile: string
    The file with the ages (default: `mc_ages.txt`)
    
  speciesfile: string
    The file with the species names (default: `mc_species.txt`)
        
  """
  if name == None:
    dirfields = directory.split("/")
    name = dirfields[len(dirfields) - 1]

  mc=Data_mc(name)  
  mc.ages=_read_ages(directory+"/"+agesfile)
  mc.species=_read_species(directory+"/"+speciesfile)
  # make ages first index, species second 
  mc.abundances=np.transpose(np.loadtxt(directory+"/"+filename))
  
  return mc

def read_umist(directory,filename="dc_molecules.dat",name=None):
  """
  Reads the results of a UMIST rate13 code model. 
  Uses the output produced by dc.pl script provided by the UMIST code distribution. 
  
  The data is provided as a :class:`prodimopy.read_mc.Data_mc` object.   
  """
  
  if name == None:
    dirfields = directory.split("/")
    name = dirfields[len(dirfields) - 1]
  
  f=open(directory+"/"+filename,"r")
  
  lines=f.readlines()
  
  header=lines[0].strip()
  header=header[2:]
  species=header.split(" ")  
    
  
  ages=list()
  abun=list()  
  
  for i in range(1,len(lines)):
    line=lines[i].strip()    
    fields=line.split(" ")    
    ages.append(float(fields[0]))
    abun.append(fields[1:])
  
  out=Data_mc(name)
  out.ages=np.array(ages,dtype='|S4').astype(np.float)
  out.species=species
  out.abundances=np.array(abun)
    
  return out

if __name__ == "__main__":
  mc=read_mc("tests/mc","mc.out")
  
  print(mc.species)
  print(mc.ages)
  print(mc.abundances[:,mc.species.index("H2")])
