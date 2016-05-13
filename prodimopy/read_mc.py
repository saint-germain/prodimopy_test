from __future__ import division 
from __future__ import print_function

'''
Holds all the output stuff from a ProDiMo Molecular Cloud model

Currently only selected data is included

@author: rab
'''
import numpy as np

class Data_mc(object):
  def __init__(self, name):  
    self.name = name
    self.ages = None
    self.species = None
    self.abundances = None
    
    
    
    
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
  '''
  Reads the final Molecular Cloud abundances.
  '''
  if name == None:
    dirfields = directory.split("/")
    name = dirfields[len(dirfields) - 1]

  mc=Data_mc(name) 
  
  f=open(directory+"/"+filename)
  
  species=list()
  abun=list()
  
  for line in f:
    fields=line.strip().split()
    species.append(fields[0])
    abun.append(float(fields[2]))
    
  mc.species=species
  mc.abundances=np.array(abun)
  
  return mc
  
'''
Reads the whol output of a molecular cloud runs
Including the ages and the species list

@author: rab
'''    
def read_mc(directory,filename,agesfile="mc_ages.txt",speciesfile="mc_species.txt",name=None):
  
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
  '''
  Reads the results of a UMIST rate13 code model
  I only can deal with the output produced by dc.pl script
  provided by the code distribution
  
  The date is provided as a Data_mc object
  TODO: use the output of the code the dc.out file 
  '''
  
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
