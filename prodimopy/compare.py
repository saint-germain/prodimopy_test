from __future__ import print_function
from __future__ import division 
from __future__ import unicode_literals

import numpy as np

from inspect import ismethod

class Compare(object):
  """
  Class for comparing to ProDiMo models.
  Initialization requires two ProDiMo_Data models, where the second one
  is the "reference" model
  
  Every compare Function returns true or false, and the relative differences
  (in case of arrays these are arrays). 
  
  Can be used in e.g. automatic testing routines or in simple command line
  tools to compare ProDiMo model results. 
  """
  def __init__(self,model,modelref):
    self.m=model
    self.mref=modelref
    # the allowed difference between the line fluxes
    self.d=1.e-2
    self.dcdnmol=1.0  # the allowed differences for the column densities (radial and vertical)
                      # chemistry is difficult and uncertain :) FIXME: this has to become better, there 
                      # a simply some columns which usually fail an require time-dependent chemistry and 
                      # the outcome seems to depend on numerical uncertainties ... 
    self.dTgas=5.e-2   # 5% Tgas is also quite sensitive, but lower would be better
    self.dLineFluxes=5.e-2  # 5% for line fluxes
    self.dAbundances=1.e-2  # chemistry is difficult and uncertain :)     
    self.lAbundances=1.e-30
    self.dZetaCR=self.d
    self.dZetaX=self.d
    self.lZetaX=1.e-25
    self.specCompare=("e-","H2","CO","H2O","N2","N2#","CO#","H2O#","H3","H3+","HCO+","HN2+","SO2","SiO",
                      "Ne+","Ne++","H+","OH","C+","S+","Si+","CN","HCN","NH3")
    
  def diffArray(self,a,aref,diff):
    '''
    Checks the relative difference between two arrays
    '''    
    da=np.absolute(a/aref-1.0)
    if da.max() >= diff:
      return False,da
    
    return True,da

  def diff(self,val,valref,diff):
    '''
    Checks the relative difference between to values
    '''    
    d=abs(val/valref-1.0)
    if d >= diff: return False,d
     
    return True,d 

  def compareLineFluxes(self): 
    '''
    Compares the line fluxes
    Currently assumes that both models include the same lines in the same order.
    '''
    if self.m.lines is None and self.mref.lines is None: return True,None
    
    if self.m.lines is not None and self.mref.lines is None: return False,None
    
    if self.m.lines is None and self.mref.lines is not None: return False,None       

    # Compare fluxes    
    for i in range(len(self.m.lines)):
        f,d=self.diff(self.m.lines[i].flux, self.mref.lines[i].flux, self.dLineFluxes)
        if f == False:
          return False,d
  
    return True,None
  
  '''
  Compares the SEDs 
  '''
  def compareSED(self):
    if self.m.sed is None and self.mref.sed is None: return True,None
    if self.m.sed is not None and self.mref.sed is None: return False,None
    if self.m.sed is None and self.mref.sed is not None: return False,None
    f,d=self.diffArray(self.m.sed.nuFnu, self.m.sed.nuFnu, self.d)
    
    if f == False:
      return False,d
    
    return True,None  
  
  def compareCdnmol(self):
    ''' 
    checks the vertical column densities
    only the outermost points are checked
    '''    
    specidxM=[self.m.spnames[idx] for idx in self.specCompare if idx in self.m.spnames]    
    #print(specidxM)    
    specidxMref=[self.mref.spnames[idx] for idx in self.specCompare if idx in self.mref.spnames] 
     
    return self.diffArray(self.m.cdnmol[:,0,specidxM],self.mref.cdnmol[:,0,specidxMref],self.dcdnmol)

  def compareRcdnmol(self):
    ''' 
    checks the radial column densities
    only the outermost points are checked
    ''' 
    
    specidxM=[self.m.spnames[idx] for idx in self.specCompare if idx in self.m.spnames]    
    #print(specidxM)    
    specidxMref=[self.mref.spnames[idx] for idx in self.specCompare if idx in self.mref.spnames] 

    return self.diffArray(self.m.rcdnmol[-1,:,specidxM],self.mref.rcdnmol[-1,:,specidxMref],self.dcdnmol)
  
  def compareDustOpac(self):
    '''
    Compares the dust opacities.
    '''
    
    f,d=self.diffArray(self.m.dust.kext, self.mref.dust.kext, self.d)
    
    if f is False:
      return False,d
  
    return True,d
  
  def compareStarSpec(self):
    '''
    Compares the input Stellar spectrum, from X-rays to mm
    '''
    f,d=self.diffArray(self.m.starSpec.Inu, self.mref.starSpec.Inu, self.d)
    
    if f is False:
      return False,d
  
    return True,d
    
  
  def compareTg(self):
    '''
    checks the gas Temperature
    '''
    return self.diffArray(self.m.tg,self.mref.tg,self.dTgas)  

  def compareTd(self):
    '''
    checks the dust Temperature
    '''
    return self.diffArray(self.m.td,self.mref.td,self.d) 
  
  def compareZetaCR(self):
    '''
    checks the cosmic ray ionisation rate
    '''
    return self.diffArray(self.m.zetaCR,self.mref.zetaCR,self.dZetaCR) 
  
  def compareZetaX(self):
    '''
    checks the Xray ionisation rate
    '''
        # set low values to zero 
    self.m.zetaX[self.m.zetaX < self.lZetaX]=self.lZetaX
    self.mref.zetaX[self.mref.zetaX < self.lZetaX]=self.lZetaX
    
    return self.diffArray(self.m.zetaX,self.mref.zetaX,self.dZetaX) 

#   def compareFlineEstimates(self):
#     '''
#     Compares the FlineEstimaes    
#     '''
#     
#     if len(self.m.lineEstimates) !=len(self.mref.lineEstimates):
#       return False,None     
#     
#     
#     self.diff(self.m.lineEstimate[i].flux, self.mref.lineEstimate[i].flux, self.dLineFluxes)
#     
#     # Compare fluxes    
#     for i in range(len(self.m.lines)):
#         f,d=
#         if f == False:
#           return False,d
#   
#     return True,None
#   
#   def compareAbundances(self):
#     if self.m.nspec != self.mref.nspec: return False,None
#     
#     # set low values to zero 
#     self.m.nmol[self.m.nmol < self.lAbundances]=self.lAbundances
#     self.mref.nmol[self.mref.nmol < self.lAbundances]=self.lAbundances
#     
#     #print(self.m.nmol[:,0,0])
#     #print(self.mref.nmol[:,0,0])
#     # only check a small selection of species
#     spec=("H2","CO","H2O","CO#","H2O#","H3+")
#     specidxM=[self.m.spnames[idx] for idx in spec if idx in self.m.spnames]    
#     #print(specidxM)    
#     specidxMref=[self.mref.spnames[idx] for idx in spec if idx in self.mref.spnames]    
#     #print(specidxMref)    
#   
#     return self.diffArray(self.m.nmol[:,:,specidxM],self.mref.nmol[:,:,specidxMref],self.dAbundances)  
  
  def doAll(self):
    '''
    Utility function to call all compare* functions and print out if the 
    test failed or not     
    '''
    
    for name in dir(self):
      if name.startswith("compare"):
        cfun = getattr(self, name)
        if ismethod(cfun):
          print("{:25s}".format(name+": "),end="")
          ok,val=cfun()
          
          if ok:
            print("{:8s}".format("OK"))
          else:
            print("{:8s}".format("FAILED"),end="")
            if val is not None:
              print("  Max/Avg rel. Error: ","{:10.2%}".format(np.max(val)),"{:10.2%}".format(np.average(val)))
            else:
              print("  Max rel. Error: ",str(val))
          