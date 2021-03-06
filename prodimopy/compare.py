from __future__ import print_function
from __future__ import division 
from __future__ import unicode_literals

import numpy as np

from inspect import ismethod
import os


class CompareAbs(object):
  """
  An "abstract" class for comparing some kind of |prodimo| model.
  
  A subclass of this class needs to implement the necessary compare routine(s).
  (example see :class:`Compare`)
  
  """
  def diffArray(self,a,aref,diff):
    """
    Checks the relative difference between two arrays. 
    
    If the arrays do not have the same shape they are considered
    as inequal (return `False`).
    
    Parameters
    ----------
    a : array_like(float,ndim=whatever)
      an array.
      
    aref: array_like(float,ndim=same as a)
      the reference array for comparison.
      
    diff : float
      if the values of the arrays differe only by <diff the are considered
      as equal.
    """
    if a.shape != aref.shape:
      return False,None
    
    da=np.absolute(a/aref-1.0)
    if da.max() >= diff:
      return False,da
    
    return True,da

  def diff(self,val,valref,diff):
    """
    Checks the relative difference between two values.
    
    Parameters
    ----------
    a : float
      a value.
      
    aref: float
      the reference value
      
    diff : float
      if the two values differe only by <diff the are considered
      as equal.
    """
    d=abs(val/valref-1.0)
    if d >= diff: return False,d
     
    return True,d 

  def doAll(self):
    """
    Utility function to call all `compare*` functions.
    
    Prints out what function failed and the errors.
    """
    for name in dir(self):
      if name.startswith("compare"):
        cfun = getattr(self, name)
        if ismethod(cfun):
          print("{:25s}".format(name+": "),end="")
          ok,val=cfun()
          
          if ok:
            print("{:8s}".format("OK"),end="")
          else:
            print("{:8s}".format("FAILED"),end="")
          
          if val is not None:
            print("  Max/Avg/Index Max rel. Error: ","{:6.3f}%".format(100.0*np.max(val)),
                  "{:6.3f}%".format(100.0*np.average(val)),
                  "{:3d}".format(np.argmax(val)))
          else:
            print("  Max rel. Error: ",str(val))


class Compare(CompareAbs):
  """
  Class for comparing to ProDiMo models of type :class:`prodimopy.read.ProDiMo_Data`       
  
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
    self.dcdnmol=0.5  # the allowed differences for the column densities (radial and vertical)
                      # chemistry is difficult and uncertain :) FIXME: this has to become better, there 
                      # a simply some columns which usually fail an require time-dependent chemistry and 
                      # the outcome seems to depend on numerical uncertainties ...
    # (e.g. 7% point in the column density goes crazy ... still oka                   
    self.fcdnmol= 0.1   
    self.dTgas=5e-1   # FIXME: 50% Tgas is also quite sensitive, but lower would be better
    self.dLineFluxes=5.e-2  # 5% for line fluxes
    self.dAbundances=1.e-2  # chemistry is difficult and uncertain :)     
    self.lAbundances=1.e-30
    self.dZetaCR=self.d
    self.dZetaX=3.e-1 # FIXME: 50% , should not be that high
    self.dHX=3.e-1 # FIXME: 50% , should not be that high
    self.lZetaX=1.e-25
#    self.specCompare=("e-","H2","CO","H2O","Ne+","Ne++","H+")
#    self.specCompare=("N2","N2#","CO#","H2O#","H3","H3+","HCO+","HN2+","SO2","SiO")
#    self.specCompare=("CO","CN")

    self.specCompare=("e-","H2","CO","H2O","N2","N2#","CO#","H2O#","H3+","HCO+","HN2+","SO2","SiO",
                      "Ne+","Ne++","H+","OH","C+","S+","Si+","N+","CN","HCN","NH3")

      
  def compareLineFluxes(self): 
    '''
    Compares the line fluxes
    Currently assumes that both models include the same lines in the same order.
    '''
    if self.m.lines is None and self.mref.lines is None: return True,None
    
    if self.m.lines is not None and self.mref.lines is None: return False,None
    
    if self.m.lines is None and self.mref.lines is not None: return False,None       

    mFluxes=np.array([x.flux for x in self.m.lines])
    mrefFluxes=np.array([x.flux for x in self.mref.lines])
    f,d=self.diffArray(mFluxes, mrefFluxes, self.dLineFluxes)
    if f == False:
      return False,d
  
    return True,d
  
  def compareLineEstimates(self): 
    '''
    Compares the fluxes from the line estimates
    Currently assumes that both models include the same lines in the same order.
    
    TODO: can actually be merged with compareLineFluxes
    '''
    if self.m.lineEstimates is None and self.mref.lineEstimates is None: return True,None
    
    if self.m.lineEstimates is not None and self.mref.lineEstimates is None: return False,None
    
    if self.m.lineEstimates is None and self.mref.lineEstimates is not None: return False,None       

    mFluxes=np.array([x.flux for x in self.m.lineEstimates])
    mrefFluxes=np.array([x.flux for x in self.mref.lineEstimates])
    f,d=self.diffArray(mFluxes, mrefFluxes, self.dLineFluxes)
    if f == False:
      return False,d
  
    return True,d
  
    
  def compareSED(self):
    """
    Compares the SEDs 
    """
    if self.m.sed is None and self.mref.sed is None: return True,None
    if self.m.sed is not None and self.mref.sed is None: return False,None
    if self.m.sed is None and self.mref.sed is not None: return False,None
    f,d=self.diffArray(self.m.sed.fnuErg, self.mref.sed.fnuErg, self.d)
    
    if f == False:
      return False,d
    
    return True,d  
  
  def compareContinuumImages(self):
    """
    Compares some of the continuum images.  
    """
    if self.m.contImages is None and self.mref.contImages is None: return True,None
    if self.m.contImages is not None and self.mref.contImages is None: return False,None
    if self.m.contImages is None and self.mref.contImages is not None: return False,None
    
    
    for wl in [1,10,100,1000]:
      imm,immwl=self.m.contImages.getImage(wl)
      imref,imrefwl=self.mref.contImages.getImage(wl)    
    
      f,d=self.diff(immwl, imrefwl, self.d)
      if f == False:
        return False,d
    
      f,d=self.diffArray(imm, imref, self.d)
    
      if f == False:
        return False,d
    
    return True,d  
  
  
  def compareCdnmol(self):
    ''' 
    checks the vertical column densities
    only the outermost points are checked
    '''    
    specidxM=[self.m.spnames[idx] for idx in self.specCompare if idx in self.m.spnames]    
    #print(specidxM)    
    specidxMref=[self.mref.spnames[idx] for idx in self.specCompare if idx in self.mref.spnames] 
    
    ok,diffarray= self.diffArray(self.m.cdnmol[:,0,specidxM],self.mref.cdnmol[:,0,specidxMref],self.dcdnmol)
    
    # if false check if it is only a certain fraction of the columns, it than can 
    # be still okay
    # is not really elegant I would say
    # TODO: also somehow return the number of failed columns
    # TODO: maybe merge rcdnmol
    if ok is False:
      ok=True # and check if any columns faild    
      for i in range(len(specidxMref)):
        faildcolumns=(diffarray[:,i]>self.dcdnmol).sum()
        if ((float(faildcolumns)/float(len(diffarray[:,i])))>self.fcdnmol):
          ok=False
                
    return ok,diffarray

  def compareRcdnmol(self):
    ''' 
    checks the radial column densities
    only the outermost points are checked
    ''' 
    
    specidxM=[self.m.spnames[idx] for idx in self.specCompare if idx in self.m.spnames]    
    #print(specidxM)    
    specidxMref=[self.mref.spnames[idx] for idx in self.specCompare if idx in self.mref.spnames] 

    ok,diffarray= self.diffArray(self.m.rcdnmol[-1,:,specidxM],self.mref.rcdnmol[-1,:,specidxMref],self.dcdnmol)

    # if false check if it is only a certain fraction of the columns, it than can 
    # be still okay
    # is not really elegant I would say
    # TODO: also somehow return the number of failed columns
    # TODO: maybe merge Cdnmol
    if ok is False:
      ok=True # and check if any columns faild    
      for i in range(len(specidxMref)):
        faildcolumns=(diffarray[:,i]>self.dcdnmol).sum()
        if ((float(faildcolumns)/float(len(diffarray[:,i])))>self.fcdnmol):
          ok=False

    return ok,diffarray 
  
  def compareDustOpac(self):
    '''
    Compares the dust opacities.
    '''    
    f,d=self.diffArray(self.m.dust.kext, self.mref.dust.kext, self.d)
    
    if f is False:
      return False,d
  
    return True,d

  def compareDustCS(self):
    '''
    Compares the dust cross-sections (from dust_opac.out).
    '''    
    f,d=self.diffArray(self.m.dust.kextcs, self.mref.dust.kextcs, self.d)
    
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

  def compareHX(self):
    '''
    checks the Xray energy deposition rate
    '''
        # set low values to zero 
#    self.m.Hx[self.m.zetaX < self.lZetaX]=self.lZetaX
#    self.mref.Hx[self.mref.zetaX < self.lZetaX]=self.lZetaX
    
    return self.diffArray(self.m.Hx,self.mref.Hx,self.dHX)  

class CompareMc(CompareAbs): 
  """
  Class for comparing two ProDiMo models of type :class:`prodimopy.read_mc.DataMc`       
  
  Every compare Function returns true or false, and the relative differences
  (in case of arrays these are arrays). 
  
  Can be used in e.g. automatic testing routines or in simple command line
  tools to compare ProDiMo model results. 
  """
 
  def __init__(self,model,modelref):
    self.m=model
    self.mref=modelref
    # the allowed relative difference between two values 
    self.d=1.e-10
    self.dabundance=1.e-5
    self.drc=1.e-10 # the difference for the rate coefficients
    
  
  def compareAbundances(self):
    """
    Compares the abundances of two molecular cloud (0D chemistry) models. 
    
    Assumes that both models used the same number of ages and species in the same order.
    """    
    # Do not consider the first age entry at it is the initial conditions 
    # that can vary from model two model and are not really a result
    
    return self.diffArray(self.m.abundances[1:,:],self.mref.abundances[1:,:],self.dabundance)

  def compareRatecoefficients(self):
    """
    Compares the rate coefficients of two molecular cloud (0D chemistry) models. 
    
    Assumes that both models have exactly the same chemical reactions in the same order.
    """    
    return self.diffArray(self.m.ratecoefficients,self.mref.ratecoefficients,self.drc)


def eval_model_type(modelDir):
  """
  Try to guess the model type (e.g. mc, full prodimo etc.).
  Default is always prodimo.
  
  Possible types: 
    `prodimo` .... full prodimo model
    `mc` ......... molecular cloud (0D model) 
  
  Returns
  -------
    str either prodmimo or mc  
  
  FIXME: this is maybe not the best place for this routine
  FIXME: provide constants for the values (what's the best way to do this in pyhton?)
  FIXME: this is just a quick hack, it would be better to use the parameters in Parameter.in
  
  """

  
  if os.path.isfile(modelDir+"/ProDiMo_01.out"):
    mtype = "prodimoTD"
  
  elif os.path.isfile(modelDir+"/Molecular_Cloud_Input.in"):
    mtype="mc"
  else:
      mtype="prodimo"
  
  return mtype



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
          