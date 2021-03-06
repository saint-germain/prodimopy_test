"""
.. module:: read 
   :synopsis: Reads the output data of a ProDiMo model.

.. moduleauthor:: Ch. Rab


"""
from __future__ import print_function
from __future__ import division 
from __future__ import unicode_literals

import numpy

from astropy import units as u
from astropy import constants as const
import numpy as np
from scipy.interpolate import interp1d
import os
from collections import OrderedDict
import math
import glob
import tarfile as tar
import io
from timeit import default_timer as timer


class Data_ProDiMo(object):
  """ 
  Data container for most of the output produced by |prodimo|. 
  
  The class also includes some convenience functions and also derives/calculates 
  some addtionialy quantities not directly included in the |prodimo| output. 

  
  .. warning:: 
  
    Not all the data included in ProDiMo.out and not all .out files 
    are considered yet. If you miss something please contact Ch. Rab, or 
    alternatively you can try to implement it yourself.      
         
  """  
  def __init__(self, name):          
    """
    Parameters
    ----------
    name : string
      The name of the model (can be empty).  

    Attributes
    ----------
          
    """
    self.name = name 
    """ string :
    The name of the model (can be empty)
    """
    self.__fpFlineEstimates = None  # The path to the FlineEstimates.out File   
    self.__tarfile=None
    self.nx = None        
    """ int : 
    The number of spatial grid points in the x (radial) direction
    """
    self.nz = None
    """ int : 
    The number of spatial grid points in the z (vertical) direction
    """
    self.nspec = None     
    """ int : 
      The number of chemical species included in the model. 
    """  
    self.nheat = None     
    """ int :
      The number of heating processes included in the model.
    """  
    self.ncool = None     
    """ int :
      The number of cooling processes included in the model.
    """  
    self.dust2gas = None
    """ float :
      The global dust to gass mass ratio (single value, given Parameter)        
    """
    self.mstar = None
    """ float :
      The stellar mass in solar units.
      is taken from ProDiMo.out        
    """
    self.x = None         
    """ array_like(float,ndim=2) :
    The x coordinates (radial direction).
    `UNIT:` au, `DIMS:` (nx,nz)    
    """  
    self.z = None         
    """ array_like(float,ndim=2) :
    The z coordinates (vertical direction).  
    `UNIT:` au, `DIMS:` (nx,nz)    
    """
    self.vol = None         
    """ array_like(float,ndim=2) :
    The volume for each grid point  
    `UNIT:` |gcm^-3|, `DIMS:` (nx,nz)    
    """
    self.rhog = None
    """ array_like(float,ndim=2) :
    The gas density.
    `UNIT:` |gcm^-3|, `DIMS:` (nx,nz)
    """
    self.rhod = None      
    """ array_like(float,ndim=2) :
    The dust density.
    `UNIT:` |gcm^-3|, `DIMS:` (nx,nz)
    """
    self.sdg = None
    """ array_like(float,ndim=2) :
    The gas vertical surface density.
    `UNIT:` |gcm^-2|, `DIMS:` (nx,nz)
    """
    self.sdd = None      
    """ array_like(float,ndim=2) :
    The dust vertical surface density.
    `UNIT:` |gcm^-2|, `DIMS:` (nx,nz)
    """
    
    self.nHtot = None 
    """ array_like(float,ndim=2) :
    The total hydrogen number density.
    `UNIT:` |cm^-3|, `DIMS:` (nx,nz)
    """
    self.muH = None
    """ float :
    The conversion constant from nHtot to rhog
    It is assume that this is constant throught the disk. It is given by 
    by `rhog/nHtot`
    `UNIT:` `g`
    """
    self.NHver = None     # 
    """ array_like(float,ndim=2) :
    Vertical total hydrogen column density. `nHtot` is integrated from the disk 
    surface to the midplane at each radial grid point. The intermediate results 
    are stored at each grid point. For example NHver[:,0] gives the total column 
    density as a function of radius.   
    `UNIT:` |cm^-2|, `DIMS:` (nx,nz)
    """
    self.NHrad = None     #   
    """ array_like(float,ndim=2) :
    Radial total hydrogen column density. Integrated along radial rays, starting 
    from the star. Otherwise same behaviour as `NHver`.
    `UNIT:` |cm^-2|, `DIMS:` (nx,nz)
    """    
    self.nd = None
    """ array_like(float,ndim=2) :
    The dust number density. 
    `UNIT:` |cm^-3|, `DIMS:` (nx,nz)
    """         
    self.tg = None        
    """ array_like(float,ndim=2) :
    The gas temperature. 
    `UNIT:` K, `DIMS:` (nx,nz)
    """     
    self.td = None         
    """ array_like(float,ndim=2) :
    The dust temperature.
    `UNIT:` K, `DIMS:` (nx,nz)
    """ 
    self.pressure = None
    """ array_like(float,ndim=2) :
    The gas pressure
    `UNIT:` |ergcm^-3|, `DIMS:` (nx,nz)
    """
    self.soundspeed = None
    """ array_like(float,ndim=2) :
    The isothermal sound speed.
    `UNIT:` |kms^-1|, `DIMS:` (nx,nz)
    """ 
    self.damean = None    
    """ array_like(float,ndim=2) :
    The mean dust particle radius
    `UNIT:` micron, `DIMS:` (nx,nz)
    """
    self.Hx=None          
    """ array_like(float,ndim=2) :
    The X-ray energy deposition rate per hydrogen nuclei
    `UNIT:` erg <H>\ :sup:`-1`, `DIMS:` (nx,nz)
    """
    self.zetaX = None     
    """ array_like(float,ndim=2) :
    X-ray ionisation rate per hydrogen nuclei.
    `UNIT:` |s^-1|, `DIMS:` (nx,nz)
    """
    self.zetaCR = None    # 
    """ array_like(float,ndim=2) :
    Cosmic-ray ionisation rate per molecular hydrogen (H2) 
    `UNIT:` |s^-1|, `DIMS:` (nx,nz)
    """
    self.zetaSTCR = None  
    """ array_like(float,ndim=2) :
    Stellar energetic particle ionisation rate per H2
    `UNIT:` |s^-1|, `DIMS:` (nx,nz)
    """
    self.tauX1 = None
    """ array_like(float,ndim=2) :
    Radial optical depth at 1 keV (for X-rays). 
    `UNIT:` , `DIMS:` (nx,nz)
    """
    self.tauX5 = None
    """ array_like(float,ndim=2) :
    Radial optical depth at 5 keV (for X-rays).
    `UNIT:` , `DIMS:` (nx,nz)
    """
    self.tauX10 = None
    """ array_like(float,ndim=2) :
    Radial optical depth at 10 keV (for X-rays).
    `UNIT:` , `DIMS:` (nx,nz)
    """
    self.AVrad = None      
    """ array_like(float,ndim=2) :
    Radial visual extinction (measerd from the star outwards). 
    `UNIT:` , `DIMS:` (nx,nz)
    """
    self.AVver = None
    """ array_like(float,ndim=2) :
    Vertical visual extinction (measured from the disk surface to the midplane).    
    `UNIT:`, `DIMS:` (nx,nz)
    """
    self.AV = None
    """ array_like(float,ndim=2) :
    Given by min([AVver[ix,iz],AVrad[ix,iz],AVrad[nx-1,iz]-AVrad[ix,iz]])
    Gives the lowest visiual extinction at a certain point. Where it is assumed radiation 
    can escape either vertically upwards, radially inwards or radially outwards. 
    `UNIT:` , `DIMS:` (nx,nz)
    """
    self.nlam = None      
    """ int : 
    The number of wavelength bands used in the continuum radiative transfer.
    """
    self.lams = None       
    """ array_like(float,ndim=1) :
    The band wavelengths considered in the radiative transfer. 
    \n`UNIT:` microns, `DIMS:` (nlam)
    """                        
    self.radFields = None
    """ array_like(float,ndim=3) :
    Radiation field (mean intensity) for each wavelength band.
    
    `UNIT:` erg |s^-1| |cm^-2| |sr^-1| |Hz^-1|, `DIMS:` (nx,nz,nlam)
    """
    self.chi = None        
    """ array_like(float,ndim=2) :
    Geometrial UV radiation field in units of the Drain field.
    `UNIT:` Draine field, `DIMS:` (nx,nz)
    """
    self.chiRT= None       
    """ array_like(float,ndim=2) :
    UV radiation field as properly calculated in the radiative transfer, in units of the Drain field.
    `UNIT:` Draine field, `DIMS:` (nx,nz)
    """
    self.kappaRos=None
    """ array_like(float,ndim=2) :
    Rosseland mean opacity. In case of gas radiative transfer for the dust plus the gas.
    `UNIT:` |cm^-1|, `DIMS:` (nx,nz)
    """            
    self.tauchem = None    
    """ array_like(float,ndim=2) :
    Chemical timescale (stead-state)
    `UNIT:` yr, `DIMS:` (nx,nz)
    """
    self.taucool = None    
    """ array_like(float,ndim=2) :
    Cooling timescale.
    `UNIT:` yr, `DIMS:` (nx,nz)
    """
    self.taudiff = None    
    """ array_like(float,ndim=2) :
    Vertical radiative diffussion timescale (using the Rosseland mean opacities).
    `UNIT:` yr, `DIMS:` (nx,nz)
    """
    self.spnames = None  # 
    """ dictionary :
    Dictionary providing the index of a particular species (e.g. spnames["CO"]). This index
    can than be used for arrays having an species dimension (like nmol). The electron is included.
    `UNIT:` , `DIMS:` (nspec)
    """
    self.spmasses = None 
    #self.spmassesTot = None # integrated species masses (over the whole model space), is an array for easier handling
    """ dictionary(float) :
    Dictionary for the masses of the individual species (same keys as spnames).
    `UNIT:` g, `DIMS:` (nspec)
    """
    self.nmol = None     
    """ array_like(float,ndim=3) :
    Number densities of all chemical species (mainly molecules)    
    `UNIT:` |cm^-3|, `DIMS:` (nx,nz,nspec)
    """
    self.cdnmol = None             
    """ array_like(float,ndim=3) :
    Vertical column number densities for each chemical species at each point in the disk. 
    Integrated from the surface to the midplane at each radial grid point.   
    `UNIT:` |cm^-2|, `DIMS:` (nx,nz,nspec)
    """
    self.rcdnmol = None             
    """ array_like(float,ndim=3) :
    Radial column number densities for each species at each point in the disk. 
    Integrated from the star outwards along fixed radial rays given by the vertical grid.
    `UNIT:` |cm^-2|, `DIMS:` (nx,nz,nspec)
    """
    self.heat = None
    """ array_like(float,ndim=3) :
    Heating rates for the various heating processes. `FIXME:` find out the unit
    `UNIT:` `DIMS:` (nx,nz,nheat)
    """
    self.cool = None
    """ array_like(float,ndim=3) :
    Cooling rates for the various coling. `FIXME:` find out the unit.
    `UNIT:` `DIMS:` (nx,nz,cool)
    """
    self.heat_names = None
    """ list (string)
    All the names of the cooling processes. 
    """
    self.cool_names = None
    """ list (string)
    All the names of the cooling processes. 
    """
    self.heat_mainidx = None
    """ array_like(float,ndim=3) :
    Index of the main heating process at the given grid point.
    `UNIT:` , `DIMS:` (nx,nz)
    """
    self.cool_mainidx = None
    """ array_like(float,ndim=3) :
    Index of the main cooling process at the given grid point.
    `UNIT:` , `DIMS:` (nx,nz)
    """
    self.lineEstimates = None  
    """ list(:class:`prodimopy.read.DataLineEstimate`) :
    All the line estimates from FlineEstimates.out. Each spectral line in FlineEstimates
    corresponds to one :class:`prodimopy.read.DataLineEstimate` object.     
    """
    self.lines = None          
    """ array_like(:class:`prodimopy.read.DataLine`) :
    Alle the spectral lines from line_flux.out (proper Linetransfer). 
    Each spectral line in line_flux.out corresponds to 
    one :class:`prodimopy.read.DataLine` object    
    """
    self.sed = None            # the spectral energy distribution (from proper ray tracing)
    """ :class:`prodimopy.read.DataSED` :
    The Spectral Energy Distribution for the model (SED) as calculated in the 
    radiative transfer with ray tracing. 
    see :class:`prodimopy.read.DataSED` for details.
    """
    self.starSpec = None
    """ :class:`prodimopy.read.DataStarSpec` :
    The (unattenuated) stellar input spectrum.  
    see :class:`prodimopy.read.DataStarSpec` for details.
    """
    self.gas = None             
    """ :class:`prodimopy.read.DataGas` :
    Holds various properties of the gas component (e.g. opacities).
    see :class:`prodimopy.read.DataGas`    
    """
    self.dust = None           
    """ :class:`prodimopy.read.DataDust` :
    Holds various properties of the dust component (e.g. opacities).
    see :class:`prodimopy.read.DataDust`    
    """
    self.env_dust = None       # dust properties for the envelope structure see :class:`prodimopy.read.DataDust`  
    """ :class:`prodimopy.read.DataDust` :
    Holds various properties of the dust component (e.g. opacities) of the envelope.
    Only relevant if |prodimo| is used in the envelope mode.
    see :class:`prodimopy.read.DataDust`    
    """
    self.elements = None
    """ :class:`prodimopy.read.DataElements` :
    Holds the initial gas phase element abundances
    """
    self.sedObs = None
    """ :class:`prodimopy.read.DataSEDobs` :
    Holds the provided SED observations (photometry, spectra, extinction etc.)
    TODO: maybe put all the observations into one object (e.g. also the lines)
    """
   
  def __str__(self):
    output = "Info ProDiMo.out: \n"
    output += "NX: " + str(self.nx) + " NZ: " + str(self.nz) + " NSPEC: " + str(self.nspec)
    output += " NLAM: " + str(self.nlam) + " NCOOL: " + str(self.ncool) + " NHEAT: " + str(self.nheat)
    output += "\n"
    output += "dust2gas: " + str(self.dust2gas)
    return output   
  
  
  def _getLineIdx(self,wl,ident=None):
    if self.lines == None: return None

    wls=numpy.array([line.wl for line in self.lines])
            
    if ident != None:      
      linestmp=[line for line in self.lines if line.ident==ident]
      if linestmp is not None and len(linestmp)>0:
        wlstmp=numpy.array([line.wl for line in linestmp])
        itmp=numpy.argmin(abs(wlstmp[:]-wl))
        # get the index of the whole line array. That should no find 
        # the exact one (e.g. used the exact same wavelength
        idx=numpy.argmin(abs(wls[:]-linestmp[itmp].wl))
        #print(self.lines[idx].ident)
        # check again
        # FIXME: causes problems with _lte lines (have the same wavelenghts)
        if self.lines[idx].ident != ident:
          print("ERROR: Something is wrong found: ident", self.lines[idx].ident," and wl ",self.lines[idx].wl,
                "for ",ident,wl)
          return None
        else:
          return idx
      else:
        print("WARN: No line found with ident", ident," and wl ",wl)
        return None
    else:      
      idx=numpy.argmin(abs(wls[:]-wl))
      
      if (abs(wls[idx]-wl)/wl) >0.01:
        print("WARN: No line found within 1% of the given wavelengeth:", wl)
        return None
      else:
        return idx 
    
      
  def getLine(self,wl,ident=None):
    '''
    Returns the spectral line closes to the given wavelentgh.
    
    Parameters
    ----------
    wl : float 
      the wavelength which is used for the search. `UNIT:` micron.
    ident : string, optional 
      A line identification string which should also be considered for the
      search. (default: `None`)
    
    Returns
    -------
    :class:`prodimopy.read.DataLine`   
      Returns `None` if no lines are included in the model.
                           
    '''
    
    idx = self._getLineIdx(wl,ident=ident)
    
    if idx is None: 
      return None
    else:
      return self.lines[idx]
  

  def getLineEstimate(self, ident, wl):
    '''
    Finds and returns the line estimate for the given species closest
    to the given wavelength.
    
    Parameters
    ----------
    ident : string 
      The line identification string.
    wl : float 
      the wavelength which is used for the search. `UNIT:` micron.
    
    Returns
    -------
    :class:`prodimopy.read.DataLineEstimate`   
      Returns `None` if the line is not found.
    
    Notes
    -----
      The parameters of this method are not consistent 
      with :meth:`~prodimopy.read.Data_ProDiMo.getLine`. This might be confusing.
      However, there is a reasong for this. The number of
      included line estimates in a |prodimo| model can be huge and just searching
      with the wavelenght might become slow. However, this probably can be improved. 
      
      To speed up the reading of `FlineEstimate.out` the extra radial info
      (see :class:`~prodimopy.read.DataLineEstimateRInfo` of all the line estimates 
      is not read. However, in this routine it is read for the single line estimate 
      found. One should use this routine to access the radial infor for a 
      single line estimate.
    '''
    found = -1
    i = 0 
    mindiff = 1.e20
    
    if self.lineEstimates is None:
      print("WARN: No lineEstimates are inluced (or not read) for this model.")
      return None
    
    # print ident, wl
    for le in self.lineEstimates: 
      if le.ident == ident:
        diff = abs(le.wl - wl)
        if diff < mindiff:
          found = i
          mindiff = diff      
      i = i + 1   
  
    if self.lineEstimates[found].rInfo == None:
      _read_lineEstimateRinfo(self, self.lineEstimates[found])
      
    return self.lineEstimates[found]
  
  def selectLineEstimates(self, ident,wlrange=None):
    """
    Returns a list of all line estimates (i.e. all transitions) for the 
    given line ident and/or in the given wavelength range.
    
    Parameters
    ----------    
    ident : string
      The line identification (species name) as defined in |prodimo|.
      The ident is not necessarely equal to the underlying chemial species 
      name (e.g. isotopologues, ortho-para, or cases like N2H+ and HN2+)
      if wlrange is set ident can be `None` in that case all line in the given 
      wlrange are returned
      
    wlrange : array_like
      All lines in the given wavelength range [start,end] and the given ident
      are returned. if the ident is `None` all lines are returned. 
      Default: `None` Units: micron
                
    Returns
    -------    
    list(:class:`prodimopy.read.DataLineEstimate`) :
      List of :class:`prodimopy.read.DataLineEstimate` objects, 
      or empty list if nothing was found.
      
    Notes
    -----
    In this routine the additional radial information of the line esitmate 
    (see :class:`~prodimopy.read.DataLineEstimateRInfo`) is not read. 
         
    """    
    lines = list()
    if wlrange is None:
      for le in self.lineEstimates:
        if le.ident == ident:
          lines.append(le)
    else:
      # TODO: this might can be done more efficient
      for le in self.lineEstimates:
        if le.wl>= wlrange[0] and le.wl<=wlrange[1]:
          if ident is None or le.ident == ident:
            lines.append(le)

    return lines

  def selectLines(self, ident):
    """
    Returns a list of all line fluxes for the 
    given line ident included in the line tranfer.
    
    Parameters
    ----------    
    ident : string
      The line identification (species name) as defined in |prodimo|.
      The ident is not necessarely equal to the underlying chemial species 
      name (e.g. isotopologues, ortho-para, or cases like N2H+ and HN2+)
                
    Returns
    -------    
    list(:class:`prodimopy.read.DataLine`) :
      List of :class:`prodimopy.read.DataLine` objects, 
      or empty list if nothing was found.
         
    """    
    if self.lines is None: return None
    
    lines = list()         
    for le in self.lines: 
      if le.ident == ident:
        lines.append(le)
  
    return lines
  
  
  def getAbun(self,spname):
    '''    
    Returns the abundance for a given species. 
        
    Parameters
    ----------
    spname : string
      The name of the chemical species as define in |prodimo|.     
                
    
    Returns
    -------
    array_like(float,ndim=2) 
      an array of dimension [nx,ny] array with species abundance or 
      `None` if the species was not found.
    
    Notes
    -----
      The abundance is given relative to the total hydrogen nuclei density 
      (i.e. ``nmol(spname)/nHtot``)
    
      A warning is printed in case the species was not found in spnames.
    '''        
    # check if it exists
    if not spname in self.spnames:
      print("WARN: getAbun: Species "+spname+" not found.")
      return None 
    
    return self.nmol[:,:,self.spnames[spname]]/self.nHtot

  def get_toomreQ(self,mstar=None):
    '''
    Returns the Toomre Q parameter as a function of radius.
    (for the midplane). 
    
    Q is given by 
    
    Q(r)=k(r)*cs(r)/(pi * G * sigma(r))
 
    for k we use the keplerian frequency, cs is the soundspeed (from the mdoel) 
    and sigma is the total gas surface density (both halfs of the disk). 
    
    Parameters
    ----------
    mstar : float
      The stellar mass in solar units (optional). 
      If `None` the value from the model is taken.
    
    Returns
    -------
    array_like(float,ndim=1) 
      an array of dimension [nx] with the Q param  
    '''
    if mstar is None:
      mstar=self.mstar
    
    mstarc=(mstar*u.M_sun).cgs.value
    grav=const.G.cgs.value
    r=(self.x[:,0]*u.au).cgs.value
    # isothermal soundspeed 
    cs=(self.soundspeed[:,0]*u.km/u.s).cgs.value
    # surface density factor two for both half of the disks
    sigma=2.0*(self.sdg[:,0]*u.g/u.cm**2).cgs.value

    # assuem keplerian rotation for Epicyclic frequency
    kappa=np.sqrt(grav*mstarc/r**3)
    Q=kappa*cs/(math.pi*grav*sigma)
    
    return Q
  
  def getSEDAnaMask(self,lam):
    '''
    Returns a numpy mask where only the grid points outside the emission 
    origin area for the given wavelength (continuum) are marked as invalid.
    
    This mask can be used in e.g. :func:`~prodimopy.Data_ProDiMo.avg_quantity`

    TODO: in case of optically thin emission only the points of one half 
          of the disk are considered. In that case an average quantity of this
          box(mask) will not be completely correct.
    
    Parameters
    ----------
    lam : float
      the wavelength for which the emission origin region should be determined:
      UNITS: micron

    Returns
    -------
    array_like(float,ndim=2)
      Numpy mask with `DIMS:` (nx,nz)
      
    '''
    sedAna = self.sed.sedAna
    x15=interp1d(sedAna.lams, sedAna.r15, bounds_error=False, fill_value=0.0, kind="linear")(lam)
    x85=interp1d(sedAna.lams, sedAna.r85, bounds_error=False, fill_value=0.0, kind="linear")(lam)
    
    mask=numpy.ones(shape=(self.nx, self.nz))
    for ix in range(self.nx):
      if self.x[ix,0]>=x15 and self.x[ix,0]<=x85:
        z85=interp1d(sedAna.lams, sedAna.z85[:,ix], bounds_error=False, fill_value=0.0, kind="linear")(lam)
        z15=interp1d(sedAna.lams, sedAna.z15[:,ix], bounds_error=False, fill_value=0.0, kind="linear")(lam)
        for iz in range(self.nz):
          if self.z[ix,iz]<=z15 and self.z[ix,iz]>=z85:
            mask[ix,iz]=0 
    
    return mask

  def getLineOriginMask(self,lineEstimate):
    '''
    Returns a numpy mask where only the grid points outside the emission 
    origin area for the given lineEstimate are marked as invalid.
    
    This mask can be used in e.g. :func:`~prodimopy.Data_ProDiMo.avg_quantity`
    
    Parameters
    ----------
    lineEstimate : `:class:`prodimopy.read.DataLineEstimate`
      a line Estimate object for which the operation should be done

    Returns
    -------
    array_like(float,ndim=2)
      Numpy mask with `DIMS:` (nx,nz) 
    '''
    fcfluxes=np.array([x.Fcolumn for x in lineEstimate.rInfo])
      
    Fcum=fcfluxes[:]
    for i in range(1,self.nx):
      Fcum[i]=Fcum[i-1]+fcfluxes[i]
    
    interx=interp1d(Fcum/Fcum[-1],self.x[:,0])
    x15=interx(0.15)
    x85=interx(0.85)
    
    mask=numpy.ones(shape=(self.nx, self.nz))
    for ix in range(self.nx):
      if (self.x[ix,0] >= x15 and self.x[ix,0] <= x85):
        for iz in range(self.nz):
          rinfo=lineEstimate.rInfo[ix]
        #print(rinfo.z15,rinfo.z85)
          if self.z[ix,iz] <= rinfo.z15 and self.z[ix,iz] >= rinfo.z85:
            mask[ix,iz]=0 
          
    return mask


  def get_KeplerOmega(self,mstar=None):
    '''
    Returns the Keplerian orbital frequency [1/s]
    (for the midplane). 
    
    omega=sqrt(G*mstar/r**3)
        
    Parameters
    ----------
    mstar : float
      The stellar mass in solar units (optional). 
      If `None` the value from the model is taken.
    
    Returns
    -------
    array_like(float,ndim=1) 
      Keplerian orbital frequency as function of radius [1/s]  
    '''
    if mstar is None:
      mstar=self.mstar
    
    mstarc=(mstar*u.M_sun).cgs.value
    grav=const.G.cgs.value
    r=(self.x[:,0]*u.au).cgs.value

    # assuem keplerian rotation for Epicyclic frequency
    omega=np.sqrt(grav*mstarc/r**3)
    
    return omega

  def avg_quantity(self,quantity,weight="gmass",weightArray=None,mask=None):
    '''
    Calculates the weighted mean value for the given field (e.g. td). 
    Different weights can be used (see `weight` parameter)
    
    TODO: option for region (e.g. a box)  
    
    Parameters
    ----------
    quantity : array_like(float,ndim=2)  
      the quantity to average. `DIMS:` (nx,nz).
    weight: str
      option for the weight. Values: `gmass` (gass mass) `dmass` (dust mass) or 
      `vol` (Volume). DEFAULT: `gmass`
  
    '''
    
    vol=np.ma.masked_array(self.vol,mask=mask)
    quantitym=np.ma.masked_array(quantity,mask=mask)
    
    if weightArray is not None:
      wA=np.ma.masked_array(weightArray,mask=mask)
      mean=np.sum(np.multiply(wA,quantitym))
      return mean/np.sum(wA)
    else: 
      if weight is "vol":
        mean=np.sum(np.multiply(vol,quantitym))
        return mean/np.sum(vol)
      else:
        if weight is "dmass":
          rho=np.ma.masked_array(self.rhod,mask=mask)
        else:
          rho=np.ma.masked_array(self.rhog,mask=mask)
        
        mass=np.multiply(vol,rho)
        mean=np.sum(np.multiply(mass,quantitym))
        mean=mean/np.sum(mass)
        return mean


class DataLineProfile():
  
  def __init__(self, nvelo):           
    self.nvelo = nvelo
    self.velo = numpy.zeros(shape=(self.nvelo))  # *u.km/u.s
    self.flux = numpy.zeros(shape=(self.nvelo))  # *u.Jansky
    self.flux_conv = numpy.zeros(shape=(self.nvelo))  # *u.Jansky
    
  def __str__(self):
    return "nvelo: " + str(self.nvelo)
        

class DataLine(object):
  """
  Data container for a single spectra line.
  """   
  def __init__(self):
    self.wl = 0.0
    self.frequency = 0.0 # is also taken from the line_flux.out
    self.prodimoInf = ""
    self.species = ""
    self.ident = ""
    # the flux in [W/m^2]
    self.flux = 0.0
    # the continuum in Jansky
    self.fcont = 0.0 
    # line profile of type data_line_profile
    self.profile = None 
  
  def __str__(self):
    text = (self.species + "/" + 
          self.prodimoInf + "/" + 
          self.ident + "/" + 
          str(self.wl) + "/" + 
          str(self.frequency) + "/" + 
          str(self.flux) + "/" + 
          str(self.fcont))
    return text 
  
  @property  
  def flux_Jy(self):    
    '''
    The flux of the line `UNIT: Jansky km s^-1`
    '''
    return _flux_Wm2toJykms(self.flux,self.frequency)


class DataLineObs(DataLine):
  '''
  Holds the observational data for one line.
  '''  
  def __init__(self, flux, flux_err, fwhm, fwhm_err, flag):
    super(DataLineObs, self).__init__()
    self.flux = flux
    self.flux_err = flux_err
    self.fwhm = fwhm
    self.fwhm_err = fwhm_err
    self.flag = flag.lower()
    self.profile=None
    self.profileErr=None 
    

class DataLineEstimate(object):
  '''
  Data container for a single line estimate.
  '''
  def __init__(self, ident, wl, jup, jlow, flux):
    """
    Parameters
    ----------
    ident : string
    wl : float
    jup : int
    jlow : int
    flux : float
        
    """""    
    self.ident = ident 
    """ string :
      The line identifications (species)
    """
    self.wl = wl
    """ float :
      The wavelength of the line.
      `UNIT: micron`
    """    
    self.jup = jup
    self.jlow = jlow
    self.flux = flux    
    self.rInfo = None  
    """ :class:`prodimopy.read.DataLineEstimateRinfo`
      The extra radial information of the line.
    """
    self.__posrInfo = None  # stores the position of the radial information to access is later on     

  @property
  def frequency(self):
    '''
    Frequency of the line in [GHz]. 
    '''
    return (self.wl * u.micrometer).to(u.GHz, equivalencies=u.spectral()).value

  @property
  def flux_Jy(self):    
    '''
    The flux of the line `UNIT: Jansky km s^-1`
    '''
    return _flux_Wm2toJykms(self.flux,self.frequency)
    
  def __str__(self):
    text = (self.ident + "/" + 
          str(self.wl) + "/" + 
          str(self.jup) + "/" + 
          str(self.jlow) + "/" + 
          str(self.flux))
    return text 


class DataLineEstimateRInfo(object):
  '''
  Data container for the extra radial info for a single line estimate.  
  '''
  def __init__(self, ix,iz, Fcolumn, tauLine, tauDust, z15, z85):   
    self.ix = ix 
    self.iz = iz 
    self.Fcolumn = Fcolumn
    self.tauLine = tauLine
    self.tauDust = tauDust
    self.z15 = z15          
    self.z85 = z85     

class DataGas(object):
  '''
  Holds the data for the gas (mainly from dust_opac.out)
  TODO: currently only the init cs is read, not the x,y data
  '''
  def __init__(self, nlam):
    self.nlam = nlam
    self.lam = numpy.zeros(shape=(nlam))
    self.energy = numpy.zeros(shape=(nlam))  # for convinence as often energy is also used for gas [eV]
    self.abs_cs = numpy.zeros(shape=(nlam))  # unit cm^2/H
    self.sca_cs = numpy.zeros(shape=(nlam))  # unit cm^2/H
  

class DataDust(object):
  '''
  Holds the data for the dust (mainly from dust_opac.out)
  TODO: dust composition is not yet read
  '''
  def __init__(self, amin, amax, apow, nsize, nlam):
    self.amin = amin
    self.amax = amax
    self.apow = apow
    self.nsize = nsize
    self.lam = numpy.zeros(shape=(nlam))
    self.energy = numpy.zeros(shape=(nlam))  # for convinence e.g. the X-ray range
    self.kext = numpy.zeros(shape=(nlam))  # in cm^2 g^-1
    self.kabs = numpy.zeros(shape=(nlam))
    self.ksca = numpy.zeros(shape=(nlam))  
    self.ksca_an = numpy.zeros(shape=(nlam))  # tweaked anisotropic scattering 
    self.kextcs = numpy.zeros(shape=(nlam))  # cross setion cm^2
    self.kabscs = numpy.zeros(shape=(nlam))  # cross setion cm^2
    self.kscacs = numpy.zeros(shape=(nlam))  # cross setion cm^2
    self.kscacs_an = numpy.zeros(shape=(nlam))  # tweaked anisotropic scattering cs   

class DataElements(object):
  '''
  Data for the Element abundances (the input). 
  
  Holds the data from `Elements.out`.
  
  The data is stored as OrderedDict with the element names as keys.
  
  Attributes
  ----------  
  '''
  def __init__(self):    
    self.abun=OrderedDict()
    """    
    OrderedDict :
      Ordered Dictionary holding the element abundaces realtive to hydrogen
    """
    self.abun12=OrderedDict()
    """
    OrderedDict :
      abundances in the +12 unit 
    """
    self.massRatio=OrderedDict()
    """
    OrderedDict :
      mass ratios 
    """
    self.amu=OrderedDict()
    """
    OrderedDict :
      atomic mass unit
    """


class DataContinuumObs(object):
  '''
  Holds the observational data for the continuum (the dust).
  
  Holds the photometric data, spectra (experimental) and radial profiles
  (experimental).
    
  '''
  def __init__(self, nlam=None):
    if nlam is not None:     
      self.nlam = nlam      
      self.lam = numpy.zeros(shape=(nlam))
      self.nu = numpy.zeros(shape=(nlam))
      self.fnuErg = numpy.zeros(shape=(nlam))
      self.fnuJy = numpy.zeros(shape=(nlam))
      self.fnuErgErr = numpy.zeros(shape=(nlam))
      self.fnuJyErr = numpy.zeros(shape=(nlam))
      self.flag=numpy.empty(nlam, dtype="U2")
    else:
      self.nlam = None      
      self.lam = None
      self.nu = None
      self.fnuErg = None
      self.fnuJy = None
      self.fnuErgErr = None
      self.fnuJyErr = None
      self.flag=None
            
    self.specs = None  # holds a list of spectra if available (wl,flux,error) 
    self.R_V = None
    self.E_BV = None
    self.A_V = None
    
    self.radprofiles=None
    

class DataSED(object):
  '''
  Holds the data for the Spectral Energy Distribution (SED).  
  '''
  def __init__(self, nlam, distance, inclination):
    self.nlam = nlam
    self.distance = distance
    self.inclination = inclination
    
    # the bolometric luminosity will be calculated by integrating over the whole
    # frequency range, considering also the distance) 
    # units are Solar luminosities    
    self.Lbol=None
    # is also calculated in read_sed
    self.Tbol=None                   
    self.lam = numpy.zeros(shape=(nlam))
    self.nu = numpy.zeros(shape=(nlam))
    self.fnuErg = numpy.zeros(shape=(nlam))
    self.nuFnuW = numpy.zeros(shape=(nlam))
    self.fnuJy = numpy.zeros(shape=(nlam))
    #self.nuFnu = numpy.zeros(shape=(nlam))
    
    # the analysis data
    self.sedAna=None

  def setLbolTbol(self):
    """
    Calculates the bolometric temperature and lumionsity for the given 
    values of the SED.
    
    quite approximate at the moment.
    """
    
    # FIXME: correctness needs to be verified
    # caluclate the bolometric luminosity
    # *-1.0 because of the ordering of nu
      #calculate Tbol see Dunham 2008 Equation 6
    # check what is here the proper value
    # in particular Tbol is very sensitive to this value
    # Lbol does not seem to change much (e.g. if 0.1 is used instead)
    # for the moment exclude the shortest wavelength as these is most likely scattering
    mask=self.lam>0.2  
    
    self.Lbol=numpy.trapz(self.fnuErg[mask],x=self.nu[mask])*-1.0
    self.Lbol=self.Lbol*4.0*math.pi*(((self.distance*u.pc).to(u.cm)).value)**2    
    self.Lbol=self.Lbol/(const.L_sun.cgs).value    
    
    #print(numpy.trapz((sed.nu*sed.fnuErg),x=sed.nu)) 
    #print(numpy.trapz(sed.fnuErg,x=sed.nu))
    self.Tbol=1.25e-11*numpy.trapz((self.nu[mask]*self.fnuErg[mask]),x=self.nu[mask])/numpy.trapz(self.fnuErg[mask],x=self.nu[mask])  


class DataSEDAna(object):
  '''
  Holds the analysis data for the Spectral Energy Distribution (SEDana.out).  
  '''
  def __init__(self,nlam,nx):
    self.lams=numpy.zeros(shape=(nlam))
    self.r15=numpy.zeros(shape=(nlam))
    self.r85=numpy.zeros(shape=(nlam))
    self.z15=numpy.zeros(shape=(nlam,nx))
    self.z85=numpy.zeros(shape=(nlam,nx))


class DataBgSpec(object):
  '''
  Backgound field input spectrum  
  '''
  def __init__(self, nlam):
    self.nlam = nlam
    self.lam = numpy.zeros(shape=(nlam))
    self.nu = numpy.zeros(shape=(nlam))    
    self.Inu = numpy.zeros(shape=(nlam))

class DataStarSpec(object):
  '''
  Stellar input spectrum  
  '''
  def __init__(self, nlam, teff,r,logg,luv):
    self.nlam = nlam
    self.teff=teff
    self.r=r
    self.logg=logg
    self.luv=luv
    self.lam = numpy.zeros(shape=(nlam))
    self.nu = numpy.zeros(shape=(nlam))    
    self.Inu = numpy.zeros(shape=(nlam))    

class Chemistry(object):
  """ 
  Data container for chemistry analysis. 
  
     
         
  """  
  def __init__(self, name):          
    """
    Parameters
    ----------
    name : string
      The name of the model (can be empty).  

    Attributes
    ----------
          
    """
    self.name = name 
    """ string :
    The name of the model (can be empty)
    """
    self.__tarfile=None
    self.farray = None        
    """ array : 
    array with the formation rate of the main formation reaction for the selected species for each grid point
    """
    self.darray = None        
    """ array : 
    array with the formation rate of the main destruction reaction for the selected species for each grid point
    """
    self.sorted_form_info = None        
    """ list : 
    list of formation processes in descending order according to their max rate value
    """
    self.sorted_dest_info = None        
    """ list : 
    list of destruction processes in descending order according to their max rate value
    """
    self.main_dest_grid = None        
    """ list : 
    list of main destruction reaction for the selected species for each grid point
    index corresponds to first index in sorted_form_info
    has the same order as grid points in chemanalysis.out
    """
    self.main_form_grid = None        
    """ list : 
    list of main formation reaction for the selected species for each grid point
    index corresponds to first index in sorted_dest_info
    has the same order as grid points in chemanalysis.out
    """
    
  def __str__(self):
    output = "Info chemistry: \n"

    return output   


def read_prodimo(directory=".", name=None, readlineEstimates=True,readObs=True, 
                 filename="ProDiMo.out", filenameLineEstimates="FlineEstimates.out", 
                 filenameLineFlux="line_flux.out",
                 td_fileIdx=None):
  '''
  Reads in all (not all yet) the output of a ProDiMo model from the given model directory.
  
  Parameters
  ----------
  directory : str 
    the directory of the model (optional).
  name : str 
    an optional name for the model (e.g. can be used in the plotting routines) 
  readlineEstimates : boolean 
    read the line estimates file (can be slow, so it is possible to deactivate it)
  readObs : boolean
    try to read observational data (e.g. SEDobs.dat ...)
  filename : str 
    the filename of the main prodimo output 
  filenameLineEstimates : str 
    the filename of the line estimates output 
  filenameLineFlux : str 
    the filename of the line flux output
  td_fileIdx : str 
    in case of time-dependent model the index of a particular output age can be provided (e.g. "03")
  
  
  Returns
  -------
  :class:`prodimopy.read.Data_ProDiMo`
    the |prodimo| model data or `None` in case of errors.
  '''
  # FIXME: tha reading in via a tar file still needs some work, disable it here for the momen
  # should be a parameter of this routine at some point
#  tarfile : string
#    the path to an e.g. tar file (also compressed, everything which works with the python tarfile). 
#    The routine tries to read the files from the archive also considering the `directory` path
  startAll=timer()
  
  tarfile=None
  
  # guess a name if not set
  if name == None:
    if directory==None or directory=="." or directory=="":
      dirfields=os.getcwd().split("/")
    else:  
      dirfields = directory.split("/")
    
    name = dirfields[-1]  
  
  # if td_fileIdx is given alter the filenames so that the timedependent
  # files can be read. However this would actually also workd with other kind
  # of indices as td_fileIdx is a strong
  if td_fileIdx != None:    
    rpstr="_"+td_fileIdx+".out"
    filename=filename.replace(".out",rpstr)
    filenameLineEstimates=filenameLineEstimates.replace(".out",rpstr)
    filenameLineFlux=filenameLineFlux.replace(".out",rpstr)  
    
  f,dummy = _getfile(filename, directory, tarfile)
  
  # FIXME: with this the calling rouinte can continue
  # The calling routine has to take care of that
  # But this is not very nice 
  if f is None: return None
  
  # read all date into the memory
  # easier to handle afterwards  
  lines = f.readlines()
  f.close()
  
  # start of the array data block
  idata = 24
   
  data = Data_ProDiMo(name)
  
  data.mstar = float(lines[0].split()[1])
  data.dust2gas = float(lines[9].split()[1])
  
  strs = lines[21].split()  
  data.nx = int(strs[1])
  data.nz = int(strs[2])
  data.nspec = int(strs[3]) + 1  # +1 for the electron
  data.nheat = int(strs[4])
  data.ncool = int(strs[5])
  data.nlam = int(strs[6])
  
  # TODO: move this to the constructure which takes nx,nz
  data.x = numpy.zeros(shape=(data.nx, data.nz))
  data.z = numpy.zeros(shape=(data.nx, data.nz))
  data.lams=numpy.zeros(shape=(data.nlam))
  data.AV = numpy.zeros(shape=(data.nx, data.nz))
  data.AVrad = numpy.zeros(shape=(data.nx, data.nz))
  data.AVver = numpy.zeros(shape=(data.nx, data.nz))  
  data.NHver = numpy.zeros(shape=(data.nx, data.nz))
  data.NHrad = numpy.zeros(shape=(data.nx, data.nz))
  data.d2g = numpy.zeros(shape=(data.nx, data.nz))
  data.tg = numpy.zeros(shape=(data.nx, data.nz))
  data.td = numpy.zeros(shape=(data.nx, data.nz))
  data.nd = numpy.zeros(shape=(data.nx, data.nz))
  data.soundspeed = numpy.zeros(shape=(data.nx, data.nz))
  data.rhog = numpy.zeros(shape=(data.nx, data.nz))
  data.pressure = numpy.zeros(shape=(data.nx, data.nz))
  data.tauchem = numpy.zeros(shape=(data.nx, data.nz))
  data.taucool = numpy.zeros(shape=(data.nx, data.nz))
  data.taudiff = numpy.zeros(shape=(data.nx, data.nz))  
  data.heat_mainidx =numpy.zeros(shape=(data.nx, data.nz),dtype=numpy.int32)
  data.cool_mainidx =numpy.zeros(shape=(data.nx, data.nz),dtype=numpy.int32)
  data.nHtot = numpy.zeros(shape=(data.nx, data.nz))
  data.damean = numpy.zeros(shape=(data.nx, data.nz))
  data.Hx=numpy.zeros(shape=(data.nx, data.nz))
  data.tauX1 = numpy.zeros(shape=(data.nx, data.nz))
  data.tauX5 = numpy.zeros(shape=(data.nx, data.nz))
  data.tauX10 = numpy.zeros(shape=(data.nx, data.nz))  
  data.zetaX = numpy.zeros(shape=(data.nx, data.nz))
  data.chi = numpy.zeros(shape=(data.nx, data.nz))
  data.chiRT = numpy.zeros(shape=(data.nx, data.nz))
  data.kappaRoss = numpy.zeros(shape=(data.nx, data.nz))  
  data.radFields=numpy.zeros(shape=(data.nx,data.nz,data.nlam))
  data.zetaCR = numpy.zeros(shape=(data.nx, data.nz))
  data.zetaSTCR = numpy.zeros(shape=(data.nx, data.nz))
  data.nmol = numpy.zeros(shape=(data.nx, data.nz, data.nspec))
  data.cdnmol = numpy.zeros(shape=(data.nx, data.nz, data.nspec)) 
  data.heat = numpy.zeros(shape=(data.nx, data.nz, data.nheat))
  data.cool = numpy.zeros(shape=(data.nx, data.nz, data.ncool))
  data.heat_names=list()
  data.cool_names=list()
  
  # Make some checks for the format 
  # new EXP format for x and z:
  newexpformat=lines[idata].find("E",0,25)>=0
  # FIXME: that is not very nice
  #        make at least some checks if the output format has changed or something  
  # number of fixed fields in ProDiMo.out (before heating and cooling rates)
  nfixFields = 21  
  # index after the J data/fields in ProDiMo 
  
  iACool = nfixFields + data.nheat + data.ncool  
  iASpec= iACool + 1 + data.nspec
  iAJJ = iASpec + data.nlam
    
  # read the species names, these are taken from the column titles
  colnames = lines[idata - 1]
  
  # FIXME: Again hardconding. Might be better to read the Species names
  # from Species.out
  # Assume that the first species is e- and search for that
  # that is more flexible in case other names change 
  #specStart = 233 + data.nheat * 13 + data.ncool * 13 + 13
  specStart=colnames.find("        e-")
  spnames = colnames[specStart:specStart + (data.nspec) * 13]
  spnames = spnames.split()
  if (len(spnames) != data.nspec): 
    print("ERROR: something is wrong with the number of Species!")
    return None  
  
  # empty dictionary
  data.spnames = OrderedDict()
  # Make a dictionary for the spnames
  for i in range(data.nspec):    
    data.spnames[spnames[i]] = i
    
  # read the heating and cooling names
  iheat=idata+data.nx*data.nz+2
  for i in range(data.nheat):    
    data.heat_names.append(lines[iheat+i][3:].strip())

  icool=idata+data.nx*data.nz+2+data.nheat+2
  for i in range(data.ncool):    
    data.cool_names.append(lines[icool+i][3:].strip())
  
  # read the band wavelenghts
  iwls=idata+data.nx*data.nz+2+data.nheat+2+data.ncool+2
  for i in range(data.nlam):      
    data.lams[i]=float((lines[iwls+i].split())[1])  
    
  
  i = 0  
  for iz in range(data.nz):
    zidx = data.nz - iz - 1
    for ix in range(data.nx):

      # stupid workaround for big disks/envelopes where x,y can be larger than 10000 au
      # there is no space between the numbers for x and z so always add one if none is there
      # this might can be removed in the future as the newest ProDiMo version use exp format now
      if (not newexpformat) and lines[idata+i][20]!= " ":
        line=lines[idata + i][:20]+" "+lines[idata + i][20:]
      else:
        line=lines[idata + i]
            
      # This line is what eats the time, but using genfromtxt for the ProDiMo.out
      # does not seem to be faster
      fields=numpy.fromstring(line,sep=" ")
      
      data.x[ix, zidx] = fields[2]
      data.z[ix, zidx] = fields[3]
      data.NHrad[ix, zidx] = fields[4]
      data.NHver[ix, zidx] = fields[5]
      data.AVrad[ix, zidx] = fields[6]
      data.AVver[ix, zidx] = fields[7]
        
      data.nd[ix, zidx] = fields[8]
      data.tg[ix, zidx] = fields[9]
      data.td[ix, zidx] = fields[10]
      data.soundspeed[ix, zidx] = fields[11]
      data.rhog[ix, zidx] = fields[12]
      data.pressure[ix, zidx] = fields[13]
      data.chi[ix, zidx] = fields[15]
      data.tauchem[ix,zidx] = fields[16]
      data.taucool[ix,zidx] = fields[17]
      data.taudiff[ix,zidx] = fields[18]
      data.heat_mainidx[ix,zidx]=fields[19]
      data.cool_mainidx[ix,zidx]=fields[20]
      data.heat[ix, zidx, :] = fields[nfixFields:nfixFields+data.nheat]
      data.cool[ix, zidx, :] = fields[(nfixFields+data.nheat):(nfixFields+data.nheat+data.ncool)]
      data.nHtot[ix, zidx] = fields[iACool]
      data.nmol[ix, zidx, :] = fields[iACool + 1:iASpec]
        
      data.radFields[ix,zidx,:]=fields[iASpec :iAJJ]
        
      data.chiRT[ix, zidx] = fields[iAJJ]
      data.kappaRoss[ix, zidx] = fields[iAJJ+1]            
      data.damean[ix, zidx] = fields[iAJJ + 3]
      data.d2g[ix, zidx] = fields[iAJJ + 4]
      data.Hx[ix, zidx] = fields[iAJJ + 5]
      data.tauX1[ix, zidx] = fields[iAJJ + 6]
      data.tauX5[ix, zidx] = fields[iAJJ + 7]
      data.tauX10[ix, zidx] = fields[iAJJ + 8]
      data.zetaX[ix, zidx] = fields[iAJJ + 9]
      data.zetaCR[ix, zidx] = fields[iAJJ + 16]      
      if len(fields) > (iAJJ + 17): data.zetaSTCR[ix, zidx] = fields[iAJJ + 17] 

      i = i + 1

  # derived quantitites
  data.rhod = data.rhog * data.d2g
  
  # AV like defined in the prodimo idl script  
  for ix in range(data.nx):
    for iz in range(data.nz):
      data.AV[ix,iz] = numpy.min([data.AVver[ix,iz],data.AVrad[ix,iz],data.AVrad[data.nx-1,iz]-data.AVrad[ix,iz]])
  
  # Read FlineEstimates.out  
  if readlineEstimates == True:
    read_lineEstimates(directory, data, filename=filenameLineEstimates,tarfile=tarfile)
  else:
    data.lineEstimates = None
    
  data.elements=read_elements(directory,tarfile=tarfile)
     
  data.dust = read_dust(directory,tarfile=tarfile)

  fileloc=directory + "/dust_opac_env.out"
  if os.path.isfile(fileloc):
    data.env_dust = read_dust(directory,filename="dust_opac_env.out",tarfile=tarfile)  

  data.starSpec = read_starSpec(directory,tarfile=tarfile)

  if os.path.exists(directory + "/gas_cs.out"):
    data.gas = read_gas(directory,tarfile=tarfile)
  
  if os.path.exists(directory + "/"+filenameLineFlux):
    data.lines = read_linefluxes(directory,filename=filenameLineFlux,tarfile=tarfile)  

  if os.path.exists(directory + "/SED.out"):
    data.sed = read_sed(directory,tarfile=tarfile)

  if os.path.exists(directory + "/Species.out"):
    # data is filled in the routine
    read_species(directory,data,tarfile=tarfile)
    #print("INFO: Calc total species masses")
    #calc_spmassesTot(data)        

  if readObs:
    data.sedObs=read_continuumObs(directory)
    
  print("INFO: Reading time: ","{:4.2f}".format(timer()-startAll)+" s")  
  
  start=timer()
  # calculate the columnd densitis
  print("INFO: Calculate column densities")
  calc_columnd(data)
  print("INFO: Calculate surface densities")
  calc_surfd(data)
  print("INFO: Calculate volumes")
  _calc_vol(data)

  # set muH
  data.muH=data.rhog[0,0]/data.nHtot[0,0]
  print("INFO: Calc time: ","{:4.2f}".format(timer()-start)+" s")  
  
  # Test for volume stuff ... total integrated mass
  #mass=np.sum(np.multiply(data.vol,data.rhog))
  #print("Total gas mass", (mass*u.g).to(u.M_sun))
  print(" ")

  return data

def read_elements(directory,filename="Elements.out",tarfile=None):
  '''
  Reads the Elements.out file. Currently only the species masses are read.
  Stores the species masses (unit g) in pdata.spmasses
  Also adds the electron "e-"
  
  Parameters
  ----------
  directory : str 
    the directory of the model
        
  filename: str 
    an alternative Filename
    
  Returns
  -------
  :class:`~prodimopy.read.DataElements`
    the Elements model data or `None` in case of errors.

  '''
  
  f,dummy = _getfile(filename, directory, tarfile)
  if f is None:  return None
  
  # skip the first line
  nelements=int(f.readline())
  
  elements=DataElements()  
  f.readline()
    
  for i in range(nelements):
    line=f.readline()
    fields=line.strip().split()
    name=fields[0].strip()
    elements.abun12[name]=float(fields[1])
    elements.abun[name]=10.0**(float(fields[1])-12.0)
    elements.amu[name]=float(fields[2])
    elements.massRatio[name]=float(fields[3])
    
  return elements


def read_species(directory,pdata,filename="Species.out",tarfile=None):
  '''
  Reads the Species.out file. Currently only the species masses are read.
  Stores the species masses (unit g) in pdata.spmasses
  Also adds the electron "e-"
  
  Parameters
  ----------
  directory : str 
    the directory of the model
    
  pdata : :class:`prodimopy.read.Data_ProDiMo` 
    the ProDiMo Model data structure, where the data is stored
    
  filename: str 
    an alternative Filename
  '''
  f,dummy = _getfile(filename, directory, tarfile)
  if f is None:  
    pdata.spmasses=None    
    return None
  
  # skip the first line
  f.readline()
  
  pdata.spmasses=OrderedDict() # empty dictonary
  for line in f:
    fields=line.strip().split()
    spname=fields[2].strip()
    mass=(float(fields[3].strip())*u.u).cgs.value
    pdata.spmasses[spname]=mass  

  pdata.spmasses["e-"]=const.m_e.cgs.value


def read_lineEstimates(directory, pdata, filename="FlineEstimates.out",tarfile=None):
  '''
  Read FlineEstimates.out Can only be done after ProDiMo.out is read.
  
  Parameters
  ----------
  directory : str 
    the directory of the model
    
  pdata : :class:`prodimopy.read.Data_ProDiMo` 
    the ProDiMo Model data structure, where the data is stored
    
  filename: str 
    an alternative Filename

  '''
  f,rfile=_getfile(filename, directory, tarfile, binary=True)
  if f is None:
    pdata.lineEstimates = None
    return None
  else:
    pdata.__fpFlineEstimates = rfile   
    pdata.__tarfile=tarfile

  # check for version currently just check if it exist
  line = f.readline().decode()  
  version = 1
  if len(line) > 68:
    # position of version 
    posver = line.find("version")
    if posver >= 0:
      version = int(line[posver + len("version"):].strip())       
    
  nlines = int(f.readline().strip())
  f.readline()
  
  pdata.lineEstimates = list()  
  nxrange = list(range(pdata.nx))
  startBytes = 0

  for i in range(nlines):
    # has to be done in fixed format
    # FIXME: it can be that nline is not really equal the nubmer of available line
    # this ir probably a bug in ProDiMo but maybe intended (see 
    # OUTPUT_FLINE_ESTIMATE in ProDiMo, Therefore add a check     
    line = f.readline()            
    if not line: break
    
    line=line.decode()    
    # FIXME: has now also a version tag!! this is for the new version
    if version == 1:
      le = DataLineEstimate((line[0:9]).strip(), float(line[10:28]), int(line[29:32]), int(line[33:37]), float(line[38:49]))
    elif version == 2:
      try:
        le = DataLineEstimate((line[0:9]).strip(), float(line[10:28]), int(line[29:34]), int(line[35:40]), float(line[41:53]))
      except ValueError as err:
        print("Conversion problems: {0}".format(err))
        print("Line (i,text): ",i,line)
        print("Field: ",line[0:9],line[10:28],line[29:34],line[35:40],line[41:53])
        raise err
    else:
      raise ValueError('Unknown version of FlineEstimates.out! version=' + str(version))          

    # Fine out the number of bytes in one radial line to use seek in the 
    # next iterations
    if i == 0:
      start = f.tell()
      le.__posrInfo = start  # store the position for reading this info if required
      for j in nxrange:      
        f.readline()
      startBytes = f.tell() - start
    else:
      le.__posrInfo = f.tell()
      f.seek(startBytes, 1)

    pdata.lineEstimates.append(le)

  f.close()

def _read_lineEstimateRinfo(pdata, lineEstimate):
  '''
  Reads the additional Rinfo data for the given lineEstimate.
  '''  
  f,dummy=_getfile(pdata.__fpFlineEstimates, None, pdata.__tarfile, binary=True)
  if f is None: return None
  
  f.seek(lineEstimate.__posrInfo, 0)
  nxrange = list(range(pdata.nx))
  lineEstimate.rInfo = list()
  for j in nxrange:            
    fieldsR = f.readline().decode().split()
    ix= int(fieldsR[0].strip())    
    iz = int(fieldsR[1].strip())
    Fcolumn = float(fieldsR[2])
    try:
      tauLine = float(fieldsR[3])
    except ValueError as e:
      # print "read_lineEstimates line: ", le.ident  ," error: ", e
      tauLine = 0.0
    tauDust = float(fieldsR[4])
    z15 = (float(fieldsR[5])*u.cm).to(u.au).value
    z85 = (float(fieldsR[6])*u.cm).to(u.au).value
    # ix -1 and iz-1 because in python we start at 0
    rInfo = DataLineEstimateRInfo(ix-1,iz-1, Fcolumn, tauLine, tauDust, z15, z85)
    lineEstimate.rInfo.append(rInfo)
      
  f.close()
  
def  read_linefluxes(directory, filename="line_flux.out",tarfile=None):
  """ Reads the line fluxes output.
  
  Parameters
  ----------
  directory : str
    the model directory.
    
  filename : str
    the filename of the output file (optional).
    
  Returns
  -------
  list(:class:`prodimopy.read.DataLine`) 
    List of lines.
  """  
  f,dummy = _getfile(filename, directory, tarfile)
  if f is None:  return None
  
  records = f.readlines()
  f.close()
  
  dist = float(records[0].split("=")[1])
  incl = float(records[1].split("=")[1])
  nlines = int(records[2].split("=")[1])
  nvelo = int(records[3].split("=")[1])
  
  # number of records in the file 
  nrecords = len(records)
  # print dist,incl,nlines,nvelo,nrecords
  
  
  lines = list()
  
  # loop over all lines
  pos = 5
  for i in range(nlines):    
    # read the data for one line      
    line = DataLine()
    rec = records[pos]      
    try:  
      line.species = (rec[10:20]).strip()
      line.ident=line.species
      line.prodimoInf = rec[21:36].strip()
      line.wl = float(rec[43:54].strip())  # *u.um
      line.frequency = float(rec[63:76].strip())  # *u.GHz
    except:
      #print("read_linefluxes: try new format") 
      line.species = (rec[10:20]).strip()
      line.ident=line.species
      line.prodimoInf = rec[21:47].strip()
      line.wl = float(rec[48:64].strip())  # *u.um
      line.frequency = float(rec[74:90].strip())  # *u.GHz
      
           
    rec = records[pos + 1]
    line.flux = float(rec.split("=")[1].strip())  # *u.Watt/u.m**2       
                    
    rec = records[pos + 2]
    line.fcont = float(rec.split("=")[1].strip())  # *u.Jansky                    
             
    # one line in the profile is for the continuum                        
    line.profile = DataLineProfile(nvelo - 1)
    
    for j in range(nvelo - 1):
      # skip the header line and the first line (continuum)?
      fields = records[pos + 5 + j].split()        
      line.profile.velo[j] = float(fields[0])  # *u.km/u.s
      line.profile.flux[j] = float(fields[1])  # *u.Jansky
      if (len(fields) > 2):
        line.profile.flux_conv[j] = float(fields[2])  # *u.Jansky
        
    lines.append(line)
    pos += (nvelo + 5)
  
  return lines

def read_lineObs(directory, nlines, filename="LINEobs.dat",tarfile=None):
  '''
  Reads the lineobs Data. the number of lines have to be known.
  '''
  f,dummy = _getfile(filename, directory, tarfile)
  if f is None:  return None
    
  records = f.readlines()
  f.close()
  
  linesobs = list()  
  versionStr=records[0].split()
  version=float(versionStr[0])
  
  for rec in records[2:2 + nlines]:  #        
    fields = rec.split()
    
    lineobs = DataLineObs(float(fields[0].strip()), \
                        float(fields[1].strip()), \
                        float(fields[2].strip()), \
                        float(fields[3].strip()), \
                        fields[4].strip())   
    
    # FIXME: not very nice
    # in case of upper limits flux might be zero in that case use sig. 
    if lineobs.flux <1.e-100:
      lineobs.flux=lineobs.flux_err
     
    linesobs.append(lineobs)
  
  # the additional data
  # check if there is actually more data
  if (len(records)>2 + nlines+1):    
    # FIXME: do this proberly (the reading etc. and different versions)  
    profile = (records[2 + nlines+1].split())[0:nlines]
    autoC = (records[2 + nlines+2].split())[0:nlines]
    vvalid = (records[2 + nlines+3].split())[0:nlines]
  
  #speres = (records[2 + nlines+4].split())[0:nlines]
     
    offset=5
    if version>2.0: offset=6
  
    # now go through the profiles  
    for i in range(nlines):       
      proffilename=records[offset+nlines+i+1].strip()
      if profile[i] == "T":       
        if "nodata"==proffilename: print("WARN: Something is wrong with line "+str(i))     
        linesobs[i].profile,linesobs[i].profileErr=read_lineObsProfile(proffilename,directory=directory)
    
  return linesobs
  
  
def read_lineObsProfile(filename,directory=".",tarfile=None):
  '''
  reads a line profile file which can be used for ProDiMo
  '''
  f,dummy = _getfile(filename, directory, tarfile)
  if f is None:  return None
  
  records = f.readlines()
  f.close()
  
  # First line number of velo points and central wavelength, which we do not
  # need here (I think this is anyway optional  
  nvelo = int(records[0].split()[0].strip())

  # skip header lines
  profile = DataLineProfile(nvelo)
  for i in range(2, nvelo + 2):
    fields = records[i].split()
    profile.velo[i - 2] = float(fields[0].strip())
    profile.flux[i - 2] = float(fields[1].strip())
    # also fill the convolved flux, which is just the same as the flux
    # for observations
    profile.flux_conv[i - 2] = profile.flux[i - 2]
    
  err=float(records[-1].split()[0].strip())
    
  return profile,err

def read_gas(directory,filename="gas_cs.out",tarfile=None):
  '''
  Reads gas_cs.out
  Returns an object of Type DataDust  
  ''' 
  f,dummy = _getfile(filename, directory, tarfile)
  if f is None:  return None
        
  nlam = int(f.readline().strip())  
  
  # skip one line  
  f.readline()  
      
  gas = DataGas(nlam)
      
  for i in range(nlam):
    fields = [float(field) for field in f.readline().split()]    
    gas.lam[i] = float(fields[0])
    gas.abs_cs[i] = float(fields[1])
    gas.sca_cs[i] = float(fields[2])
 
  gas.energy[:] = ((gas.lam[:] * u.micron).to(u.eV, equivalencies=u.spectral())).value  
 
  f.close()

  return gas  

def read_continuumObs(directory,filename="SEDobs.dat"):
  ''' 
  Reads observational continuum data (SED). 
  
  Reads the file SEDobs.dat (phtotometric points) and all files ending 
  with `*spec.dat` (e.g. the Spitzer spectrum) 
  
  FIXME: does not work yet with tar files
  '''
  contObs=None
  rfile = directory + "/"+filename
  if os.path.exists(rfile):    
    f = open(rfile, 'r')

    nlam = int(f.readline().strip())
    f.readline() # header line
    contObs = DataContinuumObs(nlam=nlam)
    for i in range(nlam):
      elems = f.readline().split()
      contObs.lam[i] = float(elems[0])
      contObs.fnuJy[i] = float(elems[1])
      contObs.fnuJyErr[i] = float(elems[2])
      contObs.flag[i] = str(elems[3])
      
    contObs.nu = (contObs.lam* u.micrometer).to(u.Hz, equivalencies=u.spectral()).value
    contObs.fnuErg = (contObs.fnuJy*u.Jy).cgs.value
    contObs.fnuErgErr = (contObs.fnuJyErr*u.Jy).cgs.value
  
  # check for the spectrum files
  fnames=glob.glob(directory+"/*spec.dat")
  if fnames is not None and len(fnames)>0:    
    if contObs is None:
      contObs=DataContinuumObs()
      
    contObs.specs=list()
    
  for fname in fnames:
    spec=numpy.loadtxt(fname, skiprows=3)    
    contObs.specs.append(spec)
    
  # check if there is some extinction data
  fn_ext= directory+"/extinct.dat" 
  if os.path.exists(fn_ext):
    if contObs is None:
      contObs=DataContinuumObs()
      
    fext= open(fn_ext,"r")
    fext.readline()
    contObs.E_BV=float(fext.readline())
    fext.readline()
    contObs.R_V=float(fext.readline())
    contObs.A_V=contObs.R_V*contObs.E_BV
    
    
#   # check if there is an image.in 
#   fn_images= directory+"/image.in"
#   if os.path.exists(fn_images):
#     if contObs is None:
#       contObs=DataContinuumObs()
#     
#     fext= open(fn_images,"r")
#     fext.readline()
#     fext.readline()
#     nimages=int(fext.readline())
#     positionAngle = float(fext.readline())
#     
#     for i in range(9):
#       line = f.readline()
    
    
  
  return contObs

def read_sed(directory,filename="SED.out",filenameAna="SEDana.out",tarfile=None):
  ''' 
  Reads the ProDiMo SED output including the analysis data.
  '''
  f,dummy = _getfile(filename, directory, tarfile)
  if f is None:  return None
  
  elems = f.readline().split()
  distance = float(elems[len(elems) - 1])
  # ignore the face-on SED
  f.readline()
  nlam = int(f.readline())
  f.readline()
  for i in range(nlam):
    f.readline()
    
  f.readline()
  elems = f.readline().split()  
  inclination = float(elems[len(elems) - 1])
  nlam = int(f.readline())
  f.readline()
  sed = DataSED(nlam, distance, inclination)
  for i in range(nlam):
    elems = f.readline().split()
    sed.lam[i] = float(elems[0])
    sed.nu[i] = float(elems[1])

    # FIXME: Workaround to catch strange values from ProDiMo .. should be fixed 
    # in ProDiMo
    try:    
      sed.fnuErg[i] = float(elems[2])
      sed.nuFnuW[i] = float(elems[3])
      sed.fnuJy[i] = float(elems[4])      
    except ValueError as err:
      print("WARN: Could not read value from SED.out: ", err)
      
  sed.setLbolTbol()
  f.close()
  
  # The analysis data, if it is not there not a big problem
  f,dummy = _getfile(filenameAna, directory, tarfile)
  if f is None:  
    sed.sedAna=None
    return sed
  
  nlam,nx=f.readline().split()
  nlam=int(nlam)
  nx=int(nx)
  
  sed.sedAna=DataSEDAna(nlam,nx)
  
  for i in range(nlam):
    elems=f.readline().split()
    sed.sedAna.lams[i]=float(elems[0])
    sed.sedAna.r15[i]=float(elems[1])
    sed.sedAna.r85[i]=float(elems[2])
    for j in range(nx):
      elems=f.readline().split()
      sed.sedAna.z15[i,j]=float(elems[0])
      sed.sedAna.z85[i,j]=float(elems[1])
      
  f.close()
    
  return sed

def read_starSpec(directory,filename="StarSpectrum.out",tarfile=None):
  ''' 
  Reads StarSpectrum.out
  '''
  f,dummy = _getfile(filename, directory, tarfile)
  if f is None:  return None
      
  teff = float((f.readline().split())[-1])
  elems=f.readline().split()  
  r=float(elems[-1])
  logg=float(elems[2])
  luv=float(f.readline().split()[-1])        
  nlam = int(f.readline())
  f.readline()

  starSpec = DataStarSpec(nlam,teff,r,logg,luv )
  for i in range(nlam):
    elems = f.readline().split()
    starSpec.lam[i] = float(elems[0])
    starSpec.nu[i] = float(elems[1])
    starSpec.Inu[i] = float(elems[2])
    
  f.close()  
  return starSpec 

def read_bgSpec(directory,filename="BgSpectrum.out",tarfile=None):
  ''' 
  Reads the BgSpectrum.out file. 
  
  Returns
  -------  
    :class:`prodimopy.read.DataBgSpec`
    the background spectra or `None` if not found.
  '''
  f,dummy = _getfile(filename, directory, tarfile)
  if f is None:  return None
             
  nlam = int(f.readline())
  f.readline()

  bgSpec = DataBgSpec(nlam)
  for i in range(nlam):
    elems = f.readline().split()
    bgSpec.lam[i] = float(elems[0])
    bgSpec.nu[i] = float(elems[1])
    bgSpec.Inu[i] = float(elems[2])
    
  return bgSpec 

def read_dust(directory,filename="dust_opac.out",tarfile=None):
  '''
  Reads dust_opac.out
  Returns an object of Type DataDust
  Dust not read the dust composition yet
  ''' 
  f,dummy = _getfile(filename, directory, tarfile)
  if f is None:  return None
          
  fields = [int(field) for field in f.readline().split()]
  ndsp = fields[0]      
  nlam = fields[1]  
  
  # skip three lines
  for i in range(ndsp):
    f.readline()
  
  # apow amax etc.
  strings = f.readline().split()  

  if len(strings) > 0:  
    amin = ((float(strings[6]) * u.cm).to(u.micron)).value
    amax = ((float(strings[7]) * u.cm).to(u.micron)).value
    apow = float(strings[8])
    nsize = int(strings[9])
    dust = DataDust(amin, amax, apow, nsize, nlam)
  else:
    dust = DataDust(-1, -1, 0.0, -1, nlam)

  
  f.readline()
      
  for i in range(nlam):
    fields = [float(field) for field in f.readline().split()]    
    dust.lam[i] = float(fields[0])
    dust.kext[i] = float(fields[1])
    dust.kabs[i] = float(fields[2])
    dust.ksca[i] = float(fields[3])      
    if len(fields) >4:  
      dust.ksca_an[i] = float(fields[5])  # skip kprn
      dust.kextcs[i] = float(fields[6])
      dust.kabscs[i] = float(fields[7])
      dust.kscacs[i] = float(fields[8])
      dust.kscacs_an[i] = float(fields[9])

  f.close()

  dust.energy[:] = ((dust.lam[:] * u.micron).to(u.eV, equivalencies=u.spectral())).value

  return dust  

def calc_NHrad_oi(data):
  '''
  Calculates the radial column density from out to in (border of model spact to star/center)
  
  TODO: move this to utils  
  '''    
  data.cdnmol = 0.0 * data.nmol
  for ix in range(data.nx):          
    for iz in range(data.nz - 2, -1, -1):  # from top to bottom        
      dz = (data.z[ix, iz + 1] - data.z[ix, iz])
      dz = dz * u.au.to(u.cm)
      nn = 0.5 * (data.nmol[ix, iz + 1, :] + data.nmol[ix, iz, :])
      data.cdnmol[ix, iz, :] = data.cdnmol[ix, iz + 1, :] + nn * dz
      
  # check if integration is correct 
  # for the total hydrogen column density the error is less than 1.e-3
  # that should be good enough for plotting
  #nHverC=data.cdnmol[:,:,data.spnames["H"]]+data.cdnmol[:,:,data.spnames["H+"]]+data.cdnmol[:,:,data.spnames["H2"]]*2.0
  #print(numpy.max(numpy.abs(1.0-nHverC[:,0]/data.NHver[:,0])))    

  NHradoi = 0.0 * data.nHtot   
  for ix in range(data.nx-2,1,-1):  # first ix point (ix=0= remains zero  
    r1= (data.x[ix+1, :]**2+data.z[ix+1, :]**2)**0.5
    r2 = (data.x[ix, :]**2+data.z[ix, :]**2)**0.5              
    dr = r1-r2
    dr = dr * u.au.to(u.cm)
    nn = 0.5 * (data.nHtot[ix+1, :] + data.nHtot[ix, :])
    NHradoi[ix, :] = NHradoi[ix+1, :] + nn * dr
      
  return NHradoi


def calc_surfd(data):
  '''
  Caluclates the gas and dust vertical surface densities at every point in the
  model.
  
  Only one half of the disk is considered. If one needs the total surface density
  simply multiply the value from the midplane (zidx=0) by two.
  '''
  data.sdg = 0.0 * data.rhog
  data.sdd = 0.0 * data.rhod
  for ix in range(data.nx):          
    for iz in range(data.nz - 2, -1, -1):  # from top to bottom        
      dz = (data.z[ix, iz + 1] - data.z[ix, iz])
      dz = dz * u.au.to(u.cm)
      nn = 0.5 * (data.rhog[ix, iz + 1] + data.rhog[ix, iz])
      data.sdg[ix, iz] = data.sdg[ix, iz + 1,] + nn * dz
      nn = 0.5 * (data.rhod[ix, iz + 1] + data.rhod[ix, iz])
      data.sdd[ix, iz] = data.sdd[ix, iz + 1,] + nn * dz
      

def _calc_vol(data):
  '''
  Inits the vol field (:attr:`~prodimpy.read.Data_ProDiMo.vol` for each individual grid point. 
  
  This is useful to estimate mean quantities which are weighted by volume
  but also by mass.
  
  The routine follows the same method as in the prodimo.pro (the IDL skript)
   
  '''
  
  tocm=(1.0*u.au).to(u.cm).value
  data.vol=numpy.zeros(shape=(data.nx, data.nz))
  fourthreepi=4.0*math.pi/3.0
  for ix in range(data.nx):
    ix1 = np.max([0,ix-1])
    ix2 = ix
    ix3 = np.min([data.nx-1,ix+1])
    x1=math.sqrt(data.x[ix1,0]*data.x[ix2,0])*tocm
    x2=math.sqrt(data.x[ix2,0]*data.x[ix3,0])*tocm      # does not depend on iz    
    for iz in range(data.nz):    
      iz1 = np.max([0,iz-1])
      iz2 = iz
      iz3 = np.min([data.nz-1,iz+1])
      tanbeta1=0.5*(data.z[ix,iz1]+data.z[ix,iz2])/data.x[ix,0]  
      tanbeta2=0.5*(data.z[ix,iz2]+data.z[ix,iz3])/data.x[ix,0]  
      data.vol[ix,iz]=fourthreepi * (x2**3-x1**3) * (tanbeta2-tanbeta1)


def calc_columnd(data):
  '''
  Calculated the vertical and radial column number densities for every species 
  at every point in the disk (from top to bottom). Very simple and rough method.
  
  Only one half of the disk is considered. If one needs the total surface density
  simply multiply the value from the midplane (zidx=0) by two.
  
  TODO: move this to utils
  '''    
  data.cdnmol = 0.0 * data.nmol
  for ix in range(data.nx):
    for iz in range(data.nz - 2, -1, -1):  # from top to bottom        
      dz = (data.z[ix, iz + 1] - data.z[ix, iz])
      dz = dz * u.au.to(u.cm)
      nn = 0.5 * (data.nmol[ix, iz + 1, :] + data.nmol[ix, iz, :])
      data.cdnmol[ix, iz, :] = data.cdnmol[ix, iz + 1, :] + nn * dz

  # check if integration is correct 
  # for the total hydrogen column density the error is less than 1.e-3
  # that should be good enough for plotting
  #nHverC=data.cdnmol[:,:,data.spnames["H"]]+data.cdnmol[:,:,data.spnames["H+"]]+data.cdnmol[:,:,data.spnames["H2"]]*2.0
  #print(numpy.max(numpy.abs(1.0-nHverC[:,0]/data.NHver[:,0])))    

  data.rcdnmol = 0.0 * data.nmol 
  for iz in range(data.nz):
    for ix in range(1,data.nx,1):  # first ix point (ix=0= remains zero  
      r1= (data.x[ix, iz]**2+data.z[ix, iz]**2)**0.5
      r2 = (data.x[ix-1, iz]**2+data.z[ix-1, iz]**2)**0.5
      dr = r1-r2
      dr = dr * u.au.to(u.cm)
      nn = 0.5 * (data.nmol[ix, iz , :] + data.nmol[ix-1, iz, :])
      data.rcdnmol[ix, iz, :] = data.rcdnmol[ix-1, iz, :] + nn * dr
#   # FIXME: test the integration error can be at most 16% ... good enough for now (most fields are better)
#   nHverC=data.rcdnmol[:,:,data.spnames["H"]]+data.rcdnmol[:,:,data.spnames["H+"]]+data.rcdnmol[:,:,data.spnames["H2"]]*2.0
#   izt=data.nz-2
#   print(nHverC[:,izt],data.NHrad[:,izt])
#   print(numpy.max(numpy.abs(1.0-nHverC[1:,:]/data.NHrad[1:,:])))

def analyse_chemistry(species,model,name=None,directory=".",filenameReactions="Reactions.out",
                      filenameChemistry='chemanalysis.out',to_txt=True,td_fileIdx=None):
    """ 
    Function that analyses the chemistry in a similar way to chemanalise.pro     
    """  
    tarfile=None

    chemistry=Chemistry(None)
    
    if name == None:
        if directory==None or directory=="." or directory=="":
            dirfields=os.getcwd().split("/")
        else:  
            dirfields = directory.split("/")

    name = dirfields[-1]  

    f,dummy = _getfile(filenameReactions, directory, tarfile) #change this

    if f is None: return None
    lines = f.readlines()
    f.close()

    fchem,dummy = _getfile(filenameChemistry, directory, tarfile)
    if fchem is None: return None
    lineschem = fchem.readlines()
    fchem.close()

    chem=np.zeros((len(lineschem)-1,5))
    
    k=0
    for i in lineschem[1:]:        
        chem[k,0]=int(i[0:4])
        chem[k,1]=int(i[4:9])
        chem[k,2]=int(i[9:14])
        chem[k,3]=int(i[14:19])
        chem[k,4]=float(i[19:30])
        k+=1

    chem=np.loadtxt(filenameChemistry,skiprows=1)


    n_rel_index = model.spnames[species]
    spchem=chem[chem[:,2]==n_rel_index]
    lreac=np.unique(spchem[:,3])

    form=[]
    dest=[]
    for i in lines:
        filt=[i[:5]=='     '[:-len(num)]+num for num in lreac.astype(int).astype(str)]
        if np.sum(filt)==1:
            if species+' ' in i[18:43]:
                dest+=[i[:82]] #82
            if species+' ' in i[47:82]:
                form+=[i[:82]]
    destix=[int(i[:5]) for i in dest]
    formix=[int(i[:5]) for i in form]

    fmrates=[]
    for i in formix:
        filt=spchem[:,3]==i
        try:
            fmrates+=[(spchem[:,4][filt]).max()]
        except:
            fmrates+=[0]

    dtrates=[]
    for i in destix:
        filt=spchem[:,3]==i
        try:
            dtrates+=[(spchem[:,4][filt]).max()]
        except:
            dtrates+=[0]     

    mainffilt=spchem[:,3]==formix[np.argmax(fmrates)]
    farix=spchem[:,0:2][mainffilt]-1
    frates=spchem[:,4][mainffilt]
    farray=np.ones((model.nx,model.nz))*frates.min()
    for i,j in zip(farix.astype(int),frates):
        farray[i[0],i[1]]=j

    maindfilt=spchem[:,3]==destix[np.argmax(dtrates)]
    darix=spchem[:,0:2][maindfilt]-1
    drates=spchem[:,4][maindfilt]
    darray=np.ones((model.nx,model.nz))*drates.min()
    for i,j in zip(darix.astype(int),drates):
        darray[i[0],i[1]]=j

    chemistry.farray=farray
    chemistry.darray=darray 

    sordix=np.argsort(dtrates)
    sordli=np.array(destix)[sordix][::-1]
    sorfix=np.argsort(fmrates)
    sorfli=np.array(formix)[sorfix][::-1]
    frg=(np.arange(len(sorfli))+1).astype(str)
    drg=(np.arange(len(sordli))+1).astype(str)

    chemistry.sorted_form_info=['     '[:-len(num)]+num+i for num,i in zip(frg,np.array(form)[np.argsort(fmrates)][::-1])]
    chemistry.sorted_dest_info=['     '[:-len(num)]+num+i for num,i in zip(drg,np.array(dest)[np.argsort(dtrates)][::-1])]

    mfrli=[]
    mdtli=[]
    for i in range(1,model.nx+1):
        for j in range(1,model.nz+1)[::-1]:
            chfil=(spchem[:,0]==i)&(spchem[:,1]==j)
            ffil=chfil&np.in1d(spchem[:,3],sorfli)
            dfil=chfil&np.in1d(spchem[:,3],sordli)
            freacs=spchem[:,3:5][ffil]
            try:
                fmax=freacs[:,0][np.argmax(freacs[:,1])]
                ifmax=frg[sorfli==fmax][0]
                mfrli+=[ifmax]
            except:
                v+=[0]
            dreacs=spchem[:,3:5][dfil]
            try:
                dmax=dreacs[:,0][np.argmax(dreacs[:,1])]
                idmax=drg[sordli==dmax][0]
                mdtli+=[idmax]
            except:
                mdtli+=[0]
    chemistry.main_form_grid=mfrli
    chemistry.main_dest_grid=mdtli
    
    
    if to_txt:
        output_chem_fname='chemistry_reactions_'+species+'.txt'
        f = open(output_chem_fname, 'w')
        f.writelines("-------------------------------------------------------\n")
        f.writelines("Formation and destruction reactions in descending order\n")
        f.writelines("according to their maximum rate across all grid points\n")
        f.writelines("species :"+species+"\n\n")
        f.writelines("Formation reactions\n")
        for i in chemistry.sorted_form_info:
            f.writelines(i+'\n')
        f.writelines("\n\n")
        f.writelines("Destruction reactions\n")
        for i in chemistry.sorted_dest_info:
            f.writelines(i+'\n')   
        f.writelines("-------------------------------------------------------\n")
        f.close()
        print("Writing results to: "+output_chem_fname)
        

    
    return chemistry




def _flux_Wm2toJykms(flux,frequency):
  '''
  Converts a flux from W m^-2 to Jy km/s
  
  Parameters
  ----------
  flux: float
    the flux in units of [W m^-2]
  
  frequency: float
    the frequency of the flux in [GHz]

  Returns
  -------
  float
    The flux of the line. `UNIT: Jansky km s^-1`
  '''
  res=flux*u.Watt/(u.m**2.0)
  ckm=const.c.to('km/s')
     
  res=(res).to(u.Jansky,equivalencies=u.spectral_density(frequency*u.GHz))
    
  return (res*ckm).value


def _getfile(filename,directory=None,tarfile=None,binary=False):
  '''
  Utility function to open a particular file from a ProDiMo model.
  '''
  pfilename = filename

  if tarfile is not None:
    if directory is not None and directory is not ".":
      pfilename = directory + "/" + filename
    
    tararchive = tar.open(tarfile,"r")
    if binary:
      f=tararchive.extractfile(pfilename)
    else:
      f=io.TextIOWrapper(tararchive.extractfile(pfilename))
  else:  
    if directory is not None:
      pfilename = directory + "/" + filename
    
    try:
      if binary:
        f = open(pfilename, 'rb')
      else:
        f = open(pfilename, 'r')
    except: 
      f=None

  if f is None:
    print(("WARN: Could not open " + pfilename + "!"))        
  else:
    print("READ: Reading File: ", pfilename, " ...")
  return f,pfilename

###############################################################################
# For testing
if __name__ == "__main__":
    # import sys
    # import phxpy.io
#     import time
#           
#     pd=read_prodimo("/home/rab/MODELS/XRTPaperNew/TEST_full")      
#     print pd  
#     print pd.nmol[pd.nx-1,0,pd.spnames["N2#"]]   
#     print pd.gas.energy/1000.0 
#     
#     start=time.time()
#     read_lineEstimates("/home/rab/MODELS/XRTPaperNew/TEST_full", pd)
#     print "Time: ",time.time()-start
#     
#     line=pd.getLineEstimate("N2H+", 1000.0)
#     line=pd.getLineEstimate("N2H+", 1000.0)
#     print line
    
    # lines=pd.selectLineEstimates("N2H+")
    # print len(lines)
            
                
  tdir = "../testdata"
  
  data = read_prodimo(tdir)  
  
                      
  linesObs = read_lineObs(tdir, len(data.lines))
  print(data.lines[0])
  print(data.lines[1])
  print(linesObs[0])
   
  profile = read_lineObsProfile(tdir + "/LineProfile_CO_21.dat")
  print(profile)
  print(profile.flux)
  print(profile.velo)
    
        
