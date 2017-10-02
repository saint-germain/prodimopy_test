"""
.. module:: read_casasim 

.. moduleauthor:: Ch. Rab


"""
from __future__ import print_function
from __future__ import division 
from __future__ import unicode_literals

import numpy

from astropy import units as u
import astropy.wcs as wcs
import astropy.io.fits as fits
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.wcs.wcs import WCS


class CasaSim():
  """
  Kind of container class for a certain set of CASA simulations (observations).
  
  The class can be used to quickly read in a CASA simulation for a single 
  spectral line. The main purpose is to read int the various different kind 
  of data which have been produced by CASA. 
    
  Can deal with the following files wiht a given PREFIX
    - `PREFIX.cube.fits` the spectral line cube
    - `PREFIX.cube.integrated.fits` the integrated line emission (zeroth moment)
    - `PREFIX.cube.integrated.radial` the azimuthally averaged radial profile (text file produce by casairing task)
    - `PREFIX.cube.specprof` the spectral line profile (text file produced by CASA)    
    - `PREFIX.cube.mom1.fits` the first moment of the line emission
    - `PREFIX.cube.pv.fits` the position velocity diagram
    - `PREFIX.cont.fits` the corrresponding continuum
    - `PREFIX.cont.radial` the azimuthally average radial continuum profile 
  
  If the class is instantiated all this files are read. However, if they do not 
  exist only a warning is printed and the according instance attributes are set 
  to `None`.
       
  """
  def __init__(self,fn_prefix,directory=".",systemic_velocity=0.0,distance=None,
               fn_cube=".cube.fits",
               fn_integrated=".cube.integrated.fits",
               fn_radprof=".cube.integrated.radial",
               fn_specprof=".cube.specprof",
               fn_mom1=".cube.mom1.fits",
               fn_pv=".cube.pv.fits",
               fn_cont=".cont.fits",
               fn_cont_radprof=".cont.radial"):
    """
    Parameters
    ----------
    
    fn_prefix : str
      the prefix used for all files the routine try to read
      
    directory : 
      the directory with the files (optional)
    
    systemic_velocity : float
      the systemic_velocity of the target (optional)
      
    distance : float
      the distance of the target (optional).
      Try to fread it from the fits file in case of synthetic ProDiMo observations

    fn* : string
      The fn* parameters can be use to change the filenames for the different files. 
    
    
    Attributes
    ----------
    systemic_velocity : float
      the systemic velocity of the target, will be passed to the other objects
      
    systemic_velocity : float
      the distance of the target, will be passed to the other objects      
    
    cube :class:`.LineCube` :
      the line cube.

    integrated : :class:`.LineIntegrated`
      the zeroth moment of the cube (integrated emission)

    radprof : :class:`.LineSpectralProfile`
      the spectra profile
      
    specprof : :class:`.RadialProfile`
      the spectra profile
      
    mom1 : :class:`.LineMom1`
      the first moment of the cube
      
    pv : :class:`.LinePV`
      the position-velocity diagram
      
    cont : :class:`.Continuum` 
      the associated continuum
      
    cont_rad : :class:`.RadialProfile` 
      the radial profile of the continuum
    
    """
    self.systemic_velocity=systemic_velocity
    self.distance=distance
    self.fn_cube=self._build_fn(directory, fn_prefix, fn_cube)
    self.fn_integrated=self._build_fn(directory, fn_prefix, fn_integrated)
    self.fn_specprof=self._build_fn(directory, fn_prefix, fn_specprof)
    self.fn_radprof=self._build_fn(directory, fn_prefix, fn_radprof)
    self.fn_mom1=self._build_fn(directory, fn_prefix, fn_mom1)
    self.fn_pv=self._build_fn(directory, fn_prefix, fn_pv)
    self.fn_cont=self._build_fn(directory, fn_prefix, fn_cont)
    self.fn_cont_radprof=self._build_fn(directory, fn_prefix, fn_cont_radprof)
    
    self.cube=LineCube(self.fn_cube,systemic_velocity=systemic_velocity)
    self.integrated=LineIntegrated(self.fn_integrated)
    
    self.radprof=RadialProfile(self.fn_radprof,bwidth=self._beam_width_arcsec(self.integrated))
    self.specprof=LineSpectralProfile(self.fn_specprof,systemic_velocity=systemic_velocity)
    self.mom1=LineMom1(self.fn_mom1,systemic_velocity=systemic_velocity)
    self.pv=LinePV(self.fn_pv,systemic_velocity=systemic_velocity)    
    self.cont=Continuum(self.fn_cont)
    self.cont_radprof=RadialProfile(self.fn_cont_radprof,bwidth=self._beam_width_arcsec(self.cont))

    return

  def _build_fn(self, directory, prefix, filename):
    return directory + "/" + prefix + filename
  
  def _beam_width_arcsec(self,image):
    if image is not None: 
        bwidth=(image.header["BMIN"]+image.header["BMAJ"])/2.0
        #FIXME: assume that it is deg to arcsec
        return bwidth*3600.0
    return None


class CASAImage(object):
  """
  Super class for an observable (CASA Simulation or real Observation) for image data including cubes.
  
  Should not be instantiated directly.
  Opens the fits file and initializes some helper attributes.
  
  """
  def __init__(self,filename=None):
    """
    Needs to be called from a subclass.
    
    Parameters
    ----------
    filename : str
      The file name of the fits file.
      
    Attributes
    ----------
    header : 
      the header of the fits file. It is just a reference to the fits header.
      
    data : array_like(float)
      the fits data. Assumes that there is only one HDU.
      
    bmaj : float
      the major axis of the beam (same units as in the fits file)

    bmaj : float
      the minor axis of the beam (same units as in the fits file)

    bPA : float
      the beam position angle
      
    centerpix : array_like(float,ndim=1)
      the pixel coordinates of the center (x and y coordinate)
    
    wcsabs : :class:`astropy.wcs.WCS` 
      an absulte astropy world coordinate system (for the spatial coordinates). 
      Can be used for plotting transformation.

    wcsrel : :class:`astropy.wcs.WCS` 
      astropy world coordinate system (for the spatial coordinates) relative the 
      the center of the image (`centerpix`) 
      Can be used for plotting transformation.
      (see :func:`~prodimopy.read_casasim.CASAImage._linear_offset_coords`)
      
    """
    if filename is not None:
      try:
        fitsdata= fits.open(filename)
        self.header=fitsdata[0].header
        self.data=fitsdata[0].data
      except FileNotFoundError:
          return None
    
    self.bmaj=self.header["BMAJ"]
    self.bmin=self.header["BMIN"]
    self.bPA=self.header["BPA"]
    # this is all a bit strange but if it is even pixels now -1 is better
    self.centerpix=[self.header["CRPIX1"],self.header["CRPIX2"]]
    self.centerpix[0]-=1
    self.centerpix[1]-=1
            
    self.wcsabs=wcs.WCS(self.header,naxis=(1,2))
    self.wcsrel=self._linear_offset_coords(self.wcsabs, self.centerpix)


  def _linear_offset_coords(self,wcsabs, center):
    """
    Returns a locally linear offset coordinate system.
    
    Given a 2-d celestial WCS object and a central coordinate, return a WCS
    that describes an 'offset' coordinate system, assuming that the
    coordinates are locally linear (that is, the grid lines of this offset
    coordinate system are always aligned with the pixel coordinates, and
    distortions from spherical projections and distortion terms are not taken
    into account)
    
    taken from: https://github.com/aplpy/aplpy/issues/8
    
    Parameters
    ----------
    wcs : `~astropy.wcs.WCS`
        The original WCS, which should be a 2-d celestial WCS
    center : array_like
        The pixel coordinates on which the offset coordinate system should be
        centered.
    """

    # Convert center to pixel coordinates
    #xp, yp = skycoord_to_pixel(center, wcs)
        
    # Set up new WCS
    new_wcs = wcs.WCS(naxis=2)
    #FIXME: do not know why +1 but otherwise the 0 points does not fit with the zero point in the pixel scale
    new_wcs.wcs.crpix = center[0]+1, center[1]+1
    new_wcs.wcs.crval = 0., 0.
    # FIXME: assumbes that the units are in degree
    scale=proj_plane_pixel_scales(wcsabs)*3600.
    scale[0]=scale[0]*-1.0 # otherwise it is the wrong direction, that does not change the orientation of the image just the labels
    new_wcs.wcs.cdelt = scale
    new_wcs.wcs.ctype = 'XOFFSET', 'YOFFSET'
    new_wcs.wcs.cunit = 'arcsec', 'arcsec'

    return new_wcs
        
  def __str__(self, *args, **kwargs):
    return("BEAM(maj,min,PA): "+str(self.bmaj)+","+str(self.bmin)+","+str(self.bPA)+"\n"
           "CENTERPIX: "+str(self.centerpix))


class LineCube(CASAImage):
  """  
  The full line cube.
  """
  def __init__(self,filename,systemic_velocity=0.0):
    self.systemic_velocity=systemic_velocity
    CASAImage.__init__(self,filename)
    self.vel=self._init_vel()
    
  def _init_vel(self):
    """
    Converts the spectral coordinate to (radio) velocities.
     
    Currently assumes that the spectral coordinate is in units of Hz.

    FIXME: read unit from fits file 

    FIXME: stolen from here ... https://github.com/astropy/astropy/issues/4529
    mabe not the best way to go
    
    FIXME: maybe could be generalized
         
    """
    wcsvel=wcs.WCS(self.header)
    xv = numpy.repeat(0, self.header['NAXIS3'])
    yv = numpy.repeat(0, self.header['NAXIS3'])
    zv = numpy.arange(self.header['NAXIS3']) 
     
    wx, wy, wz = wcsvel.wcs_pix2world(xv, yv, zv, 0)
    freq_to_vel = u.doppler_radio(self.header['RESTFRQ']*u.Hz)
    vel=(wz * u.Hz).to(u.km / u.s, equivalencies=freq_to_vel)
     
    return vel  


class LineIntegrated(CASAImage):
  """
  Zeroth moment image (Integrated emission)
  """
  def __init__(self,filename):
    CASAImage.__init__(self,filename)
  
class LineMom1(CASAImage):
  """
  First Moment image.
  """
  def __init__(self,filename,systemic_velocity=0.0):
    
    CASAImage.__init__(self,filename)
    
    self.systemic_velocity=systemic_velocity
    

class LinePV(CASAImage):
  """
  Position-Velocity diagramm.
  
  """
  def __init__(self,filename,systemic_velocity=0.0):
    
    CASAImage.__init__(self,filename)
    
    self.systemic_velocity=systemic_velocity
    
    # FIXME: CASA sets the centrpix for the velocity coordinate to a strange value
    # so if it is negative calculate from the dimension
    if self.centerpix[1] < 0:
      # FIXME: work only with odd numbers
      self.centerpix[1]=int(self.header["NAXIS2"]/2)      


class LineSpectralProfile():
  """
  The spectral line profile.
  
  Reads a spectral profile produced with CASA.
  """
  def __init__(self,filename,systemic_velocity=0.0):
    data=numpy.loadtxt(filename)
    self.vel=data[:,3]
    self.flux=data[:,4]
    self.systemic_velocity=systemic_velocity
    return
  
 
class RadialProfile():
  """
  Azimuthally averaged radial intensity profile.
  
  Reads a file produced by the casa taks casairing.
  see https://www.oso.nordic-alma.se/software-tools.php
  
  """
  def __init__(self,filename,bwidth=None):
    radprof=numpy.loadtxt(filename)
    self.arcsec=radprof[:,0]
    self.flux=radprof[:,2]
    self.flux_err=radprof[:,3]
    self.bwidth=bwidth
    return
  
class Continuum(CASAImage):
  """
  A continuum image.
  
  """  
  def __init__(self,filename): 
    CASAImage.__init__(self,filename)
