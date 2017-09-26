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


class CasaSim():
  """
  Can deal with the following files
  .cube.fits,cont.fits, cube.integrated.fits,cube.specprof
  
  TODO: provide some general parameters which can than also used for plotting e.g.
        - systemic velocity
        - supposed centeral position of target
        - velocity grid (absolute and relative) ... only radio velocity
        - threesigma value ... 
        - distance (can also be taken from fits files in case of simulations)
        - beam sizes etc.
  
  
  """
  def __init__(self,fn_prefix,directory=".",
               fn_cube=".cube.fits",
               fn_integrated=".cube.integrated.fits",
               fn_radprof=".cube.integrated.radial",
               fn_specprof=".cube.specprof",
               fn_mom1=".cube.mom1.fits",
               fn_pv=".cube.pv.fits",
               fn_cont=".cont.fits",
               fn_cont_radprof=".cont.radial"):
    
    self.fn_cube=directory+"/"+fn_prefix+fn_cube
    self.fn_integrated=directory+"/"+fn_prefix+fn_integrated
    self.fn_specprof=directory+"/"+fn_prefix+fn_specprof
    self.fn_radprof=directory+"/"+fn_prefix+fn_radprof
    self.fn_mom1=directory+"/"+fn_prefix+fn_mom1
    self.fn_pv=directory+"/"+fn_prefix+fn_pv
    self.fn_cont=directory+"/"+fn_prefix+fn_cont
    self.fn_cont_radprof=directory+"/"+fn_prefix+fn_cont_radprof
    
    self.cube=LineCube(self.fn_cube)
    self.integrated=LineIntegrated(self.fn_integrated)
    self.specprof=LineSpectralProfile(self.fn_specprof)
    self.mom1=LineMom1(self.fn_mom1)
    self.pv=LineMom1(self.fn_pv)
    self.radprof=RadialProfile(self.fn_radprof)
    self.cont=Continuum(self.fn_cont)
    self.cont_radprof=RadialProfile(self.fn_cont_radprof)
        
    self.bmaj=None
    self.bmin=None
    
    if self.cube is not None:
      self.bmaj=self.cube.header["BMAJ"]
      self.bmin=self.cube.header["BMIN"]
      self.bPA=self.cube.header["BPA"]
      print(self.bmaj,self.bmin,self.bPA,self.bmaj*self.bmin)
    else:
      print("Did not find a cube file with name"+self.fn_cube)

    
    return

class LineCube():
  """ Holds the data for a line cube produced with CASA. 
     
  """
  def __init__(self,filename):
    try:
      fitsdata= fits.open(filename)
      self.header=fitsdata[0].header
      self.data=fitsdata[0].data
    except FileNotFoundError:
      return None     
  
class LineIntegrated():
  def __init__(self,filename):
    try:
      fitsdata= fits.open(filename)
      self.header=fitsdata[0].header
      self.data=fitsdata[0].data
    except FileNotFoundError:
      return None
  
class LineMom1():
  def __init__(self,filename):
    try:
      fitsdata= fits.open(filename)
      self.header=fitsdata[0].header
      self.data=fitsdata[0].data
    except FileNotFoundError:
      return None

class LinePV():
  def __init__(self,filename):
    try:
      fitsdata= fits.open(filename)
      self.header=fitsdata[0].header
      self.data=fitsdata[0].data
    except FileNotFoundError:
      return None    

class LineSpectralProfile():
  def __init__(self,filename):
    data=numpy.loadtxt(filename)
    self.vel=data[:,3]
    self.flux=data[:,4]
    return
  
class Continuum():
  def __init__(self,filename):
    try:
      fitsdata= fits.open(filename)
      self.header=fitsdata[0].header
      self.data=fitsdata[0].data
    except FileNotFoundError:
      return None    
  
class RadialProfile():
  def __init__(self,filename):
    radprof=numpy.loadtxt(filename)
    self.arcsec=radprof[:,0]
    self.flux=radprof[:,2]
    self.flux_err=radprof[:,3]
    return
  

# class FITS_Image():
#   """ Data container for observational data (synthetic and real) in fits format.
#   
#   Notes
#   -----
#   This class is a simple container for observational data stored in fits files but also 
#   provides some additional utility functions for convenience.
#   
#   It was made mainly to read in ALAM simulation produced with casa. It can also read in several 
#   other output which can be produced with CASA.
#   these are
#   
#   - position velocity diagrams
#   - moment0 
#   
#   I can also directly read fits cubes produced by prodimo (e.g. the model output).       
#   """ 
#    
#   def __init__(self, ff,threesig,name="model",radialfile=None):
#         
#     self.name=name    
#     self.fits=ff    
#     self.header=ff[0].header
#     
#     #self.scidata = np.array(ff[0].data[:,:]*1000.0) # Convert to milli Jansky
#     self.scidata = np.array(ff[0].data[:,:]) # Convert to milli Jansky  
#     
#     if self.header['NAXIS']==3:    
#       self.scidata=self.scidata[0,:,:]
#     else:
#       self.scidata=self.scidata[:,:]  
#     # add a column so that it is squared
#     #self.scidata = np.c_[self.scidata, self.scidata[:, 0]] 
#     # make it square
#     # TODO: check if this is required
# #     nx=self.scidata.shape[0]
# #     ny=self.scidata.shape[1]
# #     
# #     print nx,ny
# #     if nx < ny: 
# #       self.scidata = np.c_[self.scidata, self.scidata[:, 0]]
# #     elif ny <nx:
# #       self.scidata = np.r_[self.scidata, self.scidata[0, :]]
#                                    
#     #print(self.scidata.shape)
#     self.object = self.header['OBJECT']    
#     #self.frequency = self.header['RESTFRQ']  # in herzt        
#     #self.pixsize = header['SCUPIXSZ']
#     # central pixel -1 as python starts to count at zero
#     #self._centerpixorig = np.array([int(header['CRPIX1']), int(header['CRPIX2'])]) - 1.0
#     self.cpix_x = int(self.header['CRPIX1'])-1
#     self.cpix_y = int(self.header['CRPIX2'])-1
#     # FIXME: make this an optional field
#     #self.date = self.header['DATE-OBS']
#     #self.ra = self.header['OBSRA']
#     #self.dec = self.header['OBSDEC']
#     #self.long = header['LONG']
#     #self.lat = header['LONG']
#     self.threesig = threesig
#     self.bmaj=(float(self.header["BMAJ"])*u.deg).to(u.arcsec).value
#     self.bmin=(float(self.header["BMIN"])*u.deg).to(u.arcsec).value
#     self.bpa=float(self.header["BPA"])
#     
#     
#     self.radprof=None
#     if radialfile is not None:
#       # FIXME: move this to separate routine
#       # FIXME: define a class for the radial profile
#       self.radprof=np.loadtxt(radialfile,usecols = (0,2,3),skiprows=1)
#       # convert to milijansky
# #      self.radprof[:,1]=self.radprof[:,1]*1000.0
# #      self.radprof[:,2]=self.radprof[:,2]*1000.0
#       self.radprof[:,1]=self.radprof[:,1]
#       self.radprof[:,2]=self.radprof[:,2]
# 
#       #print(self.radprof)
#     
#     #self.beamHPBW = beamHPBW
#     # beam solid angle for gaussian beam in sterrad
#     # 1.135 * beam HPWM **2
#     # http://docs.jach.hawaii.edu/star/sc11.htx/node25.html   
#     #
#     #self.beamSA= (1.135 * ((self.beamHPBW * u.arcsec).to(u.rad)) ** 2.0).value
#     #self.pixSA = (((self.pixsize * u.arcsec).to(u.rad)) ** 2).value
#     self.__init_wcs()    
# 
#   def __init_wcs(self):
#     """
#     Convert pixel to arcsec   
#     sets coordx (relative position to the center)    
#     """
#     # FIXME: The order of the shape parameters seems to change
#     px = np.linspace(1.0, self.scidata.shape[1], self.scidata.shape[1])
#     py = np.linspace(1.0, self.scidata.shape[0], self.scidata.shape[0])  
#   
#     self.wcs = wcs.WCS(self.header,naxis=["latitude","longitude"])
#         
# 
#     # Print out the "name" of the WCS, as defined in the FITS header
#     #print(self.wcs.wcs.name)
# 
#     # Print out all of the settings that were parsed from the header
#     #self.wcs.wcs.print_contents()
# 
#     # Some pixel coordinates of interest.
#     #pixcrd = np.array([[0, 0], [250, 250], [45, 98]], np.float_)
# 
#     # Convert pixel coordinates to world coordinates
#     # The second argument is "origin" -- in this case we're declaring we
#     # have 1-based (Fortran-like) coordinates.
#     world = self.wcs.wcs_pix2world(px,py, 1)
#     #print(world)
#     #print((world[0]*u.deg).to(u.arcsec))
#     #print((world[1]*u.deg).to(u.arcsec))    
#     
# 
#     # Convert the same coordinates back to pixel coordinates.
#     #pixcrd2 = self.wcs.wcs_world2pix(world, 1)
#     #print(pixcrd2)
# 
#     # These should be the same as the original pixel coordinates, modulo
#     # some floating-point error.
#     #assert np.max(np.abs(pixcrd - pixcrd2)) < 1e-6
#     
# #     coord = self.wcs.wcs_pix2world(px, py, 1) 
#     #print(world[0][self.cpix_x])
#     self.coordx=(world[0]*u.deg).to(u.arcsec)
#     self.coordy=(world[1]*u.deg).to(u.arcsec)
#     self.coordx=(self.coordx[self.cpix_x]-self.coordx).value
#     self.coordy=(self.coordy[self.cpix_y]-self.coordy).value
#     #print(self.coordx)    
#     #print(self.coordy)
#     # use the supixsize to convert to arcsec ... this is mainly for plotting
#      
# #     px = self.convert_pix_to_arcsec(px)
# #     py = self.convert_pix_to_arcsec(py)     
# #           
# #     px = px - px[self.cpix_x]
# #     py = py - py[self.cpix_y]     
# #          
# #     self.coordx = px
# #     self.coordy = py        
# 
#   def convert_pix_to_arcsec(self, pixels):
#     '''
#     Conversion of the pixelsize to arcsec
#     '''
#     return (pixels) * self.pixsize
#     
#   
#   def center_pixcen(self):  
#     cx=int(self.scidata.shape[1]/2)+1
#     cy=int(self.scidata.shape[1]/2)+1
#     self.cpix_x=cx
#     self.cpix_y=cy
#     
#     
#   def center_max(self):
#     """
#     Center the image to the max value of the data
#     overwrites centerpix, and coord x 
#     very simple search for the maximum value and that's it
#     """
#     maxval = np.nanmax(self.scidata)
#     centerpix = np.where(self.scidata == maxval)
#     # this order is correct!    
#     self.cpix_x = centerpix[1][0]
#     self.cpix_y = centerpix[0][0]
#     
#     print("New center: ", self.cpix_x, self.cpix_y)
# 
#     # reset also the coordinates
#     self.__init_wcs()  