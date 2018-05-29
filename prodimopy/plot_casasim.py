"""
.. module:: plot_casasim 

.. moduleauthor:: Ch. Rab

"""
from __future__ import division 
from __future__ import print_function
from __future__ import unicode_literals

from astropy import wcs
from matplotlib import patches
from matplotlib.offsetbox import AnchoredOffsetbox, AuxTransformBox
from matplotlib.ticker import MaxNLocator
import numpy
import prodimopy.plot as pplot

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np


class PlotCasasim(object):
  '''
  Plot routines for casa simulation results. 
  
  This class can be used together with :mod:`prodimoy.read_casasim`
  '''
  def __init__(self, pdf):
    """
    
    Parameters
    ----------
    pdf : class:`matplotlib.backends.backend_pdf.PdfPages` 
      this object is used to save the plots in a pdf file.
    
    """    
    self.pdf = pdf
    
  def _dokwargs(self,ax,**kwargs):
    '''
    Handles the passed kwargs elements (assumes that defaults are already set)
    TODO: make this a general function .... 
    TODO: does not work with subplots
    '''
#     if "ylim" in kwargs: 
#       ax.set_ylim(kwargs["ylim"])
#       
    if "xlim" in kwargs: 
      print(kwargs["xlim"])
      ax.set_xlim(kwargs["xlim"])

    if "ylim" in kwargs: 
      print(kwargs["ylim"])
      ax.set_ylim(kwargs["ylim"])

#               
#     if "xlog" in kwargs:
#       if kwargs["xlog"]: 
#         ax.semilogx()
#       else:
#         ax.set_xscale("linear")
#       
#     if "ylog" in kwargs:
#       if kwargs["ylog"]: 
#         ax.semilogy()
#       else:              
#         ax.set_yscale("linear")
#       
#     if "xlabel" in kwargs:
#       ax.set_xlabel(kwargs["xlabel"])  
# 
#     if "ylabel" in kwargs:
#       ax.set_ylabel(kwargs["ylabel"])
#     
#     if self.title != None:
#       if self.title.strip() != "":
#         ax.set_title(self.title.strip())
                
    if "title" in kwargs:
      if  kwargs["title"] != None and kwargs["title"].strip() != "":
        ax.set_title(kwargs["title"].strip())
      else:
        ax.set_title("")

  def add_beam(self, ax, image):
    """
    adds a beam to a plot (axis).
    
    Parameters
    ----------
    ax : class:`matplotlib.axes.Axes`
      the axis for which the beam should be added
      
    image : class:`~prodimopy.read_casasim.CASAImage`
      some kind of image (can also be a cube)
      
    
    FIXME: use the proper units from the fits file
    
    FIXME: is proberly not general yet
    
    """
    bmin = image.bmin / image.header["CDELT1"]
    bmaj = image.bmaj / image.header["CDELT1"]
    bpa = image.bPA
    
    ae = AnchoredEllipse(ax.transData, width=bmin, height=bmaj, angle=bpa,
                           loc=3, pad=0.5, borderpad=0.4, frameon=False)    
    ax.add_artist(ae)
  
  
  def plot_cube(self, cube, nrow=3, ncol=3, cvel_idx=None,step=1, 
                zlim=[None, None],rms=False,**kwargs):
    """
    Plots a spectral line cube.
    """
    
    if cube is None: return
    
    nimages = cube.data.shape[0]
    
    # wcvel=wcs.WCS(image[0].header)
    # print(wcs.find_all_wcs(image[0].header))  
    
    
    if cvel_idx is None:
      cvel_idx = int(nimages / 2)
    
    naxesh = int(nrow * ncol / 2)
    
    vmin = zlim[0]
    vmax = zlim[1]
    if vmin is None: vmin = numpy.min(cube.data)
    if vmax is None: vmax = numpy.max(cube.data)  
  
  
    # cube.header['CRVAL1'] = 0.0
    # cube.header['CRVAL2'] = 0.0
    # wcsim=wcs.WCS(cube.header,naxis=(1,2))  
    # cpix=[cube.header["CRPIX1"],cube.header["CRPIX2"]]
    # wcsrel=linear_offset_coords(wcsim, cube.centerpix)
  
    fig, axis = plt.subplots(nrow, ncol, sharex=False, sharey=False, subplot_kw=dict(projection=cube.wcsrel))
    fig.subplots_adjust(hspace=-0.1, wspace=0.0)
    # fig,axis= plt.subplots(nrow, ncol,sharex=False,sharey=False)
    # fig=plt.figure()
    # assume the widht of the image is given, scale it with the number of axis
    figsize = fig.get_size_inches()
    figsize[0] = figsize[0] * 2.0  # assumes that default size is one column in a Paper
    figsize[1] = figsize[0] / (ncol) * nrow
    fig.set_size_inches(figsize)
   
    for ir in range(nrow):
      for ic in range(ncol):
        ax = axis[ir, ic]
        iax = ir * ncol + ic         
        velidx = cvel_idx + (iax - naxesh)*step
        im = ax.imshow(cube.data[velidx, :, :], cmap="inferno", vmin=vmin, vmax=vmax,origin="lower")      
              
        # set the border of the coordinate frames to white
        ax.coords[0].frame.set_color("white")
        ax.coords[0].set_ticks(color="white", spacing=1.0 * u.arcsec)
        ax.coords[1].set_ticks(color="white", spacing=1.0 * u.arcsec)
        if not (ir == nrow - 1 and ic == 0):
          ax.coords[0].set_ticklabel_visible(False)
        else:
            ax.set_xlabel("rel. RA ['']")
            
        if ic > 0 or (not ir == nrow - 1):
          ax.coords[1].set_ticklabel_visible(False)
        else:
          ax.set_ylabel("rel. Dec. ['']")  
        
        # print the velocities relative to the systemic velocities
        props = dict(boxstyle='round', facecolor='white', edgecolor="none")
        
        ax.text(0.95, 0.95, "{:5.2f}".format(cube.vel[velidx] - cube.systemic_velocity * (u.km / u.s)), transform=ax.transAxes, fontsize=6,
          verticalalignment='top', horizontalalignment="right", bbox=props)
        
        # mark the center
        ax.plot(cube.centerpix[0], cube.centerpix[1], marker="x", color="0.6", linewidth=0.5, ms=3)
        
        # ax.axis('equal')
        self.add_beam(ax, cube)  
        
        if rms:
        # calculate the rms of the residual    
          rms=numpy.nansum((cube.data[velidx, :, :]**2.0)/cube.data[velidx, :, :].size)**0.5
          ax.text(0.05, 0.95, "rms="+"{:4.1e}".format(rms),
                 transform=ax.transAxes, fontsize=5,
                 verticalalignment='top', horizontalalignment="left", bbox=props)
        # FIXME: that would show the image like in the casaviewer
        # ax.invert_yaxis()  
                                        
    ticks = MaxNLocator(nbins=6).tick_values(vmin, vmax)
    CB = fig.colorbar(im, ax=axis.ravel().tolist(), pad=0.02,
                    format="%5.2f", fraction=0.015)  
    CB.set_ticks(ticks)
    CB.set_label("[Jy/beam]")
      
    # plt.tight_layout(pad=0.1)
    
    self._closefig(fig)

  def plot_cube_cont(self, cube, nrow=3, ncol=3, zlim=[None, None],**kwargs):
    """
    Plots a spectral line cube.
    """
    
    if cube is None: return
    
    nimages = cube.data.shape[0]
    
    # wcvel=wcs.WCS(image[0].header)
    # print(wcs.find_all_wcs(image[0].header))  
    
    
    icenter = int(nimages / 2)
    naxesh = int(nrow * ncol / 2)
    
    vmin = zlim[0]
    vmax = zlim[1]
    if vmin is None: vmin = numpy.min(cube.data)
    if vmax is None: vmax = numpy.max(cube.data)  
  
  
    # cube.header['CRVAL1'] = 0.0
    # cube.header['CRVAL2'] = 0.0
    # wcsim=wcs.WCS(cube.header,naxis=(1,2))  
    # cpix=[cube.header["CRPIX1"],cube.header["CRPIX2"]]
    # wcsrel=linear_offset_coords(wcsim, cube.centerpix)
  
    fig, axis = plt.subplots(nrow, ncol, sharex=False, sharey=False, subplot_kw=dict(projection=cube.wcsrel))
    fig.subplots_adjust(hspace=-0.1, wspace=0.0)
    # fig,axis= plt.subplots(nrow, ncol,sharex=False,sharey=False)
    # fig=plt.figure()
    # assume the widht of the image is given, scale it with the number of axis
    figsize = fig.get_size_inches()
    figsize[0] = figsize[0] * 2.0  # assumes that default size is one column in a Paper
    figsize[1] = figsize[0] / (ncol) * nrow
    fig.set_size_inches(figsize)
   
    for ir in range(nrow):
      for ic in range(ncol):
        ax = axis[ir, ic]
        iax = ir * ncol + ic
        idata = iax - naxesh
        wlidx = icenter + idata           
        im = ax.imshow(cube.data[wlidx, :, :], cmap="inferno", vmin=vmin, vmax=vmax,origin="lower")      
              
        # set the border of the coordinate frames to white
        ax.coords[0].frame.set_color("white")
        ax.coords[0].set_ticks(color="white", spacing=1.0 * u.arcsec)
        ax.coords[1].set_ticks(color="white", spacing=1.0 * u.arcsec)
        if not (ir == nrow - 1 and ic == 0):
          ax.coords[0].set_ticklabel_visible(False)
        else:
            ax.set_xlabel("rel. RA ['']")
            
        if ic > 0 or (not ir == nrow - 1):
          ax.coords[1].set_ticklabel_visible(False)
        else:
          ax.set_ylabel("rel. Dec. ['']")  
        
        # print the wavelenth 
        props = dict(boxstyle='round', facecolor='white', edgecolor="none")
        
#        ax.text(0.95, 0.95, "{:5.2f}".format(cube.wl[wlidx]* (u.micron)), transform=ax.transAxes, fontsize=6,
#          verticalalignment='top', horizontalalignment="right", bbox=props)
        
        # mark the center
        ax.plot(cube.centerpix[0], cube.centerpix[1], marker="x", color="0.6", linewidth=0.5, ms=3)
        
        # ax.axis('equal')
        self.add_beam(ax, cube)  
        
        # FIXME: that would show the image like in the casaviewer
        # ax.invert_yaxis()  
  
                                        
    ticks = MaxNLocator(nbins=6).tick_values(vmin, vmax)
    CB = fig.colorbar(im, ax=axis.ravel().tolist(), pad=0.02,
                    format="%5.2f", fraction=0.015)  
    CB.set_ticks(ticks)
    CB.set_label("[Jy/beam]")
      
    # plt.tight_layout(pad=0.1)
    
    self._closefig(fig)

  def plot_integrated_diff(self, imageObs,imageModel,imageDiff, zlim=[None, None],**kwargs):
    """
    Plots an image and its diff Model-Obs
    """
    
    if imageObs is None or imageModel is None or imageDiff is None: return
    
    vmin = zlim[0]
    vmax = zlim[1]
    if vmin is None: vmin = numpy.nanmin([imageObs.data,imageDiff.data])
    if vmax is None: vmax = numpy.nanmax([imageObs.data,imageDiff.data]) 
    # wcsim=wcs.WCS(image.header)    
    # wcsrel=linear_offset_coords(wcsim, image.centerpix)
    fig, axes = plt.subplots(1, 3, subplot_kw=dict(projection=imageObs.wcsrel),
                           figsize=pplot.scale_figs([2.1,1.0]))
    
    
    im = axes[0].imshow(imageObs.data, cmap="inferno", vmin=vmin, vmax=vmax,origin="lower")
    im2 = axes[1].imshow(imageModel.data, cmap="inferno", vmin=vmin, vmax=vmax,origin="lower")
    im3 = axes[2].imshow(imageDiff.data, cmap="inferno", vmin=vmin, vmax=vmax,origin="lower")
    
    # calculate the rms of the residual    
    rms=numpy.nansum((imageDiff.data**2.0)/imageDiff.data.size)**0.5
    
    # print the velocities relative to the systemic velocities
    props = dict(boxstyle='round', facecolor='white', edgecolor="none")
    axes[0].text(0.95, 0.95, "Observation", transform=axes[0].transAxes, 
                 fontsize=6,verticalalignment='top', horizontalalignment="right", bbox=props)
    axes[1].text(0.95, 0.95, "Model", transform=axes[1].transAxes, 
                 fontsize=6,verticalalignment='top', horizontalalignment="right", bbox=props)
    axes[2].text(0.95, 0.95, "Residual (rms="+"{:4.1e}".format(rms)+")",
                 transform=axes[2].transAxes, fontsize=6,
                 verticalalignment='top', horizontalalignment="right", bbox=props)
    # axes.coords[0].set_major_formatter('hh:mm:ss')     

    for ax in axes:
      ax.coords[1].set_ticks(color="white", spacing=1.0 * u.arcsec)
      ax.coords[0].set_ticks(color="white", spacing=1.0 * u.arcsec)
      ax.coords[0].frame.set_color("white")
      ax.set_xlabel("rel. RA ['']")
      
    # needs to be converted to pixels, that is the data unit
    self.add_beam(axes[0], imageObs)
    
    axes[0].set_ylabel("rel. Dec. ['']")
    # mark the center
    #axes[0].plot(imageObs.centerpix[0], imageObs.centerpix[1], marker="x", color="0.6", linewidth=0.5, ms=3)
  
    #axes.axvline(image.centerpix[0],color="0.6",linewidth=0.8,linestyle=":")  
    #axes.axhline(image.centerpix[1],color="0.6",linewidth=0.8,linestyle=":")
  
      # FIXME: that would show the image like in the casaviewer
    # axes.invert_yaxis()  
  
    self._dokwargs(axes[0],**kwargs)
  
    ticks = MaxNLocator(nbins=6).tick_values(vmin, vmax)
    CB = fig.colorbar(im, ax=axes.ravel().tolist(), pad=0.02,
                    format="%5.2f", fraction=0.04)  
    CB.set_ticks(ticks)
    CB.set_label("[Jy/beam km/s]")
          
    self._closefig(fig)
  
  
  def plot_integrated(self, image, zlim=[None, None],mJy=False,**kwargs):
    """
    Plots a zeroth moment image (integrated intensity) image.
    """
    scalefac=1.0
    if mJy==True:
      scalefac=1000.0
    
    if image is None: return
    
    vmin = zlim[0]
    vmax = zlim[1]
    if vmin is None: vmin = numpy.min(image.data)
    if vmax is None: vmax = numpy.max(image.data)
    vmin=vmin*scalefac
    vmax=vmax*scalefac 
  
    # wcsim=wcs.WCS(image.header)    
    # wcsrel=linear_offset_coords(wcsim, image.centerpix)
    fig, ax = plt.subplots(1, 1, subplot_kw=dict(projection=image.wcsrel))
    
    
    im = ax.imshow(image.data*scalefac, cmap="inferno", vmin=vmin, vmax=vmax,
                   origin="lower")
    # ax.coords[0].set_major_formatter('hh:mm:ss')
    ax.coords[1].set_ticks(color="white", spacing=1.0 * u.arcsec)
    ax.coords[0].set_ticks(color="white", spacing=1.0 * u.arcsec)
    ax.coords[0].frame.set_color("white")
    
    # needs to be converted to pixels, that is the data unit
    self.add_beam(ax, image)
    ax.set_xlabel("rel. RA ['']")
    ax.set_ylabel("rel. Dec. ['']")
    
    # mark the center
    ax.plot(image.centerpix[0], image.centerpix[1], marker="x", color="0.6", linewidth=0.5, ms=2)
  
    #ax.axvline(image.centerpix[0],color="0.6",linewidth=0.8,linestyle=":")  
    #ax.axhline(image.centerpix[1],color="0.6",linewidth=0.8,linestyle=":")
  
      # FIXME: that would show the image like in the casaviewer
    # ax.invert_yaxis()  
  
    self._dokwargs(ax,**kwargs)
  
    if mJy==True:
      cformat="%3.0f"
      clabel="[mJy/beam km/s]"
    else:
      cformat="%5.2f"
      clabel="[Jy/beam km/s]"
      
  
    ticks = MaxNLocator(nbins=6).tick_values(vmin, vmax)
    CB = fig.colorbar(im, ax=ax, pad=0.02,
                    format=cformat, fraction=0.04)  
    CB.set_ticks(ticks)
    CB.set_label(clabel)
          
    self._closefig(fig)


  def plot_mom1_diff(self, imageObs,imageModel,imageDiff, zlim=[None, None],**kwargs):
    """
    Plots a fits image.
    """
  
    if imageObs is None or imageModel is None: return
  
    vmin = zlim[0]
    vmax = zlim[1]
    
    veldataObs = imageObs.data - imageObs.systemic_velocity
    veldataModel = imageModel.data - imageModel.systemic_velocity    
    if vmin is None: vmin = numpy.nanmin(veldataObs)
    if vmax is None: vmax = numpy.nanmax(veldataObs)  
  
    # wcsim=wcs.WCS(image.header)    
    # wcsrel=linear_offset_coords(wcsim, image.centerpix)
    fig, axes = plt.subplots(1, 3, subplot_kw=dict(projection=imageObs.wcsrel),
                           figsize=pplot.scale_figs([2.25,1.0]))
    for ax in axes:
      ax.set_facecolor("black")
      
    im = axes[0].imshow(veldataObs, cmap="seismic", vmin=vmin, vmax=vmax,origin="lower")  
    im1 = axes[1].imshow(veldataModel, cmap="seismic", vmin=vmin, vmax=vmax,origin="lower")
    im2 = axes[2].imshow(imageDiff.data, cmap="seismic", vmin=vmin, vmax=vmax,origin="lower")
    
    
    rms=numpy.nansum((imageDiff.data**2.0)/imageDiff.data.size)**0.5
    
    # print the velocities relative to the systemic velocities
    props = dict(boxstyle='round', facecolor='white', edgecolor="none")
    axes[0].text(0.95, 0.95, "Observation", transform=axes[0].transAxes, 
                 fontsize=6,verticalalignment='top', 
                 horizontalalignment="right", bbox=props)
    axes[1].text(0.95, 0.95, "Model", transform=axes[1].transAxes, 
                 fontsize=6,verticalalignment='top', 
                 horizontalalignment="right", bbox=props)
    axes[2].text(0.95, 0.95, "Residual (rms="+"{:4.1e}".format(rms)+")", 
                 transform=axes[2].transAxes, fontsize=6,
                 verticalalignment='top', horizontalalignment="right", bbox=props)
    
    for ax in axes:
      ax.coords[1].set_ticks(color="white", spacing=1.0 * u.arcsec)
      ax.coords[0].set_ticks(color="white", spacing=1.0 * u.arcsec)
      ax.coords[0].frame.set_color("white")
      ax.set_xlabel("rel. RA ['']")
    # needs to be converted to pixels, that is the data unit
    self.add_beam(axes[0], imageObs)
    
    axes[0].set_ylabel("rel. Dec. ['']")
    
    self._dokwargs(ax,**kwargs)
  
    # FIXME: that would show the image like in the casaviewer
    # ax.invert_yaxis()  
    ticks = MaxNLocator(nbins=6).tick_values(vmin, vmax)
    CB = fig.colorbar(im, ax=axes.ravel().tolist(), pad=0.02,
                    format="%5.2f", fraction=0.04)  
    CB.set_ticks(ticks)
    CB.set_label("velocity [km/s]")
        
    self._closefig(fig)


  def plot_mom1(self, image, zlim=[None, None],**kwargs):
    """
    Plots a fits image.
    """
  
    if image is None: return
  
    vmin = zlim[0]
    vmax = zlim[1]
    
    veldata = image.data - image.systemic_velocity
    if vmin is None: vmin = numpy.nanmin(veldata)
    if vmax is None: vmax = numpy.nanmax(veldata)
  
    # wcsim=wcs.WCS(image.header)    
    # wcsrel=linear_offset_coords(wcsim, image.centerpix)
    fig, ax = plt.subplots(1, 1, subplot_kw=dict(projection=image.wcsrel))
    ax.set_facecolor("black")
      
    im = ax.imshow(veldata, cmap="seismic", vmin=vmin, vmax=vmax,origin="lower")  
    ax.coords[1].set_ticks(color="white", spacing=1.0 * u.arcsec)
    ax.coords[0].set_ticks(color="white", spacing=1.0 * u.arcsec)
    ax.coords[0].frame.set_color("white")
    
    # needs to be converted to pixels, that is the data unit
    self.add_beam(ax, image)
    ax.set_xlabel("rel. RA ['']")
    ax.set_ylabel("rel. Dec. ['']")
    
    # mark the center
    ax.plot(image.centerpix[0], image.centerpix[1], marker="x", color="0.6", linewidth=0.5, ms=3)
  
    self._dokwargs(ax,**kwargs)
  
    # FIXME: that would show the image like in the casaviewer
    # ax.invert_yaxis()  
    ticks = MaxNLocator(nbins=6).tick_values(vmin, vmax)
    CB = fig.colorbar(im, ax=ax, pad=0.02,
                    format="%5.2f", fraction=0.04)  
    CB.set_ticks(ticks)
    CB.set_label("velocity [km/s]")
        
    self._closefig(fig)


  def plot_pv(self, image, zlim=[None, None], ylim=[None, None],**kwargs):
    """
    Plots a position-velocity diagram.
    """
    
    if image is None: return
  
    vmin = zlim[0]
    vmax = zlim[1]
    if vmin is None: vmin = numpy.nanmin(image.data)
    if vmax is None: vmax = numpy.nanmax(image.data)
  
  # FIXME: stolen from here ... https://github.com/astropy/astropy/issues/4529
    # mabe not the best way to go
    # Provide the velocities already at the init step (not in the plot routines)
    
    wcsa = wcs.WCS(image.header, naxis=(1))  
    xpix = numpy.arange(image.header['NAXIS1'])  
    arcsec = wcsa.wcs_pix2world(xpix, 0)  
    xticklabels = MaxNLocator(nbins="auto", symmetric=True, prune="both").tick_values(numpy.min(arcsec), numpy.max(arcsec))
      
    xticks, = wcsa.wcs_world2pix(xticklabels, 0)
    xticks = xticks   
     
    # prepare the y coordinate, want it relative to the center (systemic velocity)
    wcsvel = wcs.WCS(image.header)  
    zv = np.arange(image.header['NAXIS2'])
    # this is necessary because of imshow, it always plots in the original pixel scale  
    # create symetric ticks in pixel space
    symticks = MaxNLocator(nbins="auto", symmetric=True, prune="both").tick_values(numpy.min(zv - image.centerpix[1]), numpy.max(zv - image.centerpix[1]))
    # this gives then the wantes ytick positions.
    yticks = symticks + image.centerpix[1]
    
    # now convert the pixel to velocities to get the ylabels
    wv = wcsvel.wcs_pix2world(yticks, yticks, 0)[1]
    freq_to_vel = u.doppler_radio(image.header['RESTFRQ'] * u.Hz)
    yticklabels = (wv * u.Hz).to(u.km / u.s, equivalencies=freq_to_vel).value
    yticklabels = yticklabels - image.systemic_velocity
  
    # convert systemtic velocity to pixel coordinates  
    sysfrequ = (image.systemic_velocity * u.km / u.s).to(u.Hz, equivalencies=freq_to_vel)
    dummy, sysvelloc = wcsvel.wcs_world2pix(sysfrequ, sysfrequ, 0)  
    
    fig, ax = plt.subplots(1, 1)
    ax.set_facecolor("black")

    # width in pixel
    # FIXME: use proper units from fits file, here it is assume to be arcsec
    bwidth = (image.bmin * 3600 / image.header["CDELT1"] + 
            image.bmaj * 3600 / image.header["CDELT1"]) / 2.0  
    
    ar = AnchoredRectangle(ax.transData, width=bwidth, height=1.0,
                           loc=3, pad=0.5, borderpad=0.4, frameon=False)    
    ax.add_artist(ar)
      
    im = ax.imshow(image.data, cmap="inferno", vmin=vmin, vmax=vmax, aspect='auto',origin="lower")
    ax.set_xticks(xticks)
    ax.set_xticklabels(map(str, xticklabels))
    ax.set_yticks(yticks)
    ax.set_yticklabels([ "{:.2f}".format(x) for x in yticklabels ])
    
    ax.set_xlabel("offset ['']")
    ax.set_ylabel(r"vel $\mathrm{[km\,s^{-1}]}$")
    ax.axvline(image.centerpix[0], color="0.6", linewidth=0.8, linestyle=":")
    ax.axhline(sysvelloc, color="0.6", linewidth=0.8, linestyle=":")
  
    ax.tick_params(colors='white', labelcolor='black')
    for spine in ax.spines.values():
      spine.set_edgecolor('white')
    
    # FIXME: why is this
    ax.invert_yaxis()
    
    self._dokwargs(ax,**kwargs)
    
    ticks = MaxNLocator(nbins=6).tick_values(vmin, vmax)
    CB = fig.colorbar(im, ax=ax, pad=0.02,
                    format="%5.2f", fraction=0.04)  
    CB.set_ticks(ticks)
    CB.set_label("[Jy/beam]")
      
    self._closefig(fig)


  def plot_specprof(self, specprof, models=None, xlim=[None, None],modelNames=None,**kwargs):
    """
    Plots a spectral profile (histogram style).
    """
    
    if specprof is None: return


    
    fig, ax = plt.subplots(1, 1)

    if models is not None:
      if modelNames is None:
        modelNames=["model"+"{:0d}".format(i+1) for i in range(len(models))]        
      
      for model,label in zip(models,modelNames):
        x, y = self.specprof_xy_hist(model)
        ax.plot(x, y, label=label)
    
    x, y = self.specprof_xy_hist(specprof)
    ax.plot(x, y, label="Observation",color="black")  
    
    #pGrayBox=0.3
    #ax.fill_between(x, y *(1.0-pGrayBox), y * (1+pGrayBox), color='0.8')
    
   
    ax.set_xlim(xlim)
    ax.set_xlabel("velocity [km/s]")
    ax.set_ylabel("flux [Jy]")
  
    handles, labels = ax.get_legend_handles_labels()
    if len(labels) > 1:
      ax.legend(handles, labels, fancybox=False)
      
    self._dokwargs(ax,**kwargs)
  
    self._closefig(fig)


  def plot_radprof(self, radprof, models=None,modelNames=None,pmGrayBox=0.25,**kwargs):
    """
    Plots a radial profile.
    """
    if radprof is None: return
    
    fig, ax = plt.subplots(1, 1)
    
    #pGrayBox=0.3
    #ax.fill_between(radprof.arcsec, radprof.flux *(1.0-pGrayBox), radprof.flux * (1+pGrayBox), color='0.8')
    
    if models is not None:
      if modelNames is None:
        modelNames=["model"+"{:0d}".format(i+1) for i in range(len(models))]        
      
      for model,label in zip(models,modelNames):
        ax.errorbar(model.arcsec, model.flux, yerr=model.flux_err, label=label)

    if pmGrayBox is not None:
      ax.fill_between(radprof.arcsec, 
                      radprof.flux - radprof.flux*pmGrayBox, 
                      radprof.flux + radprof.flux*pmGrayBox, color='0.8')
      
    # ax.plot(radprof.arcsec,radprof.flux)
    ax.errorbar(radprof.arcsec, radprof.flux, 
                yerr=radprof.flux_err, label="observation",
                color="black")
    
    # indicate the beam 
    ax.set_xlabel("radius ['']")
    ax.set_ylabel("flux [$\mathrm{Jy/beam\,km/s}$]")
    ax.set_xlim(0, None)
    
    if radprof.bwidth is not None:
      ar = AnchoredRectangle(ax.transData, width=radprof.bwidth, height=0.0,
                           loc=3, pad=0.5, borderpad=0.4, frameon=False, color="0.6")    
      ax.add_artist(ar)
          
    handles, labels = ax.get_legend_handles_labels()
    if len(labels) > 1:
      ax.legend(handles, labels, fancybox=False)
      
    self._dokwargs(ax,**kwargs)
  
    self._closefig(fig)


  def specprof_xy_hist(self, specprof):
    """
    Produce x,y coordinates to plot spectral profile in histogram style.
    
    """
    x = list()
    y = list()  
    for i in range(len(specprof.vel)):
      if i == 0:
        x.append(specprof.vel[i])
        y.append(specprof.flux[i])
      else:
        xval = (specprof.vel[i - 1] + specprof.vel[i]) / 2.0      
        x.append(xval)
        y.append(specprof.flux[i - 1])
        x.append(xval)
        y.append(specprof.flux[i])
  
    # relative to the systemic velocity
    # should be zero if not set
    x.append(specprof.vel[-1])
    y.append(specprof.flux[-1])
      
    x = np.array(x) - specprof.systemic_velocity
    return x, np.array(y)
    
    
  def _closefig(self, fig):
    '''
    save and close the plot
    
    set the transparent attribut (used rcParam savefig.transparent)
    
    FIXME: make it general together with the other plot stuff
    '''    
    self.pdf.savefig(figure=fig, transparent=False)
    plt.close(fig)
    
  
class AnchoredRectangle(AnchoredOffsetbox):
    def __init__(self, transform, width, height, loc,
                pad=0.1, borderpad=0.1, prop=None, frameon=False, color="white"):
      """
      Draw a rectangle the size in data coordinate of the give axes.
  
      pad, borderpad in fraction of the legend font size (or prop)
      adapted from :class:`AnchoredEllipse`
      """
      self._box = AuxTransformBox(transform)
      self.rectangle = patches.Rectangle((2, 1), width, height, color=color)
      self._box.add_artist(self.rectangle)
  
      AnchoredOffsetbox.__init__(self, loc, pad=pad, borderpad=borderpad,
                                 child=self._box,
                                 prop=prop,
                                 frameon=frameon)


class AnchoredEllipse(AnchoredOffsetbox):
  def __init__(self, transform, width, height, angle, loc,
              pad=0.1, borderpad=0.1, prop=None, frameon=False, color="white"):
    """
    Draw an ellipse the size in data coordinate of the give axes.

    pad, borderpad in fraction of the legend font size (or prop)
    Copied from https://matplotlib.org/mpl_toolkits/axes_grid/api/anchored_artists_api.html
    Adapted it a bit (I think)
    
    Check how to use original class properly.
    
    """
    self._box = AuxTransformBox(transform)
    self.ellipse = patches.Ellipse((0, 0), width, height, angle, color=color)
    self._box.add_artist(self.ellipse)

    AnchoredOffsetbox.__init__(self, loc, pad=pad, borderpad=borderpad,
                               child=self._box,
                               prop=prop,
                               frameon=frameon)


