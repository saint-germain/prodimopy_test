from __future__ import print_function
from __future__ import division 

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker

import numpy as np
from math import pi

import prodimopy.plot as pplot

import astropy.units as u

class PlotModels(object):

  def __init__(self, pdf,
               colors=None,
               styles=None,
               markers=None,
               fs_legend=6,  # TODO: init it with the system default legend font size
               ncol_legend=0):

    self.pdf = pdf    
    if colors == None: 
      self.colors = ["b", "g", "r", "c", "m", "y", "k", "sienna", "lime", "pink", "DarkOrange", "Olive"]
    else:
      self.colors = colors
    
    if styles == None: 
      self.styles = ["-"] * len(self.colors)
    else:
      self.styles = styles
           
    if markers == None: 
      self.markers = ["D"] * len(self.colors)
    else:
      self.markers = markers    
      
    self.fs_legend = fs_legend
    self.ncol_legend = ncol_legend              


  def _legend(self, ax):
    '''
    plots the legend, deals with multiple columns
    '''
    handles, labels = ax.get_legend_handles_labels()
    ncol = 1
    if self.ncol_legend > 1 and len(labels) > self.ncol_legend:
      ncol = int(len(labels) / self.ncol_legend)
    ax.legend(handles, labels, loc="best", fancybox=False, ncol=ncol, fontsize=self.fs_legend)

  def _dokwargs(self,ax,**kwargs):
    '''
    Handles the passed kwargs elements (assumes that defaults are already set)
    '''
    if "ylim" in kwargs: 
      ax.set_ylim(kwargs["ylim"])
      
    if "xlim" in kwargs: 
      ax.set_xlim(kwargs["xlim"])
              
    if "xlog" in kwargs:
      if kwargs["xlog"]: 
        ax.semilogx()
      else:
        ax.set_xscale("linear")
      
    if "ylog" in kwargs:
      if kwargs["ylog"]: 
        ax.semilogy()
      else:              
        ax.set_yscale("linear")
      
    if "title" in kwargs:
      if kwargs["title"].strip() != "":
        ax.set_title(kwargs["title"])      
                 
  def plot_lines(self, models, lineidents, useLineEstimate=True,jansky=False,**kwargs):
    '''
    plots a selection of FlineEstimates
    '''
    print("PLOT: plot_lines ...")
    fig, ax = plt.subplots(1, 1)  
         
    imodel = 0          
    lticks = list()
    # for some reason this dummy is require
    lticks.append("")
    for model in models:                  
      iline = 1  
      x = list()
      y = list()
      for ident in lineidents: 
        line=None             
        if useLineEstimate:
          line = model.getLineEstimate(ident[0], ident[1])          
        else:  # use the line fluxes, uses only the wavelength
          line = model.getLine(ident[1])
        x.append(iline)    
                            
        # Convert lineflux to Jansky
        
        if jansky:          
          #linetst=DataLineEstimate("HCO+", 840.38, 5, 4, 7.55e-20)
          #print(linetst)
          #print(self.si_to_jansky(linetst))
          #print(self.si_to_jansky(line))
          y.append(line.flux_Jy())        
        else:
          y.append(line.flux)
        
        if imodel == 0:
          # lticks.append(r"$\mathrm{"+pplot.spnToLatex(ident[0])+r"}$ "+r"{:.2f}".format(line.wl))
          lticks.append(r"$\mathrm{" + line.ident + r"}$ " + r"{:.2f}".format(line.wl))   
   
#         if imodel == (len(models)-1):
#           print iline, line.flux
#           ax.errorbar([iline],[line.flux],yerr=[line.flux/3.0,line.flux*3.0],fmt=".",ms=0.0,color="lightgrey",linewidth=10)    
        
        iline = iline + 1  

      mew=None
      ms=5
      if self.markers[imodel]=="+" or self.markers[imodel]=="x":
        mew=2
        ms=5      
      ax.plot(x, y, marker=self.markers[imodel], linestyle='None', ms=ms,mew=mew, color=self.colors[imodel], markeredgecolor=self.colors[imodel], label=model.name)      
         
      imodel = imodel + 1            
           
   
    ax.set_xlim(0.5, iline - 0.5)
    
    loc = plticker.MultipleLocator(base=1.0)  # this locator puts ticks at regular intervals
    ax.xaxis.set_major_locator(loc)
    
    if "ylim" in kwargs: 
      ax.set_ylim(kwargs["ylim"])
    #else:         
      #ax.set_ylim(1.e-24, 1.e-18)
      
    ax.semilogy()    
    # ax.set_xlabel(r"line")
    
    if jansky:   
      ax.set_ylabel(r" line flux [Jy km$\,$s$^{-1}$]")
    else:
      ax.set_ylabel(r" line flux [W$\,$m$^{-2}$]")
     
    ax.set_xticklabels(lticks, rotation='45', minor=False)
    zed = [tick.label.set_fontsize(7) for tick in ax.xaxis.get_major_ticks()]
     
    self._legend(ax)
    
    if "title" in kwargs and kwargs["title"] != None:
      ax.set_title(kwargs["title"])
     
    self.pdf.savefig()
    plt.close(fig)
  
  def plot_NH(self, models, **kwargs):
    '''
    Plots the total vertical column density for the given species for all the models
    as a function of radius
    '''
    print("PLOT: plot_NH ...")
    fig, ax = plt.subplots(1, 1)  
    
    xmin = 1.e100
    xmax = 0
        
    iplot = 0    
    for model in models:
      x = model.x[:, 0]
      y = model.NHver[:, 0]
      
      # print y.min() 
      ax.plot(x, y, self.styles[iplot], marker=None, color=self.colors[iplot],alpha=0.7, label=model.name)
          
      iplot = iplot + 1
      
      if min(x) < xmin: xmin = min(x)
      if max(x) > xmax: xmax = max(x)
       
  
    ax.set_xlim(xmin, xmax)
         
    ax.semilogy()     
    ax.set_xlabel(r"r [AU]")
    ax.set_ylabel(r"N$_\mathrm{<H>}$ cm$^{-2}$")
  
    self._dokwargs(ax,**kwargs)
    self._legend(ax)
    
    self.pdf.savefig(transparent=mpl.rcParams['savefig.transparent'])
    plt.close(fig)
    
  def plot_tcdspec(self, models, species, xlog=True, title=None, **kwargs):
    '''
    Plots the total vertical columndensity for the given species for all the models
    as a function of the radius
    '''
    print("PLOT: plot_tcdspec ...")
    fig, ax = plt.subplots(1, 1)  
    
    xmin = 1.e100
    xmax = 0
              
    iplot = 0    
    for model in models:
      if species in model.spnames:
        x = model.x[:, 0]
                    
        y = model.cdnmol[:, 0, model.spnames[species]]
        
        
        ax.plot(x, y, self.styles[iplot], marker=None, color=self.colors[iplot], label=model.name)
            
        iplot = iplot + 1
        
        if min(x) < xmin: xmin = min(x)
        if max(x) > xmax: xmax = max(x)      
    
    if iplot==0:
      print("WARN: Species "+species+" not found in any model!")
      plt.close(fig)
      return
    
    ax.fill_between(x, y / 3.0, y * 3.0, color='0.8')     
    ax.set_xlim(xmin, xmax)        
    ax.semilogy()
    
    ax.set_xlabel(r"r [AU]")
    # ax.set_ylabel(r"N$_\mathrm{"+pplot.spnToLatex(species)+"}$ cm$^{-2}$")
    ax.set_ylabel(r"N$_\mathrm{" + pplot.spnToLatex(species) + "}$ " + "[cm$^{-2}$]")
  
    self._dokwargs(ax,**kwargs)          
    self._legend(ax)
  
    self.pdf.savefig()
    plt.close(fig)  
    
  def plot_tauline(self, models, lineIdent, xlog=True, **kwargs):
    '''
    Plots the line optical depth as a function of radius for a given line
    for all the models
    '''
    print("PLOT: plot_tauline ...")
    fig, ax = plt.subplots(1, 1)  
    
    xmin = 1.e100
    xmax = 0
        
    iplot = 0    
    for model in models:
      x = model.x[:, 0]
      lineEstimate = model.getLineEstimate(lineIdent[0], lineIdent[1])
      y = list()
      ytauDust = list()
      for rInfo in lineEstimate.rInfo:
        y.append(rInfo.tauLine)
        ytauDust.append(rInfo.tauDust)
      
      ax.axhline(y=1.0, linestyle="-", color="black", linewidth=0.5)
      ax.plot(x, y, self.styles[iplot], marker=None, color=self.colors[iplot], label=model.name)
      # ax.plot(x,ytauDust,"--",marker=None,color="black")
      
          
      iplot = iplot + 1
      
      if min(x) < xmin: xmin = min(x)
      if max(x) > xmax: xmax = max(x)
       
  
    ax.set_xlim(xmin, xmax)
    
    if "xlim" in kwargs: ax.set_xlim(kwargs["xlim"])    
    if "xmin" in kwargs: ax.set_xlim(xmin=kwargs["xmin"])
    if "ylim" in kwargs: ax.set_ylim(kwargs["ylim"])  
    
    if xlog == True: ax.semilogx()
        
    ax.semilogy()
    
    ax.set_xlabel(r"r [AU]")    
    ax.set_ylabel(r"$\mathrm{\tau_{line}}$")
    ax.set_title(lineEstimate.ident + " " + "{:.2f}".format(lineEstimate.wl) + " $\mathrm{\mu m}$")
  
    self._legend(ax)  
  
    self.pdf.savefig()
    plt.close(fig)  
  
  def plot_avgabun(self,models,species,**kwargs):
    '''
    Plots the average abundance of the species as a function of radius
    the avergae abundance is given by NHver(species)/total nhver
    '''
    print("PLOT: plot_avgabun ...")     
          
    fig, ax = plt.subplots(1, 1)     

    iplot = 0
    for model in models:                 
      # get the species
      if (species in model.spnames):   
        y=model.cdnmol[:,0,model.spnames[species]]
        y=y/model.NHver[:,0]   
        x = model.x[:,0]            
        
        ax.plot(x, y, self.styles[iplot], marker=None, color=self.colors[iplot], label=model.name)
    
        iplot = iplot + 1              
    
    if iplot==0:
      print("Species "+species+ " not found in any model!")
      return 
            
    ax.set_xlabel(r"r [AU]")
    ax.set_ylabel(r"average $\epsilon(\mathrm{" + pplot.spnToLatex(species) + "})$")
    
    # do axis style
    ax.semilogy()     
    
    self._dokwargs(ax,**kwargs)    
    self._legend(ax)
    
    self.pdf.savefig()
    plt.close(fig) 
    
    
  def plot_abunvert(self, models, r, species, **kwargs):
    '''
    Plot vertical abundances at a certain radius.
    '''
         
    print("PLOT: plot_abunvert ...")     
          
    rstr = r"r$\approx${:.2f} AU".format(r)   
    
    fig, ax = plt.subplots(1, 1)     

    iplot = 0
    xmin = 1.e100
    xmax = 0 
    for model in models:
      # closed radial point to given radius
      ix = (np.abs(model.x[:, 0] - r)).argmin()
      
      old_settings = np.seterr(divide='ignore')     
      x = np.log10(model.NHver[ix, :])      
      np.seterr(**old_settings)  # reset to default

      y = model.nmol[ix, :, model.spnames[species]] / model.nHtot[ix, :]
      
      ax.plot(x, y, self.styles[iplot], marker=None, color=self.colors[iplot], label=model.name)
                      
      iplot = iplot + 1
      
      if min(x) < xmin: xmin = min(x)
      if max(x) > xmax: xmax = max(x)
      
    # set the limits
      
    if "xlim" in kwargs:     
      ax.set_xlim(kwargs["xlim"])
    else:
      ax.set_xlim([17.5, x.max()])
              
    if "ylim" in kwargs: ax.set_ylim(kwargs["ylim"])
     
    # print ax.get_xlim()
  
#     ax2 = ax.twiny()
#     ax2.set_xlabel("z/r")
#     ax2.set_xlim(ax.get_xlim())
#     #ax2.set_xticks(ax.get_xticks())    
#     ax2.set_xticklabels(["{:.2f}".format(x) for x in nhver_to_zr(ix, ax.get_xticks(), model)])
    
    ax.set_xlabel(r"$\log$ N$_\mathrm{H}$ [cm$^{-2}$] @" + rstr)
    ax.set_ylabel(r"$\epsilon(\mathrm{" + pplot.spnToLatex(species) + "})$")    
    
    # do axis style
    ax.semilogy()     
    
    self._legend(ax)
    # ax.text(0.025, 0.025, rstr,
    #   verticalalignment='bottom', horizontalalignment='left',
    #   transform=ax.transAxes,alpha=0.75)
      
    self.pdf.savefig()
    plt.close(fig) 
  
  
  def plot_midplane(self, models, fieldname, ylabel, **kwargs):
    '''
    Plots a quantitiy in in the midplane as a function of radius
    fieldname is any field in Data_ProDiMo
    '''
    print("PLOT: plot_midplane ...")
    fig, ax = plt.subplots(1, 1)      
    
    iplot = 0
    xmin = 1.e100
    xmax = 0
    ymin = 1.e100
    ymax = -1.e00 
    for model in models:           
      x = model.x[:, 0]
      y = getattr(model, fieldname)[:, 0]                    
      
      ax.plot(x, y, self.styles[iplot], marker=None, color=self.colors[iplot], label=model.name)
                      
      iplot = iplot + 1
      
      if min(x) < xmin: xmin = min(x)
      if max(x) > xmax: xmax = max(x)
      if min(y) < ymin: ymin = min(y)
      if max(y) > ymax: ymax = max(y)

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin, ymax)              
    ax.semilogy()
            
    ax.set_xlabel(r"r [AU]")    
    ax.set_ylabel(ylabel)    
    
    self._dokwargs(ax, **kwargs) 
    self._legend(ax)
    
    self.pdf.savefig()
    plt.close(fig)    
  
  def plot_vertical(self, models, r, fieldname, ylabel, ylog=True, **kwargs):
    '''
    Plots a quantity (fieldname) as a function of height (z/r) at a certain
    radius.    
    '''
    print("PLOT: plot_vertical ...")
    rstr = r"r$\approx${:.1f} AU".format(r) 
    
    fig, ax = plt.subplots(1, 1)      
    
    iplot = 0
    xmin = 1.e100
    xmax = -1.e100
    ymin = 1.e100
    ymax = -1.e100 
    for model in models:    
      # closed radial point to given radius
      ix = (np.abs(model.x[:, 0] - r)).argmin()
       
      x = model.z[ix, :] / model.x[ix, 0]
      y = getattr(model, fieldname)[ix, :]                    
      
      ax.plot(x, y, self.styles[iplot], marker=None, color=self.colors[iplot], label=model.name)
                      
      iplot = iplot + 1
      
      if min(x) < xmin: xmin = min(x)
      if max(x) > xmax: xmax = max(x)
      if min(y) < ymin: ymin = min(y)
      if max(y) > ymax: ymax = max(y)


    if "xlim" in kwargs: 
      ax.set_xlim(kwargs["xlim"])
    # else:
    #  ax.set_ylim(xmin,xmax)
      
    if "ylim" in kwargs: 
      ax.set_ylim(kwargs["ylim"])
    # else:
    #  ax.set_ylim(ymin,ymax)
    
    if ylog: ax.semilogy()
    
    ax.invert_xaxis()       
                
    ax.set_xlabel(r"z/r @ " + rstr)    
    ax.set_ylabel(ylabel)    
    
    self._legend(ax)
    
    self.pdf.savefig()
    plt.close(fig)    
    
    
  def plot_abunradial(self, models, species, **kwargs):
    '''
    Plots the abundance of the given species as a function of radius in the
    midplane of the disk
    '''  
    print("PLOT: plot_abunradial ...")
    fig, ax = plt.subplots(1, 1)      
    
    iplot = 0
    xmin = 1.e100
    xmax = 0 
    for model in models:           
      x = model.x[:, 0]            
      y = model.nmol[:, 0, model.spnames[species]] / model.nHtot[:, 0]
      
      ax.plot(x, y, self.styles[iplot], marker=None, color=self.colors[iplot], label=model.name)
                      
      iplot = iplot + 1
      
      if min(x) < xmin: xmin = min(x)
      if max(x) > xmax: xmax = max(x)
          
    ax.set_xlim(xmin,xmax)
    ax.semilogy()    
            
    ax.set_xlabel(r"r [AU]")    
    ax.set_ylabel(r"midplane $\epsilon(\mathrm{" + pplot.spnToLatex(species) + "})$")    
    
    self._dokwargs(ax, **kwargs)  
    self._legend(ax)
    
    self.pdf.savefig()
    plt.close(fig)    
    
  def plot_sed(self, models,plot_starSpec=True,**kwargs): 
    '''
    Plots the seds and the StarSpectrum
    '''  
    print("PLOT: plot_sed ...")
    fig, ax = plt.subplots(1, 1)      
    
    iplot = 0
    xmin = 1.e100
    xmax = 0
    ymin = 1.e100
    ymax = -1.e00 
    for model in models:       
      if model.sed == None:
        continue    
      # only use every 5 element to speed up plotting
      x = model.sed.lam
      y = model.sed.nu*model.sed.fnuErg      
      dist = ((model.sed.distance*u.pc).to(u.cm)).value                          
      ax.plot(x, y, self.styles[iplot], marker=None, color=self.colors[iplot], label=model.name)
      
      if plot_starSpec:
        # scale input Stellar Spectrum to the distance for comparison to the SED
        r= ((model.starSpec.r*u.R_sun).to(u.cm)).value      
        
        xStar = model.starSpec.lam[0::10]
        ystar= (model.starSpec.nu*model.starSpec.Inu)[0::10]
        yStar = ystar*(r**2.0*pi*dist**(-2.0))                                
        ax.plot(xStar, yStar, self.styles[iplot], marker="*",ms=2,markeredgecolor=self.colors[iplot],color=self.colors[iplot])
                      
      iplot = iplot + 1
      
      if min(x) < xmin: xmin = min(x)
      if max(x) > xmax: xmax = max(x)
      if y[-1] < ymin: ymin = y[-1]
      if max(y) > ymax: ymax = max(y)
      
    if iplot == 0: return  
      
    # set defaults, can be overwritten by the kwargs
    ax.set_xlim(xmin,xmax)
    ax.set_ylim([ymin,None])
    ax.semilogx()
    ax.semilogy()
    ax.set_xlabel(r"wavelength [$\mu$m]")    
    ax.set_ylabel(r"$\mathrm{\nu F_{\nu}\,[erg\,cm^{-2}\,s^{-1}]}$")    
      
    self._dokwargs(ax, **kwargs)                
    
    self._legend(ax)
    
    self.pdf.savefig()
    plt.close(fig)      
    
   
  def plot_line_profil(self,models,wl,**kwargs):
    '''
    Plots the line profile for the given line (id by wavelength) 
    '''  
    print("PLOT: plot_line_profile ...")
    fig, ax = plt.subplots(1, 1)      
    
    iplot = 0
    xmin = 1.e100
    xmax = -1.e100 
    ymin = 1.e100
    ymax = -1.e100 
    for model in models:
      line=model.getLine(wl)       
      if line==None: continue    
      
      # text for the title
      if iplot==0:
        lintxt=line.species+"@"+"{:.2f}".format(line.wl)+" $\mathrm{\mu m}$" 
      x = line.profile.velo           
      y = line.profile.flux-line.fcont  # remove the continuum
      
      ax.plot(x, y, self.styles[iplot], marker=None, color=self.colors[iplot], label=model.name)
                      
      iplot = iplot + 1
      
      if min(x) < xmin: xmin = min(x)
      if max(x) > xmax: xmax = max(x)
      if min(y) < ymin: ymin = min(y)
      if max(y) > ymax: ymax = max(y)
    
    if iplot==0: 
      print("WARN: No lines found: ")
      return   

    # set defaults, can be overwritten by the kwargs
    ax.set_xlim(xmin,xmax)
    ax.set_ylim([ymin,ymax*1.1])

    ax.set_xlabel(r"$\mathrm{velocity\,[km s^{-1}}$]")    
    ax.set_ylabel(r"$\mathrm{flux\,[Jy]}$")    
      
    self._dokwargs(ax, **kwargs)                
    
    ax.text(0.03, 0.97,lintxt, ha='left', va='top', transform=ax.transAxes)
    
    self._legend(ax)
               
    self.pdf.savefig()
    plt.close(fig)      
    
