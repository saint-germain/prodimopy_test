from __future__ import print_function
from __future__ import division 

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker

import numpy as np
from math import pi

import prodimopy.plot as pplot

import astropy.units as u
import astropy.constants as const
from matplotlib import patches
import types


class PlotModels(object):

  def __init__(self, pdf,
               colors=None,
               styles=None,
               markers=None,
               fs_legend=None,  # TODO: init it with the system default legend font size
               ncol_legend=0):

    self.pdf = pdf    
    if colors == None: 
      # use system default colors not fixed ones.
      # FIXME: could be done more elegant (simply set colors to None, needs changes in the plot routines)
      self.colors=[None]*20
      #self.colors = ["b", "g", "r", "c", "m", "y", "k", "sienna", "lime", "pink", "DarkOrange", "Olive","brown"]
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
      
    if fs_legend is None:
      self.fs_legend = mpl.rcParams['legend.fontsize']
    else:
      self.fs_legend=fs_legend
    self.ncol_legend = ncol_legend              

  def _sfigs(self,wfac=2.0,hfac=1.0):
    figsize=mpl.rcParams['figure.figsize']  
        
    return (figsize[0]*wfac,figsize[1]*hfac)

  def _legend(self, ax,**kwargs):
    '''
    plots the legend, deals with multiple columns
    '''
    if "slegend" in kwargs:
      if not kwargs["slegend"]: return 
    
    handles, labels = ax.get_legend_handles_labels()
    ncol = 1
    if self.ncol_legend > 1 and len(labels) > self.ncol_legend:
      ncol = int(len(labels) / self.ncol_legend)
    leg=ax.legend(handles, labels, loc="best", fancybox=False, ncol=ncol, fontsize=self.fs_legend)
    lw=mpl.rcParams['axes.linewidth']
    leg.get_frame().set_linewidth(lw)    
    
  def _set_dashes(self,line):
    '''
    Utility routine to set the dashed stuff for a line. 
    This routine should be used instead of using set_dashes directly, if the 
    default value (still hardcoded) should be used.         
    '''  
    line.set_dashes((4,4))

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
        
    if "xlabel" in kwargs:
      ax.set_xlabel(kwargs["xlabel"])  

    if "ylabel" in kwargs:
      ax.set_ylabel(kwargs["ylabel"])
      
    if "title" in kwargs:
      if kwargs["title"].strip() != "":
        ax.set_title(kwargs["title"])      
                 
  def plot_lines(self, models, lineidents, useLineEstimate=True,jansky=False,lineObs=None,**kwargs):
    '''
    Plots a selection of lines or lineEstimates
    TODO: split lines and lineestimate plots
    TODO: if not FlineEstimates than it would be also possible to plot all lines for which 
          line transfer is done
    '''
    print("PLOT: plot_lines ...")
    
    if len(lineidents)>10:
      fig, ax = plt.subplots(1, 1,figsize=self._sfigs())
    else:
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
          line = model.getLine(ident[1],ident=ident[0])
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
          # FIXME: spnToLatex does not work here with all line names ... 
          lticks.append(r"$\mathrm{" + line.ident + r"}$ " + r"{:.2f}".format(line.wl))   
           
        iline = iline + 1  

      mew=None
      ms=4
      if self.markers[imodel]=="+" or self.markers[imodel]=="x":        
        mew=1.5
        ms=5     
      elif self.markers[imodel]=="*" or self.markers[imodel]=="p":
        ms=5
              
      ax.plot(x, y, marker=self.markers[imodel], linestyle='None', ms=ms,mew=mew, 
              color=self.colors[imodel], markeredgecolor=self.colors[imodel], label=model.name,
              zorder=10)      
         
      imodel = imodel + 1  

    # boxes for factor 3 and 10 relative to the last model
    for i in range(len(x)):
      xc=x[i]-0.3
      yc=y[i]/3.0
      height=y[i]*3.0 -y[i]/3.0     
      width=0.6   
      ax.add_patch(patches.Rectangle((xc,yc),width,height,color="0.8"))

      height=y[i]*10.0 -y[i]/10.0
      yc=y[i]/10.0        
      ax.add_patch(patches.Rectangle((xc,yc),width,height,color="0.5",fill=False,linewidth=0.5))          
    
    if lineObs != None:
      nlines=len(lineObs)        
      ylinesObs=[item.flux for item in lineObs]
      ylinesObsErr=np.zeros(shape=(2,nlines))
      ylinesObsErr2=np.zeros(shape=(2,nlines))    
      ylinesObsUl=list()
      # for the errors and errorbars
      for i in range(nlines):
        ylinesObsErr[:,i]=lineObs[i].flux_err              
        ylinesObsErr2[0,i]=(lineObs[i].flux)/2.0
        ylinesObsErr2[1,i]=lineObs[i].flux      
        #print linesObs[i].flag
        if lineObs[i].flag == "ul":      
          ylinesObsUl.append(ylinesObs[i]/0.5)
          ylinesObsErr[0,i]=ylinesObs[i]*0.4
          ylinesObsErr[1,i]=0.0      
          ylinesObsErr2[0,i]=lineObs[i].flux*(1.0-1.e-10)      
        else:
          ylinesObsUl.append(0.0)  
    
      # the observations
      # takes the x frmo above
      ax.errorbar(x,ylinesObs,yerr=ylinesObsErr2,fmt=".",ms=0.0,color="lightgrey",linewidth=10,zorder=0)
      ax.errorbar(x,ylinesObs,yerr=ylinesObsErr,uplims=ylinesObsUl,fmt="o",ms=4.0,color="black",capsize=2,label="Obs.")           
   
    ax.set_xlim(0.5, iline - 0.5)
    
    loc = plticker.MultipleLocator(base=1.0)  # this locator puts ticks at regular intervals
    ax.xaxis.set_major_locator(loc)
    
    if "ylim" in kwargs: 
      ax.set_ylim(kwargs["ylim"])
      
    ax.semilogy()        
    
    if jansky:   
      ax.set_ylabel(r" line flux [Jy km$\,$s$^{-1}$]")
    else:
      ax.set_ylabel(r" line flux [W$\,$m$^{-2}$]")
      
    xgrid=np.array(x)      
    ax.vlines(xgrid-0.5,ymin=ax.get_ylim()[0],ymax=ax.get_ylim()[1],linestyle="solid",linewidth=0.5,color="grey")  
    ax.yaxis.grid(color="grey")  
  
     
    ax.set_xticklabels(lticks, rotation='70', minor=False)
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
      line, = ax.plot(x, y, self.styles[iplot], marker=None, color=self.colors[iplot],alpha=0.7, label=model.name)
      if line.is_dashed(): self._set_dashed(line)
          
      iplot = iplot + 1
      
      if min(x) < xmin: xmin = min(x)
      if max(x) > xmax: xmax = max(x)
       
  
    ax.set_xlim(xmin, xmax)
         
    ax.semilogy()     
    ax.set_xlabel(r"r [au]")
    ax.set_ylabel(r"N$_\mathrm{<H>}$ cm$^{-2}$")
  
    self._dokwargs(ax,**kwargs)
    self._legend(ax,**kwargs)
    
    self.pdf.savefig(transparent=mpl.rcParams['savefig.transparent'])
    plt.close(fig)
    
  def plot_tcdspec(self, models, species, relToH=False, **kwargs):
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
        if relToH==True:
          y=y/model.NHver[:,0]
        
        linewidth=1.5
        if self.styles[iplot]=="--": linewidth=2.5
                            
        line, = ax.plot(x, y, self.styles[iplot], marker=None,linewidth=linewidth, 
                color=self.colors[iplot], label=model.name)
        if line.is_dashed(): self._set_dashes(line)      
            
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
    
    ax.set_xlabel(r"r [au]")    
    if relToH==True:
      ax.set_ylabel(r"average $\mathrm{\epsilon(" + pplot.spnToLatex(species) + ")}$ ")
    else:
      ax.set_ylabel(r"N$_\mathrm{" + pplot.spnToLatex(species) + "}$ " + "[cm$^{-2}$]")
  
    self._dokwargs(ax,**kwargs)          
    self._legend(ax,**kwargs)
  
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
    
    ax.set_xlabel(r"r [au]")    
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
            
    ax.set_xlabel(r"r [au]")
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
          
    rstr = r"r$\approx${:.2f} au".format(r)   
    
    fig, ax = plt.subplots(1, 1)     

    iplot = 0
    for model in models:
      # closed radial point to given radius
      ix = (np.abs(model.x[:, 0] - r)).argmin()
      
      old_settings = np.seterr(divide='ignore')     
      x = np.log10(model.NHver[ix, :])      
      np.seterr(**old_settings)  # reset to default

      y = model.nmol[ix, :, model.spnames[species]] / model.nHtot[ix, :]
      
      line, = ax.plot(x, y, self.styles[iplot], marker=None, color=self.colors[iplot], label=model.name)
      if line.is_dashed(): self._set_dashes(line)
                      
      iplot = iplot + 1
      
    
    ax.set_xlim([17.5, x.max()])
    # print ax.get_xlim()
  
#     ax2 = ax.twiny()
#     ax2.set_xlabel("z/r")
#     ax2.set_xlim(ax.get_xlim())
#     #ax2.set_xticks(ax.get_xticks())    
#     ax2.set_xticklabels(["{:.2f}".format(x) for x in nhver_to_zr(ix, ax.get_xticks(), model)])
    
    ax.set_xlabel(r"$\mathrm{\log N_{<H>} [cm^{-2}]}$ @" + rstr)
    ax.set_ylabel(r"$\mathrm{\epsilon(" + pplot.spnToLatex(species) + "})$")    
    
    # do axis style
    ax.semilogy()     
    
    self._dokwargs(ax,**kwargs)
    self._legend(ax)
    # ax.text(0.025, 0.025, rstr,
    #   verticalalignment='bottom', horizontalalignment='left',
    #   transform=ax.transAxes,alpha=0.75)
      
    self.pdf.savefig()
    plt.close(fig) 
  
  
  def plot_radial(self, models, fields, ylabel,**kwargs):
    '''
    Plots a quantitiy in radial direction. Fields must have the same number
    of entries as models and must contain arrays with the dimension of nx   
    '''
    print("PLOT: radial ...")
    fig, ax = plt.subplots(1, 1)      
    
    iplot = 0
    xmin = 1.e100
    xmax = 0
    ymin = 1.e100
    ymax = -1.e00 
    for model in models:           
      x = model.x[:, 0]
      y= fields[iplot]
      
      line, = ax.plot(x, y, self.styles[iplot], marker=None, color=self.colors[iplot], label=model.name)
      if line.is_dashed(): self._set_dashes(line)
                      
      iplot = iplot + 1
      
      if min(x) < xmin: xmin = min(x)
      if max(x) > xmax: xmax = max(x)
      if min(y) < ymin: ymin = min(y)
      if max(y) > ymax: ymax = max(y)
  
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin, ymax)              
    ax.semilogy()
            
    ax.set_xlabel(r"r [au]")    
    ax.set_ylabel(ylabel)    
    
    self._dokwargs(ax, **kwargs) 
    self._legend(ax,**kwargs)
    
    self.pdf.savefig()
    plt.close(fig) 
  
  
  def plot_midplane(self, models, fieldname, ylabel, species=None,**kwargs):
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
      if species!=None:
        y = getattr(model, fieldname)[:, 0,model.spnames[species]]
      else:
        y = getattr(model, fieldname)[:, 0]                    
      
      line, = ax.plot(x, y, self.styles[iplot], marker=None, color=self.colors[iplot], label=model.name)
      if line.is_dashed(): self._set_dashes(line) 
                      
      iplot = iplot + 1
      
      if min(x) < xmin: xmin = min(x)
      if max(x) > xmax: xmax = max(x)
      if min(y) < ymin: ymin = min(y)
      if max(y) > ymax: ymax = max(y)

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin, ymax)              
    ax.semilogy()
            
    ax.set_xlabel(r"r [au]")    
    ax.set_ylabel(ylabel)    
    
    self._dokwargs(ax, **kwargs) 
    self._legend(ax,**kwargs)
    
    self.pdf.savefig()
    plt.close(fig)    
  
  def plot_vertical_nH(self, models, r, field, ylabel, species=None,**kwargs):
    '''
    Plots a quantity (field) as a function of height (column density) at a certain
    radius.    
    If field is a string it is interpreted as a field name in the ProDiMo 
    data structure. If field is a list the list is directly used. 
    The list needs to contain 2D arrays with the same shape as other ProDiMo fields
    '''
    print("PLOT: plot_vertical_nH ...")
    rstr = r"r$\approx${:.1f} au".format(r) 
    
    fig, ax = plt.subplots(1, 1)      
    
    iplot = 0
    xmin = 1.e100
    xmax = -1.e100
    ymin = 1.e100
    ymax = -1.e100 
    for model in models:    
      # closest radial point to given radius
      ix = (np.abs(model.x[:, 0] - r)).argmin()
      
      old_settings=np.seterr(divide='ignore')     
      x=np.log10(model.NHver[ix,:])    
      np.seterr(**old_settings)  # reset to default
               
      if species==None:
        
        isstr=False
        try: # this is for pyhton 2/3 compatibility
          istrs=isinstance(field, basestring)
        except NameError:
          return isinstance(field, str)
        
        if isstr:
          y = getattr(model, field)[ix, :]
        else: 
          y=(field[iplot])[ix,:]          
      else:
        y = getattr(model, field)[ix, :,model.spnames[species]]
                                  
      ax.plot(x, y, self.styles[iplot], marker=None, color=self.colors[iplot], label=model.name)
                      
      iplot = iplot + 1
      
      if min(x) < xmin: xmin = min(x)
      if max(x) > xmax: xmax = max(x)
      if min(y) < ymin: ymin = min(y)
      if max(y) > ymax: ymax = max(y)


    
    #ax.semilogx()
    ax.semilogy()
    ax.set_xlabel(r"$\mathrm{N_{H,ver}\,[cm^{-2}]}$ at "+rstr)                        
    ax.set_ylabel(ylabel)    
    
    self._dokwargs(ax,**kwargs)
    self._legend(ax,**kwargs)
    
    self.pdf.savefig()
    plt.close(fig)   
  
  
  def plot_vertical(self, models, r, fieldname, ylabel, species=None,ylog=True,zr=True,**kwargs):
    '''
    Plots a quantity (fieldname) as a function of height (z/r) at a certain
    radius.    
    '''
    print("PLOT: plot_vertical ...")
    rstr = r"r$\approx${:.1f} au".format(r) 
    
    fig, ax = plt.subplots(1, 1)      
    
    iplot = 0
    xmin = 1.e100
    xmax = -1.e100
    ymin = 1.e100
    ymax = -1.e100 
    for model in models:    
      # closest radial point to given radius
      ix = (np.abs(model.x[:, 0] - r)).argmin()
       
      if zr: 
        x = model.z[ix, :] / model.x[ix, 0]
      else: 
        x = model.z[ix, :]
        
      if species==None:
        y = getattr(model, fieldname)[ix, :]
      else:
        y = getattr(model, fieldname)[ix, :,model.spnames[species]]
                                  
      ax.plot(x, y, self.styles[iplot], marker=None, color=self.colors[iplot], label=model.name)
                      
      iplot = iplot + 1
      
      if min(x) < xmin: xmin = min(x)
      if max(x) > xmax: xmax = max(x)
      if min(y) < ymin: ymin = min(y)
      if max(y) > ymax: ymax = max(y)


    if "xlim" in kwargs: 
      ax.set_xlim(kwargs["xlim"])
    else:
      if not zr:
        ax.set_xlim(xmin,xmax)
      
    if "ylim" in kwargs: 
      ax.set_ylim(kwargs["ylim"])
    #else:
    #  ax.set_ylim(ymin,ymax)
    
    if ylog: ax.semilogy()
    
    if zr:
      ax.invert_xaxis()
      ax.set_xlabel(r"z/r @ " + rstr)
    else:
      ax.set_xlabel(r"z [au] @ " + rstr)      
      ax.semilogx()       
      ax.invert_xaxis()
                        
    ax.set_ylabel(ylabel)    
    
    # FIXME: make it general
    #self._dokwargs(ax,kwargs)
    self._legend(ax,**kwargs)
    
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
      
      if iplot==0 or iplot == (len(models)-1):       
        line, = ax.plot(x, y, self.styles[iplot], marker=None, color=self.colors[iplot], label=model.name,linewidth=2.5)
      else:
        line, = ax.plot(x, y, self.styles[iplot], marker=None, color=self.colors[iplot], label=model.name)
        if line.is_dashed(): self._set_dashes(line)
                      
      iplot = iplot + 1
      
      if min(x) < xmin: xmin = min(x)
      if max(x) > xmax: xmax = max(x)
          
    ax.set_xlim(xmin,xmax)
    ax.semilogy()    
            
    ax.set_xlabel(r"r [au]")    
    ax.set_ylabel(r"midplane $\epsilon(\mathrm{" + pplot.spnToLatex(species) + "})$")    
    
    self._dokwargs(ax, **kwargs)  
    self._legend(ax,**kwargs)
    
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
        yStar= (model.starSpec.nu*model.starSpec.Inu)[0::10]
        yStar = yStar*(r**2.0*pi*dist**(-2.0))                                
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
    
  def plot_starspec_xray(self, models,**kwargs): 
    '''
    Plots the full Stellar Spectrum
    '''  
    print("PLOT: plot_xraystarspec ...")
    fig, ax = plt.subplots(1, 1)      
    
    iplot = 0    
    xmax = 0
    ymin = 1.e100
    ymax = -1.e00 
    for model in models:              
            
      idx=np.argmin(np.abs(model.starSpec.lam-0.0124))        
      x=model.starSpec.lam[0:idx:3]
      x=x*(u.micrometer).cgs        
      x=(const.h.cgs*const.c.cgs/x).to(u.keV)                              
      x=x.value
      y= (model.starSpec.Inu)[0:idx:3]
      x=x[::-1]
      y=y[::-1]     

      ax.plot(x, y, color=self.colors[iplot],linestyle=self.styles[iplot],              
              label=model.name)

                                  
#      ax.plot(x, y, color=self.colors[iplot],linestyle=self.styles[iplot],
#              marker=self.markers[iplot],ms=2,markeredgecolor=self.colors[iplot],
#              label=model.name)
                      
      iplot = iplot + 1
            
      if max(x) > xmax: xmax = max(x)
      if np.nanmin(y) < ymin: ymin=np.nanmin(y)        
      
      
      if max(y) > ymax: ymax = max(y)          
      
    xmin=0.1 # keV  
      
    # TODO: sometimes it is just zero  
    if ymin<1.e-100: ymin=1.e-20
    # set defaults, can be overwritten by the kwargs
    ax.set_xlim(xmin,xmax)
    ax.set_ylim([ymin,ymax])
    ax.semilogx()
    ax.semilogy()
    ax.set_ylabel(r"$\mathrm{I\,[erg\,cm^{-2}\,s^{-1}\,sr^{-1}\,Hz^{-1}]}$")
    ax.set_xlabel(r"Energy [keV]")                          
      
    self._dokwargs(ax, **kwargs)                
    
    self._legend(ax)
    
    self.pdf.savefig()
    plt.close(fig)   
      
  
  def plot_starspec(self, models,**kwargs): 
    '''
    Plots the full Stellar Spectrum
    '''  
    print("PLOT: plot_starspec ...")
    fig, ax = plt.subplots(1, 1)      
    
    iplot = 0
    xmin = 1.e100
    xmax = 0
    ymin = 1.e100
    ymax = -1.e00 
    for model in models:              
      
      x = model.starSpec.lam[0::10]
      y= (model.starSpec.nu*model.starSpec.Inu)[0::10]         
                                  
      ax.plot(x, y, color=self.colors[iplot],linestyle=self.styles[iplot],
              marker=self.markers[iplot],ms=2,markeredgecolor=self.colors[iplot],
              label=model.name)
                      
      iplot = iplot + 1
      
      if min(x) < xmin: xmin = min(x)
      if max(x) > xmax: xmax = max(x)
      if y[-1] < ymin: ymin = y[-1]
      
      if max(y) > ymax: ymax = max(y)          
      
    # TODO: sometimes it is just zero  
    if ymin<1.e-100: ymin=1.e-20
    # set defaults, can be overwritten by the kwargs
    ax.set_xlim(xmin,xmax)
    ax.set_ylim([ymin,ymax])
    ax.semilogx()
    ax.semilogy()
    ax.set_ylabel(r"$\mathrm{\nu F_{\nu}\,[erg\,cm^{-2}\,s^{-1}]}$")
    ax.set_xlabel(r"wavelength [$\mathrm{\mu}$m]")    
      
    self._dokwargs(ax, **kwargs)                
    
    self._legend(ax)
    
    self.pdf.savefig()
    plt.close(fig)  
   
  def plot_line_profil(self,models,wl,ident=None,**kwargs):
    '''
    Plots the line profile for the given line (id wavelength and line ident) 
    '''  
    print("PLOT: plot_line_profile ...")
    fig, ax = plt.subplots(1, 1)      
    
    iplot = 0
    xmin = 1.e100
    xmax = -1.e100 
    ymin = 1.e100
    ymax = -1.e100 
    for model in models:      
      line=model.getLine(wl,ident=ident)       
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
    
