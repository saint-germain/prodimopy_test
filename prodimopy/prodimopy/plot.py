from __future__ import division 
from __future__ import print_function

from scipy.interpolate import interp1d

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def spnToLatex(spname):      
  name = spname        
  if spname == "HN2+": name = "N2H+" 
  if spname == "C18O": return "C^{18}O"
  if spname == "13CO": return "^13CO" 
  
  newname = ""
  for c in name:    
    if c.isdigit():
      newname += "_" + c
    elif c == "-":
      newname += "^-"
    elif c == "+":
      newname += "^+"
    elif c == "#":
      newname += "\#"       
    else:
      newname += c
      
  return newname

def nhver_to_zr(ir, nhver, model, log=True):
  zrs = model.z[ir, :] / model.x[ir, :]
    
  if log == True:   
    old_settings = np.seterr(divide='ignore')
    ipol = interp1d(np.log10(model.NHver[ir, :]), zrs, bounds_error=False, fill_value=0.0, kind="linear")
    np.seterr(**old_settings)  # reset to default
  else:
    ipol = interp1d(model.NHver[ir, :], zrs, bounds_error=False, fill_value=0.0, kind="linear")
  
  #return 0
  return ipol(nhver)  

def plog(array):      
  # ignore divide by zero in log10
  old_settings = np.seterr(divide='ignore') 
  array = np.log10(array)
  np.seterr(**old_settings)  # reset to default  
  return array

class Plot(object):
  def __init__(self, pdf, fs_legend=7.5):
    self.pdf = pdf
    self.fs_legend = fs_legend
    self.ncol_legend = 5
  
  def _legend(self, ax,loc="best"):
    '''
    plots the legend, deals with multiple columns
    '''
    handles, labels = ax.get_legend_handles_labels()
    ncol = 1
    if self.ncol_legend > 1 and len(labels) > self.ncol_legend:
      ncol = int(len(labels) / self.ncol_legend)
    
    ax.legend(handles, labels, loc=loc, fancybox=False, ncol=ncol, fontsize=self.fs_legend)
  
  
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
  
  def plot_cont(self, model, values, label="value", zlog=True, 
                zlim=[None, None],zr=True,clevels=None,contour=True,
                extend="neither",**kwargs):
    '''
    plot routine for 2D contour plots.   
    '''
  
    if zlog is True:
      pvals = plog(values)
      values[np.isnan(values)] = 0.0
      
      # TODO: that looks very complicated
      if zlim[1] == None: 
        maxval = np.log10(values.max())
      else:
        maxval = np.log10(zlim[1])
      if zlim[0] == None:    
        minval = np.log10(values[values > 0.0].min())
      else:
        minval = np.log10(zlim[0])                
    else:
      pvals = values
      if zlim[1] == None:    
        maxval = values.max()
      else:
        maxval = zlim[1]
      if zlim[0] == None:    
        minval = values.min()
      else:
        minval = zlim[0]              
      
    x = model.x  
    if zr:
      y = model.z / model.x
    else:
      y = np.copy(model.z) 
      y[:,0]=y[:,0]+0.05 
  
    levels = MaxNLocator(nbins=100).tick_values(maxval, minval)
    ticks = MaxNLocator(nbins=6, prune="both").tick_values(minval, maxval)
  
    fig, ax = plt.subplots(1, 1)   
    cmap = plt.get_cmap('jet')        
    CS = ax.contourf(x, y, pvals, levels=levels, cmap=cmap,extend=extend)
    # This is the fix for the white lines between contour levels
    for c in CS.collections:
      c.set_edgecolor("face")    

    ax.set_ylim([y.min(), y.max()])
    ax.set_xlim([x.min(), x.max()])
    ax.semilogx()        
    
    ax.set_xlabel("r [AU]")
    if zr:
      ax.set_ylabel("z/r")
    else:
      ax.set_ylabel("z [AU]")      
  
  
      #ax.text(0.27, 0.95,kwargs["title"], horizontalalignment='center',
      #   verticalalignment='center',fontsize=8,
      #   transform=ax.transAxes)
      
    self._dokwargs(ax,**kwargs)            
         
    if contour:
      if clevels != None:
        if zlog: clevels=np.log10(clevels)
        ticks=clevels
        ax.contour(CS, levels=clevels, colors='black', linestyles="solid",linewidths=1.0)
      else:
        ax.contour(CS, levels=ticks, colors='black', linestyles="dashed",linewidths=1.0)
    
    CB = fig.colorbar(CS, ax=ax,ticks=ticks,pad=0.01)
    CB.ax.tick_params(labelsize=self.fs_legend) 
    # CB.set_ticks(ticks)
    CB.set_label(label,fontsize=self.fs_legend)  

    
    
    self.pdf.savefig(transparent=mpl.rcParams['savefig.transparent'])
    plt.close(fig)  
  
  def plot_ionrates_midplane(self, model, **kwargs):                       
    
    x = np.log10(model.x[:,0])      

  #  print pdata.zetaCR[ix,:]
    y1 = model.zetaCR[:, 0]
    y2 = model.zetaX[:, 0] / 2.0  # convert to per H2 TODO: maybe do this in ProDiMo already to be consistent
    y3 = model.zetaSTCR[:, 0]  # convert to per H2 TODO: maybe do this in ProDiMo already to be consistent
  
    fig, ax = plt.subplots(1, 1)   
    ax.plot(x, y1, color="red", label="$\zeta_\mathrm{CR}$")
    ax.plot(x, y2, color="blue", label="$\zeta_\mathrm{X}$")
    ax.plot(x, y3, color="green", label="$\zeta_\mathrm{SP}$")
      
    # set the limits
      
    if "xlim" in kwargs:     
      ax.set_xlim(kwargs["xlim"])
    else:
      ax.set_xlim([x.min(),x.max()])        
    if "ylim" in kwargs: ax.set_ylim(kwargs["ylim"])
     
    # print ax.get_xlim()
    
    ax.set_xlabel(r"r [AU]")
    ax.set_ylabel("ionization rate per H$_2$ [s$^{-1}$]")
    
    # do axis style
    if "xlog" in kwargs:
      if kwargs["xlog"]:
        ax.semilogx()
    
    ax.semilogy()     
    
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="best", fancybox=True)    
      
    self.pdf.savefig(transparent=False)
    plt.close(fig)  

  
    
  def plot_ionrates(self, model, r, **kwargs):
              
    ix = (np.abs(model.x[:, 0] - r)).argmin()
    rstr = "r={:.2f} AU".format(model.x[ix, 0])   
         
    old_settings = np.seterr(divide='ignore')     
    nhver = np.log10(model.NHver[ix, :])      
    np.seterr(**old_settings)  # reset to default
  #  print pdata.zetaCR[ix,:]
    y1 = model.zetaCR[ix, :]
    y2 = model.zetaX[ix, :] / 2.0  # convert to per H2 TODO: maybe do this in ProDiMo already to be consistent
    y3 = model.zetaSTCR[ix, :]  # convert to per H2 TODO: maybe do this in ProDiMo already to be consistent
  
    fig, ax = plt.subplots(1, 1)   
    ax.plot(nhver, y1, color="red", label="$\zeta_\mathrm{CR}$")
    ax.plot(nhver, y2, color="blue", label="$\zeta_\mathrm{X}$")
    ax.plot(nhver, y3, color="green", label="$\zeta_\mathrm{SP}$")
      
    # set the limits
      
    if "xlim" in kwargs:     
      ax.set_xlim(kwargs["xlim"])
    else:
      ax.set_xlim([17.5, nhver.max()])        
    if "ylim" in kwargs: ax.set_ylim(kwargs["ylim"])
     
    # print ax.get_xlim()
  
    ax2 = ax.twiny()
    ax2.set_xlabel("z/r")
    ax2.set_xlim(ax.get_xlim())
    # ax2.set_xticks(ax.get_xticks())    
    ax2.set_xticklabels(["{:.2f}".format(x) for x in nhver_to_zr(ix, ax.get_xticks(), model)])
    
    ax.set_xlabel(r"$\log$ N$_\mathrm{H,ver}$ [cm$^{-2}$]")
    ax.set_ylabel("ionization rate per H$_2$ [s$^{-1}$]")
    
    # do axis style
    ax.semilogy()     
    
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="best", fancybox=True, framealpha=0.5)
    ax.text(0.025, 0.025, rstr,
       verticalalignment='bottom', horizontalalignment='left',
       transform=ax.transAxes, alpha=0.75)
      
    self.pdf.savefig()
    plt.close(fig)  

  def plot_avgabun(self,model,species,**kwargs):
    '''
    Plots the average abundance for the given species (can be more than one) 
    as a function of radius    
    '''
    print("PLOT: plot_avgabun ...")     
          
    fig, ax = plt.subplots(1, 1)     
  
    iplot = 0
    for spec in species:                 
    # get the species
      if (spec in model.spnames):   
        y=model.cdnmol[:,0,model.spnames[spec]]
        y=y/model.NHver[:,0]   
        x = model.x[:,0]            
        
        style="-"
        if "#" in spec: style="--"
        ax.plot(x, y, linestyle=style, marker=None, label="$\mathrm{"+spnToLatex(spec)+"}$")
    
        iplot = iplot + 1              
    
    if iplot==0:
      print("Species "+species+ " not found in any model!")
      return 
            
    ax.set_xlabel(r"r [AU]")
    ax.set_ylabel(r"average abundance")
    ax.set_xlim([x.min(),x.max()])
    
    # do axis style
    ax.semilogy()     
    
    self._dokwargs(ax,**kwargs)    
    self._legend(ax)
    
    self.pdf.savefig()
    plt.close(fig)

  def plot_midplane(self, model, fieldname, ylabel, **kwargs):
    '''
    Plots a quantitiy in in the midplane as a function of radius
    fieldname is any field in Data_ProDiMo
    '''
    print("PLOT: plot_midplane ...")
    fig, ax = plt.subplots(1, 1)      
    
    
    x = model.x[:, 0]
    y = getattr(model, fieldname)[:, 0]                    
      
    ax.plot(x, y,marker=None)
                 
    ax.set_xlim(np.min(x),np.max(x))
    ax.set_ylim(np.min(y),np.max(y))                                 
    ax.semilogy()
            
    ax.set_xlabel(r"r [AU]")    
    ax.set_ylabel(ylabel)    
    
    self._dokwargs(ax, **kwargs) 
    #self._legend(ax)
    
    self.pdf.savefig()
    plt.close(fig)     

  def plot_abun_midp(self, model,species, norm=None,styles=None,colors=None, **kwargs):
    '''
    Plots the abundances in the midplane for the given species (can be more than one)
    norm .. normalisation constant (additional to Nhto e.g. the total carbon abundance)   
    '''
    print("PLOT: plot_abun_midp ...")
    fig, ax = plt.subplots(1, 1)      
    
    iplot = 0
    xmin = 1.e100
    xmax = 0
    ymin = 1.e100
    ymax = -1.e00 
    for spec in species:           
      x = model.x[:, 0]
      y = model.nmol[:,0,model.spnames[spec]]/model.nHtot[:,0]
      if norm != None:
        y=y/norm                    
      
      # FIXME: add proper treatment for styles and colors 
      
      if styles==None:         
        style="-"
        if "#" in spec: style="--"
      else: 
        style=styles[iplot]
        
      if colors==None:
        color=None
      else:
        color=colors[iplot]  
        
      ax.plot(x, y, marker=None, linestyle=style, color=color, label="$\mathrm{"+spnToLatex(spec)+"}$")
                      
      iplot = iplot + 1
      
      if min(x) < xmin: xmin = min(x)
      if max(x) > xmax: xmax = max(x)
      if min(y) < ymin: ymin = min(y)
      if max(y) > ymax: ymax = max(y)
  
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin, ymax)              
    ax.semilogy()
    ax.set_xlabel("r [AU]")
    ax.set_ylabel(r"$\mathrm{\epsilon(X)}$")
            
    self._dokwargs(ax, **kwargs)
    self._legend(ax, loc="lower left") 
  
    self.pdf.savefig()
    plt.close(fig)

