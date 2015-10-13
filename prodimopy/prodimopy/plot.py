from __future__ import division 
from __future__ import print_function

from scipy.interpolate import interp1d

import numpy as np
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
  def __init__(self, pdf, fs_legend=6):
    self.pdf = pdf
    self.fs_legend = fs_legend
  
  
  def plot_cont(self, model, values, label="value", zlog=True, 
                zlim=[None, None],zr=True,**kwargs):
  
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
    CS = ax.contourf(x, y, pvals, levels=levels, cmap=cmap)
    # This is the fix for the white lines between contour levels
    for c in CS.collections:
      c.set_edgecolor("face")
    

    if "xlim" in kwargs:
      ax.set_xlim(kwargs["xlim"])      
    else:    
      ax.set_xlim([x.min(), x.max()])
    if "ylim" in kwargs:  
      ax.set_ylim(kwargs["ylim"])    
    else:  
      ax.set_ylim([y.min(), y.max()])

    
    if not "xlog" in kwargs:
      ax.semilogx()
    else:
      if kwargs["xlog"]: ax.semilogx()

    if "ylog" in kwargs and kwargs["ylog"]:
      ax.semilogy()
    
    
    ax.set_xlabel("r [AU]")
    if zr:
      ax.set_ylabel("z/r")
    else:
      ax.set_ylabel("z [AU]")      
    
    CB = fig.colorbar(CS, ticks=ticks)
    # CB.set_ticks(ticks)
    CB.set_label(label)  
    
    ax.contour(CS, levels=ticks, colors='black', linestyles="dashed")
    
    self.pdf.savefig()
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
    ax.legend(handles, labels, loc="best", fancybox=True, framealpha=0.5)    
      
    self.pdf.savefig()
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
