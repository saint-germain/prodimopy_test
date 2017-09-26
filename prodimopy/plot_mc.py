"""
.. module:: plot_mc 
   :synopsis: Plotting routines for molecular cloud |prodimo| models (0D chemistry).

.. moduleauthor:: Ch. Rab
"""

from __future__ import division 
from __future__ import print_function

import matplotlib.pyplot as plt
# That crashes in the bash on Mac OS X.. maybe because of installation
import prodimopy.plot as pplt 


class PlotMcModels(object):
  """
  Plot routines for a molecular cloud |prodimo| models (0D chemistry).
  """
  def __init__(self, pdf,
               colors=None,
               styles=None,
               markers=None,
               fs_legend=6,   # TODO: init it with the system default legend fontsize
               ncol_legend=0):
    """
    Parameters
    ----------
    name : pdf
      a :class:`matplotlib.backends.backend_pdf.PdfPages` object used to save the plot.
      
    colors : list
      a list of matplotlib colors used for different models. (default: None) 

    styles : list
      a list of matplotlib styles used for different models. (default: None) 

    markers : list
      a list of matplotlib markers used for different models. (default: None)
      
    Attributes
    ----------
    """       
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
      self.markers = [None] * len(self.colors)
    else:
      self.markers = markers    
       
    self.fs_legend = fs_legend
    self.ncol_legend = ncol_legend 
     
  def _dokwargs(self,ax,**kwargs):
    """    
    Handles the passed kwargs elements (assumes that defaults are already set)
    """
    if "ylim" in kwargs: 
      ax.set_ylim(kwargs["ylim"])
       
    if "xlim" in kwargs: 
      ax.set_xlim(kwargs["xlim"])
               
    if "xlog" in kwargs:
      if kwargs["xlog"]: ax.semilogx()
       
    if "ylog" in kwargs:
      if kwargs["ylog"]: ax.semilogy()  
       
  def _legend(self, ax):
    '''
    plots the legend, deals with multiple columns
    '''
    handles, labels = ax.get_legend_handles_labels()
    ncol = 1
    if self.ncol_legend > 1 and len(labels) > self.ncol_legend:
      ncol = int(len(labels) / self.ncol_legend)
    ax.legend(handles, labels, loc="best", fancybox=False, framealpha=0.75, ncol=ncol, fontsize=self.fs_legend)
     
  def plot_species(self,models,spname,**kwargs):
    fig, ax = plt.subplots(1, 1)  
     
    i = 0
    xmax=1.e-99
    xmin=1.e99
    for model in models:     
      if (spname in model.species):
        ax.plot(model.ages,model.abundances[:,model.species.index(spname)],
              self.styles[i],
              marker=self.markers[i],
              color=self.colors[i],label=model.name)
        xmax=max([max(model.ages),xmax])
        # exclude zero because of log plot        
        #nozeros=   
        xmin=min([min(model.ages[model.ages>0.0]),xmin])
 
        i = i+1            
    
    # no data to plot just return 
    if i==0:
      print("Species "+spname+" not found.")
      return
     
    ax.set_xlim(xmin, xmax)
    ax.semilogx()
    ax.semilogy()
    self._dokwargs(ax,**kwargs)
      
     
    ax.set_xlabel(r"years")         

    ax.set_ylabel(r" $\mathrm{\epsilon("+pplt.spnToLatex(spname)+")}$")
           
    self._legend(ax)
      
    #if "title" in kwargs and kwargs["title"] != None:
    #  ax.set_title(kwargs["title"])
       
    self.pdf.savefig()
    plt.close(fig) 
    
 
  def plot_abunratio(self,models,spname1,spname2,**kwargs):
    fig, ax = plt.subplots(1, 1)  
     
    i = 0
    xmax=1.e-99
    xmin=1.e99
    for model in models:     
      if (spname1 in model.species and spname2 in model.species):
        ratio=model.abundances[:,model.species.index(spname1)]/model.abundances[:,model.species.index(spname2)]
        ax.plot(model.ages,ratio,
              self.styles[i],
              marker=self.markers[i],
              color=self.colors[i],label=model.name)
        xmax=max([max(model.ages),xmax])
        # exclude zero because of log plot        
        #nozeros=   
        xmin=min([min(model.ages[model.ages>0.0]),xmin])
 
        i = i+1            
    
    if i==0:
      print("Species "+spname1+"and/or "+spname2+" not found.")
      return

    
     
    ax.set_xlim(xmin, xmax)
    ax.semilogx()
    ax.semilogy()
    self._dokwargs(ax,**kwargs)
      
     
    ax.set_xlabel(r"years")         

    ax.set_ylabel(r" $\mathrm{\epsilon("+pplt.spnToLatex(spname1)+")/"
                  +"\epsilon("+pplt.spnToLatex(spname2)+")}$")
           
    self._legend(ax)
      
    #if "title" in kwargs and kwargs["title"] != None:
    #  ax.set_title(kwargs["title"])
       
    self.pdf.savefig()
    plt.close(fig) 
    
    