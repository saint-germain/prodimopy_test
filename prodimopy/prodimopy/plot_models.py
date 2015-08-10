from __future__ import print_function
from __future__ import division 

import matplotlib.pyplot as plt

import matplotlib.ticker as plticker
import numpy as np
import prodimopy.plot as pplot

class PlotModels(object):
  def __init__(self,pdf,
               colors=None,
               styles=None,
               markers=None,
               fs_legend=6,
               ncol_legend=0):

    self.pdf=pdf    
    if colors == None: 
      self.colors=["b", "g","r", "c", "m", "y", "k","sienna","lime","pink","DarkOrange","Olive"]
    else:
      self.colors=colors
    
    if styles == None: 
      self.styles=["-"]*len(self.colors)
    else:
      self.styles=styles
           
    if markers == None: 
      self.markers=["D"]*len(self.colors)
    else:
      self.markers=markers    
      
    self.fs_legend=fs_legend
    self.ncol_legend=ncol_legend              

  # plots a selection of FlineEstimates
  def plot_lines(self,models,lineidents,**kwargs):
    fig, ax = plt.subplots(1, 1)  
         
    imodel=0          
    lticks=list()
    #for some reason this dummy is require
    lticks.append("")
    for model in models:                  
      iline = 1  
      x=list()
      y=list()
      for ident in lineidents:      
        line=model.getLineEstimate(ident[0],ident[1])
        x.append(iline)
        y.append(line.flux)
        
        if imodel == 0:
          #lticks.append(r"$\mathrm{"+pplot.spnToLatex(ident[0])+r"}$ "+r"{:.2f}".format(line.wl))
          lticks.append(r"$\mathrm{"+ident[0]+r"}$ "+r"{:.2f}".format(line.wl))
   
   
#         if imodel == (len(models)-1):
#           print iline, line.flux
#           ax.errorbar([iline],[line.flux],yerr=[line.flux/3.0,line.flux*3.0],fmt=".",ms=0.0,color="lightgrey",linewidth=10)    
        
        iline = iline + 1  
                                  
      ax.plot(x, y,marker=self.markers[imodel], linestyle='None', ms=5, color=self.colors[imodel],markeredgecolor=self.colors[imodel],label=model.name)      
         
      imodel=imodel+1            
           
   
    ax.set_xlim(0.5, iline -0.5)
    
    loc = plticker.MultipleLocator(base=1.0)  # this locator puts ticks at regular intervals
    ax.xaxis.set_major_locator(loc)
    
    if "ylim" in kwargs: 
      ax.set_ylim(kwargs["ylim"])
    else:         
      ax.set_ylim(1.e-24,1.e-18)
      
    ax.semilogy()    
    #ax.set_xlabel(r"line")
       
    ax.set_ylabel(r" line flux [W$\,$m$^{-2}$]")
     
    ax.set_xticklabels(lticks,rotation='45',minor=False)
    zed = [tick.label.set_fontsize(7) for tick in ax.xaxis.get_major_ticks()]
     
    handles, labels = ax.get_legend_handles_labels()
    ncol=1
    if self.ncol_legend >1 and len(labels)>self.ncol_legend:     
      ncol=int(len(labels)/self.ncol_legend) 
    ax.legend(handles, labels,loc="best", fancybox=False, framealpha=1.0,ncol=ncol,fontsize=self.fs_legend)     
    
    if "title" in kwargs and kwargs["title"]!=None:
      ax.set_title(kwargs["title"])
     
    self.pdf.savefig()
    plt.close(fig)
  
  def plot_NH(self,models,xlog=True,**kwargs):
    '''
    Plots the total vertical columndensity for the given species for all the models
    as a function of the radius
    '''
    fig, ax = plt.subplots(1, 1)  
    
    xmin=1.e100
    xmax=0
        
    iplot = 0    
    for model in models:
      x=model.x[:,0]
      y=model.NHver[:,0]
      
      #print y.min() 
      ax.plot(x, y, self.styles[iplot],marker=None,color=self.colors[iplot], label=model.name)
          
      iplot=iplot+1
      
      if min(x) < xmin: xmin=min(x)
      if max(x) > xmax: xmax=max(x)
       
  
    ax.set_xlim(xmin,xmax)
    
    if "xlim" in kwargs: ax.set_xlim(kwargs["xlim"])    
    if "ylim" in kwargs: ax.set_ylim(kwargs["ylim"])
    if xlog == True: ax.semilogx() 
         
    ax.semilogy()     
    ax.set_xlabel(r"r [AU]")
    ax.set_ylabel(r"N$_\mathrm{<H>}$ cm$^{-2}$")
  
    handles, labels = ax.get_legend_handles_labels()
    ncol=1
    if self.ncol_legend >1 and len(labels)>self.ncol_legend:     
      ncol=int(len(labels)/self.ncol_legend) 
    ax.legend(handles, labels, loc="best", fancybox=False, framealpha=1.0,ncol=ncol,fontsize=self.fs_legend)  
  
    self.pdf.savefig()
    plt.close(fig)
    
  def plot_tcdspec(self,models,species,xlog=True,title=None,**kwargs):
    '''
    Plots the total vertical columndensity for the given species for all the models
    as a function of the radius
    '''
    fig, ax = plt.subplots(1, 1)  
    
    xmin=1.e100
    xmax=0
              
    iplot = 0    
    for model in models:
      x=model.x[:,0]
      y=model.cdnmol[:,0,model.spnames[species]]
      
      ax.plot(x, y, self.styles[iplot],marker=None,color=self.colors[iplot], label=model.name)
          
      iplot=iplot+1
      
      if min(x) < xmin: xmin=min(x)
      if max(x) > xmax: xmax=max(x)
    
    
    ax.fill_between(x,y/3.0,y*3.0,color='0.8')     
  
    ax.set_xlim(xmin,xmax)
    
    if "xlim" in kwargs: ax.set_xlim(kwargs["xlim"])    
    if "xmin" in kwargs: ax.set_xlim(xmin=kwargs["xmin"])
    if "ylim" in kwargs: ax.set_ylim(kwargs["ylim"])  
    
    if xlog == True: ax.semilogx()
        
    ax.semilogy()
    
    ax.set_xlabel(r"r [AU]")
    #ax.set_ylabel(r"N$_\mathrm{"+pplot.spnToLatex(species)+"}$ cm$^{-2}$")
    ax.set_ylabel(r"N$_\mathrm{"+pplot.spnToLatex(species)+"}$ "+"[cm$^{-2}$]")
  
    if title != None: 
      ax.set_title(title)
  
    handles, labels = ax.get_legend_handles_labels()
    
    ncol=1
    if self.ncol_legend >1 and len(labels)>self.ncol_legend:     
      ncol=int(len(labels)/self.ncol_legend) 
    ax.legend(handles, labels, loc="best", fancybox=False, framealpha=1.0,ncol=ncol,fontsize=self.fs_legend)  
  
    self.pdf.savefig()
    plt.close(fig)  
    
  def plot_tauline(self,models,lineIdent,xlog=True,**kwargs):
    '''
    Plots the line optical depth as a function of radius for a given line
    for all the models
    '''
    fig, ax = plt.subplots(1, 1)  
    
    xmin=1.e100
    xmax=0
        
    iplot = 0    
    for model in models:
      x=model.x[:,0]
      lineEstimate=model.getLineEstimate(lineIdent[0],lineIdent[1])
      y=list()
      ytauDust=list()
      for rInfo in lineEstimate.rInfo:
        y.append(rInfo.tauLine)
        ytauDust.append(rInfo.tauDust)
      
      ax.axhline(y=1.0,linestyle="-",color="black",linewidth=0.5)
      ax.plot(x, y, self.styles[iplot],marker=None,color=self.colors[iplot], label=model.name)
      #ax.plot(x,ytauDust,"--",marker=None,color="black")
      
          
      iplot=iplot+1
      
      if min(x) < xmin: xmin=min(x)
      if max(x) > xmax: xmax=max(x)
       
  
    ax.set_xlim(xmin,xmax)
    
    if "xlim" in kwargs: ax.set_xlim(kwargs["xlim"])    
    if "xmin" in kwargs: ax.set_xlim(xmin=kwargs["xmin"])
    if "ylim" in kwargs: ax.set_ylim(kwargs["ylim"])  
    
    if xlog == True: ax.semilogx()
        
    ax.semilogy()
    
    ax.set_xlabel(r"r [AU]")    
    ax.set_ylabel(r"$\mathrm{\tau_{line}}$")
    ax.set_title(lineEstimate.ident+" "+"{:.2f}".format(lineEstimate.wl)+" $\mathrm{\mu m}$")
  
    handles, labels = ax.get_legend_handles_labels()
    
    ncol=1
    if self.ncol_legend >1 and len(labels)>self.ncol_legend:     
      ncol=int(len(labels)/self.ncol_legend) 
    ax.legend(handles, labels, loc="best", fancybox=False, framealpha=1.0,ncol=ncol,fontsize=self.fs_legend)  
  
    self.pdf.savefig()
    plt.close(fig)  
    
    
  def plot_abunvert(self,models, r, species, **kwargs):
          
    rstr = "r={:.2f} AU".format(r)   
    
    fig, ax = plt.subplots(1, 1)     

    iplot = 0
    xmin=1.e100
    xmax=0 
    for model in models:
      # TODO: this just takes the closed value to the given r
      #       that might be different from model to model
      ix=(np.abs(model.x[:, 0] - r)).argmin()
      
      old_settings = np.seterr(divide='ignore')     
      x = np.log10(model.NHver[ix, :])      
      np.seterr(**old_settings)  # reset to default

      y =model.nmol[ix,:,model.spnames[species]]/model.nHtot[ix,:]
      
      ax.plot(x, y, self.styles[iplot],marker=None,color=self.colors[iplot], label=model.name)
                      
      iplot=iplot+1
      
      if min(x) < xmin: xmin=min(x)
      if max(x) > xmax: xmax=max(x)

      
    # set the limits
      
    if "xlim" in kwargs:     
      ax.set_xlim(kwargs["xlim"])
    else:
      ax.set_xlim([17.5,x.max()])
              
    if "ylim" in kwargs: ax.set_ylim(kwargs["ylim"])
     
    #print ax.get_xlim()
  
#     ax2 = ax.twiny()
#     ax2.set_xlabel("z/r")
#     ax2.set_xlim(ax.get_xlim())
#     #ax2.set_xticks(ax.get_xticks())    
#     ax2.set_xticklabels(["{:.2f}".format(x) for x in nhver_to_zr(ix, ax.get_xticks(), model)])
    
    ax.set_xlabel(r"$\log$ N$_\mathrm{H}$ [cm$^{-2}$]")
    ax.set_ylabel(r"$\epsilon(\mathrm{"+pplot.spnToLatex(species)+"})$")
    ax.set_title(rstr)
    
    # do axis style
    ax.semilogy()     
    
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="best", fancybox=True, framealpha=0.5)
    #ax.text(0.025, 0.025, rstr,
    #   verticalalignment='bottom', horizontalalignment='left',
    #   transform=ax.transAxes,alpha=0.75)
      
    self.pdf.savefig()
    plt.close(fig) 
  
  
  def plot_midplane(self,models,fieldname,ylabel,**kwargs):
    fig, ax = plt.subplots(1, 1)      
    
    iplot = 0
    xmin=1.e100
    xmax=0
    ymin=1.e100
    ymax=-1.e00 
    for model in models:           
      x = model.x[:,0]
      y =  getattr(model,fieldname)[:,0]                    
      
      ax.plot(x, y, self.styles[iplot],marker=None,color=self.colors[iplot], label=model.name)
                      
      iplot=iplot+1
      
      if min(x) < xmin: xmin=min(x)
      if max(x) > xmax: xmax=max(x)
      if min(y) < ymin: ymin=min(y)
      if max(y) > ymax: ymax=max(y)

      
    if "ylim" in kwargs: 
      ax.set_ylim(kwargs["ylim"])
    else:
      ax.set_ylim(ymin,ymax)
              
    ax.semilogy()
            
    ax.set_xlabel(r"r [AU]")    
    ax.set_ylabel(ylabel)    
    
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels,loc="best", fancybox=False, framealpha=0.8,fontsize=self.fs_legend) 
    
    self.pdf.savefig()
    plt.close(fig)    
    
    
  def plot_abunradial(self,models,species,**kwargs):
  
    fig, ax = plt.subplots(1, 1)      
    
    iplot = 0
    xmin=1.e100
    xmax=0 
    for model in models:           
      x = model.x[:,0]            
      y =model.nmol[:,0,model.spnames[species]]/model.nHtot[:,0]
      
      ax.plot(x, y, self.styles[iplot],marker=None,color=self.colors[iplot], label=model.name)
                      
      iplot=iplot+1
      
      if min(x) < xmin: xmin=min(x)
      if max(x) > xmax: xmax=max(x)
      
    if "ylim" in kwargs: ax.set_ylim(kwargs["ylim"])
    
    ax.semilogy()
        
    
    ax.set_xlabel(r"r [AU]")    
    ax.set_ylabel(r"$\epsilon(\mathrm{"+pplot.spnToLatex(species)+"})$")    
    
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels,loc="best", fancybox=False, framealpha=0.8,fontsize=self.fs_legend) 
    
    self.pdf.savefig()
    plt.close(fig)    
    
    