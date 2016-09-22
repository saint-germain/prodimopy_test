from __future__ import division 
from __future__ import print_function

from scipy.interpolate import interp1d

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import astropy.units as u
import math


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
      
  # repair the double ionized case
  if "^+^+" in newname: newname=newname.replace("^+^+", "^{++}")
      
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
  def __init__(self, pdf, fs_legend=None,title=None):
    self.pdf = pdf
    if fs_legend is None:
      self.fs_legend = mpl.rcParams['legend.fontsize']
    else:
      self.fs_legend=fs_legend
    self.ncol_legend = 5
    self.title=title
    
    # special colors, forgot the source for it :( somewhere from the internet)
    # FIXME: make an option to aktivate them
    self.pcolors={"blue"   : "#5DA5DA",
                  "orange" : "#FAA43A",
                  "green"  : "#60BD68",
                  "pink"   : "#F17CB0",
                  "brown"  : "#B2912F",
                  "purple" : "#B276B2",
                  "yellow" : "#DECF3F",
                  "red"    : "#F15854",
                  "gray"   : "#4D4D4D"}    
  
  def _legend(self, ax,loc="best"):
    '''
    plots the legend, deals with multiple columns
    '''
    handles, labels = ax.get_legend_handles_labels()    
    
    if len(labels)>0:
      ncol = 1
      if self.ncol_legend > 1 and len(labels) > self.ncol_legend:
        ncol = int(len(labels) / self.ncol_legend)
            
      leg=ax.legend(handles, labels, loc=loc, fancybox=False, ncol=ncol, fontsize=self.fs_legend) 
      lw=mpl.rcParams['axes.linewidth']
      leg.get_frame().set_linewidth(lw)   
  
  def _dokwargs(self,ax,**kwargs):
    '''
    Handles the passed kwargs elements (assumes that defaults are already set)
    TODO: make this a general function .... 
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
    
    if self.title != None:
      if self.title.strip() != "":
        ax.set_title(self.title.strip())
      
    if "title" in kwargs:
      if  kwargs["title"] != None and kwargs["title"].strip() != "":
        ax.set_title(kwargs["title"].strip())
      else:
        ax.set_title("")
        
  def _closefig(self,fig):
    '''
    save and close the plot
    
    set the transparent attribut (used rcParam savefig.transparent)
    '''    
    
    #trans=mpl.rcParams['savefig.transparent']
    
    self.pdf.savefig(figure=fig,transparent=False)
    plt.close(fig)
    
  def _sfigs(self,**kwargs):
    '''
    Scale the figure size from matplotlibrc by the factors given in the 
    array sfigs (in kwargs) the first element is for the width the second for
    the hieght
    '''            
    figsize=mpl.rcParams['figure.figsize']
    
    if "sfigs" in kwargs:
      fac=kwargs["sfigs"]               
      return (figsize[0]*fac[0],figsize[1]*fac[1])
    else:
      return (figsize[0],figsize[1])
  
  def plot_NH(self, model, **kwargs):
    '''
    Plots the total vertical hydrogen column number density 
    as a function of radius
    '''
    print("PLOT: plot_NH ...")
    fig, ax = plt.subplots(1, 1)      
    
    x = model.x[:, 0]
    y = model.NHver[:, 0]    
    ax.plot(x, y, marker=None, color="black")
  
    ax.set_xlim(min(x), max(x))
         
    ax.semilogy()  
    ax.semilogx()   
    ax.set_xlabel(r"r [au]")
    ax.set_ylabel(r"N$_\mathrm{<H>}$ cm$^{-2}$")
  
    self._dokwargs(ax,**kwargs)
    self._legend(ax)
    self._closefig(fig)
 
 
  # FIXME: this is not really a general plot .... 
  def plot_cont_dion(self, model, zr=True,oconts=None,**kwargs):
    '''
    plot routine for 2D contour plots.
    plots the regions where either X-rays, CR or SP are the dominant H2 ionization source   
    '''
    
    values=model.zetaX[:,:]*0.0
    values[:,:]=np.nan
    values[model.zetaX*2.0>(model.zetaCR+model.zetaSTCR)]=1.0
    values[model.zetaSTCR>(model.zetaCR+model.zetaX*2.0)]=0.0
    values[model.zetaCR>(model.zetaSTCR+model.zetaX*2.0)]=-1.0
    #print(values)
    
    print("PLOT: plot_cont_dion ...")
    cX="#F15854"
    cSP="#5DA5DA"
    cCR="#4D4D4D"
        
      
    x = model.x  
    if zr:
      y = model.z / model.x
    else:
      y = np.copy(model.z) 
      y[:,0]=y[:,0]+0.05 
  
    #levels=[1.5,0.5,0.0,-0.5,-1.5]    
    #levels=MaxNLocator(nbins=4, prune="both").tick_values(-1.0,1.0)
    levels=[-1.2, -0.01,0.0, 0.01,1.2]
    ticks=[0.5,0.0,-0.5]
    
    #ticks = 
    #print(ticks)
  
    # sclae the figure size if necessary
    # TODO: maybe provide a routine for this, including scaling the figure size
    fig, ax = plt.subplots(1, 1,figsize=self._sfigs(**kwargs))                        
    
    # stupid trick to plot the masked areas
    # plot everything with one color, and than overplot the other stuff. 
    valinv=model.zetaX[:,:]*0.0
    valinv[:,:]=10
    #valinv[valinv != 10.0]=np.nan     
    CS2 = ax.contourf(x, y, valinv,levels=[9.0,10.0,11.0], colors="lightgray",hatches=["////","////","////"])
    for c in CS2.collections:
      c.set_edgecolor("face")   
    
    CS = ax.contourf(x, y, values,levels=levels,colors=(cCR,cSP,cSP,cX))            
    # This is the fix for the white lines between contour levels
    for c in CS.collections:
      c.set_edgecolor("face")   
        

    ax.set_ylim([y.min(), y.max()])
    ax.set_xlim([x.min(), x.max()])
    ax.semilogx()        
    
    ax.set_xlabel("r [au]")
    if zr:
      ax.set_ylabel("z/r")
    else:
      ax.set_ylabel("z [au]")      
         
    self._dokwargs(ax,**kwargs)            
             
    if oconts is not None:
      for cont in oconts:
        ACS=ax.contour(x, y,cont.field,levels=cont.levels, 
                       colors=cont.colors,linestyles=cont.linestyles,linewidths=cont.linewidths)
    
    CB = fig.colorbar(CS, ax=ax,ticks=ticks,pad=0.01)
    CB.ax.set_yticklabels(['X', 'SP', 'CR'])
    CB.ax.tick_params(labelsize=self.fs_legend) 
    
    
    # CB.set_ticks(ticks)
    CB.set_label("dominant ionization source",fontsize=self.fs_legend)  

    self._closefig(fig)
    
  
  def plot_cont(self, model, values, label="value", zlog=True, 
                zlim=[None, None],zr=True,clevels=None,clabels=None,contour=True,
                extend="neither",oconts=None,acont=None,acontl=None,nbins=100,
                bgcolor=None,**kwargs):
    '''
    plot routine for 2D contour plots.
    oconts needs to be an array of Countour objectes
    FIXME: acont is deprectated use oconts only 
       
    '''
    print("PLOT: plot_cont ...")
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
  
    levels = MaxNLocator(nbins=nbins).tick_values(maxval, minval)
        
    if clevels is not None:
      if zlog: clevels=np.log10(clevels)    
      ticks=clevels
    else: 
      ticks = MaxNLocator(nbins=6, prune="both").tick_values(minval, maxval)
  
    # TODO: maybe provide a routine for this, including scaling the figure size
    fig, ax = plt.subplots(1, 1,figsize=self._sfigs(**kwargs)) 

    CS = ax.contourf(x, y, pvals, levels=levels, extend=extend)
    # This is the fix for the white lines between contour levels
    for c in CS.collections:
      c.set_edgecolor("face")    

    ax.set_ylim([y.min(), y.max()])
    ax.set_xlim([x.min(), x.max()])
    ax.semilogx()        
    
    ax.set_xlabel("r [au]")
    if zr:
      ax.set_ylabel("z/r")
    else:
      ax.set_ylabel("z [au]")      
  
  
      #ax.text(0.27, 0.95,kwargs["title"], horizontalalignment='center',
      #   verticalalignment='center',fontsize=8,
      #   transform=ax.transAxes)
      
    self._dokwargs(ax,**kwargs)            
         
    if contour:
      if clevels is not None:
        #if zlog: clevels=np.log10(clevels)
        #ticks=clevels
        ax.contour(CS, levels=clevels, colors='white', linestyles="--",linewidths=1.0)
      else:
        ax.contour(CS, levels=ticks, colors='white', linestyles="--",linewidths=1.0)
    
    if oconts is not None:
      for cont in oconts:
        ACS=ax.contour(x, y,cont.field,levels=cont.levels, 
                       colors=cont.colors,linestyles=cont.linestyles,linewidths=cont.linewidths)
    
    if acont is not None:      
      print("WARN: plot_cont: please use the oconts for additional contours ...")      
      #for l in acontl:
      #  ACS=ax.contour(x, y,pvals,levels=[l], colors='black',linestyles="solid",linewidths=1.5)
      #  ax.clabel(ACS, inline=1, fontsize=8,fmt=str(l))
      ACS=ax.contour(x, y,acont,levels=acontl, colors='white',linestyles="solid",linewidths=1.5)
      # quick fix for second contour ... 
      #ACS2=ax.contour(x, y,model.nHtot,levels=[1.e6], colors='black',linestyles="solid",linewidths=2.5)
      #ax.clabel(ACS, inline=1, fontsize=8,fmt="%.0f")
      
    if bgcolor is not None:
      ax.set_axis_bgcolor(bgcolor)
    
    CB = fig.colorbar(CS, ax=ax,ticks=ticks,pad=0.01,format="%.1f")
    # FIXME: this is not very flexible and confusing
    if clabels is not None:
      CB.ax.set_yticklabels(clabels)
    #CB.ax.tick_params(labelsize=self.fs_legend) 
    # CB.set_ticks(ticks)
    CB.set_label(label)  

    self._closefig(fig)
  
  def plot_ionrates_midplane(self, model, **kwargs):                       
    
    print("PLOT: plot_ionrates_midplane ...") 
    
    cX=self.pcolors["red"]
    cSP=self.pcolors["blue"]
    cCR=self.pcolors["gray"]  
    
    x = model.x[:,0]      

  #  print pdata.zetaCR[ix,:]
    y1 = model.zetaCR[:, 0]
    y2 = model.zetaX[:, 0] * 2.0  # convert to per H2 TODO: maybe do this in ProDiMo already to be consistent
    y3 = model.zetaSTCR[:, 0]  # convert to per H2 TODO: maybe do this in ProDiMo already to be consistent
  
    fig, ax = plt.subplots(1, 1,figsize=self._sfigs(**kwargs))   
    ax.plot(x, y2, color=cX, label="$\zeta_\mathrm{X}$")
    ax.plot(x, y3, color=cSP, label="$\zeta_\mathrm{SP}$")
    ax.plot(x, y1, color=cCR, label="$\zeta_\mathrm{CR}$")
               
    # print ax.get_xlim()
    
    ax.set_xlabel(r"r [au]")
    ax.set_ylabel("$\mathrm{\zeta\,per\,H_2\,[s^{-1}]}$")
        
    ax.semilogy()
    self._dokwargs(ax,**kwargs)         
    self._legend(ax)
    #handles, labels = ax.get_legend_handles_labels()
    #ax.legend(handles, labels, loc="best", fancybox=False)    
      
    self.pdf.savefig(transparent=False)
    plt.close(fig)  
  
  
  #FIXME: this routine is also not very general (e.g. colors)  
  def plot_ionrates(self, model, r, **kwargs):
              
    cX=self.pcolors["red"]
    cSP=self.pcolors["blue"]
    cCR=self.pcolors["gray"]          
              
    ix = (np.abs(model.x[:, 0] - r)).argmin()
    rstr = "r={:.1f} au".format(model.x[ix, 0])   
         
    old_settings = np.seterr(divide='ignore')     
    nhver = np.log10(model.NHver[ix, :])      
    np.seterr(**old_settings)  # reset to default
  #  print pdata.zetaCR[ix,:]
    y1 = model.zetaCR[ix, :]
    y2 = model.zetaX[ix, :] * 2.0  # convert to per H2 TODO: maybe do this in ProDiMo already to be consistent
    y3 = model.zetaSTCR[ix, :]  
      
    fig, ax = plt.subplots(1, 1,figsize=self._sfigs(**kwargs))   
    ax.plot(nhver, y2, color=cX, label="$\zeta_\mathrm{X}$")
    ax.plot(nhver, y3, color=cSP, label="$\zeta_\mathrm{SP}$")
    ax.plot(nhver, y1, color=cCR, label="$\zeta_\mathrm{CR}$")
      
    # set the limits
      
    ax.set_xlim([17.5, nhver.max()])  
    ax.set_ylim([1.e-21,1.e-9])

    ax.set_xlabel(r"$\log$ N$_\mathrm{<H>,ver}$ [cm$^{-2}$]")
    ax.set_ylabel("$\zeta$ per H$_2$ [s$^{-1}$]")

    # do axis style
    ax.semilogy()     

    # title does not work here
    self._dokwargs(ax,title=None,**kwargs)
      
    ax2 = ax.twiny()
    ax2.set_xlabel("z/r")
    ax2.set_xlim(ax.get_xlim())
    # ax2.set_xticks(ax.get_xticks())    
    ax2.set_xticklabels(["{:.2f}".format(x) for x in nhver_to_zr(ix, ax.get_xticks(), model)])
            
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="best")
    ax.text(0.025, 0.020, rstr,
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
            
    ax.set_xlabel(r"r [au]")
    ax.set_ylabel(r"average abundance")
    ax.set_xlim([x.min(),x.max()])
    
    # do axis style
    ax.semilogy()     
    
    self._dokwargs(ax,**kwargs)    
    self._legend(ax)
    
    self.pdf.savefig(transparent=False)
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
            
    ax.set_xlabel(r"r [au]")    
    ax.set_ylabel(ylabel)    
    
    self._dokwargs(ax, **kwargs) 
    #self._legend(ax)
    
    self.pdf.savefig(transparent=False)
    plt.close(fig)   
    
  def plot_abunvert(self, model, r, species, useNH=True,
                    norm=None,styles=None,colors=None,markers=None,linewidths=None,**kwargs):
    '''
    Plot vertical abundances at a certain radius for the given species
    (can be more than one)
    '''
         
    print("PLOT: plot_abunvert ...")     
          
    rstr = r"r$\approx${:.2f} au".format(r)   
    
    fig, ax = plt.subplots(1, 1)     

    ix = (np.abs(model.x[:, 0] - r)).argmin()

    iplot=0
    ymin=1.e100
    ymax=-1.0
    
    for spec in species:           
      if useNH:
        old_settings = np.seterr(divide='ignore')     
        x = np.log10(model.NHver[ix, :])      
        np.seterr(**old_settings)  # reset to defaul
      else:
        x=model.z[ix,:]/model.x[ix,0]
      
      if spec in model.spnames:
        y = model.nmol[ix,:,model.spnames[spec]]/model.nHtot[ix,:]
        if norm != None:
          y=y/norm                    
        
        # FIXME: add proper treatment for styles and colors       
        if styles==None:         
          style="-"
          if "#" in spec: style="--"
        else: 
          style=styles[iplot]

        color=None          
        if colors!=None:
          color=colors[iplot]
          
        marker=None  
        if markers!=None:
          marker=markers[iplot]  
            
         
        lines=ax.plot(x, y, marker=marker, ms=4, markeredgecolor=color, markerfacecolor=color, 
                linestyle=style, color=color, 
                label="$\mathrm{"+spnToLatex(spec)+"}$")
              
        if linewidths != None:
          if linewidths[iplot] != None:
            lines[-1].set_linewidth(linewidths[iplot])  


                        
        iplot = iplot + 1
        if min(y) < ymin: ymin = min(y)
        if max(y) > ymax: ymax = max(y)
      
   
    if useNH:
      ax.set_xlim([17.5, x.max()])    
    ax.set_ylim(ymin,ymax)
    ax.semilogy()
                         
#     ax2 = ax.twiny()
#     ax2.set_xlabel("z/r")
#     ax2.set_xlim(ax.get_xlim())
#     #ax2.set_xticks(ax.get_xticks())    
#     ax2.set_xticklabels(["{:.2f}".format(x) for x in nhver_to_zr(ix, ax.get_xticks(), model)])
    
    if useNH:
      ax.set_xlabel(r"$\mathrm{\log\,N_{<H>}\,[cm^{-2}]}$ @" + rstr)
    else:
      ax.set_xlabel(r"z/r @" + rstr)
    ax.set_ylabel(r"$\mathrm{\epsilon(X)}$")    
    
    self._dokwargs(ax,**kwargs)
    self._legend(ax)
    self._closefig(fig)     
  
  def plot_abunrad(self, model, species, useNH=True,
                    norm=None,styles=None,colors=None,markers=None,linewidths=None,**kwargs):
    '''
    Plot vertical radial abundance in the midplane
    (can be more than one)
    Same as abunvert but radially more usefull for envelopes
    '''
         
    print("PLOT: plot_abunrad ...")                     
    
    fig, ax = plt.subplots(1, 1)     
    
    iplot=0
    ymin=1.e100
    ymax=-1.0
    
    for spec in species:           
      if useNH:
        old_settings = np.seterr(divide='ignore')     
        x = np.log10(model.NHrad[:, 0])      
        np.seterr(**old_settings)  # reset to defaul
      else:
        x=model.x[:,0]
      
      if spec in model.spnames:
        y = model.nmol[:,0,model.spnames[spec]]/model.nHtot[:,0]
        if norm != None:
          y=y/norm                    
        
        # FIXME: add proper treatment for styles and colors       
        if styles==None:         
          style="-"
          if "#" in spec: style="--"
        else: 
          style=styles[iplot]

        color=None          
        if colors!=None:
          color=colors[iplot]
          
        marker=None  
        if markers!=None:
          marker=markers[iplot]  
            
         
        lines=ax.plot(x, y, marker=marker, ms=4, markeredgecolor=color, markerfacecolor=color, 
                linestyle=style, color=color, 
                label="$\mathrm{"+spnToLatex(spec)+"}$")
              
        if linewidths != None:
          if linewidths[iplot] != None:
            lines[-1].set_linewidth(linewidths[iplot])  


                        
        iplot = iplot + 1
        if min(y) < ymin: ymin = min(y)
        if max(y) > ymax: ymax = max(y)
      
   
    if useNH:
      ax.set_xlim([17.5, x.max()])    
    ax.set_ylim(ymin,ymax)
    ax.semilogy()
                         
#     ax2 = ax.twiny()
#     ax2.set_xlabel("z/r")
#     ax2.set_xlim(ax.get_xlim())
#     #ax2.set_xticks(ax.get_xticks())    
#     ax2.set_xticklabels(["{:.2f}".format(x) for x in nhver_to_zr(ix, ax.get_xticks(), model)])
    
    if useNH:
      ax.set_xlabel(r"$\mathrm{\log\,N_{<H,rad>}\,[cm^{-2}]}$")
    else:
      ax.set_xlabel(r"r [au]")
    ax.set_ylabel(r"$\mathrm{\epsilon(X)}$")    
    
    self._dokwargs(ax,**kwargs)
    self._legend(ax)
    self._closefig(fig)       
    
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
      if spec not in model.spnames:
        print("WARN: Species "+spec+" not found") 
        continue
       
      x = model.x[:, 0]
      y = model.nmol[:,0,model.spnames[spec]]/model.nHtot[:,0]
      if norm is not None:
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
    ax.set_xlabel("r [au]")
    ax.set_ylabel(r"$\mathrm{\epsilon(X)}$")
            
    self._dokwargs(ax, **kwargs)
    self._legend(ax, loc="best") 
  
    self.pdf.savefig()
    plt.close(fig)
  
  def plot_dust_opac(self,model,dust=None,**kwargs):
    '''
    Plots the dust opacities (dust_opac.out) or the data given in the
    dust object    
    '''
    print("PLOT: dust opacities ...")
    
    fig, ax = plt.subplots(1, 1)
    
    x=model.dust.lam
    
    if dust is None: dust=model.dust
    
    ax.plot(x,dust.kabs,label="absorption")
    ax.plot(x,dust.ksca,label="scattering")
    ax.plot(x,dust.kext,label="extinction")   
    ax.set_xlabel(r"wavelength $\mathrm{[\mu m]}$")
    ax.set_ylabel(r"opacity $\mathrm{[cm^2 g(dust)^{-1}]}$")
    
    ax.set_xlim(0.09,3000.0)
    ax.set_ylim(1.e-2,None)
    
    ax.semilogx()
    ax.semilogy()
            
    self._dokwargs(ax,**kwargs)
    self._legend(ax)
    self._closefig(fig) 
  
  
  def plot_vertical(self, model, r, field, ylabel, ylog=True,zr=True,
                    xfield="zr",**kwargs):
    '''
    Plots a quantity (field) as a function of height at a certain radius
    radius.    
    '''
    print("PLOT: plot_vertical ...")
    rstr = r"r$\approx${:.1f} au".format(r) 
    
    fig, ax = plt.subplots(1, 1)      
    
    
    ix = (np.abs(model.x[:, 0] - r)).argmin()
     
    if zr and xfield=="zr": 
      x = model.z[ix, :] / model.x[ix, 0]
    elif xfield=="nH":
      old_settings = np.seterr(divide='ignore')     
      x = np.log10(model.NHver[ix, :])      
      np.seterr(**old_settings)  # reset to defaul
    elif xfield=="tg":
      x = model.tg[ix, :]
    else: 
      x = model.z[ix, :]
      
    y = field[ix, :]
                                    
    ax.plot(x, y)                    
    
    if zr:
      ax.invert_xaxis()
      ax.set_xlabel(r"z/r @ " + rstr)
    elif xfield=="nH":
      ax.set_xlabel(r"$\mathrm{\log\,N_{<H>}\,[cm^{-2}]}$ @" + rstr)  
    elif xfield=="tg":
      ax.set_xlabel(r"$\mathrm{\log\,T_{gas}\,[K]}$ @" + rstr)
      ax.invert_xaxis()  
    else:
      ax.set_xlabel(r"z [au] @ " + rstr)                   
      ax.invert_xaxis()
                        
    ax.set_ylabel(ylabel)    
    
    self._dokwargs(ax,**kwargs)
    self._legend(ax)
    
    self.pdf.savefig()
    plt.close(fig)     
    
  def plot_taus(self,model,r,**kwargs):
    '''
    Plot's taus (A_V, X-rays) as a function of vertical column density
    '''  
    ir=(np.abs(model.x[:, 0] - r)).argmin()
    rstr = "r={:.2f} au".format(model.x[ir,0])     
    
    fig, ax = plt.subplots(1, 1)
          
    old_settings = np.seterr(divide='ignore')           
          
    x=np.log10(model.NHver[ir,:])  
      
    ax.plot(x, model.tauX1[ir,:], color="blue", label=r"$\mathrm{\tau_{1\;keV}}$")
    ax.plot(x, model.tauX10[ir,:], "--",color="blue", label=r"$\mathrm{\tau_{10\;keV}}$")
    ax.plot(x, model.AVrad[ir,:], color="red", label=r"$\mathrm{A_V,rad}}$")
    ax.plot(x, model.AVver[ir,:], "--", color="red", label=r"$\mathrm{A_{V,ver}}$")
                                 
    ax.set_xlim(17.5,x.max())        
    ax.set_ylim(1.e-2,np.max([model.AVver[ir,:].max(),2.0]))
    
    np.seterr(**old_settings)  # reset to default 
      
    ax.hlines(1.0,ax.get_xlim()[0],ax.get_xlim()[1],linestyle=":") 
     
    ax2 = ax.twiny()
    ax2.set_xlabel("z/r")
    ax2.set_xlim(ax.get_xlim())
    #ax2.set_xticks(ax.get_xticks())    
    ax2.set_xticklabels(["{:.2f}".format(x) for x in nhver_to_zr(ir, ax.get_xticks(), model)])
    
    ax.set_xlabel(r"$\log$ N$_\mathrm{H}$ [cm$^{-2}$]")
    ax.set_ylabel(r"$\mathrm{A_V, \tau}$")  
    
    # do axis style
    ax.semilogy()     
    
    self._dokwargs(ax,**kwargs)
    
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="best", fancybox=False)
    ax.text(0.025, 0.025, rstr,
       verticalalignment='bottom', horizontalalignment='left',
       transform=ax.transAxes,alpha=0.75)
      
    self._closefig(fig)     

  def plot_starspec(self, model,**kwargs): 
    '''
    Plots the full Stellar Spectrum
    '''  
    print("PLOT: plot_starspec ...")
    fig, ax = plt.subplots(1, 1)      
            
    x = model.starSpec.lam[0::10]
    
    xmin=x.min()
    xmax=1000.0
    y = (model.starSpec.nu*model.starSpec.Inu)[0::10]      
         
    ax.plot(x, y, color="black")
                          
    # set defaults, can be overwritten by the kwargs
    
    ax.set_xlim([xmin,xmax])
    #ax.set_ylim([ymin,None])
    ax.semilogx()
    ax.semilogy()
    ax.set_xlabel(r"wavelength [$\mathrm{\mu}$m]")    
    ax.set_ylabel(r"$\mathrm{\nu F_{\nu}\,[erg\,cm^{-2}\,s^{-1}]}$")    
      
    self._dokwargs(ax, **kwargs)                        
    
    self.pdf.savefig()
    plt.close(fig)  



  def plot_sed(self, model,plot_starSpec=True,**kwargs): 
    '''
    Plots the seds and the StarSpectrum
    '''  
    print("PLOT: plot_sed ...")
    fig, ax = plt.subplots(1, 1)      
        
    xmin = 0.1
    ymin = 1.e-13
    if model.sed == None: return   
    # only use every 5 element to speed up plotting
    x = model.sed.lam
    y = model.sed.nu*model.sed.fnuErg      
    dist = ((model.sed.distance*u.pc).to(u.cm)).value                          
    ax.plot(x, y, marker=None, label=model.name)
    
    if plot_starSpec:
      # scale input Stellar Spectrum to the distance for comparison to the SED
      r= ((model.starSpec.r*u.R_sun).to(u.cm)).value      
      
      xStar = model.starSpec.lam[0::1]
      ystar= (model.starSpec.nu*model.starSpec.Inu)[0::1]
      yStar = ystar*(r**2.0*math.pi*dist**(-2.0))                                
      ax.plot(xStar, yStar, color="black")
                          
    # set defaults, can be overwritten by the kwargs
    
    ax.set_xlim([xmin,None])
    ax.set_ylim([ymin,None])
    ax.semilogx()
    ax.semilogy()
    ax.set_xlabel(r"wavelength [$\mathrm{\mu}$m]")    
    ax.set_ylabel(r"$\mathrm{\nu F_{\nu}\,[erg\,cm^{-2}\,s^{-1}]}$")    
      
    self._dokwargs(ax, **kwargs)                        
    
    self.pdf.savefig()
    plt.close(fig)      

class Contour(object):
  '''
  Define contourlines for one Contour for the filled contour plots.
  field needs to be an array of the same shape as the array data used for the
  filled 2D contour plots 
  '''
  def __init__(self, field,levels,colors="white",linestyles="solid",linewidths=1.5):
    self.field = field
    self.levels=levels
    self.colors=colors
    self.linestyles=linestyles
    self.linewidths=linewidths

