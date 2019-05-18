from __future__ import division 
from __future__ import print_function

from scipy.interpolate import interp1d

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import astropy.units as u
import math
import copy
import prodimopy.extinction as ext

class Plot(object):
  '''
  Plot routines for a single ProDiMo model.
  '''
  def __init__(self, pdf, fs_legend=None,title=None):
    self.pdf = pdf
    if fs_legend is None:
      self.fs_legend = mpl.rcParams['legend.fontsize']
    else:
      self.fs_legend=fs_legend
    self.ncol_legend = 5
    self.title=title
    
    # special colors, forgot the source for it :( somewhere from the internet)
    # FIXME: make an option to activate them
    self.pcolors={"blue"   : "#5DA5DA",
                  "orange" : "#FAA43A",
                  "green"  : "#60BD68",
                  "pink"   : "#F17CB0",
                  "brown"  : "#B2912F",
                  "purple" : "#B276B2",
                  "yellow" : "#DECF3F",
                  "red"    : "#F15854",
                  "gray"   : "#4D4D4D"}
  
  def _legend(self, ax,**kwargs):
    '''
    plots the legend, deals with multiple columns
    '''
    handles, labels = ax.get_legend_handles_labels()    
    
    if "loc_legend" in kwargs:
      loc=kwargs["loc_legend"]
    else:
      loc="best"
    
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
    Save and close the plot (Figure). 
    
    If self.pdf is None than nothing is done and the figure is returend
    
    #set the transparent attribut (used rcParam savefig.transparent)
    '''    
    
    #trans=mpl.rcParams['savefig.transparent']    
    if self.pdf is not None:
      self.pdf.savefig(figure=fig,transparent=False)
      plt.close(fig)
      return None
    else:
      return fig
    
  def _sfigs(self,**kwargs):
    '''
    Scale the figure size from matplotlibrc by the factors given in the 
    array sfigs (in kwargs) the first element is for the width the second for
    the heigth
    '''            
    if "sfigs" in kwargs:
      fac=kwargs["sfigs"]
      return scale_figs(fac)
    else:
      return scale_figs([1.0,1.0])
  
  def plot_NH(self, model, sdscale=False, muH=None,marker=None, **kwargs):
    '''
    Plots the total vertical hydrogen column number density 
    as a function of radius.
    
    Parameters
    ----------
    model : :class:`~prodimopy.read.Data_ProDiMo` 
      the model data
        
    sdscale : boolean 
      show additionally a scale with units in |gcm^-2|
    
    
    Returns
    -------
    :class:`~matplotlib.figure.Figure` or `None` 
      object if `self.pdf` is `None` the Figure object is reqturned, otherwise
      otherwise the plot is written directly into the pdf object(file) and 
      `None` is returned.
    '''
    print("PLOT: plot_NH ...")
    fig, ax = plt.subplots(1, 1,figsize=self._sfigs(**kwargs))      
    
    if muH is None:
      muH = model.muH
    
    x = model.x[:, 0]
    y = model.NHver[:, 0]    
    ax.plot(x, y, marker=marker,ms=3.0, color="black")
            
    ax.set_xlim(min(x), max(x))
    
         
    ax.semilogy()  
    ax.semilogx()   
    
    ax.set_xlabel(r"r [au]")
    ax.set_ylabel(r"N$_\mathrm{<H>,ver}\,\mathrm{[cm^{-2}]}$")

    self._dokwargs(ax,**kwargs)
    self._legend(ax)    
        
    # second sale on the right 
    if sdscale:
      ax2 = ax.twinx()    
      y = model.NHver[:, 0]*muH
      # just plot it again, is the easiest way (needs to be the same style etc)
      ax2.plot(x, y, color="black")
      ax2.set_ylabel(r"$\Sigma\,\mathrm{[g\,cm^{-2}]}$")
      ax2.set_xlim(min(x), max(x))
      
      # this needs to be done to get the correct scale
      ylim=np.array(ax.get_ylim())*muH
      ax2.set_ylim(ylim)
      # FIXME: check if this is required!
      ax2.semilogy()
    
    #ax.yaxis.tick_right()
    #ax.yaxis.set_label_position("right")
    #ax.yaxis.set_ticks_position('both')
        
    return self._closefig(fig)
     
 
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
    CS2 = ax.contourf(x, y, valinv,levels=[9.0,10.0,11.0], colors="0.6",hatches=["////","////","////"])
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
        if cont.filled is True:
          ACS=ax.contourf(x, y,cont.field,levels=cont.levels, 
                       colors=cont.colors,linestyles=cont.linestyles,linewidths=cont.linewidths)
        else:
          ACS=ax.contour(x, y,cont.field,levels=cont.levels, 
                       colors=cont.colors,linestyles=cont.linestyles,linewidths=cont.linewidths)
    
    CB = fig.colorbar(CS, ax=ax,ticks=ticks,pad=0.01)
    CB.ax.set_yticklabels(['X', 'SP', 'CR'])
    CB.ax.tick_params(labelsize=self.fs_legend) 
    
    
    # CB.set_ticks(ticks)
    CB.set_label("dominant ionization source",fontsize=self.fs_legend)  

    return self._closefig(fig)
    
   
  def plot_cont(self, model, values, label="value", zlog=True, 
                zlim=[None,None],zr=True,clevels=None,clabels=None,contour=True,
                extend="neither",oconts=None,acont=None,acontl=None,nbins=70,
                bgcolor=None,cb_format="%.1f",scalexy=[1,1],patches=None,
                rasterized=False,returnFig=False,**kwargs):
    '''
    Plot routine for 2D filled contour plots.
    
    Parameters
    ----------
    model : :class:`~prodimopy.read.Data_ProDiMo` 
      the model data
    
    values : array_like(float,ndim=2)
      a 2D array with numeric values for the plotting. E.g. any 2D array 
      of the :class:`~prodimopy.read.Data_ProDiMo` object.
    
    oconts : array_like(:class:`~prodimopy.plot.Contour`,ndim=1)
      list of :class:`~prodimopy.plot.Contour` objects which will be drawn 
      as additional contour levels. See also the example at the top of the page.
 
    
    scalexy: apply a scaling factor for the x and y coordinate (multiplicative)
    patches: a list of matplotlib.patches objects. For each object in the list 
    simply ax.add_patch() is called (at the very end of the routine)

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
      
    x = model.x*scalexy[0]  
    if zr:
      y = model.z / model.x
    else:
      y = np.copy(model.z)*scalexy[1] 
      y[:,0]=y[:,0]+0.05*scalexy[1] 
  
    levels = MaxNLocator(nbins=nbins).tick_values(maxval, minval)
        
    if clevels is not None:
      if zlog: clevels=np.log10(clevels)    
      ticks=clevels
    else: 
      ticks = MaxNLocator(nbins=6, prune="both").tick_values(minval, maxval)
  
    # TODO: maybe provide a routine for this, including scaling the figure size
    fig, ax = plt.subplots(1, 1,figsize=self._sfigs(**kwargs)) 

    # zorder is needed in case if rasterized is true
    CS = ax.contourf(x, y, pvals, levels=levels, extend=extend,zorder=-20)
    # This is the fix for the white lines between contour levels
    for c in CS.collections:
      c.set_edgecolor("face")    
    
    # rasterize the filled contours only, text ect. not  
    if rasterized: 
      ax.set_rasterization_zorder(-19)
  
    # axis equal needs to be done here already ... at least it seems so
    if "axequal" in kwargs:
      if kwargs["axequal"]: ax.axis('equal')

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
        if cont.showlabels:          
          ax.clabel(ACS, inline=True,inline_spacing=cont.label_inline_spacing,
                    fmt=cont.label_fmt,manual=cont.label_locations,fontsize=cont.label_fontsize)
    
    if acont is not None:      
      print("WARN: plot_cont: please use the oconts for additional contours ...")      
      #for l in acontl:
      #  ACS=ax.contour(x, y,pvals,levels=[l], colors='black',linestyles="solid",linewidths=1.5)
      #  ax.clabel(ACS, inline=1, fontsize=8,fmt=str(l))
      ACS=ax.contour(x, y,acont,levels=acontl, colors='white',linestyles="solid",linewidths=1.5)      
      # quick fix for second contour ... 
      #ACS2=ax.contour(x, y,model.nHtot,levels=[1.e6], colors='black',linestyles="solid",linewidths=2.5)
      #ax.clabel(ACS, inline=1, fontsize=8,fmt="%.0f")
      
#    ax.plot(np.sqrt(model.x[:,0]*model.x[:,0]+model.z[:,45]*model.z[:,45]),model.z[:,45],color="black")
#    ax.plot(np.sqrt(model.x[:,0]*model.x[:,0]+model.z[:,35]*model.z[:,35]),model.z[:,35],color="black")

#    ax.plot(model.x[:,0],model.z[:,48],color="black")
#    ax.plot(model.x[:,0],model.z[:,35],color="black")

      
    if bgcolor is not None:
      ax.set_axis_bgcolor(bgcolor)
    
    CB = fig.colorbar(CS, ax=ax,ticks=ticks,pad=0.01,format=cb_format)
    # FIXME: this is not very flexible and confusing
    if clabels is not None:
      CB.ax.set_yticklabels(clabels)
    #CB.ax.tick_params(labelsize=self.fs_legend) 
    # CB.set_ticks(ticks)
    CB.set_label(label)  
    
    if patches is not None:
      for patch in patches:
        ax.add_patch(patch)
    
    if returnFig:
      return fig
    else:
      return self._closefig(fig)

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
      
    return self._closefig(fig)
  
  
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
    
    return self._closefig(fig)

  def plot_avgabun(self,model,species,**kwargs):
    '''
    Plots the average abundance for the given species (can be more than one) 
    as a function of radius    
    '''
    print("PLOT: plot_avgabun ...")     
          
    fig, ax = plt.subplots(1, 1)     
  
    iplot = 0
    if(type(species)==str): species=[species]
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
    
    return self._closefig(fig)

  def plot_radial(self, model, values, ylabel, zidx=0, **kwargs):
    '''
    Plots a quantitiy along the radial grid for the given zidx (from the ProDiMo Array)
    as a function of radius. 
    
    `values` is any ProDiMo 2D array in `Data_ProDiMo` 
    or a 1D array with dim nx if `zidx` is `None` 
    
    '''
    print("PLOT: plot_radial ...")
    fig, ax = plt.subplots(1, 1)      
    
    
    
    if zidx is None:
      x = np.sqrt(model.x[:, 0]**2.0+model.z[:, 0]**2.0)
      y = values      
    else:
      x = np.sqrt(model.x[:, zidx]**2.0+model.z[:, zidx]**2.0)
      y = values[:,zidx]         
      
    ax.plot(x, y,marker=None)
                 
    ax.set_xlim(np.min(x),np.max(x))
    ax.set_ylim(np.min(y),np.max(y))
    ax.semilogy()
            
    ax.set_xlabel(r"r [au]")    
    ax.set_ylabel(ylabel)    
    
    self._dokwargs(ax, **kwargs) 
    #self._legend(ax)
    
    return self._closefig(fig)

  def plot_cdnmol(self, model, species,colors=None,styles=None, 
                  scalefacs=None,norm=None,normidx=None,ylabel="$\mathrm{N_{ver}\,[cm^{-2}}]$", 
                  patches=None,**kwargs):
    '''
    Plots the vertical column densities as a function of radius for the 
    given species.
    '''
    print("PLOT: plot_cdnmol ...")
    
    
    if colors is None:
      colors=[None]*len(species)
    
    if styles is None:
      styles=[None]*len(species)
    
    fig, ax = plt.subplots(1, 1,figsize=self._sfigs(**kwargs))      
    
    
    x = model.x[:, 0]
    
    ymin=1.e99
    ymax=1.e-99
    
    if scalefacs is None:
      scalefacs=np.ones(len(species))
      
    if normidx is not None:
      normspec=model.cdnmol[:,0,normidx]      
    
    iplot=0  
    if(type(species)==str): species=[species]
    for spec,fac in zip(species,scalefacs):
      
      ispec=-1
      try:
        ispec=model.spnames[spec]
      except:
        print("WARNING: Could not find species: ",spec)
        continue
        
      if normidx is not None:
        y=model.cdnmol[:,0,ispec]/normspec

      y=model.cdnmol[:,0,ispec]/fac
      if norm is not None:
        y=y/norm
      
      label="$"+spnToLatex(spec)+"$"
      if fac != 1.0 and norm is None:
        label+="/"+"{:3.1e}".format(fac)
              
      ax.plot(x, y,marker=None,linestyle=styles[iplot],color=colors[iplot],
              label=label)
         
      ymin=np.min([np.min(y),ymin])
      ymax=np.max([np.max(y),ymax])
      iplot=iplot+1
      
    if patches is not None:
      for patch in patches:
        ax.add_patch(copy.copy(patch))
    
    ax.set_xlim(np.min(x),np.max(x))
    ax.set_ylim(ymin,ymax)
    ax.semilogy()
            
    ax.set_xlabel(r"r [au]")
    ax.set_ylabel(ylabel)
    
    self._dokwargs(ax, **kwargs) 
    self._legend(ax,**kwargs)
    
    return self._closefig(fig)


  def plot_midplane(self, model, fieldname, ylabel, **kwargs):
    '''
    Plots a quantitiy in in the midplane as a function of radius
    fieldname is any field in Data_ProDiMo
    FIXME: remove the fieldname stuff passe  the whole array ... 
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
    
    return self._closefig(fig)

    
  def plot_abuncont(self, model, species='O', rel2H=True, label=None, zlog=True, 
                zlim=[None,None],zr=True,clevels=None,clabels=None,contour=True,
                extend="neither",oconts=None,acont=None,acontl=None,nbins=70,
                bgcolor=None,cb_format="%.1f",scalexy=[1,1],patches=None,
                rasterized=False,**kwargs):
    '''
    Plots the 2D abundance structure of a species.
    
    This is a convenience function and is simply a wrapper for 
    :func:`~prodimopy.plot.Plot.plot_cont` routine.
    
    The routine checks if the species exists, calculates the abundance and sets some 
    defaults (e.g. label) for the :func:`~prodimopy.plot.Plot.plot_cont` routine and calls it. 
    However, all the defaults can be overwritten by providing the corresponding parameter.
             
    Contributors: L. Klarmann, Ch. Rab 

    Parameters
    ----------
        
    model : :class:`prodimopy.read.Data_ProDiMo` 
      the model data 
      
    species : str 
      the name of the species as given in |prodimo|
    
    rel2H : boolean 
      plot abundances relative to the total H nuclei number density. 
      If `False` the number density of the species is plotted      
      
    label : str 
      the label for the colorbar. If None the default is plotted
    
    
    For all other parameters see :func:`~prodimopy.plot.Plot.plot_cont`

    TODO: can be improved with better and smarter default values (e.g. for the colorbar)    
    '''
    print('PLOT: plot_abuncont ...')

    # Check if species names exists
    try:
      n_rel_index = model.spnames[species]

    except KeyError:      
      print("The species "+species+ '''you want to access does not exist  
             or is spelled incorrectly. Exiting plot_abuncont routine''')
      return
    
    
    if rel2H:
      values=model.getAbun(species)      
      
      if label is None:
        label=r"$\mathrm{\epsilon("+spnToLatex(species)+")}$"
        if zlog: label="log "+label
      
      # define some default lower limit
      if zlim == [None,None]:
        zlim=[3.e-13,None]
        extend="both"
        
    else:
      values=model.nmol[:,:,n_rel_index]
      if label is None:
        label=r"$\mathrm{n("+spnToLatex(species)+") [cm^{-3}]}$"
        if zlog: label="log "+label
             
    return self.plot_cont(model, values, label=label, zlog=zlog, 
                zlim=zlim,zr=zr,clevels=clevels,clabels=clabels,contour=contour,
                extend=extend,oconts=oconts,acont=acont,acontl=acontl,nbins=nbins,
                bgcolor=bgcolor,cb_format=cb_format,scalexy=scalexy,patches=patches,
                rasterized=rasterized,**kwargs)
  
    
  def plot_abunvert(self, model, r, species, useNH=True,scaling_fac=None,
                    norm=None,styles=None,colors=None,markers=None,linewidths=None,
                    useT=False,**kwargs):
    '''
    Plot vertical abundances at a certain radius for the given species
    (can be more than one)
    '''
         
    print("PLOT: plot_abunvert ...")
          
    if colors is None:
      colors=list(self.pcolors.values())
    
    rstr = r"r$\approx${:.2f} au".format(r)   
    
    fig, ax = plt.subplots(1, 1)     

    ix = (np.abs(model.x[:, 0] - r)).argmin()

    iplot=0
    ymin=1.e100
    ymax=-1.0
    if(type(species)==str): species=[species]
    for spec in species:           
      if useT:
        x = model.td[ix, :]
      elif useNH:
        old_settings = np.seterr(divide='ignore')     
        x = np.log10(model.NHver[ix, :])      
        np.seterr(**old_settings)  # reset to defaul      
      else:
        x=model.z[ix,:]/model.x[ix,0]
      
      # check if list of names in list of names, than sum them up
      if isinstance(spec, (list, tuple, np.ndarray)):
        y=model.nHtot[ix,:]*0.0 # just to get an array
        for name in spec:
          y=y+(model.getAbun(name)[ix,:])
      
      elif spec in model.spnames:
        y = model.nmol[ix,:,model.spnames[spec]]/model.nHtot[ix,:]
      else:
        continue  
        
      if norm is not None:
        y=y/norm     
        
      if scaling_fac is not None:
        y=y*scaling_fac[iplot]
      
      # FIXME: add proper treatment for styles and colors       
      if styles==None:         
        style="-"
        if "#" in spec: style="--"
      else: 
        style=styles[iplot]

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
      
   
    if useT:
      ax.set_xlim([30, 5])
    elif useNH:
      ax.set_xlim([17.5, x.max()]) 
    else:
      ax.invert_xaxis() #(z/r=0 on the right)   
    ax.set_ylim(ymin,ymax)
    ax.semilogy()
                         
#     ax2 = ax.twiny()
#     ax2.set_xlabel("z/r")
#     ax2.set_xlim(ax.get_xlim())
#     #ax2.set_xticks(ax.get_xticks())    
#     ax2.set_xticklabels(["{:.2f}".format(x) for x in nhver_to_zr(ix, ax.get_xticks(), model)])
    
    if useT:
      ax.set_xlabel(r"$\mathrm{T_d [K]}$ @" + rstr)
    elif useNH:
      ax.set_xlabel(r"$\mathrm{\log\,N_{<H>}\,[cm^{-2}]}$ @" + rstr)
    else:
      ax.set_xlabel(r"z/r @" + rstr)
    ax.set_ylabel(r"$\mathrm{\epsilon(X)}$")    
    
    self._dokwargs(ax,**kwargs)
    self._legend(ax,**kwargs)
    return self._closefig(fig)     
  
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
    if(type(species)==str): species=[species]
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
    return self._closefig(fig)       
    
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
    if(type(species)==str): species=[species] 
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
    self._legend(ax) 
  
    return self._closefig(fig)
  
  def plot_dust_opac(self,model,dust=None,**kwargs):
    '''
    Plots the dust opacities (dust_opac.out) or the data given in the
    dust object    
    '''
    print("PLOT: dust opacities ...")
    
    fig, ax = plt.subplots(1, 1,figsize=self._sfigs(**kwargs))
    
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
    return self._closefig(fig) 
  
  
  def plot_vertical(self, model, r, field, ylabel, zr=True,
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
    
    return self._closefig(fig)
    
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
      
    return self._closefig(fig)     

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
    
    return self._closefig(fig)


  def plot_sed(self, model,plot_starSpec=True,sedObs=None,unit="erg",reddening=False,**kwargs): 
    '''
    Plots the seds and the StarSpectrum
    '''  
    print("PLOT: plot_sed ...")
    fig, ax = plt.subplots(1, 1,figsize=self._sfigs(**kwargs))
        
    xmin = 0.1
    ymin = 1.e-13
    if model.sed == None: return   
    # only use every 5 element to speed up plotting
    x = model.sed.lam
    if unit=="W":
      y = model.sed.nuFnuW
    elif unit == "Jy":
      y = model.sed.fnuJy
    else:  
      y = model.sed.nu*model.sed.fnuErg
      
    if reddening==True and sedObs is not None and sedObs.A_V is not None:
      # idx validity of extinction function
      ist=np.argmin(np.abs(x-0.0912))
      ien=np.argmin(np.abs(x-6.0))
      y[ist:ien]=y[ist:ien]/ext.reddening(x[ist:ien]*1.e4, a_v=sedObs.A_V,r_v=sedObs.R_V, model="f99")
    
    dist = ((model.sed.distance*u.pc).to(u.cm)).value
    
    if plot_starSpec:
      # scale input Stellar Spectrum to the distance for comparison to the SED
      r= ((model.starSpec.r*u.R_sun).to(u.cm)).value      
      
      xStar = model.starSpec.lam[0::1]
      yStar= (model.starSpec.nu*model.starSpec.Inu)[0::1]
      yStar = yStar*(r**2.0*math.pi*dist**(-2.0))
      
      if unit=="W":
        yStar=(yStar*u.erg/(u.s*u.cm**2)).to(u.Watt/u.m**2).value
      elif unit == "Jy":
        yStar= (model.starSpec.Inu)[0::1]
        yStar = yStar*(r**2.0*math.pi*dist**(-2.0))
        yStar=(yStar*u.erg/(u.s*u.cm**2*u.Hz)).to(u.Jy).value
      
      ax.plot(xStar, yStar, color="black")

    # plot the SED  
    ax.plot(x, y, marker=None, label=model.name)
    
    if sedObs is not None:
      okidx=np.where(sedObs.flag=="ok")

      xsedObs=sedObs.lam
      ysedObs=sedObs.nu*sedObs.fnuErg
      ysedObsErr=sedObs.nu*sedObs.fnuErgErr

      if unit=="W":
        ysedObs=sedObs.nu*((sedObs.fnuJy*u.Jy).si.value)
        ysedObsErr=sedObs.nu*((sedObs.fnuJyErr*u.Jy).si.value)
      elif unit=="Jy":
        ysedObs=sedObs.fnuJy
        ysedObsErr=sedObs.fnuJyErr

      #ax.plot(sedObs.lam[okidx],sedObs.nu[okidx]*sedObs.fnuErg[okidx],linestyle="",marker="x",color="0.5",ms=3)
      ax.errorbar(xsedObs[okidx],ysedObs[okidx], yerr=ysedObsErr[okidx],
                  fmt='o',color="0.5",ms=2,linewidth=1.0,zorder=0)

      ulidx=np.where(sedObs.flag=="ul")
      ax.plot(xsedObs[ulidx],ysedObs[ulidx],linestyle="",marker="v",color="0.5",ms=2.0)
      
      # FIXME: no proper unit treatment yet
      if sedObs.specs is not None:
        for spec in sedObs.specs:
          nu=(spec[:,0]* u.micrometer).to(u.Hz, equivalencies=u.spectral()).value
          fnuerg=(spec[:,1]* u.Jy).cgs.value
          ax.plot(spec[:,0],nu*fnuerg,linestyle="-",linewidth=0.5,color="0.5")

                          
    # set defaults, can be overwritten by the kwargs
    ax.set_xlim([xmin,None])
    ax.set_ylim([ymin,None])
    ax.semilogx()
    ax.semilogy()
    ax.set_xlabel(r"wavelength [$\mathrm{\mu}$m]")    
    if unit == "W":
      ax.set_ylabel(r"$\mathrm{\lambda F_{\lambda}\,[W\,m^{-2}]}$")
    elif unit == "Jy":
      ax.set_ylabel(r"$\mathrm{flux\,[Jy]}$")
    else:
      ax.set_ylabel(r"$\mathrm{\nu F_{\nu}\,[erg\,cm^{-2}\,s^{-1}]}$")
#    ax.yaxis.tick_right()
#    ax.yaxis.set_label_position("right")
      
    self._dokwargs(ax, **kwargs)                        
    
    return self._closefig(fig)


  def _getSEDana_boxpoints(self,lam,model,zr=True):
    '''
    Creates an array of (x,y) coordinates representing the emission origin for 
    the SEDana which can be used or the given wavelength. Those coordinates can 
    be used to draw a box on a plot (e.g. can be passed to a matplotlib Polygon)
    to draw a box (Polygon).
    
    Parameters
    ----------
    lam : float 
      the wavelength for which the emission origin should be calculated
    model : :class:`prodimopy.read.Data_ProDiMo` 
      the ProDiMo model including the SED analysis data (SEDana) 
    zr : boolean
      If `zr==True` (default) then the z coordinate of the points is returned in 
      z/r units. Optional parameter. 
  
    Returns
    -------
    array_like(float,ndim=1) :
      list of (x,y) points (in au). if zr=True the z coordinate is in z/r units.
      
      
    TODO: maybe merge somehow with :func:`~prodimopy.Data_ProDiMo.getSEDAnaMask`
    '''
    # interpolate       
    sedAna = model.sed.sedAna
    
    r15=interp1d(sedAna.lams, sedAna.r15, bounds_error=False, fill_value=0.0, kind="linear")(lam)
    r85=interp1d(sedAna.lams, sedAna.r85, bounds_error=False, fill_value=0.0, kind="linear")(lam)
    xi15=np.argmin(np.abs(model.x[:,0]-r15))
    xi85=np.argmin(np.abs(model.x[:,0]-r85))
     
    z85s=[[model.x[ix,0],interp1d(sedAna.lams, sedAna.z85[:,ix], bounds_error=False, fill_value=0.0, kind="linear")(lam)] 
           for ix in range(xi15,xi85)]
    z15s=[[model.x[ix,0],interp1d(sedAna.lams, sedAna.z15[:,ix], bounds_error=False, fill_value=0.0, kind="linear")(lam)] 
          for ix in range(xi85-1,xi15-1,-1)]
    points=z85s+z15s
    
    for point in points:
      if zr is True:
        point[1]=point[1]/point[0]

    return points


  def plot_sedAna(self, model,lams=[1.0,10,100.0,1000.0],field=None,label=None,boxcolors=None, zlog=True, 
                zlim=[None,None],zr=True,clevels=None,clabels=None,
                extend="neither",oconts=None,nbins=70,
                bgcolor=None,cb_format="%.1f",scalexy=[1,1],patches=None,
                rasterized=False,**kwargs):  
    '''
    Plots the SED analysis stuff (origin of the emission).
    '''  
    print("PLOT: plot_sedAna ...")
    
    if boxcolors is None:
      boxcolors=[self.pcolors["purple"],self.pcolors["orange"],self.pcolors["red"],"black"]
  
    if patches is None:
      patches=list()
        
    ibox=0
    for lam in lams:
      
      patch = mpl.patches.Polygon(self._getSEDana_boxpoints(lam, model, zr=True),
                                  True,fill=False,color=boxcolors[ibox],zorder=100,linewidth=2.0)
        
      patches.append(patch)
      ibox+=1
  
    if field is None:
      field=model.nHtot
      label=r"log $n_\mathrm{<H>}\,\mathrm{[cm^{-3}]}$"
      oconts=[Contour(model.AV,[1],linestyles="--",colors=self.pcolors["gray"])]
      
    fig =self.plot_cont(model, field, label=label, zlog=zlog, 
                zlim=zlim,zr=zr,clevels=clevels,clabels=clabels,contour=False,
                extend=extend,oconts=oconts,nbins=nbins,
                bgcolor=bgcolor,cb_format=cb_format,scalexy=scalexy,patches=patches,
                rasterized=rasterized,returnFig=True,**kwargs)
    

    ax=fig.axes[0]    
    
    ibox=0    
    for lam in lams:
      ax.text(0.02, 0.92-ibox/18.0,"$"+"{:5.1f}".format(lam)+r"\,\mathrm{\mu m}$",
              horizontalalignment='left',
              verticalalignment='bottom',fontsize=6,
              transform = ax.transAxes,color=boxcolors[ibox],
              bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
      ibox+=1    

    self._dokwargs(ax, **kwargs)                        
    
    return self._closefig(fig)

  

  def plot_taulines(self, model, lineIdents, **kwargs):
    '''
    Plots the line optical depth as a function of radius for the given lines.
    The lines are identified via a list of lineIdents containt of an array with 
    ident and wavelength of the line e.g. ["CO",1300.0]. 
    It searches for the closest lines.
    
    TODO: there is no options for linestyles and colors yet (the defaults are used). 
    '''
    print("PLOT: plot_taulines ...")
    fig, ax = plt.subplots(1, 1,figsize=self._sfigs(**kwargs))  
    
    xmin = 1.e100
    xmax = 0
        
    iplot = 0   
    for lineIdent in lineIdents: 
      x = model.x[:, 0]
      lineEstimate = model.getLineEstimate(lineIdent[0], lineIdent[1])      
      y = list()
      # FIXME: why is there a loop
      for rInfo in lineEstimate.rInfo:
        y.append(rInfo.tauLine)
      
      ax.axhline(y=1.0, linestyle="-", color="black", linewidth=0.5)
      label=r"$\mathrm{"+spnToLatex(lineEstimate.ident) + "}$ " + "{:.2f}".format(lineEstimate.wl) + " $\mathrm{\mu m}$"      
      ax.plot(x, y, marker=None, label=label)     
          
      iplot = iplot + 1
      
      if min(x) < xmin: xmin = min(x)
      if max(x) > xmax: xmax = max(x)
         
    ax.set_xlim(xmin, xmax)
       
    ax.semilogx()       
    ax.semilogy()
    
    ax.set_xlabel(r"r [au]")    
    ax.set_ylabel(r"$\mathrm{\tau_{line}}$")
    
  
    self._dokwargs(ax,**kwargs)
    self._legend(ax,**kwargs)  
  
    return self._closefig(fig)

  
  def plot_line_origin(self,model,ids,field, label="value", boxcolors=None, 
                       boxlinestyles=None,showLineLabels=True, zlog=True, 
                zlim=[None,None],zr=True,clevels=None,clabels=None,contour=False,
                extend="neither",oconts=None,nbins=70,
                bgcolor=None,cb_format="%.1f",scalexy=[1,1],patches=None,
                rasterized=False,showContOrigin=False,**kwargs): 
    '''
    Plots the line origins for a list of lineestimates given by their ids
    ([lineid,wavelength])
    
    Does not give the exact same results as the corresponding idl routines
    as we do not use interpolation here. We rather use the same method for 
    calculating the averaged values over the emission area and for plotting
    this area.
    '''
    if boxcolors is None:
      boxcolors=[self.pcolors["gray"],self.pcolors["orange"],self.pcolors["brown"],self.pcolors["red"],
                 self.pcolors["purple"]]
    
    if boxlinestyles is None:
      boxlinestyles=["-"]*5
  
    lestimates=list()
    for id in ids:
      lestimates.append(model.getLineEstimate(id[0],id[1]))

    if patches is None:
      patches=list()
    
    ibox=0
    for lesti in lestimates:
      # to be consistent we use the LineOriginMask to determine the box
      # as a result the plotted region is not necessarely the same as in idl
      # as we do not interpolate here
      xmasked=np.ma.masked_array(model.x,mask=model.getLineOriginMask(lesti))
      x15=np.min(xmasked)
      x85=np.max(xmasked)      
      xi15=np.argmin(np.abs(model.x[:,0]-x15))
      xi85=np.argmin(np.abs(model.x[:,0]-x85))

      z85s=[[model.x[rp.ix,0],rp.z85] for rp in lesti.rInfo[xi15:xi85+1]]
      z15s=[[model.x[rp.ix,0],rp.z15] for rp in lesti.rInfo[xi85:xi15-1:-1]]
      points=z85s+z15s

      if zr is True:
        for point in points:
          point[1]=point[1]/point[0]
      
      patch = mpl.patches.Polygon(points,True,fill=False,color=boxcolors[ibox],
                                  linestyle=boxlinestyles[ibox],
                                  zorder=100,linewidth=2.0)
        
      patches.append(patch)
      
      if showContOrigin is True:
        pointsc=self._getSEDana_boxpoints(lesti.wl, model,zr)
        patchc = mpl.patches.Polygon(pointsc,True,fill=False,color=boxcolors[ibox],zorder=100,linewidth=1.0,linestyle="--")
        patches.append(patchc)
      
      ibox+=1
  
    fig =self.plot_cont(model, field, label=label, zlog=zlog, 
                zlim=zlim,zr=zr,clevels=clevels,clabels=clabels,contour=False,
                extend=extend,oconts=oconts,acont=None,acontl=None,nbins=nbins,
                bgcolor=bgcolor,cb_format=cb_format,scalexy=scalexy,patches=patches,
                rasterized=rasterized,returnFig=True,**kwargs)
    
    ax=fig.axes[0]
    
    if showLineLabels:
      ibox=0    
      for idline in ids:
        ax.text(0.02, 0.92-ibox/18.0,"$"+spnToLatex(idline[0])+"$ "+str(idline[1]),
                horizontalalignment='left',
                verticalalignment='bottom',fontsize=6,
                transform = ax.transAxes,color=boxcolors[ibox],
                bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
        ibox+=1
      
    return self._closefig(fig)

  def plot_lineprofile(self,model,wl,ident=None,linetxt=None,lineObs=None,**kwargs):
    '''
    Plots the line profile for the given line (id wavelength and line ident) 
    '''  
    print("PLOT: plot_line_profile ...")
    fig, ax = plt.subplots(1, 1,figsize=self._sfigs(**kwargs))      
    
    line=model.getLine(wl,ident=ident)
    if line==None: 
      print("WARN: line "+ str(ident)+" at "+str(line.wl)+" micron not found")  
    
    # text for the title
 
    x = line.profile.velo           
    y = line.profile.flux-line.fcont  # remove the continuum

    if linetxt is None:
      linetxt=line.species+"@"+"{:.2f}".format(line.wl)+" $\mathrm{\mu m}$"    
    ax.plot(x, y,marker=None, label=linetxt)
    
    # plot the line profile if it exists
    if lineObs is not None:
      # FIXME: this is not very nice
      # make a warning if lineObs and line Data are not consistent
      # it could be that they are for one model       
      lineIdx=model._getLineIdx(wl,ident=ident)

      line=lineObs[lineIdx]
      if line.profile is not None:
        x=line.profile.velo
        y = line.profile.flux  # remove the continuum      
        if line.profileErr is not None:
          ax.fill_between(x, y-line.profileErr , y+line.profileErr, color='0.8',zorder=0)
        ax.plot(x, y, marker=None, color="black", label="Obs.",zorder=0)


    ax.set_xlabel(r"$\mathrm{velocity\,[km\;s^{-1}}$]")    
    ax.set_ylabel(r"$\mathrm{flux\,[Jy]}$")    
      
    self._dokwargs(ax, **kwargs)                
    self._legend(ax,**kwargs)
               
    return self._closefig(fig)
  

  def plot_heat_cool(self,model,zr=True,oconts=None,**kwargs):
    """
    Plots the dominant heating and cooling processes. 
    
    The initial python code for this routine is from Frank Backs 
    
    TODO: check if too many heating/cooling processes (e.g. out of bound for colors)
    TODO: check for insignificant dominant processes (e.g. only one pixel)
    TODO: possibility to have different oconts for the heating and cooling figures 
    """
    colors = np.array([(230, 25, 75), (60, 180, 75), (255, 225, 25), (0, 130, 200),
                        (245, 130, 48), (145, 30, 180), (70, 240, 240), (240, 50, 230),
                        (210, 245, 60), (250, 190, 190), (0, 128, 128), (230, 190, 255),
                        (170, 110, 40), (255, 250, 200), (128, 0, 0), (170, 255,95),
                        (128, 128, 0), (255, 215, 180), (0, 0, 128), (128, 128, 128),
                        (0,0,0),(255,255,255)], dtype=float)
    colors/=255  
 
    if zr:
      z=model.z/model.x
    else:
      z=model.z
 
    # used for plotting
    max_idx=np.zeros(shape=(model.nx, model.nz), dtype='int16')
 
    # list of all the dominant heating processes
    idxlisth=np.unique(model.heat_mainidx)    
    idxlistc=np.unique(model.cool_mainidx)    
    
    fig, axarr = plt.subplots(1, 2, figsize=self._sfigs(sfigs=[2.0,1.3])) 
    plt.subplots_adjust(bottom=0.3)
    axh=axarr[0]
    axc=axarr[1]

    # this if for the labels, and also maps the colors to the names  
    for i in range(len(idxlisth)):      
      # -1 because python starts at zeror
      axh.scatter(0, 0, marker="s", color=colors[i], label=model.heat_names[idxlisth[i]-1])
      
      # this is necessary to have the fields with increasing number without 
      # gaps, otherwhise the colormapping in pcolormesh does not work
      max_idx[model.heat_mainidx==idxlisth[i]]=i
      
    cMap = mpl.colors.ListedColormap(colors[0:len(idxlisth)-1])   
    axh.pcolormesh(model.x, z, max_idx, linewidth=0,cmap=cMap)
    axh.legend(loc='upper center', bbox_to_anchor=(0.5, -0.175), ncol=2, 
               frameon=False,fontsize=5.5)
    axh.set_title("dominant heating processes")
    

    # now for the cooling
    # this if for the labels, and also maps the colors to the names  
    for i in range(len(idxlistc)):      
      # -1 because python starts at zeror
      axc.scatter(0, 0, marker="s", color=colors[i], label=model.cool_names[idxlistc[i]-1])
      
      # this is necessary to have the fields with increasing number without 
      # gaps, otherwhise the colormapping in pcolormesh does not work
      max_idx[model.cool_mainidx==idxlistc[i]]=i
      
    cMap = mpl.colors.ListedColormap(colors[0:len(idxlistc)-1])   
    axc.pcolormesh(model.x, z, max_idx, linewidth=0,cmap=cMap)
    axc.legend(loc='upper center', bbox_to_anchor=(0.5, -0.175), ncol=2, 
               frameon=False,fontsize=5.5)
    axc.set_title("dominant cooling processes")


    for ax in [axh,axc]:
      ax.set_xlim(np.min(model.x),None)             
      ax.semilogx()
      ax.set_xlabel("r [au]")
      if zr:
        ax.set_ylim(0,None)
        ax.set_ylabel("z/r")
      else:
        ax.set_ylabel("z [au]")
      
      # Additional Contours, plot for both plots at the moment   
      if oconts is not None:
        for cont in oconts:
          ACS=ax.contour(model.x, z,cont.field,levels=cont.levels, 
                         colors=cont.colors,linestyles=cont.linestyles,linewidths=cont.linewidths)
          if cont.showlabels:          
            ax.clabel(ACS, inline=True,inline_spacing=cont.label_inline_spacing,
                      fmt=cont.label_fmt,manual=cont.label_locations,fontsize=cont.label_fontsize)
  
      
    self._dokwargs(axh,**kwargs)
    self._dokwargs(axc,**kwargs)
  
    return self._closefig(fig)


class Contour(object):
  '''
  Define an additional contour which can be used in the contour plotting
  routines.
  
  Object of this class can be passed to e.g. the :func:`~prodimopy.plot.Plot.plot_cont` routine and will be drawn their. 
  
  TODO: provide a field for label strings (arbitrary values) need to be the same size as levels
  '''
  def __init__(self, field,levels,colors="white",linestyles="solid",linewidths=1.5,
               showlabels=False,label_locations=None,label_fmt="%.1f",
               label_fontsize=7,
               label_inline_spacing=5,
               filled=False):
    '''
    Attributes
    ----------  
      
    '''
    self.field = field
    """ array_like(float,ndim=2) :
      A 2D array of values used for the Contours. Needs to have the same 
      dimensions as the array used for the contour plotting routine. So any 
      2D array of the :class:`prodimopy.read.Data_ProDiMo` will do.
    """ 
    self.levels=levels
    """ array_like(float,ndim=1) :
      list of values for which contour lines should be drawn.
    """    
    self.colors=colors
    """ array_like(ndim=1) : 
      list of colors for the idividual contours. If only a single value is 
      provided (i.e. no array) this value is applied to all contours.
      The values of colors can be given in the same way as it is done for 
      matplotlib.
    """
    self.linestyles=linestyles
    """ array_like(ndim=1) :
      linestyles for the contours. Works like the `colors` parameter. 
      Any style that matplotlib understands will work.
    """ 
    self.linewidths=linewidths    
    """ array_like(ndim=1) :
      linewidths for the individual contour levels. Works the same as the 
      `colors` parameter.
    """    
    self.showlabels=showlabels
    """ boolean : `False`
      show text label for each level or not (default: False)
      Still kind of experimental
    """      
    self.label_locations=label_locations
    self.label_fmt=label_fmt
    self.label_fontsize=label_fontsize
    self.label_inline_spacing=label_inline_spacing
    self.filled=filled


def spnToLatex(spname):
  """
  Utilitiy function to convert species names to proper latex math strings.
  
  The returned string can directly be embedded in a latex $ $ statement. 
  """
  # use string in case it is a binary format (python 3 comaptibility)      
  name = str(spname)      
  # TODO: make this a bit smarter    
  if str(spname) == "HN2+": name = "N2H+" 
  if str(spname) == "C18O": return "C^{18}O"
  if str(spname) == "13CO": return "^{13}CO"
  if str(spname) == "H13CO+": return "H^{13}CO^+"
  
  newname = ""  
  previous_char=None
  for c in name:    
    if c.isdigit():
      newname += "_" + c
    # deal with ortho and para (o-, p-) species   
    elif c == "-" and not (previous_char=="o" or previous_char=="p"):
      newname += "^-"
    elif c == "+":
      newname += "^+"
    elif c == "#":
      newname += "\#"       
    else:
      newname += c
      
    previous_char=c
      
  # for line names (species)    
  if "_H" in newname:
    newname=newname.replace("_H","\_H")
      
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


def scale_figs(scale):
  '''
  Scale the figure size from matplotlibrc by the factors given in the 
  array scale the first element is for the width the second for
  the heigth.
  '''            
  figsize=mpl.rcParams['figure.figsize']
  
  return (figsize[0]*scale[0],figsize[1]*scale[1])
