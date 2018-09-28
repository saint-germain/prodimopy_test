"""
Default script for the plotting the results of a single ProDiMo model. 

A script called `pplot`, which can directly by called from the command line, 
will be installed automatically during the installation process.  
"""
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals

# the above two statement are for phyton2 pyhton3 compatibility.
# With this statmens you can write python3 code and it should also work
# in python2 (depends on your code) 

# this is for the argument parser included in python
import argparse

# this is for the PDF output
from matplotlib.backends.backend_pdf import PdfPages

# Thats the prodimpy modules 
# The first one is for reading a ProDiMo Model 
# The second one is for plotting
import prodimopy.read as pread
import prodimopy.plot as pplot
import numpy

# The main routine is require to have an entry point. 
# It is not necessary if you want to write your own script.
def main(args=None):
  ###############################################################################
  # Command line parsing
  # this is optional you do not have to do it this way. 
  # You can use the prodimopy module in any way you want
  parser = argparse.ArgumentParser(description='prodimopy simple example for plotting ')
  parser.add_argument('-dir', required=False, default=".", help='The directory of the input files DEFAULT: "."')
  parser.add_argument('-output', required=False, default="./out.pdf", help='The output filename, DEFAULT: "./out.pdf"')
  parser.add_argument('-td_fileIdx', required=False, default=None, help='The index for the time dependent output file e.g. "01", DEFAULT: None')
  args = parser.parse_args()
  
  print("-dir: ", args.dir)
  print("-output: ", args.output)
  print("-td_fileIdx: ", args.td_fileIdx)
  
  # thats for time dependent models, also optional
  outfile=args.output
  if args.td_fileIdx != None:
    outfile=outfile.replace(".pdf","_"+args.td_fileIdx+".pdf")
  
  # This reads the output of a ProDiMo model
  # there are more optional arguments 
  pd = pread.read_prodimo(args.dir,td_fileIdx=args.td_fileIdx)
  
  # Here the plotting happens. This produces one pdf file
  with PdfPages(outfile) as pdf:
    # Init the prodimo plotting module for one model and an optional title
    # it is required to pass a PdfPages object  
    pp=pplot.Plot(pdf,title=pd.name)
      
    # This plots the Stellar spectrum  
    pp.plot_starspec(pd,ylim=[1.e-1,1.e11])  
    
    if pd.sed is not None: 
      pp.plot_sed(pd,sedObs=pd.sedObs)
      
    pp.plot_dust_opac(pd,ylim=[1.e-1,None])
      
    pp.plot_NH(pd)
    
    # here the default contour plots are used not very fancy currently. 
    # you can use latex for the labels!
    # note most of the parameters are optional
    pp.plot_cont(pd, pd.nHtot, r"$\mathrm{n_{<H>} [cm^{-3}]}$")  
  
    pp.plot_cont(pd, pd.nHtot, r"$\mathrm{n_{<H>} [cm^{-3}]}$",contour=True,zr=False,xlog=False,ylog=False,
                   zlim=[1.e4,None], extend="both")  
    pp.plot_cont(pd, pd.rhod, r"$\mathsf{\rho_{dust} [g\;cm^{-3}]}$",zr=True,xlog=True,ylog=False,
                    zlim=[1.e-25,None],extend="both")
    pp.plot_cont(pd, pd.nd, r"$\mathrm{n_{dust} [cm^{-3}]}$",contour=True,zr=True,xlog=True,ylog=False,
                    zlim=[1.e-8,None],extend="both")
    
    # Plot radiation field at certain wavelengths (here simply done with the index)
    pp.plot_cont(pd,pd.radFields[:,:,0],zr=False,xlog=False,label="lam="+str(pd.lams[0]))
    pp.plot_cont(pd,pd.radFields[:,:,5],zr=False,xlog=False,label="lam="+str(pd.lams[5]))

    pp.plot_cont(pd, pd.tg, r"$\mathrm{T_{gas} [K]}$",contour=True,
                    zlim=[5,5000],extend="both")

    pp.plot_cont(pd, pd.tg, r"$\mathrm{T_{dust} [K]}$",contour=True,
                    zlim=[5,1500],extend="both")
    
    
  
    # this plots the vertical abundances at different radii for different species
    rs=[1,10,100]
    species=["CO","HCO+","e-","S+"]
    for r in rs:
      pp.plot_abunvert(pd, r, species,ylim=[3.e-15,3.e-3])
      
    # plots the average abundances of the above species
    pp.plot_avgabun(pd, species,ylim=[3.e-15,3.e-3])
    
    # contour plots for a couple of species
    species=["CO","H2O","HCO+","HN2+"]
    for spname in species:
      pp.plot_abuncont(pd, spname, zr=True,xlog=True,ylog=False, zlim=[1.e-4,3.e-15],extend="both")
      
  
    cth2lines=pd.selectLineEstimates("C2H2_H")
    fluxes=[lest.flux for lest in cth2lines]
    imax=numpy.argmax(fluxes)
    print(cth2lines[imax],cth2lines[imax].wl,cth2lines[imax].ident)    
  
    pp.plot_line_origin(pd,[[cth2lines[imax].ident,cth2lines[imax].wl]],
                        pd.nmol[:,:,pd.spnames["C2H2"]],xlim=[None,30],ylim=[None,0.4],
                        zlim=[1.e-3,1.e8],extend="both")
  
  