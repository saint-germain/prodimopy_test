.. _sec_plot:

Plotting routines for a single ProDiMo model
============================================

Collections of plotting routines for ProDiMo model output. All the routines use matplotlib.
Typically the output is a pdf file. 

Usage example
-------------

.. code-block:: python

   # the prodimopy modules for reading and plotting
   import prodimopy.read as pread
   import prodimopy.plot as pplot
   # this is for the PDF output
   from matplotlib.backends.backend_pdf import PdfPages

   # read the model from the current directory
   model=pread.read_prodimo()

   # Create an out.pdf where all the various plots are stored
   with PdfPages("out.pdf") as pdf:
      # create a prodimo Plot object for the various plotting routines
      pp=pplot.Plot(pdf)
      # vertical hydrogen number density 
      pp.plot_NH(model,ylim=[1.e20,None])
      # a generic contour plot routine with many options
      pp.plot_cont(model, model.nHtot, r"$\mathsf{n_{<H>} [cm^{-3}]}$",
                   zlim=[1.e4,None],extend="both")


Source documentation
--------------------  
.. automodule:: prodimopy.plot

