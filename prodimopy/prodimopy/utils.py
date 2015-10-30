from __future__ import print_function
from __future__ import division 
from __future__ import unicode_literals

from astropy import units as u

def calc_dustcolumnd(model):
  '''
  Calculated the vertical column density for every species at every point 
  in the disk (from top to bottom). Very simple and rough method.
  
  :param model: the ProDiMo model  
  :return cd: mass column density for the dust at every point               
              
  TODO: maye put this back into the data structure as util function or both=
  '''    
  cd = 0.0 * model.rhod
  for ix in range(model.nx):          
    for iz in range(model.nz - 2, -1, -1):  # from top to bottom        
      dz = (model.z[ix, iz + 1] - model.z[ix, iz])
      dz = dz * u.au.to(u.cm)      
      nn = 0.5 * (model.rhod[ix, iz + 1] + model.rhod[ix, iz])
      cd[ix, iz] = cd[ix, iz + 1] + nn * dz
  
  return cd