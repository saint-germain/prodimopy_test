from __future__ import print_function
from __future__ import division 
from __future__ import unicode_literals

import numpy 
from astropy import units as u
import os

'''
Holds all the outputput stuff from a ProDiMo Model

Currently only selected data is include

@author: rab
'''
class Data_ProDiMo(object):
  def __init__(self, name):  
    self.name = name
    self.__fpFlineEstimates = None  # The path to the FlineEstimates.out File
    self.nx = None
    self.nz = None 
    self.nspec = None
    self.nheat = None
    self.ncool = None
    self.nlam = None
    self.dust2gas = None
    # 2D quantities 
    self.x = None
    self.z = None  
    self.tg = None
    self.td = None 
    self.nd = None
    self.rho = None
    self.rhod = None   
    self.nHtot = None    
    self.damean = None  # mean dust radius
    self.zetaX = None  # the X-ray ionisation rate at every point
    self.zetaCR = None
    self.zetaSTCR = None
    self.tauX1 = None
    self.tauX5 = None
    self.tauX10 = None
    self.NHver = None
    self.NHrad = None
    self.AVrad = None
    self.AVver = None
    self.dummyH2 = None
    self.spnames = None  # is a dictionary to access the indices for nmol
    self.nmol = None          
    self.cdnmol = None  # vertical columnd densities
        
    self.lineEstimates = None  # all the line estimate results
    self.lines = None  # all the lines from the line Transfer
    self.sed = None  # the SED (only the "real" one)
    self.starSpec = None
    self.dust = None  # dust properites (mainly dust opacities)     
  
  def __str__(self):
    output = "Info ProDiMo.out: \n"
    output += "NX: " + str(self.nx) + " NZ: " + str(self.nz) + " NSPEC: " + str(self.nspec)
    output += " NLAM: " + str(self.nlam) + " NCOOL: " + str(self.ncool) + " NHEAT: " + str(self.nheat)
    output += "\n"
    output += "dust2gas: " + str(self.dust2gas)
    return output   
  
  
  
  def selectLineEstimates(self, ident):
    '''
    Returns a list with all lineEstimates with given ident
    '''
    lines = list()     
    for le in self.lineEstimates: 
      if le.ident == ident:
        lines.append(le)
  
    return lines
  
  
  def getLine(self,wl):
    '''
    Finds the line closest to the given wl
    '''
    if self.lines == None: return None
    wls=numpy.array([line.wl for line in self.lines])    
    idx=numpy.argmin(abs(wls[:]-wl))
    return self.lines[idx]
    
  
  def getLineEstimate(self, ident, wl):
    '''
    Finds one particular line from lineEstimates array
    is probably slow
    '''
    found = -1
    i = 0 
    mindiff = 1.e20     
    # print ident, wl
    for le in self.lineEstimates: 
      if le.ident == ident:
        diff = abs(le.wl - wl)
        if diff < mindiff:
          found = i
          mindiff = diff      
      i = i + 1   
  
    if self.lineEstimates[found].rInfo == None:
      # TODO: read the radial data if required
      _read_lineEstimateRinfo(self, self.lineEstimates[found])
  
    return self.lineEstimates[found]
   

class DataLineProfile():
  
  def __init__(self, nvelo):
    '''
    Constructor
    '''
    self.nvelo = nvelo
  
    self.velo = numpy.zeros(shape=(self.nvelo))  # *u.km/u.s
    self.flux = numpy.zeros(shape=(self.nvelo))  # *u.Jansky
    self.flux_conv = numpy.zeros(shape=(self.nvelo))  # *u.Jansky
    
  def __str__(self):
    return "nvelo: " + str(self.nvelo)
        

class DataLine(object):
  '''
  holds the data for one line
  '''
  # 
  def __init__(self):
    '''
    Constructor
    '''
    self.wl = 0.0
    self.frequency = 0.0
    self.prodimoInf = ""
    self.species = ""
    self.ident = ""
    # the flux in [W/m^2]
    self.flux = 0.0
    # the continuum in Jansky
    self.fcont = 0.0 
    # line profile of type data_line_profile
    self.profile = None 
  
  def __str__(self):
    text = (self.species + "/" + 
          self.prodimoInf + "/" + 
          self.ident + "/" + 
          str(self.wl) + "/" + 
          str(self.frequency) + "/" + 
          str(self.flux) + "/" + 
          str(self.fcont))
    return text 


class DataLineObs(DataLine):
  '''
  Holds the observational data for one line
  '''  
  def __init__(self, flux, flux_err, fwhm, fwhm_err, flag):
    super(DataLineObs, self).__init__()
    self.flux = flux * u.W / u.m ** 2
    self.flux_err = flux_err * u.W / u.m ** 2
    self.fwhm = fwhm * u.km / u.s
    self.fwhm_err = fwhm_err * u.km / u.s
    self.flag = flag


class DataLineEstimateRInfo(object):
  '''
  Holds the data for the radial information for one line and on radial position
  '''
  def __init__(self, iz, Fcolumn, tauLine, tauDust, z15, z85):    
    self.iz = iz 
    self.Fcolumn = Fcolumn
    self.tauLine = tauLine
    self.tauDust = tauDust
    self.z15 = z15          
    self.z85 = z85     

class DataLineEstimate(object):
  '''
  Holds the data for one line in FlineEstimates
  '''
  def __init__(self, ident, wl, jup, jlow, flux):
    self.ident = ident 
    self.wl = wl
    self.frequency = None
    self.jup = jup
    self.jlow = jlow
    self.flux = flux
    self.__posrInfo = None  # stores the position of the radial information to access is later on
    self.rInfo = None  # The radial information holds class DataLineEstimateRinfo 
    
  def __str__(self):
    text = (self.ident + "/" + 
          str(self.wl) + "/" + 
          str(self.jup) + "/" + 
          str(self.jlow) + "/" + 
          str(self.flux))
    return text 

class DataGas(object):
  '''
  Holds the data for the gas (mainly from dust_opac.out)
  TODO: currently only the init cs a read, not the x,y data
  '''
  def __init__(self, nlam):
    self.nlam = nlam
    self.lam = numpy.zeros(shape=(nlam))
    self.energy = numpy.zeros(shape=(nlam))  # for convinence as often energy is also used for gas [eV]
    self.abs_cs = numpy.zeros(shape=(nlam))  # unit cm^2/H
    self.sca_cs = numpy.zeros(shape=(nlam))  # unit cm^2/H
  

class DataDust(object):
  '''
  Holds the data for the dust (mainly from dust_opac.out)
  TODO: dust composition is not yet read
  '''
  def __init__(self, amin, amax, apow, nsize, nlam):
    self.amin = amin
    self.amax = amax
    self.apow = apow
    self.nsize = nsize
    self.lam = numpy.zeros(shape=(nlam))
    self.kext = numpy.zeros(shape=(nlam))  # in cm^2 g^-1
    self.kabs = numpy.zeros(shape=(nlam))
    self.ksca = numpy.zeros(shape=(nlam))  
    self.ksca_an = numpy.zeros(shape=(nlam))  # tweaked anisotropic scattering 
    self.kextcs = numpy.zeros(shape=(nlam))  # cross setion cm^2
    self.kabscs = numpy.zeros(shape=(nlam))  # cross setion cm^2
    self.kscacs = numpy.zeros(shape=(nlam))  # cross setion cm^2
    self.kscacs_an = numpy.zeros(shape=(nlam))  # tweaked anisotropic scattering cs   

class DataSED(object):
  '''
  Holds the data for the SED  
  '''
  def __init__(self, nlam, distance, inclination):
    self.nlam = nlam
    self.distance = distance
    self.inclination = inclination
    self.lam = numpy.zeros(shape=(nlam))
    self.nu = numpy.zeros(shape=(nlam))
    self.fnuErg = numpy.zeros(shape=(nlam))
    self.fnuJy = numpy.zeros(shape=(nlam))
    self.nuFnu = numpy.zeros(shape=(nlam))

class DataStarSpec(object):
  '''
  Stellar input spectrum  
  '''
  def __init__(self, nlam, teff,r,logg,luv):
    self.nlam = nlam
    self.teff=teff
    self.r=r
    self.logg=logg
    self.luv=luv
    self.lam = numpy.zeros(shape=(nlam))
    self.nu = numpy.zeros(shape=(nlam))    
    self.Inu = numpy.zeros(shape=(nlam))    

def read_lineEstimates(directory, pdata, filename="FlineEstimates.out"):
  '''
  Read FlineEstimates.out Can only be done after ProDiMo.out is read
  '''
  # do not read the whole file as it is huge
  try: 
    # nedd binary read to be compatible to python 3
    f = open(directory + "/" + filename, 'rb')
    pdata.__fpFlineEstimates = directory + "/" + filename   
  except:
    print(("WARN: Could not open " + directory + "/" + filename + "!"))
    pdata.lineEstimates = None
    return
  
  # check for version currently just check if it exist
  line = f.readline().decode()  
  version = 1
  if len(line) > 68:
    # position of version 
    posver = line.find("version")
    if posver >= 0:
      version = int(line[posver + len("version"):].strip())       
    
  nlines = int(f.readline().strip())
  f.readline()
  
  pdata.lineEstimates = list()  
  nxrange = list(range(pdata.nx))
  startBytes = 0  
  for i in range(nlines):
    # has to be done in fixed format
    # FIXME: it can be that nline is not really equal the nubmer of available line
    # this ir probably a But in ProDiMo but maybe intended (see 
    # OUTPUT_FLINE_ESTIMATE in ProDiMo, Therefore add a check     
    line = f.readline()            
    if not line: break
    
    line=line.decode()    
    # FIXME: has now also a version tag!! this is for the new version
    if version == 1:
      le = DataLineEstimate((line[0:9]).strip(), float(line[10:28]), int(line[29:32]), int(line[33:37]), float(line[38:49]))
    elif version == 2:
      try:
        le = DataLineEstimate((line[0:9]).strip(), float(line[10:28]), int(line[29:34]), int(line[35:40]), float(line[41:53]))
      except ValueError as err:
        print("Conversion problems: {0}".format(err))
        print("Line (i,text): ",i,line)
        print("Field: ",line[0:9],line[10:28],line[29:34],line[35:40],line[41:53])
        raise err
    else:
      raise ValueError('Unknown version of FlineEstimates.out! version=' + str(version))          
            
    le.frequency = (le.wl * u.micrometer).to(u.GHz, equivalencies=u.spectral()).value
          
    # Fine out the number of bytes in one radial line to use seek in the 
    # next iterations    
    if i == 0:
      start = f.tell()
      le.__posrInfo = start  # store the position for reading this info if required
      for j in nxrange:      
        f.readline()
      startBytes = f.tell() - start
    else:
      le.__posrInfo = f.tell()
      f.seek(startBytes, 1)

    pdata.lineEstimates.append(le)

  f.close()

def _read_lineEstimateRinfo(pdata, lineEstimate):
  '''
  Reads the additional Rinfo data for the given lineEstimate  
  '''  
  try:
    f = open(pdata.__fpFlineEstimates, 'rb')
  except: 
    print(("WARN: Could not open " + pdata.__fpFlineEstimates))     
    return None
  
  f.seek(lineEstimate.__posrInfo, 0)
  nxrange = list(range(pdata.nx))
  lineEstimate.rInfo = list()
  for j in nxrange:            
    fieldsR = f.readline().decode().split()     
    iz = int(fieldsR[1].strip())
    Fcolumn = float(fieldsR[2])
    try:
      tauLine = float(fieldsR[3])
    except ValueError as e:
      # print "read_lineEstimates line: ", le.ident  ," error: ", e
      tauLine = 0.0
    tauDust = float(fieldsR[4])
    z15 = float(fieldsR[5])
    z85 = float(fieldsR[6])               
    rInfo = DataLineEstimateRInfo(iz, Fcolumn, tauLine, tauDust, z15, z85)
    lineEstimate.rInfo.append(rInfo)
      
  f.close()
  
def  read_linefluxes(directory, filename="line_flux.out"):      
  '''
  Reads the lines_fluxes.out and puts the data into Data_line
  '''  
  try:
    f = open(directory + "/" + filename, 'r')
  except: 
    print(("WARN: Could not open " + directory + "/" + filename + "!"))     
    return None
    
  
  records = f.readlines()
  f.close()
  
  dist = float(records[0].split("=")[1])
  incl = float(records[1].split("=")[1])
  nlines = int(records[2].split("=")[1])
  nvelo = int(records[3].split("=")[1])
  
  # number of records in the file 
  nrecords = len(records)
  # print dist,incl,nlines,nvelo,nrecords
  
  
  lines = list()
  
  # loop over all lines
  pos = 5
  for i in range(nlines):    
    # read the data for one line      
    line = DataLine()
    rec = records[pos]        
    line.species = (rec[10:20]).strip()
    line.prodimoInf = rec[21:36].strip()
    line.wl = float(rec[43:54].strip())  # *u.um
    line.frequency = float(rec[63:76].strip())  # *u.GHz
           
    rec = records[pos + 1]
    line.flux = float(rec.split("=")[1].strip())  # *u.Watt/u.m**2       
                    
    rec = records[pos + 2]
    line.fcont = float(rec.split("=")[1].strip())  # *u.Jansky                    
             
    # one line in the profile is for the continuum                        
    line.profile = DataLineProfile(nvelo - 1)
    
    for j in range(nvelo - 1):
      # skip the header line and the first line (continuum)?
      fields = records[pos + 5 + j].split()        
      line.profile.velo[j] = float(fields[0])  # *u.km/u.s
      line.profile.flux[j] = float(fields[1])  # *u.Jansky
      if (len(fields) > 2):
        line.profile.flux_conv[j] = float(fields[2])  # *u.Jansky
        
    lines.append(line)
    pos += (nvelo + 5)
  
  return lines

def read_lineObs(directory, nlines, filename="LINEobs.dat"):
  '''
  Reads the lineobs Data. the number of lines have to be known
  this is for version 1 LINEObs.dat
  '''
  try:
    f = open(directory + "/" + filename, 'r')
  except: 
    print(("WARN: Could not open " + directory + "/" + filename + "!"))     
    return None
    
  records = f.readlines()
  f.close()
  
  linesobs = list()  
  for rec in records[2:2 + nlines]:  #        
    fields = rec.split()
    
    lineobs = DataLineObs(float(fields[0].strip()), \
                        float(fields[1].strip()), \
                        float(fields[2].strip()), \
                        float(fields[3].strip()), \
                        fields[4].strip())    
    linesobs.append(lineobs)
    
  return linesobs
  
  
def read_lineObsProfile(filename):
  '''
  reads a line profile file which can be used for ProDiMo
  '''
  f = open(filename, 'r')
  records = f.readlines()
  f.close()
  
  # First line number of velo points and central wavelength, which we do not
  # need here (I think this is anyway optional  
  nvelo = int(records[0].split()[0].strip())

  # skip header lines
  profile = DataLineProfile(nvelo)
  for i in range(2, nvelo + 2):
    fields = records[i].split()
    profile.velo[i - 2] = float(fields[0].strip())
    profile.flux[i - 2] = float(fields[1].strip())
    # also fill the convolved flux, which is just the same as the flux
    # for observations
    profile.flux_conv[i - 2] = profile.flux[i - 2]
    
  return profile

def read_gas(directory):
  '''
  Reads gas_cs.out
  Returns an object of Type DataDust  
  ''' 
  try:
    f = open(directory + "/gas_cs.out", 'r')
  except:
    print(("WARN: Could not open " + directory + "/gas_cs.out" + "!"))
    return None
        
  nlam = int(f.readline().strip())  
  
  # skip one line  
  f.readline()  
      
  gas = DataGas(nlam)
      
  for i in range(nlam):
    fields = [float(field) for field in f.readline().split()]    
    gas.lam[i] = float(fields[0])
    gas.abs_cs[i] = float(fields[1])
    gas.sca_cs[i] = float(fields[2])
 
  gas.energy[:] = ((gas.lam[:] * u.micron).to(u.eV, equivalencies=u.spectral())).value  
 
  f.close()

  return gas  

def read_sed(directory):
  ''' 
  Reads SED.out
  '''
  rfile = directory + "/SED.out"    
  try:
    f = open(rfile, 'r')
  except:
    print("WARN: Could not read " + rfile + "!")
    return None  
  
  
  elems = f.readline().split()
  distance = float(elems[len(elems) - 1])
  # ignore the face-on SED
  f.readline()
  nlam = int(f.readline())
  f.readline()
  for i in range(nlam):
    f.readline()
    
  f.readline()
  elems = f.readline().split()  
  inclination = float(elems[len(elems) - 1])
  nlam = int(f.readline())
  f.readline()
  sed = DataSED(nlam, distance, inclination)
  for i in range(nlam):
    elems = f.readline().split()
    sed.lam[i] = float(elems[0])
    sed.nu[i] = float(elems[1])
    sed.fnuErg[i] = float(elems[2])
    sed.nuFnu[i] = float(elems[3])
    sed.fnuJy[i] = float(elems[4])
    
  return sed      

def read_starSpec(directory):
  ''' 
  Reads StarSpectrum.out
  '''
  rfile = directory + "/StarSpectrum.out"    
  try:
    f = open(rfile, 'r')
  except:
    print("WARN: Could not read " + rfile + "!")
    return None  
      
  teff = float((f.readline().split())[-1])
  elems=f.readline().split()  
  r=float(elems[-1])
  logg=float(elems[2])
  luv=float(f.readline().split()[-1])        
  nlam = int(f.readline())
  f.readline()

  starSpec = DataStarSpec(nlam,teff,r,logg,luv )
  for i in range(nlam):
    elems = f.readline().split()
    starSpec.lam[i] = float(elems[0])
    starSpec.nu[i] = float(elems[1])
    starSpec.Inu[i] = float(elems[2])
    
  return starSpec 

def read_dust(directory):
  '''
  Reads dust_opac.out
  Returns an object of Type DataDust
  Dust not read the dust composition yet
  ''' 
  try:
    f = open(directory + "/dust_opac.out", 'r')
  except:
    print(("WARN: Could not open " + directory + "/dust_opac.out" + "!"))
    return None
        
  fields = [int(field) for field in f.readline().split()]
  ndsp = fields[0]      
  nlam = fields[1]  
  
  # skip three lines
  for i in range(ndsp):
    f.readline()
  
  # apow amax etc.
  strings = f.readline().split()  
  
  amin = ((float(strings[6]) * u.cm).to(u.micron)).value
  amax = ((float(strings[7]) * u.cm).to(u.micron)).value
  apow = float(strings[8])
  nsize = int(strings[9])
  
  dust = DataDust(amin, amax, apow, nsize, nlam)
  
  f.readline()
      
  for i in range(nlam):
    fields = [float(field) for field in f.readline().split()]    
    dust.lam[i] = float(fields[0])
    dust.kext[i] = float(fields[1])
    dust.kabs[i] = float(fields[2])
    dust.ksca[i] = float(fields[3])    
    dust.ksca_an[i] = float(fields[5])  # skip kprn
    dust.kextcs[i] = float(fields[6])
    dust.kabscs[i] = float(fields[7])
    dust.kscacs[i] = float(fields[7])
    dust.kscacs_an[i] = float(fields[9])

  f.close()

  return dust  

def read_prodimo(directory, name=None, readlineEstimates=True, filename="ProDiMo.out", filenameLineEstimates="FlineEstimates.out"):
  '''
  Reads the ProDiMo model.
  Reads ProDiMo.out FlineEstimates.out and dust_opac.out
  '''
  
  if name == None:
    dirfields = directory.split("/")
    name = dirfields[len(dirfields) - 1]  
    
  pfilename = directory + "/" + filename
  f = open(pfilename, 'r')
  print("READ: Reading File: ", pfilename, " ...")
  # read all date into the memory
  # easier to handle afterwards
  lines = f.readlines()  
  f.close()
  idata = 24
   
  data = Data_ProDiMo(name)
  
  data.dust2gas = float(lines[9].split()[1])
  
  strs = lines[21].split()  
  data.nx = int(strs[1])
  data.nz = int(strs[2])
  data.nspec = int(strs[3]) + 1  # +1 for the electron
  data.nheat = int(strs[4])
  data.ncool = int(strs[5])
  data.nlam = int(strs[6])
  
  # TODO: move this to the constructure which takes nx,nz
  data.x = numpy.zeros(shape=(data.nx, data.nz))
  data.z = numpy.zeros(shape=(data.nx, data.nz))
  data.AVrad = numpy.zeros(shape=(data.nx, data.nz))
  data.AVver = numpy.zeros(shape=(data.nx, data.nz))  
  data.NHver = numpy.zeros(shape=(data.nx, data.nz))
  data.NHrad = numpy.zeros(shape=(data.nx, data.nz))
  data.tg = numpy.zeros(shape=(data.nx, data.nz))
  data.td = numpy.zeros(shape=(data.nx, data.nz))
  data.nd = numpy.zeros(shape=(data.nx, data.nz))
  data.rho = numpy.zeros(shape=(data.nx, data.nz))
  data.nHtot = numpy.zeros(shape=(data.nx, data.nz))
  data.damean = numpy.zeros(shape=(data.nx, data.nz))
  data.tauX1 = numpy.zeros(shape=(data.nx, data.nz))
  data.tauX5 = numpy.zeros(shape=(data.nx, data.nz))
  data.tauX10 = numpy.zeros(shape=(data.nx, data.nz))  
  data.zetaX = numpy.zeros(shape=(data.nx, data.nz))
  data.dummyH2 = numpy.zeros(shape=(data.nx, data.nz))
  data.zetaCR = numpy.zeros(shape=(data.nx, data.nz))
  data.zetaSTCR = numpy.zeros(shape=(data.nx, data.nz))
  data.nmol = numpy.zeros(shape=(data.nx, data.nz, data.nspec))
  data.cdnmol = numpy.zeros(shape=(data.nx, data.nz, data.nspec)) 
    
    # number of fixed fields in ProDiMo.out (before heating and cooling rates)
  nfixFields = 21  
  # index after the J data/fields in ProDiMo 
  iAJJ = nfixFields + data.nheat + data.ncool + 1 + data.nspec + data.nlam
  iACool = nfixFields + data.nheat + data.ncool  
  
  # read the species names, these are taken from the column titles
  colnames = lines[idata - 1]
  specStart = 233 + data.nheat * 13 + data.ncool * 13 + 13  
  spnames = colnames[specStart:specStart + (data.nspec) * 13]
  spnames = spnames.split()
  if (len(spnames) != data.nspec): 
    print("ERROR: something is wrong with the number of Species!")
    return None  
  
  # empty dictionary
  data.spnames = {}
  # Make a dictionary for the spans
  for i in range(data.nspec):    
    data.spnames[spnames[i]] = i  
  
  i = 0            
  for iz in range(data.nz):
    for ix in range(data.nx): 
      fields = lines[idata + i].split()      
      zidx = data.nz - iz - 1            
      data.x[ix, zidx] = float(fields[2])
      data.z[ix, zidx] = float(fields[3])
      data.NHrad[ix, zidx] = float(fields[4])
      data.NHver[ix, zidx] = float(fields[5])
      data.AVrad[ix, zidx] = float(fields[6])
      data.AVver[ix, zidx] = float(fields[7])      
      data.nd[ix, zidx] = float(fields[8])
      data.tg[ix, zidx] = float(fields[9])
      data.td[ix, zidx] = float(fields[10])
      data.rho[ix, zidx] = float(fields[12])
      data.nHtot[ix, zidx] = float(fields[iACool])
      data.damean[ix, zidx] = float(fields[iAJJ + 3])
      data.tauX1[ix, zidx] = float(fields[iAJJ + 6])
      data.tauX5[ix, zidx] = float(fields[iAJJ + 7])
      data.tauX10[ix, zidx] = float(fields[iAJJ + 8])      
      data.zetaX[ix, zidx] = float(fields[iAJJ + 9])
      data.zetaCR[ix, zidx] = float(fields[iAJJ + 16])      
      if len(fields) > (iAJJ + 17): data.zetaSTCR[ix, zidx] = float(fields[iAJJ + 17]) 
      data.dummyH2[ix, zidx] = float(fields[iACool + 5])      
      data.nmol[ix, zidx, :] = numpy.array(list(map(float, fields[iACool + 1:iACool + 1 + data.nspec])))
      
      i = i + 1

  data.rhod = data.rho * data.dust2gas

  # Read FlineEstimates.out
  if readlineEstimates == True:
    print(("READ: " + directory + "/" + filenameLineEstimates))
    read_lineEstimates(directory, data, filename=filenameLineEstimates)
  else:
    data.lineEstimates = False  
   
  print("READ: " + directory + "/dust_opac.out")
  data.dust = read_dust(directory)

  print("READ: " + directory + "/StarSpectrum.out")
  data.starSpec = read_starSpec(directory)

  if os.path.exists(directory + "/gas_cs.out"):
    print("READ: " + directory + "/gas_cs.out")
    data.gas = read_gas(directory)

  if os.path.exists(directory + "/line_flux.out"):  
    print("READ: " + directory + "/line_flux.out")
    data.lines = read_linefluxes(directory)  

  if os.path.exists(directory + "/SED.out"):
    print("READ: " + directory + "/SED.out")
    data.sed = read_sed(directory)
  
  # calculate the columnd densitis
  print("INFO: Calculate column densities")
  calc_columnd(data)
  print(" ")

  return data

def calc_columnd(data):
  '''
  Calculated the vertical column density for every species at every point 
  in the disk (from top to bottom). Very simple and rough method
  '''    
  data.cdnmol = 0.0 * data.nmol
  for ix in range(data.nx):          
    for iz in range(data.nz - 2, -1, -1):  # from top to bottom        
      dz = (data.z[ix, iz + 1] - data.z[ix, iz])
      dz = dz * u.au.to(u.cm)
      nn = 0.5 * (data.nmol[ix, iz + 1, :] + data.nmol[ix, iz, :])
      data.cdnmol[ix, iz, :] = data.cdnmol[ix, iz + 1, :] + nn * dz

###############################################################################
# For testing
if __name__ == "__main__":
    # import sys
    # import phxpy.io
#     import time
#           
#     pd=read_prodimo("/home/rab/MODELS/XRTPaperNew/TEST_full")      
#     print pd  
#     print pd.nmol[pd.nx-1,0,pd.spnames["N2#"]]   
#     print pd.gas.energy/1000.0 
#     
#     start=time.time()
#     read_lineEstimates("/home/rab/MODELS/XRTPaperNew/TEST_full", pd)
#     print "Time: ",time.time()-start
#     
#     line=pd.getLineEstimate("N2H+", 1000.0)
#     line=pd.getLineEstimate("N2H+", 1000.0)
#     print line
    
    # lines=pd.selectLineEstimates("N2H+")
    # print len(lines)
            
                
  tdir = "../testdata"
  
  data = read_prodimo(tdir)  
  
                      
  linesObs = read_lineObs(tdir, len(data.lines))
  print(data.lines[0])
  print(data.lines[1])
  print(linesObs[0])
   
  profile = read_lineObsProfile(tdir + "/LineProfile_CO_21.dat")
  print(profile)
  print(profile.flux)
  print(profile.velo)
    
        
