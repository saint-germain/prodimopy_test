from __future__ import print_function
from __future__ import division 

import unittest
import prodimopy.read as pread
from prodimopy.read import Data_ProDiMo, read_linefluxes

# class TestDataProDiMo(unittest.TestCase):
#     def test___init__(self):
#         # data__pro_di_mo = Data_ProDiMo(name)
#         assert False  # TODO: implement your test here
# 
#     def test___str__(self):
#         # data__pro_di_mo = Data_ProDiMo(name)
#         # self.assertEqual(expected, data__pro_di_mo.__str__())
#         assert False  # TODO: implement your test here
# 
#     def test_getLineEstimate(self):
#         # data__pro_di_mo = Data_ProDiMo(name)
#         # self.assertEqual(expected, data__pro_di_mo.getLineEstimate(ident, wl))
#         assert False  # TODO: implement your test here
# 
#     def test_selectLineEstimates(self):
#         # data__pro_di_mo = Data_ProDiMo(name)
#         # self.assertEqual(expected, data__pro_di_mo.selectLineEstimates(ident))
#         assert False  # TODO: implement your test here

# class TestDataLine(unittest.TestCase):
#     def test___init__(self):
#         # data_line = DataLine()
#         assert False  # TODO: implement your test here
# 
#     def test___str__(self):
#         # data_line = DataLine()
#         # self.assertEqual(expected, data_line.__str__())
#         assert False  # TODO: implement your test here
# 
# class TestDataLineObs(unittest.TestCase):
#     def test___init__(self):
#         # data_line_obs = DataLineObs(flux, flux_err, fwhm, fwhm_err, flag)
#         assert False  # TODO: implement your test here
# 
# class TestDataLineEstimateRInfo(unittest.TestCase):
#     def test___init__(self):
#         # data_line_estimate_r_info = DataLineEstimateRInfo(iz, Fcolumn, tauLine, tauDust, z15, z85)
#         assert False  # TODO: implement your test here
# 
# class TestDataLineEstimate(unittest.TestCase):
#     def test___init__(self):
#         # data_line_estimate = DataLineEstimate(ident, wl, jup, jlow, flux)
#         assert False  # TODO: implement your test here
# 
#     def test___str__(self):
#         # data_line_estimate = DataLineEstimate(ident, wl, jup, jlow, flux)
#         # self.assertEqual(expected, data_line_estimate.__str__())
#         assert False  # TODO: implement your test here
# 
# class TestDataGas(unittest.TestCase):
#     def test___init__(self):
#         # data_gas = DataGas(nlam)
#         assert False  # TODO: implement your test here
# 
# class TestDataDust(unittest.TestCase):
#     def test___init__(self):
#         # data_dust = DataDust(amin, amax, apow, nsize, nlam)
#         assert False  # TODO: implement your test here
# 
# class TestReadLineEstimates(unittest.TestCase):
#     def test_read_line_estimates(self):
#         # self.assertEqual(expected, read_lineEstimates(directory, pdata, filename))
#         assert False  # TODO: implement your test here
# 
# class TestReadLinefluxes(unittest.TestCase):
#   def test_read_linefluxes(self):
#     # self.assertEqual(expected, read_linefluxes(directory, filename))
#     assert False  # TODO: implement your test here
# 
# class TestReadLineObs(unittest.TestCase):
#   def test_read_line_obs(self):   
#     assert False  # TODO: implement your test here
# 
# class TestReadLineObsProfile(unittest.TestCase):
#     def test_read_line_obs_profile(self):
#         # self.assertEqual(expected, read_lineObsProfile(filename))
#         assert False  # TODO: implement your test here
# 
# class TestReadGas(unittest.TestCase):
#     def test_read_gas(self):
#         # self.assertEqual(expected, read_gas(directory))
#         assert False  # TODO: implement your test here

class TestReadDust(unittest.TestCase):
  def test_read_dust(self):        
      self.assertNotEqual(pread.read_dust("."), None)

class TestReadSed(unittest.TestCase):
  def test_read_sed(self):
    sed = pread.read_sed(".")
    self.assertNotEqual(sed, None)
    self.assertEqual(sed.distance, 140.0)
    self.assertEqual(sed.inclination, 45.0)
    self.assertEqual(sed.fnuErg[0], 1.723472e-27)

class TestReadStarSpec(unittest.TestCase):
  def test_read_starSpec(self):
    starSpec = pread.read_starSpec(".")
    self.assertNotEqual(starSpec, None)
    self.assertEqual(starSpec.r, 2.0862466372196695)    
    self.assertEqual(starSpec.Inu[0], 8.2394407e-030)
    
class TestGetLine(unittest.TestCase):
  def test_get_line(self):
    # dummy Object
    data =Data_ProDiMo("dummy")
    lines=read_linefluxes(".")
    data.lines=lines
    line=data.getLine(63.0)
    self.assertNotEqual(line, None)
    self.assertEqual(line.species, "OI")
    self.assertEqual(line.wl, 63.18367)
    
    line=data.getLine(867.0)
    self.assertNotEqual(line, None)
    self.assertEqual(line.species, "CO")
    self.assertEqual(line.wl, 866.96325)

class TestReadProdimo(unittest.TestCase):
  def test_read_prodimo(self):
    data = pread.read_prodimo(".")
    self.assertNotEqual(data, None)
    self.assertNotEqual(data.lineEstimates, None)
    self.assertNotEqual(data.dust, None)
    self.assertNotEqual(data.sed, None)
    self.assertNotEqual(data.starSpec, None)
    self.assertEqual(data.x[0, data.nz - 1], 0.07)
            
      

# class TestCalcColumnd(unittest.TestCase):
#     def test_calc_columnd(self):
#         # self.assertEqual(expected, calc_columnd(data))
#         assert False  # TODO: implement your test here

if __name__ == '__main__':
    unittest.main()
