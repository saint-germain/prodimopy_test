'''
  Simple script to compare to ProDiMo models. 
  Uses the doAll() method of the compare module.

  This script only provides text output. 
  To call it simply type 
    pcompare model modelref
    
  Type pcompare --help for help. 
'''
from __future__ import print_function
from __future__ import division 
from __future__ import unicode_literals
import argparse


import prodimopy.read as pread
import prodimopy.compare as pcomp

def main(args=None): 
  parser = argparse.ArgumentParser(description='Plot Line observations comparison')
  parser.add_argument('model', help='The directory/path of the particular model used for the comparison')
  parser.add_argument('modelref', help='The directory/path of the reference model used for comparison')
  args = parser.parse_args()
  
  print("Compare model "+args.model+" to "+args.modelref)
  
  model=pread.read_prodimo(args.model,readlineEstimates=False)
  modelref=pread.read_prodimo(args.modelref,readlineEstimates=False)
  
  compare=pcomp.Compare(model, modelref)
  
  compare.doAll()
