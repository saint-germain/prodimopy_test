'''
  Simple script to compare to ProDiMo models. 
  Uses the doAll() method of the compare module.

  This script only provides text output. 
  To call it, simply type: 
    pcompare modeldir1 modeldir2
    
  Type pcompare --help for help. 
'''
from __future__ import print_function
from __future__ import division 
from __future__ import unicode_literals
import argparse


import prodimopy.read as pread
import prodimopy.compare as pcomp

def main(args=None): 
  parser = argparse.ArgumentParser(description='Compares to ProDiMo models')
  parser.add_argument('model1', help='The directory/path of the first model used for the comparison.')
  parser.add_argument('model2', help='The directory/path of the second/reference model used for comparison.')
  parser.add_argument('-tdIdx1', required=False, default=None,help='The index for the time dependent resulst for model1 (e.g. 02). Default: None')
  parser.add_argument('-tdIdx2', required=False, default=None,help='The index for the time dependent results for model2 (e.g. 02). Default: None')
  args = parser.parse_args()
  
  print("Compare model "+args.model1+" to "+args.model2)
  
  model1=pread.read_prodimo(args.model1,readlineEstimates=False,td_fileIdx=args.tdIdx1)
  model2=pread.read_prodimo(args.model2,readlineEstimates=False,td_fileIdx=args.tdIdx1)
  
  compare=pcomp.Compare(model1, model2)
  
  compare.doAll()
