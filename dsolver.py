#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Written by Rui Xu, Jan. 2018
# Email: ruix@princeton.edu
import sys
import time
import disp
import utils
import logging 
import argparse
import numpy as np
from scipy import interp
from scipy.optimize import root
__author__ = 'ruix'

def main(args):
  """ parse command line arguments"""
  tstart = time.clock()
  if not args.log :
    logging.basicConfig(level=args.loglevel or logging.INFO)
  else:
    logging.basicConfig(filename='log', filemode='w', level=logging.DEBUG)
  logger = logging.getLogger(__name__)

  """ read plasma parameters """
  param = utils.read_param(args)
  
  """ iterate through wavenumber  """
  dk     = (param['kend'][0]-param['kstart'][0])/(param['ksteps'][0]-1)
  fzeta  = np.empty(int(param['ksteps'][0]),dtype=complex)
  wave_k = np.empty(int(param['ksteps'][0]))
  zeta_guess = complex(param['omega_r'][0],param['omega_i'][0])    

  for n in range(int(param['ksteps'][0])):
    logger.info('%d th iteration in %d ksteps \n' ,n,param['ksteps'][0])
    wave_k[n] = param['kstart'][0]+n*dk
    """ find dispersion relation root """
    data = (args, param, wave_k[n])
    try:
        sol = root(disp.det,(zeta_guess.real,zeta_guess.imag), \
            args=data,method='lm',tol=param['sol_err'][0]) 
        fzeta[n] = complex(sol.x[0],sol.x[1])
        logger.info("solution: zeta= %1.1e +%1.1e i",fzeta[n].real,fzeta[n].imag)
    except ValueError:
      logger.info('ERROR in root finding: wave_k =%f',wave_k[n])
      sys.exit()
    zeta_guess = fzeta[n]

  """ save results to output file """
  param['wave_k'] = wave_k
  param['fzeta'] = fzeta
  if args.output :
    np.save(args.output+'.npy',param)
  else:
    np.save('output.npy', param)

  tend = time.clock()
  logger.info('\n Total time lapse: %1.2f\n',tend-tstart)
  return  0

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Specify command line arguments')
  parser.add_argument('-i','--input', help='Input file name',required=False)
  parser.add_argument('-o','--output',help='Output file name')
  parser.add_argument('-l','--log', action='store',help='log file')
  parser.add_argument('-v', '--verbose', help='Verbose (debug) logging',
          action='store_const', const=logging.DEBUG,dest='loglevel')
  args = parser.parse_args()
  #print (parser.parse_args())
  sys.exit(main(args))

 
 
