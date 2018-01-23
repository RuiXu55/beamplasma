#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
mpl.rc('font', family='sans-serif')

parser = argparse.ArgumentParser(description='Specify command line arguments')
parser.add_argument('-i','--input', help='Input file name',required=False)
args = parser.parse_args()

if __name__ == '__main__':
  if not args.input :
    data = np.load('output.npy').item()
  else:
    data = np.load(args.input+'.npy').item()
  x = data['wave_k']/np.sqrt(data['bperp'][0]*data['mu'][1])
  #plt.plot(data['wave_k'],data['fzeta'].imag/abs(data['Omega'][1]),lw=2,c='r',label='Imag')
  #plt.plot(data['wave_k'],data['fzeta'].real/abs(data['Omega'][1]),lw=2,c='k',label='Real')
  plt.plot(data['wave_k'],data['fzeta'].imag,lw=2,c='r',label='Imag')
  plt.plot(data['wave_k'],data['fzeta'].real,lw=2,c='k',label='Real')

  plt.xlabel(r'$k\rho_i$')
  #plt.xlabel(r'$kc/\omega_{pe}$')
  plt.ylabel('$\omega/\Omega_i$')
  plt.rcParams['font.size'] = 16
  #plt.ylim([0,6e-4])
  #plt.savefig('sample.pdf',format='pdf')
  plt.legend()
  plt.show()
