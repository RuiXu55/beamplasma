#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.size'] = 16
mpl.rc('font', family='sans-serif')

parser = argparse.ArgumentParser(description='Specify command line arguments')
parser.add_argument('-i','--input', help='Input file name',required=False)
args = parser.parse_args()

col = ['k','r','b','g']
if __name__ == '__main__':
  if not args.input :
    data = np.load('output.npy').item()
  else:
    data = np.load(args.input+'.npy').item()
  x = data['wave_k']#/np.sqrt(data['bperp'][0]*data['mu'][0]/data['mu'][2])  # x = kd_e
  fig, ax = plt.subplots(2,1)
  #ax[0].plot(x,data['fzeta'].imag*data['mu'][2],lw=2,c=col[0],linestyle='-',label='imag')
  #ax[1].plot(x,data['fzeta'].real*data['mu'][2],lw=2,c=col[1],linestyle='-',label='real')
  ax[0].plot(x,data['fzeta'].imag,lw=2,c=col[0],linestyle='-',label='imag')
  ax[1].plot(x,data['fzeta'].real,lw=2,c=col[1],linestyle='-',label='real')
  #plt.xlabel(r'$2\pi k c/\omega_{pe}$')
  #ax[1].set_xlabel(r'$ k c/\omega_{pe}$')
  ax[1].set_xlabel(r'$ k \rho_{i}$')
  ax[0].set_ylabel(r'$\gamma/\Omega_{pi}$')
  ax[1].set_ylabel(r'$\gamma/\Omega_{pi}$')
  ax[0].legend()
  ax[1].legend()

  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.tight_layout()
  #plt.savefig('test_gary3.png', format='png', dpi=400, bbox_inches='tight')
  plt.show()
