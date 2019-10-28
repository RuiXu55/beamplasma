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

col = ['k','r','b','g']
if __name__ == '__main__':
  if not args.input :
    data = np.load('output.npy').item()
    x = data['wave_k']/np.sqrt(data['bperp'][0])
    plt.plot(x, data['fzeta'].imag,linestyle='-', label='imag')
    plt.plot(x, data['fzeta'].real,linestyle='--', label='real')
  else:
    data = np.load(args.input+'.npy').item()
    x = data['wave_k']/np.sqrt(data['bperp'][0])
    plt.plot(x,data['fzeta'].imag,lw=2,c=col[0],linestyle='-',label='imag')
    plt.plot(x,data['fzeta'].real,lw=2,c=col[0],linestyle='--', label='real')
  plt.xlabel(r'$k \rho_i$')
  plt.ylabel(r'$\gamma/\Omega_{pi}$')
  plt.legend()
  plt.show()
