#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import numpy as np
import matplotlib as mpl
from matplotlib.colors import LogNorm
mpl.rcParams['font.size'] = 15
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Specify command line arguments')
parser.add_argument('-i','--input', help='Input file name',required=False)
args = parser.parse_args()

if __name__ == '__main__':
  if not args.input :
    data = np.load('output.npy').item()
  else:
    data = np.load(args.input+'.npy').item()

  wave_k  = data['wave_k']/np.sqrt(data['bperp'][0]*data['mu'][0]/data['mu'][1])  # x = kd_e
  #plt.plot(x,data['fzeta'].imag*data['va'][0]/np.sqrt(data['mu'][0]/data['mu'][1]),lw=2,c=col[0],linestyle='-',label='img')
  X, Y = np.meshgrid(wave_k, data['theta'])
  gamma = data['fzeta'].imag.clip(0)*1e4
  plt.contourf(X, Y,
          gamma*data['va'][0]/np.sqrt(data['mu'][0]/data['mu'][1]),
          cmap=plt.cm.jet)
  plt.axis([0.04, 0.14, 60, 68])
  plt.xlabel('$kd_e$')
  plt.ylabel(r'$\theta$')
  cbar = plt.colorbar()
  cbar.ax.set_ylabel(r'$\Gamma/\omega_{pe}\times10^4$')
  plt.show()
