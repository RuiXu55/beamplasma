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
    x = data['wave_k']/np.sqrt(data['bperp'][0]*data['mu'][0]/data['mu'][1])  # x = kd_e
    plt.plot(x,data['fzeta'].imag*data['va'][0]/np.sqrt(data['mu'][0]/data['mu'][1]),lw=2,c=col[0],linestyle='-',label='img')
    #plt.plot(x,data['fzeta'].imag, lw=2)
    plt.xlabel(r'$k c/\omega_{pe}$')
    plt.ylabel(r'$\gamma/\omega_{pe}$')
  else:
    data = np.load(args.input+'.npy').item()
    x = data['wave_k']/np.sqrt(data['bperp'][0]*data['mu'][1])
    plt.plot(x,data['fzeta'].imag/abs(data['Omega'][1]),lw=2,c=col[0],linestyle='-')
    plt.plot(x,data['fzeta'].real/abs(data['Omega'][1]),lw=2,c=col[0],linestyle='--')
    #plt.xlabel(r'$k\rho_i$')
    plt.xlabel(r'$kc/\omega_{pe}$')
    plt.ylabel('$\omega/\Omega_e$')
    plt.rcParams['font.size'] = 16
    #plt.ylim([0,1])
    #plt.savefig('emic.pdf',format='pdf')
    plt.legend()
  plt.show()
