#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rc('font', family='sans-serif')

parser = argparse.ArgumentParser(description='Specify command line arguments')
parser.add_argument('-i','--input', help='Input file name',required=False)
args = parser.parse_args()

if __name__ == '__main__':
  if not args.input :
    data = np.load('output.npy').item()
  else:
    data = np.load(args.input).item()
  plt.semilogx(data['wave_k'],data['fzeta'].imag,lw=2,label='Imag')
  plt.semilogx(data['wave_k'],data['fzeta'].real,lw=2,label='Real')
  plt.xlabel('$k$')
  plt.ylabel('$\gamma$')
  #plt.ylim([1e-6,1])
  plt.rcParams['font.size'] = 16
  #plt.savefig('sample.pdf',format='pdf')
  plt.legend()
  plt.show()
