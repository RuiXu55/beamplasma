#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
mpl.rc('font', family='sans-serif')

ang = [60, 63, 65, 68, 70]
col = ['k','r','b','g','c','m']
if __name__ == '__main__':
    for i in range(len(ang)):
        data = np.load('output%2d.npy' %ang[i]).item()
        x = data['wave_k']/np.sqrt(data['bperp'][0]*data['mu'][1])  # x = kd_e
        '''
        plt.plot(x,data['fzeta'].imag/abs(data['Omega'][1]),lw=2,c=col[i],linestyle='-',label='ang=%d'%ang[i])
        plt.xlabel(r'$k c/\omega_{pe}$')
        plt.ylabel(r'$\gamma/\Omega_{ce}$')
        '''
        #plt.plot(x,data['fzeta'].imag*data['va'][0]*np.sqrt(data['mu'][0]/data['mu'][1]),lw=2,c=col[i],linestyle='-',label='ang=%d'%ang[i])
        plt.plot(x,data['fzeta'].imag,lw=2,c=col[i],linestyle='-',label='ang=%d'%ang[i])
        plt.xlabel(r'$k c/\omega_{pe}$')
        plt.ylabel(r'$\gamma/\omega_{ci}$')
        plt.legend()
    plt.show()
