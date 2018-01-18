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

''' integrate complex function along real-axis '''
def intgrand(x,*args):
    return (x-0.1)**2/(x-0.1)

if __name__ == '__main__':
    from scipy.integrate import quad
    sol = quad(intgrand,-2,2,points=[0.1])
    print (sol)

