import sys
import scipy
import logging
import numpy as np
import scipy.special as sp
from numpy import linalg as LA
from scipy.integrate import quad

def complex_quadrature(intgrand,args,p,wave_k,i):
  data = (args,p,wave_k,i)
  def real_func(x,*args):
    return scipy.real(intgrand(x,*args))
  def imag_func(x,*args):
    return scipy.imag(intgrand(x,*args))
  int_real = quad(real_func, 1, np.inf, args=data)
  int_imag = quad(imag_func, 1, np.inf, args=data)
  return (int_real[0] + 1j*int_imag[0])

def phasecut(ph,ls,le):
  if ph<ls:
    ph += 2*np.pi
  elif ph>le:
    ph -= 2*np.pi
  return ph

def intgranda(u,*arg):
    args, p, wave_k,i= arg
    theta = p['theta'][0]*np.pi/180.
    rp   = np.abs(u+1)
    rm   = np.abs(u-1)
    phip = np.angle(u+1)
    phim = np.angle(u-1)
    phip = phasecut(phip,-np.pi,np.pi)
    phim = phasecut(phim,0,2*np.pi)

    gamma = 1j*(rp*rm)**(-0.5)*np.exp(-1j*(phip+phim)/2.) 
    alphap= np.cos(theta)/p['beta'][0] +
    1j*np.sin(theta)*np.sqrt(1./p['beta'][0]**2-1.)
    hu    = 
    rho = p['mu'][i]*(1.-p['beta'][i]*u*np.cos(theta))*gamma
    return rho

def intgrandb(u,*arg):
    args, p, wave_k= arg
    theta = p['theta'][0]*np.pi/180.
    for i in range(2):
      gamma = 1./np.sqrt(1.- p['beta'][i]**2)
      rho = p['mu'][i]*(1.-p['beta'][i]*u*np.cos(theta))*gamma
    return rho

def intgrandc(u,*arg):
    args, p, wave_k= arg
    theta = p['theta'][0]*np.pi/180.
    for i in range(2):
      gamma = 1./np.sqrt(1.- p['beta'][i]**2)
      rho = p['mu'][i]*(1.-p['beta'][i]*u*np.cos(theta))*gamma
    return rho

def det(z0,*arg):
  args, p, wave_k= arg
  logging.basicConfig(level=args.loglevel or logging.INFO)
  logger  = logging.getLogger(__name__)

  sol = complex_quadrature(intgranda,args,p,wave_k,0)
  A = 2*np.pi*sol
  print ('A=',A)
  disp_det  = sol
  #logger.debug('for omega=',omega,"disp_det = %e+%ei \n",disp_det.real,disp_det.imag)
  return (disp_det.real,disp_det.imag)

if __name__ == '__main__':
  print ('dispersion relation dielectric tensor.')
