import sys
import scipy
import logging
import numpy as np
import scipy.special as sp
from numpy import linalg as LA
from scipy.integrate import quad

def complex_quadrature(intgrand,args,p,wave_k,omega,i):
  data = (args,p,wave_k,omega,i)
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
    args, p, wave_k,omega,i= arg
    theta = p['theta'][0]*np.pi/180.
    zeta = omega/wave_k 
    for i in range(2):
        if i==1:
            u = zeta
        rp   = np.abs(u+1)
        rm   = np.abs(u-1)
        phip = np.angle(u+1)
        phim = np.angle(u-1)
        phip = phasecut(phip,-np.pi,np.pi)
        phim = phasecut(phim,0,2*np.pi)

        gamma = 1j*(rp*rm)**(-0.5)*np.exp(-1j*(phip+phim)/2.) 
        alphap= np.cos(theta)/p['beta'][0] + 1j*np.sin(theta)*np.sqrt(1./p['beta'][0]**2-1.)
        alpham= np.cos(theta)/p['beta'][0] - 1j*np.sin(theta)*np.sqrt(1./p['beta'][0]**2-1.)
        sp    = np.abs(u-alphap)
        xip   = np.angle(u-alphap)
        xip   = phasecut(xip,-1.5*np.pi,np.pi/2.)
        sm    = np.abs(u-alpham)
        xim   = np.angle(u-alpham)
        xim   = phasecut(xim,-.5*np.pi,1.5*np.pi)

        hu    = p['mu'][i]*np.abs(p['beta'][i])*gamma*np.sqrt(sp*sm)*np.exp(1j*(xip+xim)/2.) 
        rho   = p['mu'][i]*(1.-p['beta'][i]*u*np.cos(theta))*gamma
        nu    = p['mu'][i]*p['beta'][i]

        if i==0:
            gu    = -2*np.pi*gamma*(rho**2+2*nu**2+hu*(rho**2+2*nu**2)+\
                    nu**2*hu**2)/hu**5*np.exp(-hu)/wave_k
        else:
            gzeta = -2*np.pi*gamma*(rho**2+2*nu**2+hu*(rho**2+2*nu**2)+\
                    nu**2*hu**2)/hu**5*np.exp(-hu)/wave_k

    return (gu-gzeta)/(u-zeta) + gzeta*np.log((zeta-1)/(zeta+1)) 

def det(z0,*arg):
  args, p, wave_k= arg
  logging.basicConfig(level=args.loglevel or logging.INFO)
  logger  = logging.getLogger(__name__)

  omega = z0[0] + 1j*z0[1]

  A = complex_quadrature(intgranda,args,p,wave_k,omega,0)
  print ('A=',A)
  sys.exit()
  disp_det  = A
  #logger.debug('for omega=',omega,"disp_det = %e+%ei \n",disp_det.real,disp_det.imag)
  return (disp_det.real,disp_det.imag)

if __name__ == '__main__':
  print ('dispersion relation dielectric tensor.')
