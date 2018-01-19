import sys
import math
import mpmath as mp
import scipy
import logging
import numpy as np
from scipy.integrate import quad

'''
integrate complex function over real axis, 
integrate real & imaginary part seperately
'''
def complex_quadrature(intgrand,data):
  def real_func(x,*args):
    return scipy.real(intgrand(x,*args))
  def imag_func(x,*args):
    return scipy.imag(intgrand(x,*args))
  int_real = quad(real_func, -1., 1., args=data,epsrel=1e-8)
  int_imag = quad(imag_func, -1., 1., args=data,epsrel=1e-8)
  return (int_real[0] + 1j*int_imag[0])

'''
phase constraint for a complex parameter: ls<ph<le
'''
def phase(param,ls,le):
  ph = np.angle(param)
  if ph<ls:
    ph += 2*np.pi
  elif ph>le:
    ph -= 2*np.pi
  return ph

'''
calculate gamma =1/sqrt(1-u**2) with phase info.
'''
def calgamma(u):
    phip = phase(u + 1., -np.pi,np.pi)
    phim = phase(u - 1., 0, 2*np.pi)
    return 1j/np.sqrt(np.abs(u+1)*np.abs(u-1))*np.exp(-1j*(phip+phim)/2.) 

def calhu(gamma,theta,u,p,i):
    alphap= np.cos(theta)/p['beta'][i] + \
            1j*np.sin(theta)*np.sqrt(1./p['beta'][i]**2-1.)
    alpham= np.cos(theta)/p['beta'][i] - \
            1j*np.sin(theta)*np.sqrt(1./p['beta'][i]**2-1.)
    xip   = phase(u-alphap,-1.5*np.pi,0.5*np.pi)
    xim   = phase(u-alpham,-0.5*np.pi,1.5*np.pi)

    hu    = p['mu'][i]*np.abs(p['beta'][i])*gamma*np.sqrt(np.abs(u-alphap)*\
            np.abs(u-alpham))*np.exp(1j*(xip+xim)/2.) 
    return hu

'''
intgrand for Aj. See Appendix A2 in Bret. A, 2010 PRE
'''
def intgranda(u,*arg):
    args, p, wave_k,omega,i= arg
    theta = p['theta'][0]*np.pi/180.
    gamma = calgamma(u)
    hu    = calhu(gamma,theta,u,p,i)

    rho   = p['mu'][i]*(1.-p['beta'][i]*u*np.cos(theta))*gamma
    nu    = p['mu'][i]*p['beta'][i]*np.sin(theta)
    # hu = np.sqrt(rho**2-nu**2)
    gu    = 2*np.pi*gamma*(rho**2+2*nu**2+hu*(rho**2+2*nu**2)+\
            nu**2*hu**2)/hu**5*np.exp(-hu)
    return gu/(omega-u*wave_k) 

'''
calculate det(T) for kinetic beam plasma system
'''
def det(z0,*arg):
  args, p, wave_k= arg
  
  theta = p['theta'][0]*np.pi/180.
  logging.basicConfig(level=args.loglevel or logging.INFO)
  logger  = logging.getLogger(__name__)

  omega = z0[0] + 1j*z0[1]
  epxx = 1.
  epzz = 1.
  epxz = 0

  for i in range(1):
    A = complex_quadrature(intgranda,(args,p,wave_k,omega,i))
    B = 0
    C = 0
    print ('A',A)
    sys.exit()
    gammaj = 1./np.sqrt(1-p['beta'][i]**2)
    Lam   = p['mu'][i]/(4*np.pi*gammaj**2*sp.yv(0,p['mu'][i]/gammaj))
    epxx += -p['den'][i]/omega**2*(omega-wave_k*np.cos(theta)*p['beta'][i])*Lam*\
            (A*np.cos(theta)**2+B*np.sin(theta)**2+2*C*np.sin(theta)*np.cos(theta))
    epzz += -p['den'][i]/omega**2*p['mu'][i]*p['beta'][i]**2-\
            p['den'][i]/omega**2*p['mu'][i]*Lam*(A*np.sin(theta)**2 +\
            B*np.cos(theta)**2 - 2*C*np.sin(theta)*np.cos(theta))
    epxz += -p['den'][i]/omega**2*(omega-wave_k*np.cos(theta)*p['beta'][i])*p['mu'][i]*\
            Lam*((B-A)*np.sin(theta)*np.cos(theta)+C*(np.cos(theta)**2-np.sin(theta)**2))
  det = (omega**2*epxx-(wave_k*np.cos(theta))**2)*(omega**2*epzz-\
        (wave_k*np.sin(theta))**2)-(omega**2*epxz+wave_k**2*np.sin(theta)*np.cos(theta))**2
  return (det.real,det.imag)

if __name__ == '__main__':
  print ('dielectric tensor for beam-plasma system!')
