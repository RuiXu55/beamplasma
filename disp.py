import sys
import math
import mpmath as mp
import scipy as sp
from scipy.integrate import quad
import logging
import numpy as np

'''
integrate complex function over real axis, 
integrate real & imaginary part seperately
'''
def complex_quadrature(intgrand,arg):
  def real_func(x,*args):
    return sp.real(intgrand(x,*args))
  def imag_func(x,*args):
    return sp.imag(intgrand(x,*args))
  sol_real = quad(real_func, -1., 1.,arg)
  sol_imag = quad(imag_func, -1., 1.,arg)
  return (sol_real[0] + 1j*sol_imag[0])

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
calculate gamma =1/sqrt(1-u**2) and h(u) with phase info.
'''
def calgamma(u):
    phip = phase(u+1., -np.pi,np.pi)
    phim = phase(u-1., 0., 2.*np.pi)
    return 1j/np.sqrt(np.abs(u+1.)*np.abs(u-1.))*np.exp(-1j*(phip+phim)/2.) 

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
    args, p, wave_k,omega,i,k= arg
    theta = p['theta'][0]*np.pi/180.
    gamma = calgamma(u)
    hu    = calhu(gamma,theta,u,p,i)

    rho   = p['mu'][i]*(1.-p['beta'][i]*u*np.cos(theta))*gamma
    nu    = p['mu'][i]*p['beta'][i]*np.sin(theta)
    if k == 1:
        gu = 2*np.pi*gamma*(rho**2+2*nu**2+hu*(rho**2+2*nu**2)+\
             nu**2*hu**2)/hu**5*np.exp(-hu)
    elif k == 2:
        gu = 2*np.pi*gamma**3*u**2*(2*rho**2+nu**2+hu*(2*rho**2+nu**2)+\
             rho**2*hu**2)/hu**5*np.exp(-hu)
    elif k == 3:
        gu = -2*np.pi*nu*gamma**2*rho*u*(3.+3*hu+hu**2)/hu**5*np.exp(-hu)
    return gu/(omega-u*wave_k) 


'''
calculate det(T) for kinetic beam plasma system
'''
def det(z0,*arg):
  args, p, wave_k= arg
  logging.basicConfig(level=args.loglevel or logging.INFO)
  logger  = logging.getLogger(__name__)

  theta = p['theta'][0]*np.pi/180.
  omega = z0[0] + 1j*z0[1]
  epxx = 1.
  epzz = 1.
  epxz = 0

  for i in range(2):
    A = complex_quadrature(intgranda,(args,p,wave_k,omega,i,1))
    B = complex_quadrature(intgranda,(args,p,wave_k,omega,i,2))
    C = complex_quadrature(intgranda,(args,p,wave_k,omega,i,3))
    print ('A',A,'B',B,'C',C)
    sys.exit()
    gamj  = 1./np.sqrt(1-p['beta'][i]**2)
    Lam   = p['mu'][i]/(4*np.pi*gamj**2*sp.special.yv(0,p['mu'][i]/gamj))
    epxx += -p['den'][i]/omega**2*(omega-wave_k*np.cos(theta)*\
            p['beta'][i])*p['mu'][i]*Lam*(A*np.cos(theta)**2+\
            B*np.sin(theta)**2+2*C*np.sin(theta)*np.cos(theta))
    epzz += p['den'][i]/omega**2*p['mu'][i]*p['beta'][i]**2-\
            p['den'][i]/omega**2*(omega-wave_k*np.cos(theta)*\
            p['beta'][i])*p['mu'][i]*Lam*(A*np.sin(theta)**2 +\
            B*np.cos(theta)**2 - 2*C*np.sin(theta)*np.cos(theta))
    epxz += -p['den'][i]/omega**2*(omega-wave_k*np.cos(theta)*p['beta'][i])*\
            p['mu'][i]*Lam*((B-A)*np.sin(theta)*np.cos(theta)+\
            C*(np.cos(theta)**2-np.sin(theta)**2))
  det = (omega**2*epxx-(wave_k*np.cos(theta))**2)*(omega**2*epzz-\
        (wave_k*np.sin(theta))**2)-(omega**2*epxz+\
        wave_k**2*np.sin(theta)*np.cos(theta))**2
  return (det.real,det.imag)

def ing(x,args):
    return (np.sqrt(x)+1j)/(x-1-1j)
if __name__ == '__main__':
  print ('dielectric tensor for beam-plasma system!')
  sol = complex_quadrature(ing,(1))
  print (sol)
