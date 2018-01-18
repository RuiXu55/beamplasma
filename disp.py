import sys
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
  int_real = quad(real_func, -1, 1, args=data)
  int_imag = quad(imag_func, -1, 1, args=data)
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
calculate gu function. See Bret. A, 2010 eqn (7)
'''
def calgu(u,arg):
    args, p, wave_k,omega,i= arg
    theta = p['theta'][0]*np.pi/180.
    rp   = np.abs(u+1)
    rm   = np.abs(u-1)
    phip = phase(u+1,-np.pi,np.pi)
    phim = phase(u-1,0,2*np.pi)
    gamma = 1j*(rp*rm)**(-0.5)*np.exp(-1j*(phip+phim)/2.) 

    alphap= np.cos(theta)/p['beta'][i] + \
            1j*np.sin(theta)*np.sqrt(1./p['beta'][i]**2-1.)
    alpham= np.cos(theta)/p['beta'][i] - \
            1j*np.sin(theta)*np.sqrt(1./p['beta'][i]**2-1.)
    Sp    = np.abs(u-alphap)
    xip   = phase(u-alphap,-1.5*np.pi,np.pi/2.)
    Sm    = np.abs(u-alpham)
    xim   = phase(u-alpham,-.5*np.pi,1.5*np.pi)

    hu    = p['mu'][i]*np.abs(p['beta'][i])*gamma*\
            np.sqrt(Sp*Sm)*np.exp(1j*(xip+xim)/2.) 
    rho   = p['mu'][i]*(1.-p['beta'][i]*u*np.cos(theta))*gamma
    nu    = p['mu'][i]*p['beta'][i]

    gu    = -2*np.pi*gamma*(rho**2+2*nu**2+hu*(rho**2+2*nu**2)+\
            nu**2*hu**2)/hu**5*np.exp(-hu)/wave_k
    return gu

'''
intgrand for Aj. See Appendix A2 in Bret. A, 2010 PRE
'''
def intgranda(u,*arg):
    args, p, wave_k,omega,i= arg
    gu    = calgu(u,arg)
    gzeta = calgu(omega/wave_k,arg)
    return (gu-gzeta)/(u-omega/wave_k) 

'''
calculate det(T) for kinetic beam plasma system
'''
def det(z0,*arg):
  args, p, wave_k= arg
  logging.basicConfig(level=args.loglevel or logging.INFO)
  logger  = logging.getLogger(__name__)

  omega = z0[0] + 1j*z0[1]
  i = 0
  data = (args,p,wave_k,omega,i)
  gzeta = calgu(omega/wave_k,data)
  sol = complex_quadrature(intgranda,data)
  A = sol +  gzeta*np.log((omega/wave_k-1.)/(omega/wave_k+1.))
  print ('gzeta=',gzeta, 'A=',A)

  sys.exit()
  disp_det  = A
  #logger.debug('for omega=',omega,"disp_det = %e+%ei \n",disp_det.real,disp_det.imag)
  return (disp_det.real,disp_det.imag)

if __name__ == '__main__':
  print ('dielectric tensor for beam-plasma system!')
