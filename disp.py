import sys
import logging
import zetaf as f
import numpy as np
import scipy.special as sp
from numpy import linalg as LA

'''
calculate det(T) for kinetic beam plasma system
'''
def det(u0,*arg):
  args, p, wave_k= arg
  logging.basicConfig(level=args.loglevel or logging.INFO)
  logger= logging.getLogger(__name__)
  theta = p['theta'][0]*np.pi/180.
  omega = u0[0]+1j*u0[1]
    
  D     = np.zeros((3,3),dtype=complex)
  for i in range(int(p['Nsp'][0])):
    al    = (wave_k*np.sin(theta))**2/2.0/p['Omega'][i]**2*\
            p['bperp'][i]/p['bperp'][0]*p['mu'][i]
    delta = p['bperp'][i]/p['bpara'][i]
    kVOmp = wave_k*np.cos(theta)*p['V'][i]/np.sqrt(p['bperp'][0])
    kpwp  = wave_k*np.cos(theta)*np.sqrt(p['bpara'][i]/p['bperp'][0]*p['mu'][i])
    pref  = p['mu'][i]*p['den'][i]/omega
    
    D[2,2]+= pref*2*p['V'][i]/wave_k/np.cos(theta)*\
             np.sqrt(p['bperp'][0]/p['mu'][i])/p['bperp'][i]*\
             np.sqrt(p['den'][i]/p['mu'][i])

    for n in range(-int(p['m'][0]),int(p['m'][0])+1):
      ze   = (omega-kVOmp- n*p['Omega'][i])/kpwp
      An   = (delta-1)/omega+((omega-kVOmp-n*p['Omega'][i])*delta\
             + n*p['Omega'][i])/kpwp/omega*f.Z(ze)
      Bn   = 1.-kVOmp/omega + (omega-n*p['Omega'][i])*An

      D[0,0]+= n**2*sp.ive(n,al)*An/al*pref
      D[0,1]+= 1j*n*f.dive(n,al,1)*An*pref
      D[1,0]+= -1j*n*f.dive(n,al,1)*An*pref
      D[0,2]+= n*sp.ive(n,al)*Bn/np.tan(theta)/al/p['Omega'][i]*pref
      D[2,0]+= n*sp.ive(n,al)*Bn/np.tan(theta)/al/p['Omega'][i]*pref
      D[1,1]+= An*pref*(n**2/al*sp.ive(n,al)-2*al*f.dive(n,al,1))
      D[1,2]+= -1j*f.dive(n,al,1)*Bn/np.tan(theta)/p['Omega'][i]*pref
      D[2,1]+= 1j*f.dive(n,al,1)*Bn/np.tan(theta)/p['Omega'][i]*pref
      D[2,2]+= 2.*sp.ive(n,al)*Bn*(omega-n*p['Omega'][i])/(wave_k**2*\
                np.cos(theta)**2*p['bperp'][i]/p['bperp'][0])/p['mu'][i]*pref
        
  # add (va_c term comes from displacement current)
  rn      = wave_k/np.sqrt(p['bperp'][0])/omega
  idty    = p['va'][0]**2
  D[0,0] += idty - (rn*np.cos(theta))**2
  D[0,2] += rn**2*np.cos(theta)*np.sin(theta)
  D[1,1] += idty - rn**2
  D[2,0] += rn**2*np.cos(theta)*np.sin(theta)
  D[2,2] += idty - (rn*np.sin(theta))**2
  res = LA.det(D)
  return (res.real,res.imag)

if __name__ == '__main__':
  print ('dielectric tensor for beam-plasma system!')
