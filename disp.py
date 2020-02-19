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
  omega = u0[0]+1j*u0[1]  # omega = omega/Omega_ci
  D = np.identity(3,dtype=complex)
  if(p['theta'][0]>0.001): # oblique propagation
      # nn = kc/omega
      kn = wave_k/omega/p['va'][0]/np.sqrt(p['bperp'][0])
      D[0,0] += -kn**2*np.cos(theta)**2
      D[0,2] += kn**2*np.cos(theta)*np.sin(theta)
      D[1,1] += -kn**2
      D[2,2] += -kn**2*np.sin(theta)**2
      # loop for different species
      for i in range(int(p['Nsp'][0])):
        if p['den'][i]<1e-10: continue  # neglect zero density species/zero mass
        Delta = p['bperp'][i]/p['bpara'][i]  # temperature anisotropy
        Lambda = (wave_k*np.sin(theta))**2/2.0/p['Omega'][i]**2*\
                p['bperp'][i]/p['bperp'][0]/p['mu'][i]
        kV = wave_k*np.cos(theta)*p['V'][i]/np.sqrt(p['bperp'][0])
        kpwp  = wave_k*np.cos(theta)*np.sqrt(p['bpara'][i]/p['bperp'][0]/p['mu'][i]) #k_para*Vs_parall
        kpwperp  = kpwp*np.sqrt(Delta) #k_para*Vs_parall
        omp_omg2 = p['den'][i]/p['mu'][i]/p['va'][0]**2/omega**2  #(omega_ps/omega)^2

        for n in range(-int(p['m'][0]),int(p['m'][0])+1):
          ze = (omega-kV-n*p['Omega'][i])/kpwp
          An = Delta - 1. + (Delta*ze+n*p['Omega'][i]/kpwp)*f.Z(ze)  # An = omega*An
          Bn = omega - kV + (omega-n*p['Omega'][i])*An  # Bn=k_para*omega*Bn/Omega_ci
          Cn = sp.ive(n,Lambda) - f.dIne(n,Lambda)
          Fn = n**2/Lambda*sp.ive(n,Lambda) + 2*Lambda*sp.ive(n,Lambda) \
                  - 2*Lambda*f.dIne(n,Lambda)

          D[0,0] += n**2*sp.ive(n,Lambda)*An/Lambda*omp_omg2
          D[0,1] += -1j*n*Cn*An*omp_omg2
          D[0,2] += n*sp.ive(n,Lambda)/Lambda*Bn*np.tan(theta)/p['Omega'][i]*omp_omg2
          D[1,1] += An*Fn*omp_omg2
          D[1,2] += 1j*Cn*Bn*np.tan(theta)/p['Omega'][i]*omp_omg2
          D[2,2] += 2.*sp.ive(n,Lambda)*Bn*(omega-n*p['Omega'][i])/kpwperp**2*omp_omg2
        D[2,2] += 2*p['den'][i]**2*np.sqrt(p['bperp'][0])*\
                p['V'][i]/(wave_k*np.cos(theta)*p['va'][0]**2*omega*p['bperp'][i])

      D[1,0] = -D[0,1]
      D[2,0] = D[0,2]
      D[2,1] = -D[1,2]

      res = LA.det(D/omega**p['expo'][0])
  else: # parallel propagation
      pol = +1  # +1 for whistler mode, -1 for firehose mode
      res = p['va'][0]**2*omega**2 - wave_k**2/p['bperp'][0]
      for i in range(3):
        if p['den'][i]<1e-10: continue
        delta = 1. - p['bperp'][i]/p['bpara'][i]
        kV = wave_k*p['V'][i]/np.sqrt(p['bperp'][0])
        kpwp  = wave_k*np.sqrt(p['bpara'][i]/p['bperp'][0]/p['mu'][i])
        ze    = (omega-kV+pol*p['Omega'][i])/kpwp
        ze0   = (omega-kV)/kpwp
        res  += (ze0*f.Z(ze) + delta*f.dp(ze,0)/2.)*p['den'][i]/p['mu'][i]
      res *= omega**p['expo'][0]
  return (res.real,res.imag)

if __name__ == '__main__':
  print ('dielectric tensor for beam-plasma system!')
