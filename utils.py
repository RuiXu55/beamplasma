import numpy as np
def read_param(args):
  param = {}
  if not args.input :
    args.input = 'input'
  with open(args.input,'r') as f:
    for line in f:
      if not line.isspace():
        data = line.split()
        if data[0] != '<':   
          key = data[0]
          val = float(data[2])
          if key in param.keys():
            param[key].append(val)
          else:
            param[key] = [] 
            param[key].append(val)
  # calculate some params
  param['Omega'] = []
  param['Nsp']   = []
  param['Nsp'].append(np.array(param['mu']).shape[0])
  for i in range(param['Nsp'][0]):
    param['Omega'].append(param['mu'][i]*param['q'][i])

  # now only work for electron beam
  if param['den'][2]>0:
    param['V'].append(-param['V'][1]*param['den'][1]/param['den'][2])
  else:
    param['V'].append(0)

  return param
