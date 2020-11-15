#!/usr/bin/python3
# this script will calculate sensitivites: (1) from local perturbation, (2) using GSA

from SALib.sample import saltelli
from SALib.analyze import sobol
import matplotlib.pyplot as plot
import numpy as np, pandas as pd
from numpy import sin,cos, exp, pi

def Rastigrin(x):
	h = cos(2*x[0]-0.5)/2 + cos(x[0]) + 2
	g = (sin(x[1]) - sin(2*x[1]) + sin(3*x[1]))/3 - ((x[1]**2)/(x[1]+1))*sin(4*x[1])/4
	return g*h
	#A = 10
	#x0 = x[0]**2 - A*cos(2*pi*x[0])
	#x1 = x[1]**2 - A*cos(2*pi*x[1])
	#return A*2 + x0 + x1

print(Rastigrin([0.4,1]))

problem = {'num_vars' : 2,
	    'names' : ['x0','x1'],
	    'bounds' : [[-20, 20],[-20,20]]}
print('generating parameter samples')
param_values = saltelli.sample(problem,10000)
print('parameter matrix is (%d,%d)' % (param_values.shape[0],param_values.shape[1]))

# initialize an output array:
Y = np.zeros([param_values.shape[0]])

for i,X in enumerate(param_values):
	Y[i] = Rastigrin(X)

Si = sobol.analyze(problem,Y)
print('S_x1 :: %0.2f' % Si['ST'][0])
print('S_x2 :: %0.2f' % Si['ST'][1])

# now lets compute using local perturbations
xo,yo=1.0,1.0
x1,y1 = xo+1,yo+1
S_x_loc = np.abs((Rastigrin([xo,yo])-Rastigrin([x1,yo]))/(xo-x1))
S_y_loc = np.abs((Rastigrin([xo,yo])-Rastigrin([xo,y1]))/(yo-y1))
norm = S_x_loc + S_y_loc
S_x_loc = S_x_loc/norm
S_y_loc = S_y_loc/norm

print('S_x1 using local perturbation :: %0.2f' % S_x_loc)
print('S_x2 using local perturbation :: %0.2f' % S_y_loc)

