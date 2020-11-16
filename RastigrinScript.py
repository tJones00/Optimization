#
import numpy as np, matplotlib.pyplot as plot, pandas as pd
from Rastigrin import *
from scipy.optimize import minimize, differential_evolution
from SALib.sample import saltelli
from SALib.analyze import sobol

# step 1: calculate sensitivities (run this about 1,000 times)

# sobol's method

X_loc = np.random.uniform(low = -5, high = 5, size = 10)
Y_loc = np.random.uniform(low = -5, high = 5, size = 10)

# initialize storage arrays:
S1,S2,x1,x2 = [np.zeros((len(X_loc),)) for _ in range(4)]
M_Global,M_local = [np.zeros((len(X_loc),2)) for _ in range(2)]

for n in range(len(X_loc)):
	print('Realization %d' % n)
	problem = {'num_vars' : 2, 'names' : ['x','y'],
		   'bounds' : [[-5,5],[-5,5]]}
	#print('generating saltelli samples')
	param_values = saltelli.sample(problem,1000)
	Y = np.zeros([param_values.shape[0]])
	for i,X in enumerate(param_values):	
		Y[i] = Rastigrin(X)	
	Si = sobol.analyze(problem,Y)
	S1[n] = Si['ST'][0]
	S2[n] = Si['ST'][1]
	
	# now lets compute a local derivative:
	xi,yi = X_loc[n],Y_loc[n]
	x1[n],x2[n] = dRastigrin([xi,yi],np.random.uniform(-1,1,2))
	#print(x1[n]+x2[n])
	
	# lets now compute the function minimum:
	result = differential_evolution(Rastigrin,[(-5,5),(-5,5)],strategy='best1bin')
	local = minimize(Rastigrin,[xi,yi],method='CG')
	M_Global[n,:] = result.x
	M_local[n,:] = local.x
	
print([S1,S2])

plot.figure(1)
plot.plot(S1,'bx')
plot.title('S1')
plot.ylim(0,1)

plot.figure(2)
plot.plot(S2,'bx')
plot.title('S2')
plot.ylim(0,1)

plot.figure(3)
plot.plot(x1,'rx')
plot.title('Sensitivity of X1 using Perturbation Method')
plot.ylim(0,1)

plot.figure(4)
plot.plot(x2,'rx')	
plot.title('Sensitivity of X2 using Perturbation Method')
plot.ylim(0,1)

plot.figure(5)
plot.plot(M_Global[:,0],'gx')
plot.title('X1 at minimum (Global Search)')
plot.ylim(-5,5)

plot.figure(6)
plot.plot(M_Global[:,1],'gx')
plot.title('X2 at minimum (Global Search)')
plot.ylim(-5,5)

plot.figure(7)
plot.plot(M_local[:,0],'kx')
plot.title('X1 at minimum (Local Search)')
plot.ylim(-5,5)

plot.figure(8)
plot.plot(M_local[:,1],'kx')
plot.title('X2 at minimum (Local Search)')
plot.ylim(-5,5)
plot.show()
