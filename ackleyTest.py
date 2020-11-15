#!/usr/bin/python3
# this script will run a newton solver and differential evolution optimizer usinon the atchely problem. The atchley function has a well-defined global minimum, but contains a series of lcal minima -- presumably the newton solver will find alocal minimum, while the DE solver will find the global solution. 


from scipy.optimize import newton,fsolve
from scipy.optimize import differential_evolution
from scipy.optimize import minimize
import numpy as np, matplotlib.pyplot as plot
from numpy import sin, cos, exp, pi,sqrt
from mpl_toolkits.mplot3d import Axes3D 

def Ackley(x):
	return -20*np.exp(-0.2*np.sqrt(0.5*(x[0]**2 + x[1]**2))) - np.exp(0.5*(cos(2*np.pi*x[0]) + cos(2*np.pi*x[1]))) + np.exp(1) + 20

def Rastigrin(x):

	A = 10
	x0 = x[0]**2 - A*cos(2*pi*x[0])
	x1 = x[1]**2 - A*cos(2*pi*(x[1]))
	x2 = x[2]**2 - A*cos(2*pi*(x[2]))
	return A*len(x) + x0 + x1 + x2

x = np.linspace(-4,4)
y = np.linspace(-4,4)
XX,YY = np.meshgrid(x,y)
	
#Z = Rastigrin([XX,YY])

#fig = plot.figure()
#ax = fig.gca(projection='3d')
#surf = ax.plot_surface(XX,YY,Z, cmap='jet')


# now we will move into the optimization scheme:
bounds = [(-2,2),(-2,2),(-2,2)]
de = differential_evolution(Rastigrin,bounds)
print('Differential Evolution')
print(de.x,de.fun)

newt = minimize(Rastigrin,[1.0,1.0,1.0],method = 'CG')
#newt = newton(Rastigrin,[1.0,1.0,1.0],maxiter=10000)
print('Conjugate-Gradient method\n')
print(newt)

plot.show()



