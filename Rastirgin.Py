# this script will contain a set of functions for calculating derivatives and function evalutions of the Rastigrin benchmarking function

from numpy import sin, cos, exp, pi, sqrt

def Rastigrin(x):
	A = 10
	y = 0
	for i,X in enumerate(x):
		y = y + X**2 - A*cos(2*pi*X)
	
	return A*len(x) + y

def dRastigrin(x,perturb):
	if len(x) != perturb:
		raise Exception('Please define perturbation for each input parameter')
	Ro = Rastigrin(x)
	R1 = Rastigrin([x[0],perturb[1]])
	R2 = Rastigrin([perturb[0],x[1]])
	dx1 = x[0]-perturb[0]
	dx2 = x[1]-perturb[1]

	norm = R1 + R2 # total system response
	dR1dx1 = ((R2-Ro)/dx1)/norm
	dR2dx2 = ((R1-Ro)/dx2)/norm
	return dR1dx1, dR2dx2	
