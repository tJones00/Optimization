# this script will contain a set of functions for calculating derivatives and function evalutions of the Rastigrin benchmarking function

from numpy import sin, cos, exp, pi, sqrt,abs

def Rastigrin(x):
	A = 10
	y = 0
	for i,X in enumerate(x):
		y = y + X**2 - A*cos(2*pi*X)
	
	return A*len(x) + y

def dRastigrin(x,perturb):
	A = 10.0
	#if len(x) != perturb:
	#	raise Exception('Please define perturbation for each input parameter')
	xo,x1 = x[0],perturb[0]
	yo,y1 = x[1],perturb[1]

	Ro = Rastigrin([xo,yo])
	R1 = Rastigrin([x1,yo])
	R2 = Rastigrin([xo,y1])
	dx1 = xo-x1
	dx2 = yo-y1
	


	#norm = R1 + R2 # total system response
	dR1dx1 = abs(((R2-Ro)/dx1))
	dR2dx2 = abs(((R1-Ro)/dx2))
	
	#dR1dx1 = 2*xo + 2*pi*A*sin(2*pi*xo) 
	#dR2dx2 = 2*x1 + 2*pi*A*sin(2*pi*x1)
	norm = dR1dx1 + dR2dx2

	return dR1dx1/norm, dR2dx2/norm	
