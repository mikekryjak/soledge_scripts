from math import ceil
import numpy as np
from tkinter				import messagebox
try:
	from scipy.misc 			import comb
except:
	from scipy.special 			import comb

def bernstein_poly(i, n, t):

#     The Bernstein polynomial of n, i as a function of t

	return comb(n, i) * ( t**(n-i) ) * (1 - t)**i

def bezier(xPoints, yPoints, nOut=500):

#       Given a set of control points, return the
#       bezier curve defined by the control points.

#       points should be to numpy arrays
#        nOut is the number of time steps, defaults to 1000

#       See http://processingjs.nihongoresources.com/bezierinfo/

	nPoints = len(xPoints)

	t = np.linspace(0.0, 1.0, nOut)

	polynomial_array = np.array([ bernstein_poly(i, nPoints-1, t) for i in range(0, nPoints)   ])

	xVals = np.dot(xPoints, polynomial_array)
	yVals = np.dot(yPoints, polynomial_array)

	return xVals[::-1], yVals[::-1]

def bezier3(xControl, yControl, nOut=500):

# Piece-wise cubic Bezier splines
# "xControl/yControl" is a 2D array of shape [Nctrlpts,3] with the following dimensions:
#   - dim 1 = index of control point triplets (2 minimum)
#   - dim 2 = 0 for left point, 1 for central point, 2 for right point
# nOut is the total number of points wanted in the whole interval of definition of the function
# This number will be automatically increased by the function if 2 control points are too close to each other

# Define empty output by default
	xVals = np.empty((0), dtype='f8')
	yVals = np.empty((0), dtype='f8')

# Check inputs
	nctrl = xControl.shape[0]
	if(nctrl < 2):
		messagebox.showerror("Bezier3", "function called with insufficient number of control points triplets.\nShould have at least 2.\nReturning empty curve.")
		return xVals, yVals

	intervalWidth = xXontrol[-1,1] - xControl[0,1]				# Width of total interval along x coordinate

#	Loop on interval between control triplets
	for ictrl in range(nctrl-1):
		npts = max(2,ceil(nOut*(xControl[ictrl+1,1] - xControl[ictrl,1])/intervalWidth))
		t = np.linspace(0., 1., npts)								# Parametric variable to define cubic Bezier spline

#		Remove first point of t unless this is the last interval
		if(ictrl > 0): t = t[1:]

#		Compute spline interpolation for local interval
		polt1 = (1-t)**3
		polt2 = 3*(1-t)**2.*t
		polt3 = 3*(1-t)*t**2
		polt4 = t**3
		newX = xControl[ictrl,1]*polt1 + xControl[ictrl,2]*polt2 + xControl[ictrl+1,0]*polt3 + xControl[ictrl+1,1]*polt4
		newY = yControl[ictrl,1]*polt1 + yControl[ictrl,2]*polt2 + yControl[ictrl+1,0]*polt3 + yControl[ictrl+1,1]*polt4

#		Agregate
		xVals = np.append(xVals, newX)
		yVals = np.append(yVals, newX)

	return xVals, yVals

