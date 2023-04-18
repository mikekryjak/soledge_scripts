import numpy as np
from math import sqrt

def find_closest_segment(xSegStart, ySegStart, xSegEnd, ySegEnd, xPts, yPts):
	vx	= xSegEnd - xSegStart
	vy	= ySegEnd - ySegStart
	lenght2 = vx**2+vy**2
	try:
		nSegs = len(lenght2)
	except:
		nSegs = 1

	if(nSegs > 1):
		iTruSeg = np.where(lenght2 > 0)[0]								#Remove collapsed segments
		if(len(iTruSeg) < len(vx)):
			vx		= vx[iTruSeg]
			vy		= vx[iTruSeg]
			lenght2	= lenght2[iTruSeg]

		vx1	= xPts - xSegStart[iTruSeg]
		vy1	= yPts - ySegStart[iTruSeg]
		vx2	= xPts - xSegEnd[iTruSeg]
		vy2	= yPts - ySegEnd[iTruSeg]

	elif(lenght2 == 0.):												# one segment collapsed to a point 
		return 1e20, 0
	else:
		vx1	= xPts - xSegStart
		vy1	= yPts - ySegStart
		vx2	= xPts - xSegEnd
		vy2	= yPts - ySegEnd


	p	= (vx*vx1+vy*vy1)/lenght2
	d	= np.where((p > 0) & (p<1), \
								np.abs(vx*vy1-vy*vx1)/np.sqrt(lenght2), \
								np.where(np.sqrt(vx1**2+vy1**2) < np.sqrt(vx2**2+vy2**2), \
										 np.sqrt(vx1**2+vy1**2),  np.sqrt(vx2**2+vy2**2)))
	i	= np.argmin(d)
	d	= d[i]

	return d, iTruSeg[i]

def get_dist_from_segment(xSegStart, ySegStart, xSegEnd, ySegEnd, xPts, yPts):
	vx		= xSegEnd - xSegStart
	vy		= ySegEnd - ySegStart
	lenght2 = vx**2+vy**2
	try:
		nSegs = len(lenght2)
	except:
		nSegs = 1

	if(nSegs > 1):
		iTruSeg = np.where(lenght2 > 0)[0]
		if(len(iTruSeg) < nSegs):
			vx		= vx[iTruSeg]
			vy		= vx[iTruSeg]
			lenght2	= lenght2[iTruSeg]

		vx1	= xPts - xSegStart[iTruSeg]
		vy1	= yPts - ySegStart[iTruSeg]
		vx2	= xPts - xSegEnd[iTruSeg]
		vy2	= yPts - ySegEnd[iTruSeg]

	elif(lenght2 == 0.):											# One segment collapsed to a point 
		return 1e20
	else:
		vx1	= xPts - xSegStart
		vy1	= yPts - ySegStart
		vx2	= xPts - xSegEnd
		vy2	= yPts - ySegEnd


	p	= (vx*vx1+vy*vy1)/lenght2
	d	= np.where((p > 0) & (p<1), \
								np.abs(vx*vy1-vy*vx1)/np.sqrt(lenght2), \
								np.where(np.sqrt(vx1**2+vy1**2) < np.sqrt(vx2**2+vy2**2), \
										 np.sqrt(vx1**2+vy1**2),  np.sqrt(vx2**2+vy2**2)))

	return d