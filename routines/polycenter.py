import numpy						as np

#	Computes the centroid of a polygon
#	If "filled" is non zero, returns the centroid of the surface delimited by the polygon
#	If "filled" is 0, returns the centroid of the line defined by the polygon

def polycenter(Rpoly, Zpoly, filled=0):

	if(filled != 0):														# Mean of points
		Rmean = np.mean(Rpoly)
		Zmean = np.mean(Zpoly)

		surfTri = 0.5*(Rpoly[:-1]*Zpoly[1:]-Zpoly[:-1]*Rpoly[1:])				# Surface of sub-triangles

		RTri = (Rpoly[:-1]+Rpoly[1:])/3.										# Centroid of sub-triangles
		ZTri = (Rpoly[:-1]+Rpoly[1:])/3.

		surf	= np.sum(surfTri)												# Centroid of filled polygon
		Rcenter = np.sum(Rtri*surfTri)/surf
		Zcenter = np.sum(Ztri*surfTri)/surf
	else:																		# Centers of each segment
		Rseg = 0.5*((Rpoly[:-1]+Rpoly[1:]))
		Zseg = 0.5*((Zpoly[:-1]+Zpoly[1:]))

		lenSeg = np.sqrt((Rpoly[1:]-Rpoly[:-1])**2+(Zpoly[1:]-Zpoly[:-1])**2)	# Length of each segment

		len = np.sum(lenSeg)													# Centroid of empty polygon
		Rcenter = np.sum(Rseg*lenSeg)/len
		Zcenter = np.sum(Zseg*lenSeg)/len

	return Rcenter,Zcenter