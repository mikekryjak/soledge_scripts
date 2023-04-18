import numpy					as np
from routines.polycenter		import polycenter

def find_core(Config):

	r2D		= Config.r2D						#z2D (same z on the row)
	z2D		= Config.z2D						#z2D (same z on the row)
	Bphi2D	= Config.Bphi2D						#Bphi2D[z,r]
	flux2D	= Config.flux2D
	Wall	= Config.Walls[Config.iExtWalls[0]]
	
#	Find local mimima of flux function
	dPsi_dZ, dPsi_dR = np.gradient(flux2D)
	iZmin, iRmin = np.where((np.diff(np.sign(dPsi_dR), axis=1)[:-1,:] > 0) & (np.diff(np.sign(dPsi_dZ), axis=0)[:,:-1] > 0))

	Rmin = r2D[iZmin, iRmin]
	Zmin = z2D[iZmin, iRmin]

#	Determine barycenter of the main chamber
#	It is assumed that the first wall of type 1 is the main one

#	Compute barycenter
	Rcenter, Zcenter = polycenter(Wall.Rwall, Wall.Zwall, 0)

#	We keep the minimum the closest to the barycenter of the chamber
	dist2 = (Rmin-Rcenter)**2+(Zmin-Zcenter)**2
	imin  = np.argmin(dist2)
	Rcore = Rmin[imin]
	Zcore = Zmin[imin]
	Bcore = Bphi2D[iZmin[imin], iRmin[imin]] 

	return Rcore,Zcore,Bcore