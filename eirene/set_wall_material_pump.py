import numpy as np
from math 	 import sqrt

def set_wall_material_pump(Eirene, IsPump, IndMatPump, R, Z):
	
#	find closest segments

	nWalls = len(Eirene.Wall.R12)
	dMin   = 1.e30
	for iWall in range(nWalls):
		WallRc = 0.5*(Eirene.Wall.R12[iWall][:,0]+Eirene.Wall.R12[iWall][:,1])
		WallZc = 0.5*(Eirene.Wall.Z12[iWall][:,0]+Eirene.Wall.Z12[iWall][:,1])
		
		iPt = np.empty(3, dtype = 'i4')
		dTot   = 0. 
		for k in range(3):
			iPt[k] = np.argmin((WallRc-R[k])**2 + (WallZc-Z[k])**2)
			dTot += sqrt((WallRc[iPt[k]]-R[k])**2 + (WallZc[iPt[k]]-Z[k])**2)
		if(dTot < dMin):
			dMin = dTot
			iPtSet = np.copy(iPt)
			iWallSet = iWall
	
	iPtMin = min(iPtSet[0], iPtSet[2])
	iPtMax = max(iPtSet[0], iPtSet[2])
	if((iPtMin <=iPtSet[1]) and (iPtMax >=iPtSet[1])):
		Eirene.Wall.TypeMat[iWallSet][iPtMin:iPtMax+1] = IndMatPump
		Eirene.Wall.IsPump[iWallSet][iPtMin:iPtMax+1]  = IsPump
	else:
		Eirene.Wall.TypeMat[iWallSet][iPtMax:]		= IndMatPump
		Eirene.Wall.TypeMat[iWallSet][:iPtMin+1]	= IndMatPump
		Eirene.Wall.IsPump[iWallSet][iPtMax:]		= IsPump
		Eirene.Wall.IsPump[iWallSet][:iPtMin+1]		= IsPump
	Eirene.Wall.iWall  					= iWallSet
