from tkinter 				import messagebox
import types
import numpy as np
from math 	 				import sqrt
from routines.globals		import BC_CORE, BC_WALL

def set_wall_puff(Eirene, PuffInd, RPuff, ZPuff):
	
	nWalls = len(Eirene.Wall.R12)
	dMin   = 1.e30
	for iWall in range(nWalls):
		WallRc = 0.5*(Eirene.Wall.R12[iWall][:,0]+Eirene.Wall.R12[iWall][:,1])
		WallZc = 0.5*(Eirene.Wall.Z12[iWall][:,0]+Eirene.Wall.Z12[iWall][:,1])

		iPt = np.argmin((WallRc-RPuff)**2 + (WallZc-ZPuff)**2)
		d		= sqrt((WallRc[iPt]-RPuff)**2 + (WallZc[iPt]-ZPuff)**2)
	
		if(d < dMin):
			dMin 	= d
			iPtPuff = iPt
			iWallSet = iWall

	nTri   = Eirene.WallTriangles.ntri[Eirene.Wall.TriSequences[iWallSet][iPtPuff]]

	if  (Eirene.Triangles.BC1[nTri] > BC_CORE): side = 0
	elif(Eirene.Triangles.BC2[nTri] > BC_CORE): side = 1
	elif(Eirene.Triangles.BC3[nTri] > BC_CORE): side = 2
	else:
		print("\t\tnTri=",nTri)
		print("\t\tBC1, BC2, BC3=",Eirene.Triangles.BC1[nTri], Eirene.Triangles.BC2[nTri], Eirene.Triangles.BC3[nTri])
		messagebox.showwarning("Invalid wall side check code")
		return

	if(PuffInd == -1): Eirene.Wall.Puff.append(types.SimpleNamespace())
	Eirene.Wall.Puff[PuffInd].ntri  = nTri
	Eirene.Wall.Puff[PuffInd].side  = side
	Eirene.Wall.Puff[PuffInd].R	    = 0.5*(Eirene.Wall.R12[iWallSet][iPtPuff,0]+Eirene.Wall.R12[iWallSet][iPtPuff,1])
	Eirene.Wall.Puff[PuffInd].Z	    = 0.5*(Eirene.Wall.Z12[iWallSet][iPtPuff,0]+Eirene.Wall.Z12[iWallSet][iPtPuff,1])
	Eirene.Wall.Puff[PuffInd].iWall	= iWallSet

	return