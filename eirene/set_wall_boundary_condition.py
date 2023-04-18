import numpy as np
from routines.globals					import BC_CORE, DEBUG

#=================================================================================
# This routine fix set eirebe Boundary Condition at wall with materials and pumps
#=================================================================================

def set_wall_boundary_condition(Eirene):

	if(DEBUG > 0):	print("set_wall_boundary_condition")
		
	Eirene.Wall.Pump[0].Number = np.empty((len(Eirene.Wall.Pump)), dtype='i4')
	for k in range(len(Eirene.Wall.Pump)):
			Eirene.Wall.Pump[0].Number[k] = 2 + k									#Python index

	Eirene.Wall.Material[0].Number = np.empty((len(Eirene.Wall.Material)), dtype='i4')
	for k in range(len(Eirene.Wall.Material)):
		Eirene.Wall.Material[0].Number[k] = 2 + len(Eirene.Wall.Pump) + k			#Python index

	TriSequences = Eirene.Wall.TriSequences
	for iWall in range(len(TriSequences)):
		TriSeq  = TriSequences[iWall]
		IsPump  = Eirene.Wall.IsPump[iWall]
		TypeMat = Eirene.Wall.TypeMat[iWall]
		for n in range(len(TriSeq)):
			iTriPump = np.where(IsPump)[0]
				
			if(len(iTriPump) > 0):
				iPump	 = TypeMat[iTriPump]
				nTriPump = Eirene.WallTriangles.ntri[TriSeq[iTriPump]]

				iBC	 = np.where(Eirene.Triangles.BC1[nTriPump] > BC_CORE)[0]
				if(len(iBC) > 0): Eirene.Triangles.BC1[nTriPump[iBC]] = Eirene.Wall.Pump[0].Number[iPump[iBC]]

				iBC	 = np.where(Eirene.Triangles.BC2[nTriPump] > BC_CORE)[0]
				if(len(iBC) > 0): Eirene.Triangles.BC2[nTriPump[iBC]] = Eirene.Wall.Pump[0].Number[iPump[iBC]]

				iBC	 = np.where(Eirene.Triangles.BC3[nTriPump] > BC_CORE)[0]
				if(len(iBC) > 0): Eirene.Triangles.BC3[nTriPump[iBC]] = Eirene.Wall.Pump[0].Number[iPump[iBC]]

			iTriMat = np.where(~IsPump)[0]
				
			if(len(iTriMat) > 0):
				iMat	 = TypeMat[iTriMat]
				nTriMat = Eirene.WallTriangles.ntri[TriSeq[iTriMat]]

				iBC	 = np.where(Eirene.Triangles.BC1[nTriMat] > BC_CORE)[0]
				if(len(iBC) > 0): Eirene.Triangles.BC1[nTriMat[iBC]] = Eirene.Wall.Material[0].Number[iMat[iBC]]

				iBC	 = np.where(Eirene.Triangles.BC2[nTriMat] > BC_CORE)[0]
				if(len(iBC) > 0): Eirene.Triangles.BC2[nTriMat[iBC]] = Eirene.Wall.Material[0].Number[iMat[iBC]]

				iBC	 = np.where(Eirene.Triangles.BC3[nTriMat] > BC_CORE)[0]
				if(len(iBC) > 0): Eirene.Triangles.BC3[nTriMat[iBC]] = Eirene.Wall.Material[0].Number[iMat[iBC]]

	if(DEBUG > 0):	print("set_wall_boundary_condition: Completed")
