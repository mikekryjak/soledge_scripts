
import numpy as np
import os
from files.save_eirene_triangles		import save_eirene_triangles
from files.save_neighbors				import save_neighbors
from routines.globals					import BC_CORE, DEBUG, FILE_TYPE_SOLEDGE2D, FILE_TYPE_SOLEDGE3X
#=================================================
# This routine write .....
#=================================================

def save_puffs(Dir, Eirene, FileType=FILE_TYPE_SOLEDGE2D):

	if(DEBUG > 0):	print("save_neighbors_knots_puff")
		
	Eirene.Wall.Pump[0].Number = np.empty((len(Eirene.Wall.Pump)), dtype='i4')
	for k in range(len(Eirene.Wall.Pump)):
			Eirene.Wall.Pump[0].Number[k] = 2 + k									#Python index

	Eirene.Wall.Material[0].Number = np.empty((len(Eirene.Wall.Material)), dtype='i4')
	for k in range(len(Eirene.Wall.Material)):
		Eirene.Wall.Material[0].Number[k] = 2 + len(Eirene.Wall.Pump) + k			#Python index

	TriSequences = Eirene.Wall.TriSequences
	for iSeq in range(len(TriSequences)):
		TriSeq   = TriSequences[iSeq]
		IsPump   = Eirene.Wall.IsPump[iSeq]
		TypeMat  = Eirene.Wall.TypeMat[iSeq]
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

	save_neighbors(Dir, Eirene, FileType=FILE_TYPE_SOLEDGE2D)
	save_eirene_triangles(Dir, Eirene, FileType=FILE_TYPE_SOLEDGE2D)

	try:
		os.mkdir(Dir)
	except OSError:
		pass

	if(DEBUG > 1):	print("\tsaving to: "+ Dir + "puffs")

	fid = open(Dir + 'puffs','w')
		
	for n in range(len(Eirene.Wall.Puff)):
		fid.write('{:d}\t'.format(n+1))									#python index
		fid.write('{:s}\t'.format(Eirene.Wall.Puff[n].Name))
		fid.write('{:d}\t'.format(Eirene.Wall.Puff[n].ntri+1))				#python index
		fid.write('{:d}\t'.format(Eirene.Wall.Puff[n].side+1))				#python index
		fid.write('{:f}\t'.format(Eirene.Wall.Puff[n].R))
		fid.write('{:f}\n'.format(Eirene.Wall.Puff[n].Z))

	fid.close()

	if(DEBUG > 0):	print("save_neighbors_knots_puff: Completed")
