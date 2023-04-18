# Function definition is here

import os
import numpy							as np
from routines.globals					import DEBUG, BC_CORE, BC_WALL, FILE_TYPE_SOLEDGE2D, FILE_TYPE_SOLEDGE3X

#=================================================
# This routine 
#=================================================

def save_neighbors(Dir, Eirene, FileType=FILE_TYPE_SOLEDGE2D):

	if(DEBUG > 0):	print("save_neighbors")

	try:
		os.mkdir(Dir)
	except OSError:
		pass

	if(FileType == FILE_TYPE_SOLEDGE2D):	prefix = "soledge2D"
	elif(FileType == FILE_TYPE_SOLEDGE3X):	prefix = "raptorX"
	else:
		print("\tERROR: I am sorry invalid FileType=",FileType," Check code!")
		exit()

	if(DEBUG > 1):	print("\tsaving to: " + Dir + prefix + ".neighbors")

	fid = open(Dir + prefix + ".neighbors",'w')
	
	Triangles  = Eirene.Triangles
	nTriangles = len(Triangles)
	fid.write('{:d}\n'.format(nTriangles))
	for n in range(nTriangles):
		fid.write('{:d}\t'.format(n + 1))										#Python index to Matalab/Fortran indexes
		fid.write('{:d}\t'.format(Triangles.neigh1[n] + 1))					#Python index to Matalab/Fortran indexes
		fid.write('{:d}\t'.format(Triangles.typeneigh1[n] + 1))				#Python index to Matalab/Fortran indexes
		fid.write('{:d}\t'.format(Triangles.BC1[n]))
		fid.write('{:d}\t'.format(Triangles.neigh2[n] + 1))
		fid.write('{:d}\t'.format(Triangles.typeneigh2[n] + 1))
		fid.write('{:d}\t'.format(Triangles.BC2[n]))
		fid.write('{:d}\t'.format(Triangles.neigh3[n] + 1 ))
		fid.write('{:d}\t'.format(Triangles.typeneigh3[n] + 1))
		fid.write('{:d}\t'.format(Triangles.BC3[n]))
		fid.write('{:d}\t'.format(0))
		fid.write('{:d}\n'.format(0))

	fid.close()

	if(FileType == FILE_TYPE_SOLEDGE3X):	return

	if(DEBUG > 1):	print("\tsaving to: " + Dir + "wall_sequence_properties")

	fid = open(Dir + "wall_sequence_properties",'w')
	WallTriangles	= Eirene.WallTriangles
	for iSeq in range(len(Eirene.Wall.TriSequences)):
		TypeMat		= Eirene.Wall.TypeMat[iSeq]
		IsPump		= Eirene.Wall.IsPump[iSeq]
		TriSeq		= Eirene.Wall.TriSequences[iSeq]

		data = np.zeros((len(TriSeq),5), dtype='i4')

		data[:,0]	= WallTriangles[TriSeq].ntri + 1				#Python index to Matalab/Fortran indexes
		data[:,1]	= WallTriangles[TriSeq].side + 1				#Python index to Matalab/Fortran indexes
		data[:,2]	= BC_WALL

		iPumps = np.where(IsPump)[0]
		data[iPumps,3] = 0
		data[iPumps,4] = 1

		iMats = np.where(~IsPump)[0]
		data[iMats,3] = TypeMat[iMats] + 1							#Python index to Matalab/Fortran indexes
		data[iMats,4] = 0

		for n in range(len(TriSeq)):
			fid.write('{:d} \t {:d} \t {:d} \t {:d} \t {:d} \n'.format(data[n,0],data[n,1],data[n,2],data[n,3],data[n,4]))

	fid.close()

	if(DEBUG > 0):	print("save_neighbors: Completed")
