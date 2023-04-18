import numpy as np
import os
from routines.globals import BC_CORE, DEBUG, FILE_TYPE_SOLEDGE2D, FILE_TYPE_SOLEDGE3X

#=================================================
# This routine soledge2D.npco_char file
#=================================================

def save_wall_back_interp(Dir, Eirene, FileType=FILE_TYPE_SOLEDGE2D):

	if(DEBUG > 0):	print("save_wall_back_interp")

	Triangles = Eirene.Triangles
	WallTri   = np.where((Triangles.BC1 == BC_WALL) | (Triangles.BC2 == BC_WALL) | (Triangles.BC3 == BC_WALL))[0]
	nWallTri  = len(WallTri)

	try:
		os.mkdir(Dir)
	except OSError:
		pass

	if(FileType == FILE_TYPE_SOLEDGE2D):	prefix = "soledge2D"
	elif(FileType == FILE_TYPE_SOLEDGE2D):	prefix = "raptorX"
	else:
		print("\tERROR: I am sorry invalid FileType=",FileType," Check code!")
		exit()

	if(DEBUG > 1):	print("\tsaving to: " +Dir + prefix + ".backinterp")

	fid = open(Dir + prefix + ".backinterp",'w')
	fid.write('{:d}\n'.format(nWallTri))
	for n in range(nWallTri):
		fid.write('{:d}\t'.format(WallTri[n]+1))						#Python to Matlab/Fortran indexes
		fid.write('{:d}\t'.format(Triangles.k[WallTri[n]]+1))			#Python to Matlab/Fortran indexes
		fid.write('{:d}\t'.format(Triangles.i[WallTri[n]]+1))
		fid.write('{:d}\n'.format(Triangles.j[WallTri[n]]+1))

	fid.close()


	nTriangles = len(Triangles)

	fid = open(Dir + prefix + ".backinterp_full",'w')
	for n in range(nTriangles):
		fid.write('{:d}\t'.format(n+1))								#Python to Matlab/Fortran indexes
		fid.write('{:d}\t'.format(Triangles.k[n]+1))				#Python to Matlab/Fortran indexes
		fid.write('{:d}\t'.format(Triangles.i[n]+1))
		fid.write('{:d}\n'.format(Triangles.j[n]+1))

	fid.close()

	if(DEBUG > 0):	print("save_wall_back_interp: Completed")
