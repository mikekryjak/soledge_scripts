import numpy				as np
import os
from routines.globals		import DEBUG, FILE_TYPE_SOLEDGE2D, FILE_TYPE_SOLEDGE3X

#=================================================
# This routine write soledge2D.elemente file
#=================================================

def save_elemente(Dir, Eirene, FileType=FILE_TYPE_SOLEDGE2D):

	if(DEBUG > 0):	print("save_elemente")
	try:
		os.mkdir(Dir)
	except OSError:
		pass

	if(FileType == FILE_TYPE_SOLEDGE2D):	prefix = "soledge2D"
	else:									prefix = "raptorX"

	if(DEBUG > 1):	print("\tsaving to: "+ Dir + prefix + ".elemente")
	
	fid = open(Dir + prefix + ".elemente",'w')

	Triangles  = Eirene.Triangles
	nTriangles = len(Triangles)
	fid.write('{:d}\n'.format(nTriangles))
	for n in range(nTriangles):
		fid.write('{:d}\t'.format(n + 1))					#Python index to Matalab/Fortran indexes
		fid.write('{:d}\t'.format(Triangles.p1[n]+1))
		fid.write('{:d}\t'.format(Triangles.p2[n]+1))
		fid.write('{:d}\n'.format(Triangles.p3[n]+1))

	fid.close()

	if(DEBUG > 0):	print("save_elemente: Completed")


