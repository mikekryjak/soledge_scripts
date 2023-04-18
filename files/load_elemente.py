import numpy as np
import os
from routines.input_routines import scan_array

#=================================================
# This routine read soledge2D.elemente file
#=================================================

def load_elemente(Dir, FileType=FILE_TYPE_SOLEDGE2D):

	if(FileType == FILE_TYPE_SOLEDGE2D):	prefix = "soledge2D"
	else:									prefix = "raptorX"

	fid = open(Dir + prefix + ".elemente",'r')

	nTriangles	= scan_array(fid.readline(),"int")
	Triangles = np.empty((nTriangles), dtype=[('ntri','i4'), ('p1','i4'), ('p2','i4'), ('p3','i4')])

	for n in range(nTriangles):
		Data = scan_array(fid.readline(),"int", 4)
		Triangles.ntri[n] = Data[0] - 1							#Matlab/Fortran to Python ndex
		Triangles.p1[n]	  = Data[1] - 1
		Triangles.p2[n]	  = Data[2] - 1
		Triangles.p3[n]	  = Data[3] - 1

	fid.close()

	return Triangles

