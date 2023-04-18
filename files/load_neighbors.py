# Function definition is here

import numpy as np
import os
from routines.globals			import DEBUG

#=================================================
# This routine 
#=================================================

def load_neighbors(Dir="./"):

	try:
		Neigh	 = np.loadtxt(Dir+'soledge2D.neighbors', skiprows=1, dtype='i4')
	except:
		print("\tError: Not found "+ Dir+"soledge2D.neighbors")
		exit()

	nTriangles	= Neigh.shape[0]
	Triangles = np.empty((nTriangles), dtype=[('BC1','i4'),		    ('BC2','i4'),		     ('BC3','i4'),
											  ('neigh1','i4'),		('neigh2','i4'),		 ('neigh3','i4'),
											  ('typeneigh1','i4'),	('typeneigh2','i4'), ('typeneigh3','i4')])
	Triangles = Triangles.view(np.recarray)

	Triangles.neigh1	 = Neigh[:,1] - 1							#Matlab/Fortran to Python ndex
	Triangles.typeneigh1 = Neigh[:,2] - 1
	Triangles.BC1		 = Neigh[:,3]

	Triangles.neigh2	 = Neigh[:,4] - 1
	Triangles.typeneigh2 = Neigh[:,5] - 1
	Triangles.BC2		 = Neigh[:,6]

	Triangles.neigh3	 = Neigh[:,7] - 1
	Triangles.typeneigh3 = Neigh[:,8] - 1
	Triangles.BC3		 = Neigh[:,9]

	return Triangles
	