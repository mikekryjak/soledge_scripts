
# Function definition is here

import h5py
import numpy						as np
from routines.intersection			import intersection
from eirene.get_wall_triangle		import get_wall_triangle
from files.load_neighbors			import load_neighbors
from routines.h5_routines			import h5_read

#==============================================================================
# This routine plots ne/ni and Te/Ti ionization and gas pressure on eirene mesh
#==============================================================================

def get_triangles_on_line(path, RLine, ZLine):

	print("get_triangles_on_line")
#	print("\dget_triangles_on_line: RLine, ZLine=",RLine, ZLine)

	OtherSides  = np.array([[1,2],[0,2],[0,1]])
	Vertex	    = np.array([[0,1],[1,2],[2,0]])

	theta_line  = np.arctan2(np.array([ZLine[1]-ZLine[0]]), np.array([RLine[1]-RLine[0]]))*180./np.pi
	theta_line	= float(theta_line[0])
	OldTriangle, OldSide, Ri, Zi = get_wall_triangle(path=path, rz0_line=[RLine[0],ZLine[0]], theta_line=theta_line, no_plot=1, no_print=1)
#	print("\dget_triangles_on_line: OldTriangle, OldSide, Ri, Zi=",OldTriangle, OldSide, Ri, Zi)

	IntTriangles = np.array([OldTriangle])
	IntSides	 = np.array([OldSide])
	IntR		 = np.array([Ri])
	IntZ		 = np.array([Zi])
	Dist		 = np.array([0.])
	nFounds	 	 = 1

	Triangles		= load_neighbors(Dir=path)
	TriNeigh		= np.array([Triangles.neigh1,		Triangles.neigh2,		Triangles.neigh3]).T
	TriTypeneigh	= np.array([Triangles.typeneigh1,	Triangles.typeneigh2,	Triangles.typeneigh3]).T
	TriTypeFace		= np.array([Triangles.BC1,			Triangles.BC2,			Triangles.BC3]).T
	Triangles		= 1

#	Read Eirene mesh

	if_tri	 = h5py.File(path+"triangles.h5", "r")

	TriKnots = h5_read(if_tri,"triangles/tri_knots") - 1						#Matlab/Fortan to python indexes
	R		 = h5_read(if_tri,"knots/R")*0.01
	Z		 = h5_read(if_tri,"knots/Z")*0.01

	if_tri.close()

#	Find all triangles along line

	Found	 = True
	StopNext = False
	while (Found):
		Found  = True
		for iSide in OtherSides[OldSide,:]:
			RVertex = np.array([R[TriKnots[OldTriangle,Vertex[iSide,0]]],R[TriKnots[OldTriangle,Vertex[iSide,1]]]])
			ZVertex = np.array([Z[TriKnots[OldTriangle,Vertex[iSide,0]]],Z[TriKnots[OldTriangle,Vertex[iSide,1]]]])
			Ri,Zi,Found = intersection(RVertex, ZVertex, RLine, ZLine)
#			print("\tdget_triangles_on_line: iSide, RVertex, ZVertex, Ri, Zi, Found=",iSide, RVertex, ZVertex, Ri, Zi, Found)
#			print("\tget_triangles_on_line: Ri, Zi, Found=",Ri, Zi, Found)
			
			if(Found):
				nFounds		+= 1
				IntR		 = np.append(IntR, Ri)
				IntZ		 = np.append(IntZ, Zi)
				if(StopNext):
					IntSides	 = np.append(IntSides, iSide)
					IntTriangles = np.append(IntTriangles, OldTriangle)
					Found = False
				else:
					OldSide      = TriTypeneigh[OldTriangle,iSide]
					OldTriangle	 = TriNeigh[OldTriangle,iSide]
					IntSides	 = np.append(IntSides, OldSide)
					IntTriangles = np.append(IntTriangles, OldTriangle)
					if(TriTypeFace[OldTriangle,iSide] != 0): StopNext   = True	# Found boundary 
				break

	Dist	 		= np.cumsum(np.append(0., np.sqrt((IntR[1:] - IntR[0:nFounds-1])**2 + (IntZ[1:] - IntZ[0:nFounds-1])**2)))
	IntKnots		= np.empty((nFounds,2), dtype = 'i4')
	IntWKnots		= np.empty((nFounds,2), dtype = 'f8')
	IntKnots[:,0]	= TriKnots[IntTriangles,Vertex[IntSides,0]]
	IntKnots[:,1]	= TriKnots[IntTriangles,Vertex[IntSides,1]]
	IntWKnots[:,1]	= np.sqrt((IntR - R[IntKnots[:,0]])**2			 + (IntZ - Z[IntKnots[:,0]])**2) / \
					  np.sqrt((R[IntKnots[:,1]]-R[IntKnots[:,0]])**2 + (Z[IntKnots[:,1]]-Z[IntKnots[:,0]])**2)	
	IntWKnots[:,0]	= 1. - IntWKnots[:,1]

	print("get_triangles_on_line: completed")

	return IntTriangles, IntSides, Dist, IntR, IntZ, IntKnots, IntWKnots