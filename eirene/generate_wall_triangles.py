import os
import numpy				as np
import numpy.matlib			as mat
from routines.globals		import *

# Function definition is here
#=======================================================================================================================================

#=========================================================
#  This (with save_wall_back_interp.py) replace identify_wall_triangles.m routine
#=========================================================

def generate_wall_triangles(Triangles):

	if(DEBUG > 0):	print("generate_wall_triangles")

	triNum1		= np.where(Triangles.BC1 == BC_WALL)[0]
	triNum2		= np.where(Triangles.BC2 == BC_WALL)[0]
	triNum3		= np.where(Triangles.BC3 == BC_WALL)[0]
	triFace1	= np.zeros_like(triNum1)
	triFace2	= np.ones_like(triNum2)
	triFace3 	= np.ones_like(triNum3)*2
	triNum  	= np.concatenate((triNum1, triNum2, triNum3))
	triFace		= np.concatenate((triFace1, triFace2, triFace3))

	iSort		= np.argsort(triNum)										#this is not needed
	triNum		= triNum[iSort]
	triFace		= triFace[iSort]

#	WallTriangles	= Is like triangles in MatLab versione
#	k, i, j		= Indexes of quadrangle to which the triangle belongs
#	p1, p2, p3	= Number of knots at three vertexes
#	neigh#		= number of neighbor triangle on side #
#	typeneigh#	= side number of neighbor triangle on side # (0,1,2)
#	BC#			= boundary condiction (BC_UKNOWN=-1, BC_TRIANGLE = 0, BC_WALL = 3, BC_CORE = 1)
#	step			= triangle type
#	ntri			= number of global triangle (not WallTriangle)
#	side			= number triangle side to wall (0,1,2)
#	surf			= toroidal surface of side to the wall
#	list_#		= (k,i,j) indexes of last uncutted quadrangle (East, West, North, South)
#	weight_#		= fraction of power last uncutted quadrangle to this triangle (East, West, North, South)
	
	nWallTriangles = len(triNum)
	WallTriangles  = np.empty((nWallTriangles), dtype=[('k','i4'), ('i','i4'), ('j','i4'),('p1','i4'), ('p2','i4'), ('p3','i4'), ('Area','i4'),
												('neigh1','i4'), ('typeneigh1','i4'), ('BC1','i4'),
												('neigh2','i4'), ('typeneigh2','i4'), ('BC2','i4'),
												('neigh3','i4'), ('typeneigh3','i4'), ('BC3','i4'),
												('step','i4'), ('ntri','i4'), ('side','i4'), ('surf','f8'),
												('list_e', 'i4',(3)), ('weight_e','f8'),
												('list_w', 'i4',(3)), ('weight_w','f8'),
												('list_n', 'i4',(3)), ('weight_n','f8'),
												('list_s', 'i4',(3)), ('weight_s','f8')])
	WallTriangles  = WallTriangles.view(np.recarray)

	WallTriangles.k				= np.copy(Triangles.k[triNum])
	WallTriangles.i				= np.copy(Triangles.i[triNum])
	WallTriangles.j				= np.copy(Triangles.j[triNum])
	WallTriangles.p1			= np.copy(Triangles.p1[triNum])
	WallTriangles.p2			= np.copy(Triangles.p2[triNum])
	WallTriangles.p3			= np.copy(Triangles.p3[triNum])
	WallTriangles.neigh1		= np.copy(Triangles.neigh1[triNum])
	WallTriangles.typeneigh1	= np.copy(Triangles.typeneigh1[triNum])
	WallTriangles.BC1			= np.copy(Triangles.BC1[triNum])
	WallTriangles.neigh2		= np.copy(Triangles.neigh2[triNum])
	WallTriangles.typeneigh2	= np.copy(Triangles.typeneigh2[triNum])
	WallTriangles.BC2			= np.copy(Triangles.BC2[triNum])
	WallTriangles.neigh3		= np.copy(Triangles.neigh3[triNum])
	WallTriangles.typeneigh3	= np.copy(Triangles.typeneigh3[triNum])
	WallTriangles.BC3			= np.copy(Triangles.BC3[triNum])
	WallTriangles.step			= np.copy(Triangles.step[triNum])
	WallTriangles.ntri			= np.copy(triNum)
	WallTriangles.side			= np.copy(triFace)

	WallTriangles.list_e		= -np.ones(WallTriangles.list_e.shape, dtype='i4') 
	WallTriangles.list_w		= -np.ones(WallTriangles.list_w.shape, dtype='i4')
	WallTriangles.list_n		= -np.ones(WallTriangles.list_n.shape, dtype='i4')
	WallTriangles.list_s		= -np.ones(WallTriangles.list_s.shape, dtype='i4')

	WallTriangles.weight_e		= np.zeros(WallTriangles.weight_e.shape, dtype='f8')
	WallTriangles.weight_w		= np.zeros(WallTriangles.weight_w.shape, dtype='f8')
	WallTriangles.weight_n		= np.zeros(WallTriangles.weight_n.shape, dtype='f8')
	WallTriangles.weight_s		= np.zeros(WallTriangles.weight_s.shape, dtype='f8')


	iTri = np.where((WallTriangles.typeneigh1 == -1) | (WallTriangles.typeneigh2 == -1) |(WallTriangles.typeneigh3 == -1))[0]
	if(len(iTri) != nWallTriangles):
		print("\tlen(iTri)=",len(iTri)," nWallTriangles=",nWallTriangles)
		print("\tError in findind wall sides!!!")
		exit()

	if(DEBUG > 0):	print("generate_wall_triangles: Completed")

	return WallTriangles