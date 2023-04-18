
# Function definition is here

import types
import h5py
import numpy							as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot 				as pyp
from routines.h5_routines 				import h5_read
from routines.intersect_contour			import intersect_contour
from routines.find_closest_segment		import find_closest_segment
from eirene.get_wall_triangles			import get_wall_triangles


#=========================================================
# This routine plot triangles
#=========================================================

def get_wall_triangle(Eirene, rz0_line = [2.,0.], theta_line=0., no_plot=0, no_print=0, no_triangles=1):

	print("get_wall_triangle")

	RWallTriangles, ZWallTriangles, iWallKnots, iWallTriangles, iWallSide = get_wall_triangles(Eirene)					#Get sorted wall triangles
	nWallTriangles	= len(iWallTriangles)
	
	cWall			= types.SimpleNamespace()
	cWall.arc		= [types.SimpleNamespace()]
	cWall.arc[0].x	= np.append(RWallTriangles, RWallTriangles[0])
	cWall.arc[0].y	= np.append(ZWallTriangles, ZWallTriangles[0])

	rMax			= np.max(np.sqrt((cWall.arc[0].x - rz0_line[0])**2 + (cWall.arc[0].y - rz0_line[1])**2))*2
	theta_line		= theta_line*np.pi/180.

	cLine			= types.SimpleNamespace()
	cLine.arc		= [types.SimpleNamespace()]
	cLine.arc[0].x	= np.array([rz0_line[0], rz0_line[0]+rMax*np.cos(theta_line)])
	cLine.arc[0].y	= np.array([rz0_line[1], rz0_line[1]+rMax*np.sin(theta_line)])
	
	X				= intersect_contour(cWall, cLine)
	RIntersect		= X[0].x
	ZIntersect		= X[0].y
	d, iTri			= find_closest_segment(RIntersect, ZIntersect, cWall.arc[0].x, cWall.arc[0].y)

	nTriIntersect  = iWallTriangles[iTri]
	nSideIntersect = iWallSide[iTri]
	if(no_print == 0):
		print("\tTriangle n.",nTriIntersect+1)
		print("\tSide     n.",nSideIntersect+1)

	if(no_plot == 0):
#		Plot triangles and wall

		Fig   = pyp.figure()
		Ax   = Fig.add_subplot(1,1,1)
		Ax.set_xlabel("R")
		Ax.set_ylabel("Z")
		Ax.set_aspect(1.)
		
		Ax.triplot(RKnots, ZKnots, WallTriKnots, 'b-')
		try:
			if_mesh = h5py.File(path+"/mesh.h5", "r")
			RWall	= h5_read(if_mesh, "wall/R")
			ZWall	= h5_read(if_mesh, "wall/Z")
			if_mesh.close()

			Ax.plot(RWall,ZWall,'k-')
		except:
			if_mesh = 0

#		Ax.plot(RIntersect,ZIntersect,'ro')
		Ax.plot(RLine[0,:],ZLine[0,:],'r-')
		Ax.plot([RKnots[TriKnots[nTriIntersect, nSideIntersect]],RKnots[TriKnots[nTriIntersect, ind2[nSideIntersect]]]],
				[ZKnots[TriKnots[nTriIntersect, nSideIntersect]],ZKnots[TriKnots[nTriIntersect, ind2[nSideIntersect]]]],'go')
					
		pyp.show()

	print("get_wall_triangle: Completed")

	if(no_triangles == 0):
		iOrder = np.arange(len(RWallTriangles))
		iOrder = np.append(iOrder[iTri:], iOrder[:iTri]) 
		return nTriIntersect, nSideIntersect, RIntersect, ZIntersect, RWallTriangles[iOrder], ZWallTriangles[iOrder], iWallTriangles[iOrder], iWallSide[iOrder], iWallKnots[iOrder] 
	else:	
		return nTriIntersect, nSideIntersect, RIntersect, ZIntersect
				
