# Function definition is here

import numpy as np

#=========================================================
# This create Triangles structure
#=========================================================

def triangles_new(nTriangles):

	Triangles = np.empty((nTriangles), dtype=[  ('k','i4'), ('i','i4'), ('j','i4'),
												('np_k','i4'), ('np_i','i4'), ('np_j','i4'),
												('p1','i4'), ('p2','i4'), ('p3','i4'), ('Area','i4'),
												('neigh1','i4'), ('typeneigh1','i4'), ('BC1','i4'),
												('neigh2','i4'), ('typeneigh2','i4'), ('BC2','i4'),
												('neigh3','i4'), ('typeneigh3','i4'), ('BC3','i4'),
												('step','i4'),('PlasmaVacuum','i4')])
	Triangles		= Triangles.view(np.recarray)
	Triangles.BC1   = -np.ones(nTriangles, dtype='i4')			#Boudary condition to -1 as unknown
	Triangles.BC2   = np.copy(Triangles.BC1)
	Triangles.BC3   = np.copy(Triangles.BC1)

	Triangles.neigh1 = -np.ones(nTriangles, dtype='i4')		#near triangle non exist
	Triangles.neigh2 = np.copy(Triangles.neigh1)
	Triangles.neigh3 = np.copy(Triangles.neigh1)

	Triangles.typeneigh1 = -np.ones(nTriangles, dtype='i4')	#type near triangle non exist
	Triangles.typeneigh2 = np.copy(Triangles.typeneigh1)
	Triangles.typeneigh3 = np.copy(Triangles.typeneigh1)
	Triangles.PlasmaVacuum = np.zeros(nTriangles, dtype='i4')
	
	return Triangles

#=========================================================
# This copy Triangles structure
#=========================================================

def triangles_copy(OldTriangles):

	Triangles = np.copy(OldTriangles)
	Triangles = Triangles.view(np.recarray)

	return Triangles

#=========================================================
# This change Triangles structure
#=========================================================

def triangles_reshape(OldTriangles, nTriangles=0):

	nOldTriangles = OldTriangles.k.shape[0]

	if(nTriangles == 0) : nTriangles = nOldTriangles + 1
	Triangles 	  = triangles_new(nTriangles)
		
	if(nOldTriangles > 0): Triangles[:nOldTriangles] = OldTriangles
		
	Triangles.BC1[nOldTriangles:]		 = -1						#Boudary condition to -1 as unknown
	Triangles.BC2[nOldTriangles:]		 = -1
	Triangles.BC3[nOldTriangles:]		 = -1

	Triangles.neigh1[nOldTriangles:]	 = -1						#near triangle non exist
	Triangles.neigh2[nOldTriangles:]	 = -1
	Triangles.neigh3[nOldTriangles:]	 = -1

	Triangles.typeneigh1[nOldTriangles:] = -1						#type near triangle non exist
	Triangles.typeneigh2[nOldTriangles:] = -1
	Triangles.typeneigh3[nOldTriangles:] = -1

	return Triangles

#=========================================================
# This create TriKnots structure
#=========================================================

def triknots_new(nTriKnots):
	return np.empty((nTriKnots,3), dtype='i')

#=========================================================
# This change TriKnots structure
#=========================================================

def triknots_reshape(OldTriKnots, nTriKnots=0):

	nOldTriKnots = OldTriKnots.shape[0]

	if nTriKnots == 0 : nTriKnots = nOldTriKnots + 1
	TriKnots = np.empty((nTriKnots,3), dtype='i')
		
	if nOldTriKnots > 0: TriKnots[0:nOldTriKnots,:] = OldTriKnots

	return TriKnots

#==========================================
# This routine set triangles and Nodes
#==========================================

def triangles_points_add(k, i, j, p1, p2, p3, Area, Direct, OldTriangles, OldTriKnots):

	"""
#	Intructions to debug problems with triangles generation

	if((np.min(p1) < 0) or (np.min(p2) < 0) or (np.min(p3) < 0)):
		raise ValueError("triangles_points_add: Invalid p1/p2/p3 values")

	if((np.min(np.abs(p1-p2)) == 0) or (np.min(np.abs(p1-p3)) == 0) or (np.min(np.abs(p2-p3)) == 0)):
		pi1 = [p1,p1,p2]
		pi2 = [p2,p3,p3]
		for ki in range(3):
			ii=np.where((pi1[ki]-pi2[ki]) == 0)[0]
			if(len(ii) > 0):
				print("\ttriangles_points_add: two identical points for:")
				print("k = ",k)
				print("i = ",i[ii])
				print("j = ",j[ii])
				print("p = ",pi1[ki][ii])
				print("p = ",pi2[ki][ii])
		raise ValueError("triangles_points_add: Invalid p1=p2=p3")
	"""

	nPoints       = len(i)
	nTrianglesOld = len(OldTriangles)
	nTriangles    = nTrianglesOld+nPoints
	Triangles  	  = triangles_reshape(OldTriangles, nTriangles)
	TriKnots   	  = triknots_reshape(OldTriKnots, nTriangles)

	Triangles.k[nTrianglesOld:nTriangles]	= k
	Triangles.i[nTrianglesOld:nTriangles]	= i
	Triangles.j[nTrianglesOld:nTriangles]	= j
	Triangles.p1[nTrianglesOld:nTriangles]	= p1

	if Direct != 0:
		Triangles.p2[nTrianglesOld:nTriangles]	= p2
		Triangles.p3[nTrianglesOld:nTriangles]	= p3
	else:
		Triangles.p2[nTrianglesOld:nTriangles]	= p3
		Triangles.p3[nTrianglesOld:nTriangles]	= p2
	Triangles.Area[nTrianglesOld:nTriangles]	= Area

	TriKnots[nTrianglesOld:nTriangles, 0] = Triangles.p1[nTrianglesOld:nTriangles]
	TriKnots[nTrianglesOld:nTriangles, 1] = Triangles.p2[nTrianglesOld:nTriangles]
	TriKnots[nTrianglesOld:nTriangles, 2] = Triangles.p3[nTrianglesOld:nTriangles]

	return Triangles, TriKnots
	

#==========================================
# This routine set triangles and Nodes
#==========================================

def triangles_point_add(k, i, j, p1, p2, p3, Area, Direct, OldTriangles, OldTriKnots):

	nTriangles = len(OldTriangles)
	Triangles  = triangles_reshape(OldTriangles)
	TriKnots   = triknots_reshape(OldTriKnots)

	Triangles.k[nTriangles]		= k
	Triangles.i[nTriangles]		= i
	Triangles.j[nTriangles]		= j
	Triangles.p1[nTriangles]	= p1

	if(Direct != 0):
		Triangles.p2[nTriangles]	= p2
		Triangles.p3[nTriangles]	= p3
	else:
		Triangles.p2[nTriangles]	= p3
		Triangles.p3[nTriangles]	= p2

	Triangles.Area[nTriangles]	= Area

	TriKnots[nTriangles, 0] = Triangles.p1[nTriangles]
	TriKnots[nTriangles, 1] = Triangles.p2[nTriangles]
	TriKnots[nTriangles, 2] = Triangles.p3[nTriangles]

	return Triangles, TriKnots


#=========================================================
# This create KnotsInterp structure
#=========================================================

def KnotsInterp_new(nKnots):


	KnotsInterp = np.zeros((nKnots), dtype=[('assp','i4'), ('nsol','i4'), ('neir', 'i4'), ('sol', 'i4', (8,3)), ('eir', 'i4', (3))])
	KnotsInterp = KnotsInterp.view(np.recarray)

	KnotsInterp.assp = np.zeros(nKnots, dtype='i4')
	KnotsInterp.nsol = np.zeros(nKnots, dtype='i4')
	KnotsInterp.neir = np.zeros(nKnots, dtype='i4')
	KnotsInterp.sol  = -np.ones((nKnots,8,3), dtype='i4')							#set to -1 non defined index (in python 0 is a good index)
	KnotsInterp.eir  = -np.ones((nKnots,3), dtype='i4')
	
	return KnotsInterp


#=========================================================
# This copy KnotsInterp structure
#=========================================================

def KnotsInterp_copy(OldKnotsInterp):

	KnotsInterp = np.copy(OldKnotsInterp)
	KnotsInterp = KnotsInterp.view(np.recarray)

	return KnotsInterp

#=========================================================
# This change KnotsInterp structure
#=========================================================

def KnotsInterp_reshape(OldKnotsInterp, nKnotsInterp=0):

	nOldKnotsInterp = OldKnotsInterp.assp.shape[0]

	if(nKnotsInterp == 0): nKnotsInterp = nOldKnotsInterp + 1
	KnotsInterp = KnotsInterp_new(nKnotsInterp)
		
	if(nOldKnotsInterp > 0): KnotsInterp[0:nOldKnotsInterp] = OldKnotsInterp

	return KnotsInterp


