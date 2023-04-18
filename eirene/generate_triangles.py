import numpy as np
from eirene.triangles_routines			import triangles_points_add
from eirene.triangles_routines			import triangles_point_add
from eirene.triangles_routines			import triangles_new
from eirene.triangles_routines			import triangles_copy
from eirene.triangles_routines			import triknots_new
from routines.globals					import DEBUG, CORE_NEIGHBOUR

# Function definition is here
#=======================================================================================================================================

#=================================================
# This routine generate triangles in plasma
#=================================================

def generate_triangles_in_the_plasma(Zones, Direct):

	if(DEBUG > 0):	print("generate_triangles_in_the_plasma")

	nZones	= len(Zones)
	Triangles = triangles_new(0)
	TriKnots  = triknots_new(0)
	
	for k in range(nZones):

		for iPwall in range(Zones[k].IsCrossed.shape[3]):
			IsCrossed  = np.sum(Zones[k].IsCrossed[:,:,:,iPwall], axis=2)
			ii, jj	   = np.where(((IsCrossed == 0) & (Zones[k].InPlasma == 1)) | ((Zones[k].IsAligned == 1) & (Zones[k].InPlasmaRZ == 1)))
			nPoints = len(ii)
			if(nPoints > 0):

#				Like the matlab find_triangle_area.m

				if(Zones[k].Neighbour.south == CORE_NEIGHBOUR):
					Area = np.where(ii == 0, np.ones(nPoints, dtype='i4'), np.zeros(nPoints, dtype='i4'))			#set core condition for triangles in first core quadrangles
				else:
					Area = np.zeros(nPoints, dtype='i4')

#				premier triangle
				Triangles, TriKnots = triangles_points_add(k, ii, jj, Zones[k].KnotA[ii,jj], Zones[k].KnotB[ii,jj], Zones[k].KnotC[ii,jj], Area, Direct, Triangles, TriKnots)

#				deuxieme triangle
				Triangles, TriKnots = triangles_points_add(k, ii, jj, Zones[k].KnotA[ii,jj], Zones[k].KnotC[ii,jj], Zones[k].KnotD[ii,jj], Area, Direct, Triangles, TriKnots)


	if(DEBUG > 1):	print("\tgenerated ", len(Triangles)," triangles in plasma")

	if(DEBUG > 0):	print("generate_triangles_in_the_plasma: Completed")

	return Triangles, TriKnots

#=======================================================================================================================================

#==========================================
# This routine generate triangles on wall
#==========================================

def generate_triangles_on_wall(Zones, OldTriangles, OldTriKnots, Direct):

	if(DEBUG > 0):	print("generate_triangles_on_wall")
	
	nZones	= len(Zones)
	Triangles = triangles_copy(OldTriangles)
	TriKnots  = np.copy(OldTriKnots)

	Area   = -1
	for k in range(nZones):
		for iPwall in range(Zones[k].IsCrossed.shape[3]):
			IsCrossed  = np.sum(Zones[k].IsCrossed[:,:,:,iPwall], axis=2)
			ii, jj = np.where((Zones[k].IsAligned == 0) & (IsCrossed > 0) & (Zones[k].InPlasma == 1))
			if(len(ii) > 0):
#				if(DEBUG > 1):	print("\tCondition 1")
			
#				Condition 1  (H Corner)
#				Condition 1a create 1 triangle (AEH)
				ll = np.where((Zones[k].KnotE[ii,jj] != -1) & (Zones[k].KnotH[ii,jj] != -1) & (Zones[k].KnotA[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotA[ii[ll], jj[ll]], Zones[k].KnotE[ii[ll], jj[ll]], Zones[k].KnotH[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
#				create 3 triangles (EBH, BCH, HCD)
				ll = np.where((Zones[k].KnotE[ii,jj] != -1) & (Zones[k].KnotH[ii,jj] != -1) & (Zones[k].KnotA[ii,jj] == -1) & (Zones[k].KnotB[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotE[ii[ll], jj[ll]], Zones[k].KnotB[ii[ll], jj[ll]], Zones[k].KnotH[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
				ll = np.where((Zones[k].KnotE[ii,jj] != -1) & (Zones[k].KnotH[ii,jj] != -1) & (Zones[k].KnotA[ii,jj] == -1) & (Zones[k].KnotB[ii,jj] != -1) & (Zones[k].KnotC[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotB[ii[ll], jj[ll]], Zones[k].KnotC[ii[ll], jj[ll]], Zones[k].KnotH[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
				ll = np.where((Zones[k].KnotE[ii,jj] != -1) & (Zones[k].KnotH[ii,jj] != -1) & (Zones[k].KnotA[ii,jj] == -1) & (Zones[k].KnotC[ii,jj] != -1) & (Zones[k].KnotD[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotH[ii[ll], jj[ll]], Zones[k].KnotC[ii[ll], jj[ll]], Zones[k].KnotD[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
			
#				if(DEBUG > 1):	print("\tCondition 2")
			
#				Condition 2  (F corner)
#				Condition 2a create 1 triangle (BFE)
				ll = np.where((Zones[k].KnotE[ii,jj] != -1) & (Zones[k].KnotF[ii,jj] != -1) & (Zones[k].KnotB[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotB[ii[ll], jj[ll]], Zones[k].KnotF[ii[ll], jj[ll]], Zones[k].KnotE[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
#				Condition 2b create 3 triangles (EFA, AFD, FCD)
				ll = np.where((Zones[k].KnotE[ii,jj] != -1) & (Zones[k].KnotF[ii,jj] != -1) & (Zones[k].KnotB[ii,jj] == -1) & (Zones[k].KnotA[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotE[ii[ll], jj[ll]], Zones[k].KnotF[ii[ll], jj[ll]], Zones[k].KnotA[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
				ll = np.where((Zones[k].KnotE[ii,jj] != -1) & (Zones[k].KnotF[ii,jj] != -1) & (Zones[k].KnotB[ii,jj] == -1) & (Zones[k].KnotA[ii,jj] != -1) & (Zones[k].KnotD[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotA[ii[ll], jj[ll]], Zones[k].KnotF[ii[ll], jj[ll]], Zones[k].KnotD[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
				ll = np.where((Zones[k].KnotE[ii,jj] != -1) & (Zones[k].KnotF[ii,jj] != -1) & (Zones[k].KnotB[ii,jj] == -1) & (Zones[k].KnotC[ii,jj] != -1) & (Zones[k].KnotD[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotF[ii[ll], jj[ll]], Zones[k].KnotC[ii[ll], jj[ll]], Zones[k].KnotD[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
			
#				if(DEBUG > 1):	print("\tCondition 3")
			
#				Condition 3  (B corner)
#				Condition 3a on cr 2 triangles (BGE, BCG)
				ll = np.where((Zones[k].KnotE[ii,jj] != -1) & (Zones[k].KnotG[ii,jj] != -1) & (Zones[k].KnotB[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotB[ii[ll], jj[ll]], Zones[k].KnotG[ii[ll], jj[ll]], Zones[k].KnotE[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
			
				ll = np.where((Zones[k].KnotE[ii,jj] != -1) & (Zones[k].KnotG[ii,jj] != -1) & (Zones[k].KnotB[ii,jj] != -1) & (Zones[k].KnotC[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotB[ii[ll], jj[ll]], Zones[k].KnotC[ii[ll], jj[ll]], Zones[k].KnotG[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
#				Condition 3b on create 2 triangles (EDA, EGD)
				ll = np.where((Zones[k].KnotE[ii,jj] != -1) & (Zones[k].KnotG[ii,jj] != -1) & (Zones[k].KnotB[ii,jj] == -1) & (Zones[k].KnotA[ii,jj] != -1) & (Zones[k].KnotD[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotE[ii[ll], jj[ll]], Zones[k].KnotD[ii[ll], jj[ll]], Zones[k].KnotA[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
				ll = np.where((Zones[k].KnotE[ii,jj] != -1) & (Zones[k].KnotG[ii,jj] != -1) & (Zones[k].KnotB[ii,jj] == -1) & (Zones[k].KnotD[ii,jj] != -1))[0]
				if len(ll) > 0:
						Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotE[ii[ll], jj[ll]], Zones[k].KnotG[ii[ll], jj[ll]], Zones[k].KnotD[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
#				if(DEBUG > 1):	print("\tCondition 4")
			
#				Condition 4a on create 1 triangle (FCG)
				ll = np.where((Zones[k].KnotG[ii,jj] != -1) & (Zones[k].KnotF[ii,jj] != -1) & (Zones[k].KnotC[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotF[ii[ll], jj[ll]], Zones[k].KnotC[ii[ll], jj[ll]], Zones[k].KnotG[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
#				Condition 4b create 3 triangles (BFA, AFD, FGD)
				ll = np.where((Zones[k].KnotG[ii,jj] != -1) & (Zones[k].KnotF[ii,jj] != -1) & (Zones[k].KnotC[ii,jj] == -1) & (Zones[k].KnotA[ii,jj] != -1) & (Zones[k].KnotB[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotB[ii[ll], jj[ll]], Zones[k].KnotF[ii[ll], jj[ll]], Zones[k].KnotA[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
				ll = np.where((Zones[k].KnotG[ii,jj] != -1) & (Zones[k].KnotF[ii,jj] != -1) & (Zones[k].KnotC[ii,jj] == -1) & (Zones[k].KnotA[ii,jj] != -1) & (Zones[k].KnotD[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotA[ii[ll], jj[ll]], Zones[k].KnotF[ii[ll], jj[ll]], Zones[k].KnotD[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
				ll = np.where((Zones[k].KnotG[ii,jj] != -1) & (Zones[k].KnotF[ii,jj] != -1) & (Zones[k].KnotC[ii,jj] == -1) & (Zones[k].KnotD[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotF[ii[ll], jj[ll]], Zones[k].KnotG[ii[ll], jj[ll]], Zones[k].KnotD[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
#				if(DEBUG > 1):	print("\tCondition 5")
			
#				Condition 5  (H Corner)
#				Condition 5a create triangle (HGD)
				ll = np.where((Zones[k].KnotG[ii,jj] != -1) & (Zones[k].KnotH[ii,jj] != -1) & (Zones[k].KnotD[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotH[ii[ll], jj[ll]], Zones[k].KnotG[ii[ll], jj[ll]], Zones[k].KnotD[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
	#			Condition 5b create 3 triangles (ABH, BCH,HCG)
				ll = np.where((Zones[k].KnotG[ii,jj] != -1) & (Zones[k].KnotH[ii,jj] != -1) & (Zones[k].KnotD[ii,jj] == -1) & (Zones[k].KnotA[ii,jj] != -1) & (Zones[k].KnotB[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotA[ii[ll], jj[ll]], Zones[k].KnotB[ii[ll], jj[ll]], Zones[k].KnotH[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
				ll = np.where((Zones[k].KnotG[ii,jj] != -1) & (Zones[k].KnotH[ii,jj] != -1) & (Zones[k].KnotD[ii,jj] == -1) & (Zones[k].KnotC[ii,jj] != -1) & (Zones[k].KnotB[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotB[ii[ll], jj[ll]], Zones[k].KnotC[ii[ll], jj[ll]], Zones[k].KnotH[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
				ll = np.where((Zones[k].KnotG[ii,jj] != -1) & (Zones[k].KnotH[ii,jj] != -1) & (Zones[k].KnotD[ii,jj] == -1) & (Zones[k].KnotC[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotH[ii[ll], jj[ll]], Zones[k].KnotC[ii[ll], jj[ll]], Zones[k].KnotG[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
#				if(DEBUG > 1):	print("\tCondition 6")
			
#				Condition 6
#				cas 6a create 2 triangles (BHA, BFH)
				ll = np.where((Zones[k].KnotF[ii,jj] != -1) & (Zones[k].KnotH[ii,jj] != -1) & (Zones[k].KnotB[ii,jj] != -1) & (Zones[k].KnotA[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotB[ii[ll], jj[ll]], Zones[k].KnotH[ii[ll], jj[ll]], Zones[k].KnotA[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
				ll = np.where((Zones[k].KnotF[ii,jj] != -1) & (Zones[k].KnotH[ii,jj] != -1) & (Zones[k].KnotB[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotB[ii[ll], jj[ll]], Zones[k].KnotF[ii[ll], jj[ll]], Zones[k].KnotH[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
#				Condition 6b create 2 triangles (FDH, FCD)
				ll = np.where((Zones[k].KnotF[ii,jj] != -1) & (Zones[k].KnotH[ii,jj] != -1) & (Zones[k].KnotB[ii,jj] == -1) & (Zones[k].KnotD[ii,jj] != -1))[0]
				if len(ll) > 0:
						Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotF[ii[ll], jj[ll]], Zones[k].KnotD[ii[ll], jj[ll]], Zones[k].KnotH[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
			
				ll = np.where((Zones[k].KnotF[ii,jj] != -1) & (Zones[k].KnotH[ii,jj] != -1) & (Zones[k].KnotB[ii,jj] == -1) & (Zones[k].KnotC[ii,jj] != -1) & (Zones[k].KnotD[ii,jj] != -1))[0]
				if len(ll) > 0:
					Triangles, TriKnots = triangles_points_add(k, ii[ll], jj[ll], Zones[k].KnotF[ii[ll], jj[ll]], Zones[k].KnotC[ii[ll], jj[ll]], Zones[k].KnotD[ii[ll], jj[ll]], Area, Direct, Triangles, TriKnots)
		

	if(DEBUG > 1):	print("\tgenerated ", len(Triangles)-len(OldTriangles)," triangles in wall")
	if(DEBUG > 0):	print("generate_triangles_on_wall: Completed")

	return Triangles, TriKnots

#=======================================================================================================================================

#==========================================
# This routine generate triangles in center
#==========================================

def generate_triangles_center(Zones, nKnotsCenter, OldTriangles, OldTriKnots, Direct):

	if(DEBUG > 0):	print("generate_triangles_center")

	Triangles = triangles_copy(OldTriangles)
	TriKnots  = np.copy(OldTriKnots)

	nZones = len(Zones)
	Area   = 1
	for k in range(nZones):
		if(Zones[k].Neighbour.south == CORE_NEIGHBOUR): 							#center
			jj	 = np.arange((Zones[k].Nz), dtype='i4')
			mone = -np.ones((Zones[k].Nz), dtype='i4')
			Triangles, TriKnots = triangles_points_add(-1, mone, mone, Zones[k].KnotA[0,jj], nKnotsCenter, Zones[k].KnotD[0,jj], Area, Direct, Triangles, TriKnots)

	if(DEBUG > 1):	print("\tgenerated ", len(Triangles)-len(OldTriangles)," triangles in center")
	if(DEBUG > 0):	print("generate_triangles_center: Completed")

	return Triangles, TriKnots

