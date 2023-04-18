import numpy as np
from interfaces.plot_and_ask	import plot_and_ask
from eirene.triangles_routines	import triangles_new
from eirene.triangles_routines	import triangles_copy
from eirene.triangles_routines	import triknots_new
from routines.globals			import DEBUG, BC_WALL

# This function replace MATLab function remove_non_connex_triangles
# Function definition is here
#=======================================================================================================================================

def remove_unused_triangles(Zones, RKnots, ZKnots, Triangles, TriKnots):

	if(DEBUG > 0):	print("remove_unused_triangles")

	iRemove = np.where(Triangles.step == 0)[0]
	
	if(len(iRemove) == 0): 
		Knots  		= np.arange(len(RKnots))
		if(DEBUG > 1):	print("remove_unused_triangles: Completed")
		return Knots, RKnots, ZKnots, Triangles, TriKnots

	nTriangles		= len(Triangles)
	iKeep			= np.where(Triangles.step == 1)[0]
	
	nRemoveTriangles	= len(iRemove)

#	Define keep triangles
	
	nKeepTriangles		= nTriangles - nRemoveTriangles
	KeepTriangles		= triangles_new(nKeepTriangles)
	KeepTriangles		= Triangles[iKeep]

#	Define and transfer keep triangles in TriKnots
	
	KeepTriKnots		= triknots_new(nKeepTriangles)

	print("\tnRemoveTriangles=", nRemoveTriangles)

#	Renumber triangles in TrinAtSide
	
	one				= np.ones(nTriangles, dtype='i4')
	zero			= np.zeros(nTriangles, dtype='i4')
	nNewTri			= np.where(Triangles.step == 1, one, zero)
	nNewTri			= np.cumsum(nNewTri) - 1
	iTriAtSide		= np.where(KeepTriangles.neigh1 > -1); iTriAtSide = iTriAtSide[0]
	if(len(iTriAtSide) > 0): KeepTriangles.neigh1[iTriAtSide] = nNewTri[KeepTriangles.neigh1[iTriAtSide]]
	
	iTriAtSide		= np.where(KeepTriangles.neigh2 > -1); iTriAtSide = iTriAtSide[0]
	if(len(iTriAtSide) > 0): KeepTriangles.neigh2[iTriAtSide] = nNewTri[KeepTriangles.neigh2[iTriAtSide]]
	
	iTriAtSide		= np.where(KeepTriangles.neigh3 > -1); iTriAtSide = iTriAtSide[0]
	if(len(iTriAtSide) > 0): KeepTriangles.neigh3[iTriAtSide] = nNewTri[KeepTriangles.neigh3[iTriAtSide]]

#	Find unconected nodes after removing triangles
	
	nKeepKnots			= np.union1d(np.union1d(KeepTriangles.p1,      KeepTriangles.p2),     KeepTriangles.p3)
	nRemoveKnots		= np.union1d(np.union1d(Triangles.p1[iRemove], Triangles.p2[iRemove]),Triangles.p3[iRemove])
	nRemoveKnots		= np.setdiff1d(nRemoveKnots,nKeepKnots)

#	Manage removing of non connected nodes

	if(len(nRemoveKnots) > 0):
		
		nKnots			 		= len(RKnots)
		zero				 	= np.zeros(nKnots, dtype='i4')
		one				 		= np.ones(nKnots,  dtype='i4')
		KeepKnots			 	= np.copy(zero)
		KeepKnots[nKeepKnots]	= one[nKeepKnots]	
		KeepKnots			 	= np.cumsum(KeepKnots)
		KeepKnots[nRemoveKnots]	= zero[nRemoveKnots]						
		KeepKnots				= KeepKnots - 1					#index python style

#		Renumber Knot in Triangles and TriKnots
	
		KeepTriangles.p1 = KeepKnots[KeepTriangles.p1]
		KeepTriangles.p2 = KeepKnots[KeepTriangles.p2]
		KeepTriangles.p3 = KeepKnots[KeepTriangles.p3]
		
		KeepTriKnots[:,0] = np.copy(KeepTriangles.p1)
		KeepTriKnots[:,1] = np.copy(KeepTriangles.p2)
		KeepTriKnots[:,2] = np.copy(KeepTriangles.p3)
	
		nZones = len(Zones)
		
#		Renumber Knot in Zones

		for k in range(nZones):		
			muno  = -np.ones(Zones[k].KnotA.shape, dtype = 'i4')
			
			ii,jj = np.where(Zones[k].KnotA >= 0)
			if(len(ii)> 0):
				KnotTmp				  = np.copy(Zones[k].KnotA)
				Zones[k].KnotA		  = np.copy(muno)
				Zones[k].KnotA[ii,jj] = np.copy(KeepKnots[KnotTmp[ii,jj]])
				
			ii,jj = np.where(Zones[k].KnotB >= 0)
			if(len(ii)> 0):
				KnotTmp				  = np.copy(Zones[k].KnotB)
				Zones[k].KnotB		  = np.copy(muno)
				Zones[k].KnotB[ii,jj] = np.copy(KeepKnots[KnotTmp[ii,jj]])
				
			ii,jj = np.where(Zones[k].KnotC >= 0)
			if(len(ii)> 0):
				KnotTmp				  = np.copy(Zones[k].KnotC)
				Zones[k].KnotC		  = np.copy(muno)
				Zones[k].KnotC[ii,jj] = np.copy(KeepKnots[KnotTmp[ii,jj]])
				
			ii,jj = np.where(Zones[k].KnotD >= 0)
			if(len(ii)> 0):
				KnotTmp				  = np.copy(Zones[k].KnotD)
				Zones[k].KnotD		  = np.copy(muno)
				Zones[k].KnotD[ii,jj] = np.copy(KeepKnots[KnotTmp[ii,jj]])
				
			ii,jj = np.where(Zones[k].KnotE >= 0)
			if(len(ii)> 0):
				KnotTmp				  = np.copy(Zones[k].KnotE)
				Zones[k].KnotE		  = np.copy(muno)
				Zones[k].KnotE[ii,jj] = np.copy(KeepKnots[KnotTmp[ii,jj]])
				
			ii,jj = np.where(Zones[k].KnotF >= 0)
			if(len(ii)> 0):
				KnotTmp				  = np.copy(Zones[k].KnotF)
				Zones[k].KnotF		  = np.copy(muno)
				Zones[k].KnotF[ii,jj] = np.copy(KeepKnots[KnotTmp[ii,jj]])
				
			ii,jj = np.where(Zones[k].KnotG >= 0)
			if(len(ii)> 0):
				KnotTmp				  = np.copy(Zones[k].KnotG)
				Zones[k].KnotG		  = np.copy(muno)
				Zones[k].KnotG[ii,jj] = np.copy(KeepKnots[KnotTmp[ii,jj]])
				
			ii,jj = np.where(Zones[k].KnotH >= 0)
			if(len(ii)> 0):
				KnotTmp				  = np.copy(Zones[k].KnotH)
				Zones[k].KnotH		  = np.copy(muno)
				Zones[k].KnotH[ii,jj] = np.copy(KeepKnots[KnotTmp[ii,jj]])
	
#		Remove non connected knots from knodes R and Z arrays
	
		ii = np.where(KeepKnots >= 0); ii = ii[0]
		
		KeepRKnots = np.copy(RKnots[ii])
		KeepZKnots = np.copy(ZKnots[ii])
		KeepKnots  = np.arange(len(KeepRKnots))

		if(DEBUG > 0):	print("remove_unused_triangles: Completed")

		return KeepKnots, KeepRKnots, KeepZKnots, KeepTriangles, KeepTriKnots
		
	else:
		KeepTriKnots	= np.copy(TriKnots[iKeep, :])
		Knots  			= np.arange(len(RKnots))
		
		if(DEBUG > 0):	print("remove_unused_triangles: Completed")

		return Knots, RKnots, ZKnots, KeepTriangles, KeepTriKnots, 

