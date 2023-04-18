import numpy as np
from routines.globals import *

# Function definition is here
#=======================================================================================================================================


def find_neighbors(OldTriangles, TriKnots, RKnots, ZKnots):

	if(DEBUG > 0):	print("eirene_find_neighbors")

	Triangles  = np.copy(OldTriangles)
	Triangles  = Triangles.view(np.recarray)
	nTriangles = len(Triangles)

	if(DEBUG > 1):	print("\tprocessing ",nTriangles," triangles")
	
	TriBC		 =  np.empty((TriKnots.shape[0],3), dtype='i4')
	TriNeigh	 = -np.ones((TriKnots.shape[0],3), dtype='i4')
	TriTypeNeigh = -np.ones((TriKnots.shape[0],3), dtype='i4')

	TriBC[:,0]		= Triangles.BC1
	TriBC[:,1]		= Triangles.BC2
	TriBC[:,2]		= Triangles.BC3
	
	SideIndex = np.array([[0,1],[1,2],[2,0]])
	for k1 in range(TriKnots.shape[0]):					#Loop on triangles
		for iSide1 in range(3):							#For each triangle neighbor on all sides
			if(TriNeigh[k1,iSide1] == -1):
				Found = False
				for iSide2 in range(3):
					k2 = np.where(((TriKnots[k1,SideIndex[iSide1,0]] == TriKnots[k1+1:,SideIndex[iSide2,0]]) & 
								    (TriKnots[k1,SideIndex[iSide1,1]] == TriKnots[k1+1:,SideIndex[iSide2,1]])) |
								   ((TriKnots[k1,SideIndex[iSide1,0]] == TriKnots[k1+1:,SideIndex[iSide2,1]]) &
								    (TriKnots[k1,SideIndex[iSide1,1]] == TriKnots[k1+1:,SideIndex[iSide2,0]])))[0]
					if(len(k2) > 0):
						if(len(k2)>1):
							k2 = k2 + k1 + 1
							print("\t\tfind_neighbors: Error in finding neighbors, more than one match")
							print("\t\tFor Triangle    k1=",k1," and  side=", iSide1)
							print("\t\t             k,i,j=",Triangles.k[k1],Triangles.i[k1],Triangles.j[k1])
							print("\t\t          p1,p2,p3=",Triangles.p1[k1],Triangles.p2[k1],Triangles.p3[k1])
							print("\t\t R(p1),R(p2),R(p3)=",RKnots[Triangles.p1[k1]],RKnots[Triangles.p2[k1]],RKnots[Triangles.p3[k1]])
							print("\t\t Z(p1),Z(p2),Z(p3)=",ZKnots[Triangles.p1[k1]],ZKnots[Triangles.p2[k1]],ZKnots[Triangles.p3[k1]])
							print("\t\tFound Triangles k2=",k2," with side=", iSide2)
							print("\t\t             k,i,j=",Triangles.k[k2],Triangles.i[k2],Triangles.j[k2])
							print("\t\t          p1,p2,p3=",Triangles.p1[k2],Triangles.p2[k2],Triangles.p3[k2])
							print("\t\t R(p1),R(p2),R(p3)=",RKnots[Triangles.p1[k2]],RKnots[Triangles.p2[k2]],RKnots[Triangles.p3[k2]])
							print("\t\t Z(p1),Z(p2),Z(p3)=",ZKnots[Triangles.p1[k2]],ZKnots[Triangles.p2[k2]],ZKnots[Triangles.p3[k2]])
							exit()
				
						k2 = k2[0] + k1 + 1	
						TriBC[k1,iSide1]		= BC_TRIANGLE
						TriNeigh[k1,iSide1]		= k2 
						TriTypeNeigh[k1,iSide1]	= iSide2
						
						TriBC[k2,iSide2]		= BC_TRIANGLE
						TriNeigh[k2,iSide2]		= k1 
						TriTypeNeigh[k2,iSide2]	= iSide1
						
						Found = True
						break
					
				if(not Found):
					if(Triangles.Area[k1] == 1): TriBC[k1,iSide1] = BC_CORE		
					else:						 TriBC[k1,iSide1] = BC_WALL	

	Triangles.BC1			= TriBC[:,0]
	Triangles.neigh1		= TriNeigh[:,0]
	Triangles.typeneigh1	= TriTypeNeigh[:,0]

	Triangles.BC2			= TriBC[:,1]
	Triangles.neigh2		= TriNeigh[:,1]
	Triangles.typeneigh2	= TriTypeNeigh[:,1]

	Triangles.BC3			= TriBC[:,2]
	Triangles.neigh3		= TriNeigh[:,2]
	Triangles.typeneigh3	= TriTypeNeigh[:,2]
	
	ii= np.where((Triangles.typeneigh1 < 0) | (Triangles.typeneigh2 < 0) | (Triangles.typeneigh3 < 0))
	if(DEBUG > 1):	print("\tneighbors not found for ",len(ii[0])," triangles")
	ii= np.where((Triangles.BC1 == BC_WALL) | (Triangles.BC2 == BC_WALL) | (Triangles.BC3 == BC_WALL))
	if(DEBUG > 1):	print("\tneighbors found wall ",len(ii[0])," sides")
	ii= np.where((Triangles.BC1 == BC_CORE) | (Triangles.BC2 == BC_CORE) | (Triangles.BC3 == BC_CORE))
	if(DEBUG > 1):	print("\tneighbors found core ",len(ii[0])," sides")

	if(DEBUG > 0):	print("eirene_find_neighbors: Completed")

	return Triangles
