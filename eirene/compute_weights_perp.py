import numpy				as np
#import matplotlib.pyplot	as pyp
from routines.globals		import DEBUG

# Function definition is here
#=======================================================================================================================================

#=========================================================
# This routine compute triangles weights
#=========================================================

#def compute_weights_perp(Zones, WallTriangles, RKnots, ZKnots, Wall, PlotCheck=0):
def compute_weights_perp(Zones, WallTriangles, RKnots, ZKnots):

	if(DEBUG > 0):	print("compute_weights_perp")

	"""
	if PlotCheck > 0:
		fig  = pyp.figure()
		ax   = fig.add_subplot(1,1,1)
		ax.set_xlabel("R")
		ax.set_ylabel("Z")
		ax.set_aspect(1.)
		ax.set_title("Weights perpendicular flux surface")
	"""

	nZones = len(Zones)

	RKnotsTri  = np.array([RKnots[WallTriangles.p1],RKnots[WallTriangles.p2],RKnots[WallTriangles.p3],RKnots[WallTriangles.p1]])
	ZKnotsTri  = np.array([ZKnots[WallTriangles.p1],ZKnots[WallTriangles.p2],ZKnots[WallTriangles.p3],ZKnots[WallTriangles.p1]])
	SideIndex1 = np.array([1,2,0])
	SideIndex2 = np.array([0,1,2])

	if(DEBUG > 1):	print("\tNorth case")

	for k in range(nZones):
		Nx = Zones[k].Nx
		Nz = Zones[k].Nz

		nListMax 					=  max(1,np.max(Zones[k].list_tri_perp_nums))
		Zones[k].list_tri_n_nnums	=  np.zeros((Nx, Nz), dtype='i4')
		Zones[k].list_tri_n_nums	= -np.ones((Nx, Nz, nListMax), dtype='i4') 
		Zones[k].list_tri_n_weights	=  np.zeros((Nx, Nz, nListMax), dtype='f8') 

		ii,jj = np.where((Zones[k].Chi2[1:-1, 1:-1] == 0) & (Zones[k].Chi2[2:, 1:-1] == 1))		#Case 01
		for l in range(len(ii)):
			i = ii[l]
			j = jj[l]
			nList = Zones[k].list_tri_perp_nnums[i,j]
			List  = Zones[k].list_tri_perp_nums[i,j,0:nList] 
				
#			bR = Zones[k].R2[i+2,j+1] - Zones[k].R2[i+1,j+1]				#matlab method
#			bZ = Zones[k].Z2[i+2,j+1] - Zones[k].Z2[i+1,j+1]

			bR	  = np.zeros((len(List)), dtype='f8')
			bZ	  = np.zeros((len(List)), dtype='f8')
			for lTri in range(len(List)):
				kTri  = WallTriangles.k[List[lTri]]
				iTri  = WallTriangles.i[List[lTri]]
				jTri  = WallTriangles.j[List[lTri]]
				bR[lTri] = 0.5*((Zones[kTri].gridR[iTri+1,jTri] + Zones[kTri].gridR[iTri+1,jTri+1]) - (Zones[kTri].gridR[iTri,jTri] + Zones[kTri].gridR[iTri,jTri+1]))
				bZ[lTri] = 0.5*((Zones[kTri].gridZ[iTri+1,jTri] + Zones[kTri].gridZ[iTri+1,jTri+1]) - (Zones[kTri].gridZ[iTri,jTri] + Zones[kTri].gridZ[iTri,jTri+1]))
			bMod = np.sqrt(bR**2 + bZ**2)
			bR  /= bMod
			bZ  /= bMod

			i1 = SideIndex1[WallTriangles.side[List]]
			i2 = SideIndex2[WallTriangles.side[List]]

			cR		= (RKnotsTri[i1, List] + RKnotsTri[i2, List])*0.5
			cZ		= (ZKnotsTri[i1, List] + ZKnotsTri[i2, List])*0.5
			dSR		= RKnotsTri[i1, List] - RKnotsTri[i2, List]
			dSZ		= ZKnotsTri[i1, List] - ZKnotsTri[i2, List]

			WallTriangles.surf[0:nList] = np.sqrt(dSR**2+dSZ**2)*2*np.pi*cR

			"""
			if PlotCheck > 0:
				ax.plot(RKnotsTri[:,List],ZKnotsTri[:,List],'g-')
				ax.plot([cR,cR+dSR], [cZ,cZ+dSZ],   'b-')
				ax.plot([cR,cR+bR],  [cZ,cZ+bZ],    'r-')
			"""
				
			weights =bR*dSZ-bZ*dSR
			weights = np.where(weights > 0.,weights, np.zeros(nList,dtype='f8')) 

			sum_weights = np.sum(weights)
			if(sum_weights != 0.): weights	 /= sum_weights
			else:				   weights[:] = 0.

			Zones[k].list_tri_n_nnums[i,j]			 = nList
			Zones[k].list_tri_n_nums[i,j,0:nList]	 = List
			Zones[k].list_tri_n_weights[i,j,0:nList] = weights


	if(DEBUG > 1):	print("\tSouth case")

	for k in range(nZones):
		Nx = Zones[k].Nx
		Nz = Zones[k].Nz

		nListMax 					=  max(1,np.max(Zones[k].list_tri_perp_nums))
		Zones[k].list_tri_s_nnums	=  np.zeros((Nx, Nz), dtype='i4')
		Zones[k].list_tri_s_nums	= -np.ones((Nx, Nz, nListMax), dtype='i4') 
		Zones[k].list_tri_s_weights	=  np.zeros((Nx, Nz, nListMax), dtype='f8') 

		ii,jj = np.where((Zones[k].Chi2[1:-1, 1:-1] == 0) & (Zones[k].Chi2[:-2, 1:-1] == 1))		#Case West
		for l in range(len(ii)):
			i = ii[l]
			j = jj[l]
			nList = Zones[k].list_tri_perp_nnums[i,j]
			List  = Zones[k].list_tri_perp_nums[i,j,0:nList] 

#			bR = Zones[k].R2[i,j+1] - Zones[k].R2[i+1,j+1]							#matalb method
#			bZ = Zones[k].Z2[i,j+1] - Zones[k].Z2[i+1,j+1]

			bR	  = np.zeros((len(List)), dtype='f8')
			bZ	  = np.zeros((len(List)), dtype='f8')
			for lTri in range(len(List)):
				kTri  = WallTriangles.k[List[lTri]]
				iTri  = WallTriangles.i[List[lTri]]
				jTri  = WallTriangles.j[List[lTri]]
				bR[lTri] = 0.5*((Zones[kTri].gridR[iTri,jTri] + Zones[kTri].gridR[iTri,jTri+1]) - (Zones[kTri].gridR[iTri+1,jTri] + Zones[kTri].gridR[iTri+1,jTri+1]))
				bZ[lTri] = 0.5*((Zones[kTri].gridZ[iTri,jTri] + Zones[kTri].gridZ[iTri,jTri+1]) - (Zones[kTri].gridZ[iTri+1,jTri] + Zones[kTri].gridZ[iTri+1,jTri+1]))

			bMod = np.sqrt(bR**2 + bZ**2)
			bR  /= bMod
			bZ  /= bMod

			i1 = SideIndex1[WallTriangles.side[List]]
			i2 = SideIndex2[WallTriangles.side[List]]

			cR		= (RKnotsTri[i1, List] + RKnotsTri[i2, List])*0.5
			cZ		= (ZKnotsTri[i1, List] + ZKnotsTri[i2, List])*0.5
			dSR		= RKnotsTri[i1, List] - RKnotsTri[i2, List]
			dSZ		= ZKnotsTri[i1, List] - ZKnotsTri[i2, List]

			WallTriangles.surf[0:nList] = np.sqrt(dSR**2+dSZ**2)*2*np.pi*cR

			"""
			if PlotCheck > 0:
				ax.plot(RKnotsTri[:,List],ZKnotsTri[:,List],'g-')
				ax.plot([cR,cR+dSR], [cZ,cZ+dSZ],   'b-')
				ax.plot([cR,cR+bR],  [cZ,cZ+bZ],    'r-')
			"""

			weights = bR*dSZ-bZ*dSR
				
			weights = np.where(weights > 0.,weights, np.zeros(nList,dtype='f8')) 

			sum_weights = np.sum(weights)
			if(sum_weights != 0.): weights	 /= sum_weights
			else:				   weights[:] = 0.

			Zones[k].list_tri_s_nnums[i,j]			 = nList
			Zones[k].list_tri_s_nums[i,j,0:nList]	 = List
			Zones[k].list_tri_s_weights[i,j,0:nList] = weights

	"""
	if PlotCheck > 0:
		print("\tplot edge triangles and mesh & triangle distance")
		print("\t\tgreen:   triangles")
		print("\t\tred:     mesh distance perpendicular flux surface")
		print("\t\blue:     triangles distance")
		ax.plot(Wall[:,0], Wall[:,1], 'k-')
		pyp.show()
	"""
	
	if(DEBUG > 0):	print("compute_weights_perp: Completed")




         
