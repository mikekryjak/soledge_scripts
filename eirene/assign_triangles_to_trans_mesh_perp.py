import os
import numpy					as np
from interfaces.plot_and_ask	import plot_and_ask
from routines.utils_routines	import append_RZ
from routines.utils_walls		import set_XYlines_walls
from routines.globals			import DEBUG, TRANS_MESH_DMAX

# Function definition is here
#=======================================================================================================================================

#=========================================================
# This routine assigne trianges to mesh
#=========================================================

def assign_triangles_to_trans_mesh_perp(Root, WallTriangles, RKnots, ZKnots):

	if(DEBUG > 0):	print("assign_triangles_to_trans_mesh_perp")

	Zones = Root.Config.Zones

	if(DEBUG > 2):
		PointsIR = np.empty(0,dtype='f8')
		PointsIZ = np.empty(0,dtype='f8')
		PointsWR = []; 	PointsWZ = []
		Points1R = []; 	Points1Z = []
		Points2R = [];	Points2Z = []
		Points3R = [];	Points3Z = []
		PointsnR = [];	PointsnZ = []
		RKnotsTri	   = np.array([RKnots[WallTriangles.p1],RKnots[WallTriangles.p2],RKnots[WallTriangles.p3],RKnots[WallTriangles.p1]])
		ZKnotsTri	   = np.array([ZKnots[WallTriangles.p1],ZKnots[WallTriangles.p2],ZKnots[WallTriangles.p3],ZKnots[WallTriangles.p1]])

	nZones   = len(Zones)


	for k in range(nZones):
		Nx = Zones[k].Nx
		Nz = Zones[k].Nz

		Zones[k].list_tri_perp_nnums = np.zeros((Nx, Nz), dtype='i4')
		Zones[k].list_tri_perp_nums	 = np.zeros((Nx, Nz, 16), dtype='i4')


		ii,jj = np.where(((Zones[k].Chi2[1:Nx+1, 1:Nz+1] == 0) & (Zones[k].Chi2[2:Nx+2, 1:Nz+1] == 1)) |
						 ((Zones[k].Chi2[1:Nx+1, 1:Nz+1] == 0) & (Zones[k].Chi2[0:Nx,   1:Nz+1] == 1)))		#cas 01 et 10 parallel

		if len(ii) > 0:
			for l in range(len(ii)):
				ncase   = 0									#%mesh a considerer: %essayer avec n_case maille autour
				PointOK = False
				while(not PointOK):
					ncase  += 3
					List	= np.empty((2*ncase+1,3), dtype='i4')
					ListTri = np.empty(4*ncase+2, dtype='i4')
	
					k_ = k
					i_ = ii[l]
					j_ = jj[l]
					iListL     = ncase
					List[iListL,:] = [k_, i_, j_]
					for n in range(ncase):
						k_2 = k_
						i_2 = i_
						j_2 = j_
						k_  = Zones[k_2].point_north_k[i_2,j_2]
						i_  = Zones[k_2].point_north_i[i_2,j_2]
						j_  = Zones[k_2].point_north_j[i_2,j_2]
						if k_ == -1: break									#test proximity
	
						d = np.sqrt((Zones[k_2].gridRc[i_2,j_2] - Zones[k_].gridRc[i_,j_])**2 + (Zones[k_2].gridZc[i_2,j_2] - Zones[k_].gridZc[i_,j_])**2)
						if d > TRANS_MESH_DMAX: break
	
						iListL		 -= 1
						List[iListL,:] = [k_, i_, j_]
	
					k_ = k
					i_ = ii[l]
					j_ = jj[l]
					iListH     = ncase
					for n in range(ncase):
						k_2 = k_
						i_2 = i_
						j_2 = j_
						k_  = Zones[k_2].point_south_k[i_2,j_2]
						i_  = Zones[k_2].point_south_i[i_2,j_2]
						j_  = Zones[k_2].point_south_j[i_2,j_2]
						if k_ == -1: break									#test proximity
	
						d = np.sqrt((Zones[k_2].gridRc[i_2,j_2] - Zones[k_].gridRc[i_,j_])**2 + (Zones[k_2].gridZc[i_2,j_2] - Zones[k_].gridZc[i_,j_])**2)
						if d > TRANS_MESH_DMAX: break
	
						iListH		 += 1
						List[iListH,:] = [k_, i_, j_]
	
					List = List[iListL:iListH+1,:]
	
					if(DEBUG > 2):
						for kk in range(len(List[:,0])):
							PointsWR = np.append(PointsWR, Zones[List[kk,0]].gridRc[List[kk,1],List[kk,2]])
							PointsWZ = np.append(PointsWZ, Zones[List[kk,0]].gridZc[List[kk,1],List[kk,2]])
					
					nListTri = 0
					for n in range(List.shape[0]):
						iTri = np.where((WallTriangles.k==List[n,0]) & (WallTriangles.i==List[n,1]) & (WallTriangles.j==List[n,2])); iTri=iTri[0]
						if len(iTri) > 0:
							nTri = len(iTri)
							ListTri[nListTri: nListTri+nTri] = iTri
							nListTri += nTri

					if((nListTri == 0) and (ncase > 50)):
#						look at 2D vicinity last chance
						k_ = k
						i_ = ii[l]
						j_ = jj[l]
						List2D=np.array([[k_,i_,j_],[k_,i_-1,j_],[k_,i_+1,j_],[k_,i_,j_-1],[k_,i_,j_+1],[k_,i_-1,j_-1],[k_,i_+1,j_-1],[k_,i_-1,j_+1],[k_,i_+1,j_+1]])
						ListTri = np.empty(2*9, dtype='i4')
						for n in range(List2D.shape[0]):
							iTri = np.where((WallTriangles.k == List2D[n,0]) & (WallTriangles.i == List2D[n,1]) & (WallTriangles.j == List2D[n,2]))[0]
							if(len(iTri) > 0):
								nTri = len(iTri)
								ListTri[nListTri: nListTri+nTri] = iTri
								nListTri += nTri

					if(nListTri > 0):
						PointOK = True
						ListTri = ListTri[0:nListTri]

						if(DEBUG > 2):
							if(nListTri == 1):
								Points1R, Points1Z = append_RZ(Points1R, Points1Z, RKnotsTri[:, ListTri[0]], ZKnotsTri[:, ListTri[0]])
							elif(nListTri == 2):
								Points2R, Points2Z = append_RZ(Points2R, Points2Z, RKnotsTri[:, ListTri[0]], ZKnotsTri[:, ListTri[0]])
							elif(nListTri == 3):
								Points3R, Points3Z = append_RZ(Points3R, Points3Z, RKnotsTri[:, ListTri[0]], ZKnotsTri[:, ListTri[0]])
							else:
								PointsnR, PointsnZ = append_RZ(PointsnR, PointsnZ, RKnotsTri[:, ListTri[0]], ZKnotsTri[:, ListTri[0]])
		
						Zones[k].list_tri_perp_nnums[ii[l],jj[l]] 			= nListTri
						Zones[k].list_tri_perp_nums[ii[l],jj[l],:] 		 	= -1
						Zones[k].list_tri_perp_nums[ii[l],jj[l],0:nListTri] = ListTri
	
					elif(ncase > 50):
						print("\tError no triangle found on south-north direction!")
		
						print("l, k,ii[l], jj[l]=",l, k,ii[l], jj[l])
						RZKnots  = np.array([RKnots, ZKnots])
						TriNodes = np.array([WallTriangles.p1, WallTriangles.p2, WallTriangles.p3]).T

						k_ = k
						i_ = ii[l]
						j_ = jj[l]
						WallCellR = np.array([Zones[k].gridR[i_,j_], Zones[k].gridR[i_+1,j_], Zones[k].gridR[i_+1,j_+1], Zones[k].gridR[i_,j_+1], Zones[k].gridR[i_,j_]])
						WallCellZ = np.array([Zones[k].gridZ[i_,j_], Zones[k].gridZ[i_+1,j_], Zones[k].gridZ[i_+1,j_+1], Zones[k].gridZ[i_,j_+1], Zones[k].gridZ[i_,j_]])

						if(List.shape[0] > 1):
							NearCellsR = []
							NearCellsZ = []
							for n in range(List.shape[0]-1):
								k_ = List[n,0]
								i_ = List[n,1]
								j_ = List[n,2]
								if((i_ < Zones[k_].gridR.shape[0] -1) and (j_ < Zones[k_].gridR.shape[1] -1)):
									NearCellsR, NearCellsZ = append_RZ(NearCellsR, NearCellsZ, \
																	np.array([Zones[k_].gridR[i_,j_], Zones[k_].gridR[i_+1,j_], Zones[k_].gridR[i_+1,j_+1], Zones[k_].gridR[i_,j_+1], Zones[k_].gridR[i_,j_]]), \
																	np.array([Zones[k_].gridZ[i_,j_], Zones[k_].gridZ[i_+1,j_], Zones[k_].gridZ[i_+1,j_+1], Zones[k_].gridZ[i_,j_+1], Zones[k_].gridZ[i_,j_]]))

							Xarrs,  Yarrs, lines = set_XYlines_walls(Root.Config, wall_line="c--")

							Xarrs.extend([RZKnots,  NearCellsR, WallCellR])
							Yarrs.extend([TriNodes, NearCellsZ, WallCellZ])
							lines.extend(["b.-", "g-", "r-"])
						else:
							Xarrs,  Yarrs, lines = set_XYlines_walls(Root.Config, wall_line="c--")

							Xarrs.extend([RZKnots,  WallCellR])
							Yarrs.extend([TriNodes, WallCellZ])
							lines.extend(["b-", "r-"])

						LinesData = [Xarrs, Yarrs, lines]
						choice = plot_and_ask(Root, LinesData=LinesData, title="Red quadrangles without triangle")

						return False

			if(DEBUG > 2):
				PointsIR = np.append(PointsIR, Zones[k].gridRc[ii,jj])
				PointsIZ = np.append(PointsIZ, Zones[k].gridZc[ii,jj])
				
	if(DEBUG > 2):
		print("\tplot nodes & first found triangle in quadrangles crossing wall flux directon")
		print("\t\tblue_o:  centers of soledge cells cutted by the walll")
		print("\t\treed:    centers of nearest soledge cells in poloidal directon")
		print("\t\tgreen: 	one triangles in nearest cells")
		print("\t\tcyan: 	two triangles in nearest cells")
		print("\t\tmagenta: three triangles in nearest cells")
		print("\t\tyellow: 	more than three triangles in nearest cells")

		Xarrs, Yarrs, lines = set_XYlines_walls(Root.Config)
		Xarrs.extend([PointsIR, PointsWR, Points1R, Points2R, Points3R, PointsnR])
		Yarrs.extend([PointsIZ, PointsWZ, Points1Z, Points2Z, Points3Z, PointsnZ])
		lines.extend(['bo', 'r.', 'g-', 'c-', 'm-', 'y-'])
		LinesData = [Xarrs, Yarrs, lines]
		choice = plot_and_ask(Root, LinesData=LinesData, title="Mesh & Edges points")

	if(DEBUG > 0):	print("assign_triangles_to_trans_mesh: Completed")

	return True
