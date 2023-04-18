import numpy				as np
from routines.globals		import DEBUG

#=========================================================
# This remapping parameters
#=========================================================

def remapping_para(Zones, WallTriangles, RKnots, ZKnots):

	if(DEBUG > 0):	print("remapping_para")

	nWallTriangles = len(WallTriangles)

	nZones = len(Zones)

	if(DEBUG > 1):	print("\tremapping East")

	nFound	= np.zeros(nWallTriangles, dtype='i4')
	for k in range(nZones):
		ii,jj = np.where(((Zones[k].Chi2[1:-1, 1:-1] == 0) & (Zones[k].Chi2[1:-1, 2:  ] == 1)) |
						 ((Zones[k].Chi2[1:-1, 1:-1] == 0) & (Zones[k].Chi2[1:-1, :-2] == 1)))		#cas 01 et 10 parallel

		for l in range(len(ii)):
			i = ii[l]
			j = jj[l]
			nListTri = Zones[k].list_tri_e_nnums[i,j]
			for iListTri in range(nListTri):
				n = Zones[k].list_tri_e_nums[i,j,iListTri]
				if(WallTriangles.weight_e[n] > 0.):
					print("\tERROR in remapping East")
					print("k=",k+1," i=",i+1," j=",j+1)
					print("list_tri_e_nums[k,i,j] = ",Zones[k].list_tri_e_nums[i,j,0:nListTri]+1)
					print("n                      = ",n+1)
					print("list_e[n,:]            = ",WallTriangles.list_e[n, :]+1)
					print("weight_e[n]            = ",WallTriangles.weight_e[n]+1)
				else:
					WallTriangles.list_e[n, :] = [i,j,k]
					WallTriangles.weight_e[n]  = Zones[k].list_tri_e_weights[i,j,iListTri]

	if(DEBUG > 1):	print("\tremapping West")

	nFound	= np.zeros(nWallTriangles, dtype='i4')
	for k in range(nZones):
		ii,jj = np.where(((Zones[k].Chi2[1:-1, 1:-1] == 0) & (Zones[k].Chi2[1:-1, 2:  ] == 1)) |
						 ((Zones[k].Chi2[1:-1, 1:-1] == 0) & (Zones[k].Chi2[1:-1,  :-2] == 1)))		#cas 01 et 10 parallel

		for l in range(len(ii)):
			i = ii[l]
			j = jj[l]
			nListTri = Zones[k].list_tri_w_nnums[i,j]
			for iListTri in range(nListTri):
				n = Zones[k].list_tri_w_nums[i,j,iListTri]
				if(WallTriangles.weight_w[n] > 0.):
					print("\tERROR in remapping West")
					print("k=",k+1," i=",i+1," j=",j+1)
					print("list_tri_w_nums[k,i,j] = ",Zones[k].list_tri_w_nums[i,j,0:nListTri]+1)
					print("n                      = ",n+1)
					print("list_w[n,:]            = ",WallTriangles.list_w[n, :]+1)
					print("weight_w[n]            = ",WallTriangles.weight_w[n]+1)
				else:
					WallTriangles.list_w[n, :] = [i,j,k]
					WallTriangles.weight_w[n]  = Zones[k].list_tri_w_weights[i,j,iListTri]

	if(DEBUG > 0):	print("remapping_para: Completed")

