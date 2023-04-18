import numpy				as np
from routines.globals		import DEBUG

#=========================================================
# This remapping parameters
#=========================================================

def remapping_perp(Zones, WallTriangles, RKnots, ZKnots):

	if(DEBUG > 0):	print("remapping_perp")

	nWallTriangles = len(WallTriangles)
	nZones = len(Zones)

	if(DEBUG > 1):	print("\tremapping South")

	for k in range(nZones):
		ii,jj = np.where(((Zones[k].Chi2[1:-1, 1:-1] == 0) & (Zones[k].Chi2[2:,  1:-1] == 1)) |
						 ((Zones[k].Chi2[1:-1, 1:-1] == 0) & (Zones[k].Chi2[:-2, 1:-1] == 1)))		#cas 01 et 10 parallel
		for l in range(len(ii)):
			i = ii[l]
			j = jj[l]
			nListTri = Zones[k].list_tri_s_nnums[i,j]
			for iListTri in range(nListTri):
				n = Zones[k].list_tri_s_nums[i,j,iListTri]
				if(WallTriangles.weight_s[n] > 0.):
					print("\tWARNING in remapping South")
					print("k=",k+1," i=",i+1," j=",j+1)
					print("list_tri_s_nums[k,i,j]       = ",Zones[k].list_tri_s_nums[i,j,0:nListTri]+1)
					print("n                            = ",n+1)
					print("list_s[n,:]                  = ",WallTriangles.list_s[n, :]+1)
					print("weight_s[n]                  = ",WallTriangles.weight_s[n])
				else:
					WallTriangles.list_s[n, :] = [i,j,k]
					WallTriangles.weight_s[n] = Zones[k].list_tri_s_weights[i,j,iListTri]

	if(DEBUG > 1):	print("\tremapping North")

	for k in range(nZones):
		ii,jj = np.where(((Zones[k].Chi2[1:-1, 1:-1] == 0) & (Zones[k].Chi2[2:, 1:-1] == 1)) |
						 ((Zones[k].Chi2[1:-1, 1:-1] == 0) & (Zones[k].Chi2[:-2,1:-1] == 1)))		#cas 01 et 10 parallel
		for l in range(len(ii)):
			i = ii[l]
			j = jj[l]
			nListTri = Zones[k].list_tri_n_nnums[i,j]
			for iListTri in range(nListTri):
				n = Zones[k].list_tri_n_nums[i,j,iListTri]
				if(WallTriangles.weight_n[n] > 0.):
					print("\tWARNING in remapping North")
					print("k=",k+1," i=",i+1," j=",j+1)
					print("list_tri_n_nums[k,i,j]       = ",Zones[k].list_tri_n_nums[i,j,0:nListTri]+1)
					print("n                            = ",n+1)
					print("list_n[n,:]                  = ",WallTriangles.list_n[n, :]+1)
					print("weight_n[n]                  = ",WallTriangles.weight_n[n]+1)
				else:
					WallTriangles.list_n[n, :] = [i,j,k]
					WallTriangles.weight_n[n] = Zones[k].list_tri_n_weights[i,j,iListTri]

	if(DEBUG > 0):	print("remapping_perp: Completed")
