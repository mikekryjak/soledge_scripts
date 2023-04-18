import numpy as np
from eirene.generate_wall_triangles	import generate_wall_triangles

def generate_vacuum_wall_triangles(Triangles, OldWallTriangles):

	print("generate_vacuum_wall_triangles")

	WallTriangles	 = generate_wall_triangles(Triangles)

	TriNum			= WallTriangles.ntri
	oldTriNum		= OldWallTriangles.ntri
	nWallTriangles	= len(WallTriangles.ntri)

	for iTri in range(nWallTriangles):
		ntri	= TriNum[iTri]
		iTriOld	= np.where(oldTriNum == ntri)[0]
		if(len(iTriOld) == 1):
			iTriOld = iTriOld[0]
			WallTriangles.list_e[iTri,:] = OldWallTriangles.list_e[iTriOld,:] 
			WallTriangles.list_w[iTri,:] = OldWallTriangles.list_w[iTriOld,:]
			WallTriangles.list_n[iTri,:] = OldWallTriangles.list_n[iTriOld,:]
			WallTriangles.list_s[iTri,:] = OldWallTriangles.list_s[iTriOld,:]

			WallTriangles.weight_e[iTri] = OldWallTriangles.weight_e[iTriOld]
			WallTriangles.weight_w[iTri] = OldWallTriangles.weight_w[iTriOld]
			WallTriangles.weight_n[iTri] = OldWallTriangles.weight_n[iTriOld]
			WallTriangles.weight_s[iTri] = OldWallTriangles.weight_s[iTriOld]
		elif(len(iTriOld) > 1):
			print("\tMore than one old Plasma triangle for one new triangle!!!")
			print("\tlen(iTriOld)=",len(iTriOld))
			exit()

	print("generate_vacuum_wall_triangles: completed")

	return WallTriangles
