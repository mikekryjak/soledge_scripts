import numpy					as np
from interfaces.plot_and_ask	import plot_and_ask
from routines.globals			import DEBUG, MIN_AREA_TRIANGLES, MIN_ORDER_TRIANGLES, BC_WALL, BC_CORE


def check_remapping_weights(Root, Zones, WallTriangles):

	if(DEBUG > 0):	print("check_remapping_weights")

	nWallTriangles = len(WallTriangles)
	nZones = len(Zones)

	TypeFace = np.array([WallTriangles.BC1, WallTriangles.BC2, WallTriangles.BC3]).T
	TypeFace = np.where(TypeFace > 0, 1, 0)
	sType	 = np.sum(TypeFace)
	if(sType != nWallTriangles):
		print("ERROR: invalid boudary values")

	if(DEBUG > 1):	print("\tcheck mesh North")

	"""
	for i in range(nWallTriangles):
		print("iTri = {:d}, weight: n={:0.3f}, s={:0.3f}, e={:0.3f}, w={:0.3f}".format(i+1,WallTriangles.weight_n[i],WallTriangles.weight_s[i],WallTriangles.weight_e[i],WallTriangles.weight_w[i]))
	"""

	Error_list_n = np.empty(0,dtype='i4')
	Error_diff_n = np.empty(0,dtype='f8')
	Used_tri	 = np.empty(0,dtype='i4')
	for k in range(nZones):
		ii,jj = np.where((Zones[k].Chi2[1:-1, 1:-1] == 0) & (Zones[k].Chi2[2:,  1:-1] == 1))
		for l in range(len(ii)):
			tList = np.where((WallTriangles.list_n[:, 0] == ii[l]) & (WallTriangles.list_n[:, 1] == jj[l]) & (WallTriangles.list_n[:, 2] == k))[0]	
			sTypeFace = np.sum(TypeFace[tList,WallTriangles.side[tList]])
			if(sTypeFace != len(tList)):
				print("\t\tERROR in side")

			if(len(tList) > 0):
				sum_weights = np.sum(WallTriangles.weight_n[tList])
				Used_tri = np.append(Used_tri, tList)
			else:				sum_weights = 0.
			if(abs(sum_weights -1.) > 1e-6):
				if(Error_list_n.shape[0] == 0): Error_list_n = np.array([k,ii[l],jj[l]]).reshape(1,3)
				else:							Error_list_n = np.append(Error_list_n,  np.array([k,ii[l],jj[l]]).reshape(1,3), axis=0)
				Error_diff_n = np.append(Error_diff_n, sum_weights - 1.)
	if(len(Used_tri) < len(np.where(WallTriangles.list_n[:, 0] > -1)[0])):
		print("\tERROR in remapping at north at not all triangles used")
		

	if(len(Error_diff_n) > 0):
		print("\tWARNING in remapping at north at following grid points!")
		for l in range(Error_list_n.shape[0]):
			k = Error_list_n[l,0]
			i = Error_list_n[l,1]
			j = Error_list_n[l,2]
			print("\t\tAt (R,Z)           = ({:0.3f},{:0.3f})".format(Zones[k].gridRc[i,j],Zones[k].gridZc[i,j]))
			print("\t\t[k,i,j]            = [{:d},{:d},{:d}])".format(k+1,i+1,j+1))
			print("\t\tlist_tri_n_nnums   = ",Zones[k].list_tri_n_nums[i,j,:Zones[k].list_tri_n_nnums[i,j]]+1)
			print("\t\tlist_tri_n_weights = ",Zones[k].list_tri_n_weights[i,j,:Zones[k].list_tri_n_nnums[i,j]])
			print("\t\tError_diff_n       = ",Error_diff_n[l])
			
		


	if(DEBUG > 1):	print("\tcheck mesh South")

	Error_list_s = np.empty(0,dtype='i4')
	Error_diff_s = np.empty(0,dtype='f8')
	Used_tri	 = np.empty(0,dtype='i4')
	for k in range(nZones):
		ii,jj = np.where((Zones[k].Chi2[1:-1, 1:-1] == 0) & (Zones[k].Chi2[:-2,  1:-1] == 1))
		for l in range(len(ii)):
			tList = np.where((WallTriangles.list_s[:, 0] == ii[l]) & (WallTriangles.list_s[:, 1] == jj[l]) & (WallTriangles.list_s[:, 2] == k))[0]	
			sTypeFace = np.sum(TypeFace[tList,WallTriangles.side[tList]])
			if(sTypeFace != len(tList)):
				print("\t\tERROR in side")
			if(len(tList) > 0):
				sum_weights = np.sum(WallTriangles.weight_s[tList])
				Used_tri = np.append(Used_tri, tList)
			else:				sum_weights = 0.
			if(abs(sum_weights -1.) > 1e-6):
				if(Error_list_s.shape[0] == 0): Error_list_s = np.array([k,ii[l],jj[l]]).reshape(1,3)
				else:							Error_list_s = np.append(Error_list_s,  np.array([k,ii[l],jj[l]]).reshape(1,3), axis=0)
				Error_diff_s = np.append(Error_diff_s, sum_weights - 1.)

	if(len(Used_tri) < len(np.where(WallTriangles.list_s[:, 0] > -1)[0])):
		print("\tERROR in remapping at south not all triangles used")

	if(len(Error_diff_s) > 0):
		print("\tWARNING in remapping at south at following grid points")
		for l in range(Error_list_s.shape[0]):
			k = Error_list_s[l,0]
			i = Error_list_s[l,1]
			j = Error_list_s[l,2]
			print("\t\tAt (R,Z)           = ({:0.3f},{:0.3f})".format(Zones[k].gridRc[i,j],Zones[k].gridZc[i,j]))
			print("\t\t[k,i,j]            = [{:d},{:d},{:d}])".format(k+1,i+1,j+1))
			print("\t\tlist_tri_s_nnums   = ",Zones[k].list_tri_s_nums[i,j,:Zones[k].list_tri_s_nnums[i,j]]+1)
			print("\t\tlist_tri_s_weights = ",Zones[k].list_tri_s_weights[i,j,:Zones[k].list_tri_s_nnums[i,j]])
			print("\t\tError_diff_s       = ",Error_diff_s[l])

	if(DEBUG > 1):	print("\tcheck mesh East")

	Error_list_e = np.empty(0,dtype='i4')
	Error_diff_e = np.empty(0,dtype='f8')
	Used_tri	 = np.empty(0,dtype='i4')
	for k in range(nZones):
		ii,jj = np.where((Zones[k].Chi2[1:-1, 1:-1] == 0) & (Zones[k].Chi2[1:-1,  2:] == 1))
		for l in range(len(ii)):
			tList = np.where((WallTriangles.list_e[:, 0] == ii[l]) & (WallTriangles.list_e[:, 1] == jj[l]) & (WallTriangles.list_e[:, 2] == k))[0]	
			sTypeFace = np.sum(TypeFace[tList,WallTriangles.side[tList]])
			if(sTypeFace != len(tList)):
				print("\t\tERROR in side")
			if(len(tList) > 0):
				sum_weights = np.sum(WallTriangles.weight_e[tList])
				Used_tri = np.append(Used_tri, tList)
			else:				sum_weights = 0.
			if(abs(sum_weights -1.) > 1e-6):
				if(Error_list_e.shape[0] == 0): Error_list_e = np.array([k,ii[l],jj[l]]).reshape(1,3)
				else:							Error_list_e = np.append(Error_list_e,  np.array([k,ii[l],jj[l]]).reshape(1,3), axis=0)
				Error_diff_e = np.append(Error_diff_e, sum_weights - 1.)

	if(len(Used_tri) < len(np.where(WallTriangles.list_e[:, 0] > -1)[0])):
		print("\tERROR in remapping at east not all triangles used")

	if(len(Error_diff_e) > 0):
		print("\tWARNING in remapping at east at following grid points")
		for l in range(Error_list_e.shape[0]):
			k = Error_list_e[l,0]
			i = Error_list_e[l,1]
			j = Error_list_e[l,2]
			print("\t\tAt (R,Z)           = ({:0.3f},{:0.3f})".format(Zones[k].gridRc[i,j],Zones[k].gridZc[i,j]))
			print("\t\t[k,i,j]            = [{:d},{:d},{:d}])".format(k+1,i+1,j+1))
			print("\t\tlist_tri_e_nnums   = ",Zones[k].list_tri_e_nums[i,j,:Zones[k].list_tri_e_nnums[i,j]]+1)
			print("\t\tlist_tri_e_weights = ",Zones[k].list_tri_e_weights[i,j,:Zones[k].list_tri_e_nnums[i,j]])
			print("\t\tError_diff_e       = ",Error_diff_e[l])

	if(DEBUG > 1):	print("\tcheck mesh West")

	Error_list_w = np.empty(0,dtype='i4')
	Error_diff_w = np.empty(0,dtype='f8')
	Used_tri	 = np.empty(0,dtype='i4')
	for k in range(nZones):
		ii,jj = np.where((Zones[k].Chi2[1:-1, 1:-1] == 0) & (Zones[k].Chi2[1:-1, :-2] == 1))
		for l in range(len(ii)):
			tList = np.where((WallTriangles.list_w[:, 0] == ii[l]) & (WallTriangles.list_w[:, 1] == jj[l]) & (WallTriangles.list_w[:, 2] == k))[0]	
			sTypeFace = np.where(TypeFace[tList,WallTriangles.side[tList]] == 0)[0]
			if(sTypeFace != len(tList)):
				print("\t\tERROR in side")
			if(len(tList) > 0):
				sum_weights = np.sum(WallTriangles.weight_w[tList])
				Used_tri = np.append(Used_tri, tList)
			else:				sum_weights = 0.
			if(abs(sum_weights -1.) > 1e-6):
				if(Error_list_w.shape[0] == 0): Error_list_w = np.array([k,ii[l],jj[l]]).reshape(1,3)
				else:							Error_list_w = np.append(Error_list_w,  np.array([k,ii[l],jj[l]]).reshape(1,3), axis=0)
				Error_diff_w = np.append(Error_diff_w, sum_weights - 1.)

	if(len(Used_tri) < len(np.where(WallTriangles.list_w[:, 0] > -1)[0])):
		print("\tERROR in remapping at west not all triangles used")

	if(len(Error_diff_w) > 0):
		print("\tWARNING in remapping at west at following grid points")
		for l in range(Error_list_w.shape[0]):
			k = Error_list_w[l,0]
			i = Error_list_w[l,1]
			j = Error_list_w[l,2]
			print("\t\tAt (R,Z)           = ({:0.3f},{:0.3f})".format(Zones[k].gridRc[i,j],Zones[k].gridZc[i,j]))
			print("\t\t[k,i,j]            = [{:d},{:d},{:d}])".format(k+1,i+1,j+1))
			print("\t\tlist_tri_w_nnums   = ",Zones[k].list_tri_w_nums[i,j,:Zones[k].list_tri_w_nnums[i,j]]+1)
			print("\t\tlist_tri_w_weights = ",Zones[k].list_tri_w_weights[i,j,:Zones[k].list_tri_w_nnums[i,j]])
			print("\t\tError_diff_w       = ",Error_diff_w[l])

		
	if(DEBUG > 0):	print("check_remapping_weights: Completed")
