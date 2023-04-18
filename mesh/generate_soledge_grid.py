import numpy						as np
import scipy.interpolate			as interp
from matplotlib.path				import Path
from mesh.clean_chis				import clean_chis
from routines.find_closest_segment	import find_closest_segment
from routines.globals				import DEBUG, EXTERNAL_PLASMA_WALL

def generate_soledge_grid(Config, CenterPointChi=True):

	if(DEBUG > 0):	print("generate_soledge_grid")

#	define Wall path

	Zones = Config.Zones
	
	Config.Bphi2D = np.abs(Config.Bphi2D)
	f_Br	  	  = interp.RectBivariateSpline(Config.r2D[0,:], Config.z2D[:,0], Config.Br2D.T)										#RectBivariateSpline wants [x,y] order
	f_Bz	  	  = interp.RectBivariateSpline(Config.r2D[0,:], Config.z2D[:,0], Config.Bz2D.T)
	f_Bphi	  	  = interp.RectBivariateSpline(Config.r2D[0,:], Config.z2D[:,0], Config.Bphi2D.T)

	for k in range(len(Zones)):
		Zones[k].gridRc = 0.25*(Zones[k].gridR[:-1,:-1] + Zones[k].gridR[:-1,1:] + Zones[k].gridR[1:,1:] + Zones[k].gridR[1:,:-1])
		Zones[k].gridZc = 0.25*(Zones[k].gridZ[:-1,:-1] + Zones[k].gridZ[:-1,1:] + Zones[k].gridZ[1:,1:] + Zones[k].gridZ[1:,:-1])

		Zones[k].Nx		= Zones[k].gridRc.shape[0]
		Zones[k].Nz		= Zones[k].gridRc.shape[1]

		"""
		nTot	 	 = Zones[k].gridRc.size
		RZ1		 	 = np.empty((nTot,2), dtype='f8')
		RZ1[:,0]	 = Zones[k].gridRc.reshape(nTot)
		RZ1[:,1]	 = Zones[k].gridZc.reshape(nTot)

		IsIn		 = Config.WallPath.contains_points(RZ1)								#find point iside wall contour
		Zones[k].Chi = np.where(IsIn, np.zeros(nTot, dtype='f8'), np.ones(nTot, dtype='f8'))
		Zones[k].Chi = Zones[k].Chi.reshape(Zones[k].gridRc.shape)
		"""

		Zones[k].Br	  = np.reshape(f_Br.ev(Zones[k].gridRc.reshape(-1),   Zones[k].gridZc.reshape(-1)), Zones[k].gridRc.shape)	#RectBivariateSpline want 1D array in evaluation
		Zones[k].Bz	  = np.reshape(f_Bz.ev(Zones[k].gridRc.reshape(-1),   Zones[k].gridZc.reshape(-1)), Zones[k].gridRc.shape)
		Zones[k].Bphi = np.reshape(f_Bphi.ev(Zones[k].gridRc.reshape(-1), Zones[k].gridZc.reshape(-1)), Zones[k].gridRc.shape)

	define_chi(Config, CenterPointChi=CenterPointChi)
	clean_chis(Config)
	
	if(DEBUG > 0):	print("generate_soledge_grid: Completed")


#	This routines se Chi=1 if all mesh corner points are on the wall or outside the wall
#	===============================================

def define_chi(Config, CenterPointChi=True):
	
	if(DEBUG > 0):	print("define_chi")

	Zones = Config.Zones

	if(not CenterPointChi):
		Index_i = [0,1,1,0]
		Index_j = [0,0,1,1]
		nTotNearWall = 0

	nTotNearWall = 0
	for k in range(len(Zones)):
		if(CenterPointChi):
			Chi  = np.zeros(Zones[k].gridRc.size, dtype='f8')
		else:
			Nx	 = Zones[k].Nx
			Nz	 = Zones[k].Nz
			Chi  = np.zeros((Zones[k].gridRc.size,4), dtype='f8')

		iPwalls = Config.iPwalls
		for iPwall in range(len(iPwalls)):
			iWall= iPwalls[iPwall]
			Wall = Config.Walls[iWall]
			WallPath = Wall.WallPath

			if(CenterPointChi):																				#Like Matlab chi based on center mesh
				RZ1  = np.array([Zones[k].gridRc.reshape(-1), Zones[k].gridZc.reshape(-1)]).T
				IsIn = WallPath.contains_points(RZ1)														#find points iside wall contour

				if(Wall.Type == EXTERNAL_PLASMA_WALL):
					Chi = np.where((Chi == 0.) & IsIn , 0., 1.)					#To be in point must be inside the wall 
				else:
					Chi = np.where((Chi == 1.) | IsIn , 1., 0.)

			else:																							#Just one point inside the wall (not on the wall)
				for iNode in range(4):
					RZ1  = np.array([Zones[k].gridR[Index_i[iNode]:Nx+Index_i[iNode], Index_j[iNode]:Nz+Index_j[iNode]].reshape(-1),
									 Zones[k].gridZ[Index_i[iNode]:Nx+Index_i[iNode], Index_j[iNode]:Nz+Index_j[iNode]].reshape(-1)]).T
				
					IsIn = WallPath.contains_points(RZ1)														#find points iside wall contour

					Radius = 1.e-4
					if(Wall.Type == EXTERNAL_PLASMA_WALL):
						if(WallPath.contains_point([Wall.Rwall[0],Wall.Zwall[0]], radius= -Radius)): Radius = -Radius				#reverse radius if needed
					if(Wall.Type == INTERNAL_PLASMA_WALL):
						if(WallPath.contains_point([Wall.Rwall[0],Wall.Zwall[0]], radius= Radius)): Radius = -Radius				#reverse radius if needed
					IsNear = np.where((np.where(WallPath.contains_points(RZ1, radius= Radius), 1, 0) - \
									   np.where(WallPath.contains_points(RZ1, radius=-Radius), 1, 0)) > 0)[0]	#find points close to wall contour

					nTotNearWall += len(IsNear)
					nOnWall  = 0																				#find points on wall contour
					if(len(IsNear) > 0):
						OnWall	  = np.zeros(len(IsIn), dtype = 'i4')
						for i in range(len(IsNear)):
							iKnot = IsNear[i]
							d, iSeg = find_closest_segment(RZ1[iKnot,0], RZ1[iKnot,1], Rwall, Zwall)
		#					print("\tcheck as near R,Z=",RZ1[iKnot,0], RZ1[iKnot,1])
							if(d < 1e-6):
		#						print("\t\tIt is on wall")
								OnWall[iKnot] = 1
								nOnWall		 += 1

					if(Wall.Type == EXTERNAL_PLASMA_WALL):
						if(nOnWall == 0): Chi[:,iNode] = np.where((Chi[:,iNode] == 0.) & IsIn , 0., 1.)								#To be in point must be inside the wall 
						else:			  Chi[:,iNode] = np.where((Chi[:,iNode] == 0.) & (IsIn & (OnWall == 0)) , 0., 1.)			#but not on the wall
					else:
						if(nOnWall == 0): Chi[:,iNode] = np.where((Chi[:,iNode] == 1.) | IsIn , 1., 0.)
						else:			  Chi[:,iNode] = np.where((Chi[:,iNode] == 1.) | IsIn | (OnWall == 0) , 1., 0.)

				Chi = np.where(np.sum(Chi, axis=2) < 4., 0., 1.)							#To be in at least one point must be in

			Zones[k].Chi = Chi.reshape(Zones[k].gridRc.shape)

	if((not CenterPointChi) and (nTotNearWall == 0)):
		print("\tATTENTION!! COULD BE SOMETHING WRONG IN SEARCHING NODES CLOSE TO THE WALL")
		print("\tATTENTION!! COULD BE SOMETHING WRONG IN SEARCHING NODES CLOSE TO THE WALL")
		print("\tATTENTION!! COULD BE SOMETHING WRONG IN SEARCHING NODES CLOSE TO THE WALL")

	if(DEBUG > 0):	print("define_chi: Completed")
