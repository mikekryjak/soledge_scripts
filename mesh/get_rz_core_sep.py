import numpy				as np
import scipy.interpolate	as interp
from routines.utils_walls	import get_spline_min_walls
from routines.globals		import DEBUG, CORE_NEIGHBOUR, WALL_NEIGHBOUR

def get_rz_core_sep(Config, core_and_sep=True, use_mag_zones = True):

	if(DEBUG > 1): print("get_rz_core_sep")

	if(use_mag_zones):
		Zones	  = Config.MagZones
		Megazones = Config.MagMegazones
	else:
		Zones	  = Config.Zones
		Megazones = Config.Megazones

#	Search Core zone

	k = 0
	while (k < len(Zones)):
		South = Zones[k].Neighbour.south
		if(South == CORE_NEIGHBOUR):
			CoreZone = k
			Rcore	 = Zones[CoreZone].gridR[0,:]
			Zcore	 = Zones[CoreZone].gridZ[0,:]
			East 	 = Zones[CoreZone].Neighbour.east
			break
		k += 1

	while (East != CoreZone):
		Rcore	 = np.append(Rcore, Zones[East].gridR[0, 1:])
		Zcore	 = np.append(Zcore, Zones[East].gridZ[0, 1:])
		East 	 = Zones[East].Neighbour.east

	CoreMegazone = Zones[CoreZone].mz

	if(not core_and_sep): 
		if(DEBUG > 1): print("get_rz_core_sep: completed")
		return Rcore, Zcore, CoreMegazone

	if(len(Config.X_points) > 0):																					#Manage X configuration
		SepMegazone  = CoreMegazone																					# Search separatrix Megazone

		north = Zones[CoreZone].Neighbour.north
		while ((north >= 0) and Megazones[Zones[north].mz].isperiodic):
			SepMegazone  = Zones[north].mz
			north = Zones[north].Neighbour.north

		if(Megazones[Zones[north].mz].isperiodic):
			print("\tError in get_rz_core_sep, not found non-periodic zone in X configuration")
			raise

		SepZone  = Megazones[SepMegazone].list[0]
		jSep	 = -1

	else:																											#Manage limiter configuration
		f_psi		 = interp.RectBivariateSpline(Config.r2D[0,:], Config.z2D[:,0], Config.flux2D.T)				#RectBivariateSpline want [x,y] order
		psi_wall_min = get_spline_min_walls(Config, f_psi)
		RMesh		 = Zones[CoreZone].gridRc[:,0]
		ZMesh		 = Zones[CoreZone].gridZc[:,0]
		iZones		 = np.array([CoreZone])

		north		= Zones[CoreZone].Neighbour.north
		while (north >= 0):
			RMesh  = np.append(RMesh, Zones[north].gridRc[:,0])
			ZMesh  = np.append(ZMesh, Zones[north].gridZc[:,0])
			iZones = np.append(iZones, north)
			north = Zones[north].Neighbour.north

		psi_mesh = f_psi.ev(RMesh, ZMesh)
		f_psi	 = 0

		jSep	 = np.argmin(np.abs(psi_mesh - psi_wall_min))

		SepZone = -1
		while (jSep > 0):
			SepZone += 1
			jSep   -= Zones[SepZone].gridZc.shape[0]

		jSep += Zones[SepZone].gridZc.shape[0]
		SepMegazone = Zones[SepZone].mz	

	Rsep	 = Zones[SepZone].gridR[jSep,:]
	Zsep	 = Zones[SepZone].gridZ[jSep,:]
	East 	 = Zones[SepZone].Neighbour.east
	while (East != SepZone):
		Rsep	 = np.append(Rsep, Zones[East].gridR[jSep, 1:])
		Zsep	 = np.append(Zsep, Zones[East].gridZ[jSep, 1:])
		East 	 = Zones[East].Neighbour.east

	RZcore = np.array([Rcore, Zcore]).T
	RZsep  = np.array([Rsep, Zsep]).T

	if(DEBUG > 1): print("get_rz_core_sep: completed")

	return RZcore, RZsep, CoreMegazone, SepMegazone, jSep