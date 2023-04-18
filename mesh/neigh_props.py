import types
import numpy				as np
import scipy.interpolate	as interp

def neigh_props(Config):
	
	Zones = Config.Zones

	f_Br	= interp.RectBivariateSpline(Config.r2D[0,:], Config.z2D[:,0], Config.Br2D.T)										#RectBivariateSpline wants [x,y] order
	f_Bz	= interp.RectBivariateSpline(Config.r2D[0,:], Config.z2D[:,0], Config.Bz2D.T)
	f_Bphi	= interp.RectBivariateSpline(Config.r2D[0,:], Config.z2D[:,0], Config.Bphi2D.T)

	nZones = len(Zones)
	for k in range(nZones):
		Zones[k].NorthP			= types.SimpleNamespace()
		Zones[k].NorthP.R		= (Zones[k].gridR[-1, 1:] + Zones[k].gridR[-1, :-1])*0.5
		Zones[k].NorthP.Z		= (Zones[k].gridZ[-1, 1:] + Zones[k].gridZ[-1, :-1])*0.5
		Zones[k].NorthP.Br		= f_Br.ev(Zones[k].NorthP.R,   Zones[k].NorthP.Z)
		Zones[k].NorthP.Bz		= f_Bz.ev(Zones[k].NorthP.R,   Zones[k].NorthP.Z)
		Zones[k].NorthP.Bphi	= f_Bphi.ev(Zones[k].NorthP.R, Zones[k].NorthP.Z)
		Zones[k].NorthP.B		= np.sqrt(Zones[k].NorthP.Br**2 + Zones[k].NorthP.Bz**2 + Zones[k].NorthP.Bphi**2)

		Zones[k].SouthP			= types.SimpleNamespace()
		Zones[k].SouthP.R		= (Zones[k].gridR[0, 1:] + Zones[k].gridR[0, :-1])*0.5
		Zones[k].SouthP.Z		= (Zones[k].gridZ[0, 1:] + Zones[k].gridZ[0, :-1])*0.5
		Zones[k].SouthP.Br		= f_Br.ev(Zones[k].SouthP.R,   Zones[k].SouthP.Z)
		Zones[k].SouthP.Bz		= f_Bz.ev(Zones[k].SouthP.R,   Zones[k].SouthP.Z)
		Zones[k].SouthP.Bphi	= f_Bphi.ev(Zones[k].SouthP.R, Zones[k].SouthP.Z)
		Zones[k].SouthP.B		= np.sqrt(Zones[k].SouthP.Br**2 + Zones[k].SouthP.Bz**2 + Zones[k].SouthP.Bphi**2)

		Zones[k].EastP			= types.SimpleNamespace()
		Zones[k].EastP.R		= (Zones[k].gridR[1:, -1] + Zones[k].gridR[:-1, -1])*0.5
		Zones[k].EastP.Z		= (Zones[k].gridZ[1:, -1] + Zones[k].gridZ[:-1, -1])*0.5
		Zones[k].EastP.Br		= f_Br.ev(Zones[k].EastP.R,   Zones[k].EastP.Z)
		Zones[k].EastP.Bz		= f_Bz.ev(Zones[k].EastP.R,   Zones[k].EastP.Z)
		Zones[k].EastP.Bphi		= f_Bphi.ev(Zones[k].EastP.R, Zones[k].EastP.Z)
		Zones[k].EastP.B		= np.sqrt(Zones[k].EastP.Br**2 + Zones[k].EastP.Bz**2 + Zones[k].EastP.Bphi**2)

		Zones[k].WestP			= types.SimpleNamespace()
		Zones[k].WestP.R		= (Zones[k].gridR[1:, 0] + Zones[k].gridR[:-1, 0])*0.5
		Zones[k].WestP.Z		= (Zones[k].gridZ[1:, 0] + Zones[k].gridZ[:-1, 0])*0.5
		Zones[k].WestP.Br		= f_Br.ev(Zones[k].WestP.R,   Zones[k].WestP.Z)
		Zones[k].WestP.Bz		= f_Bz.ev(Zones[k].WestP.R,   Zones[k].WestP.Z)
		Zones[k].WestP.Bphi		= f_Bphi.ev(Zones[k].WestP.R, Zones[k].WestP.Z)
		Zones[k].WestP.B		= np.sqrt(Zones[k].WestP.Br**2 + Zones[k].WestP.Bz**2 + Zones[k].WestP.Bphi**2)

