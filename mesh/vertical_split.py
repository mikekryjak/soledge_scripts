import types
import numpy			as np
from math 				import floor
from mesh.P_routines 	import copyP

def vertical_split(Zones, ZonesToSplit, nSplits):
	
	nz		  = ZonesToSplit
	ns		  = len(ZonesToSplit)
	nZones	  = len(Zones)
	nZonesOld = nZones
	
	if(Zones[ZonesToSplit[0]].Neighbour.south < -1):
		southToNorth = True
	else:
		southToNorth = False

	for k in range(ns):
		zb		= np.copy(Zones[nz[k]].zb)
		xb		= np.copy(Zones[nz[k]].xb)
		z		= np.copy(Zones[nz[k]].z)
		x		= np.copy(Zones[nz[k]].x)
		
		Nx		= Zones[nz[k]].Nx
		Nz		= Zones[nz[k]].Nz
		Chi		= np.copy(Zones[nz[k]].Chi)
		R		= np.copy(Zones[nz[k]].gridRc)
		Z		= np.copy(Zones[nz[k]].gridZc)
		Rr		= np.copy(Zones[nz[k]].gridR)
		Zr		= np.copy(Zones[nz[k]].gridZ)
		Br		= np.copy(Zones[nz[k]].Br)
		Bz		= np.copy(Zones[nz[k]].Bz)
		Bphi	= np.copy(Zones[nz[k]].Bphi)
		
		West	= Zones[nz[k]].Neighbour.west
		East	= Zones[nz[k]].Neighbour.east
		North	= Zones[nz[k]].Neighbour.north
		South	= Zones[nz[k]].Neighbour.south
		SouthP	= copyP(Zones[nz[k]].SouthP)
		NorthP	= copyP(Zones[nz[k]].NorthP)
		EastP	= copyP(Zones[nz[k]].EastP)
		WestP	= copyP(Zones[nz[k]].WestP)
		magz	= np.copy(Zones[nz[k]].magz)

		Nznew	= np.zeros(nSplits, dtype='i4')
		Nxnew	= np.zeros(nSplits, dtype='i4')

		iZnew	= np.linspace(0,Nz,nSplits+1).astype(np.int)
		for k2 in range(nSplits):
			Nznew[k2]		 = iZnew[k2+1]-iZnew[k2]
			Nznew[nSplits-1] = Nz-np.sum(Nznew[:nSplits-1])
			if(k2 ==  0):
				Zones[nz[k]].zb			= zb[:, :Nznew[k2]+1]
				Zones[nz[k]].xb			= xb[:, :Nznew[k2]+1]
				Zones[nz[k]].z			= z[:, :Nznew[k2]]
				Zones[nz[k]].x			= x[:, :Nznew[k2]]
				Zones[nz[k]].gridRc		= R[:, :Nznew[k2]]
				Zones[nz[k]].gridZc		= Z[:, :Nznew[k2]]
				Zones[nz[k]].gridR		= Rr[:, :Nznew[k2]+1]
				Zones[nz[k]].gridZ		= Zr[:, :Nznew[k2]+1]
				Zones[nz[k]].Br			= Br[:, :Nznew[k2]]
				Zones[nz[k]].Bz			= Bz[:, :Nznew[k2]]
				Zones[nz[k]].Bphi		= Bphi[:, :Nznew[k2]]
				Zones[nz[k]].Chi		= Chi[:, :Nznew[k2]]
				
				Zones[nz[k]].WestP		= copyP(WestP)

				Zones[nz[k]].NorthP.Br	= NorthP.Br[:Nznew[k2]]
				Zones[nz[k]].NorthP.Bz	= NorthP.Bz[:Nznew[k2]]
				Zones[nz[k]].NorthP.Bphi= NorthP.Bphi[:Nznew[k2]]
				Zones[nz[k]].NorthP.R	= NorthP.R[:Nznew[k2]]
				Zones[nz[k]].NorthP.Z	= NorthP.Z[:Nznew[k2]]
				
				Zones[nz[k]].SouthP.Br	= SouthP.Br[:Nznew[k2]]
				Zones[nz[k]].SouthP.Bz	= SouthP.Bz[:Nznew[k2]]
				Zones[nz[k]].SouthP.Bphi= SouthP.Bphi[:Nznew[k2]]
				Zones[nz[k]].SouthP.R	= SouthP.R[:Nznew[k2]]
				Zones[nz[k]].SouthP.Z	= SouthP.Z[:Nznew[k2]]
				if(k2 == nSplits-1):
					Zones[nz[k]].EastP		= copyP(EastP)
					Zones[nz[k]].Neighbour.east	= East
				else:
					Zones[nz[k]].EastP		= types.SimpleNamespace()
					Zones[nz[k]].EastP.Br	= Br[:, Nznew[k2]]
					Zones[nz[k]].EastP.Bz	= Bz[:, Nznew[k2]]
					Zones[nz[k]].EastP.Bphi	= Bphi[:, Nznew[k2]]
					Zones[nz[k]].EastP.R	= R[:, Nznew[k2]]
					Zones[nz[k]].EastP.Z	= Z[:, Nznew[k2]]
					Zones[nz[k]].Neighbour.east= nZones

				Zones[nz[k]].magz	= np.copy(magz)
				Zones[nz[k]].Nx 	= Nx
				Zones[nz[k]].Nz 	= Nznew[k2]
				curb = Nznew[k2]+1
				cur  = Nznew[k2]
			else:
				nZones = nZones+1
				Zones.append(types.SimpleNamespace())
				Zones[-1].zb			= zb[:, curb-1:curb+Nznew[k2]]
				Zones[-1].xb			= xb[:, curb-1:curb+Nznew[k2]]
				Zones[-1].z				= z[:,cur:cur+Nznew[k2]]
				Zones[-1].x				= x[:,cur:cur+Nznew[k2]]
				Zones[-1].gridRc		= R[:,cur:cur+Nznew[k2]]
				Zones[-1].gridZc		= Z[:,cur:cur+Nznew[k2]]
				Zones[-1].gridR			= Rr[:,cur:cur+Nznew[k2]+1]
				Zones[-1].gridZ			= Zr[:,cur:cur+Nznew[k2]+1]
				Zones[-1].Br			= Br[:,cur:cur+Nznew[k2]]
				Zones[-1].Bz			= Bz[:,cur:cur+Nznew[k2]]
				Zones[-1].Bphi			= Bphi[:,cur:cur+Nznew[k2]]
				Zones[-1].Chi			= Chi[:,cur:cur+Nznew[k2]]
				
				Zones[-1].WestP			= types.SimpleNamespace()
				Zones[-1].WestP.Br 		= Br[:,cur-1]
				Zones[-1].WestP.Bz 		= Bz[:,cur-1]
				Zones[-1].WestP.Bphi	= Bphi[:,cur-1]
				Zones[-1].WestP.R		= R[:,cur-1]
				Zones[-1].WestP.Z		= Z[:,cur-1]
				
				Zones[-1].NorthP		= types.SimpleNamespace()
				Zones[-1].NorthP.Br 	= NorthP.Br[cur:cur+Nznew[k2]]
				Zones[-1].NorthP.Bz 	= NorthP.Bz[cur:cur+Nznew[k2]]
				Zones[-1].NorthP.Bphi	= NorthP.Bphi[cur:cur+Nznew[k2]]
				Zones[-1].NorthP.R		= NorthP.R[cur:cur+Nznew[k2]]
				Zones[-1].NorthP.Z		= NorthP.Z[cur:cur+Nznew[k2]]
				
				Zones[-1].SouthP		= types.SimpleNamespace()
				Zones[-1].SouthP.Br		= SouthP.Br[cur:cur+Nznew[k2]]
				Zones[-1].SouthP.Bz		= SouthP.Bz[cur:cur+Nznew[k2]]
				Zones[-1].SouthP.Bphi	= SouthP.Bphi[cur:cur+Nznew[k2]]
				Zones[-1].SouthP.R		= SouthP.R[cur:cur+Nznew[k2]]
				Zones[-1].SouthP.Z		= SouthP.Z[cur:cur+Nznew[k2]]
				
				Zones[-1].Neighbour		= types.SimpleNamespace()
				if(k2 == nSplits-1):
					Zones[-1].EastP = copyP(EastP)
					Zones[-1].Neighbour.east= East
				else:
					Zones[-1].EastP			= types.SimpleNamespace()
					Zones[-1].EastP.Br 		= Br[:,cur+Nznew[k2]]
					Zones[-1].EastP.Bz 		= Bz[:,cur+Nznew[k2]]
					Zones[-1].EastP.Bphi	= Bphi[:,cur+Nznew[k2]]
					Zones[-1].EastP.R		= R[:,cur+Nznew[k2]]
					Zones[-1].EastP.Z		= Z[:,cur+Nznew[k2]]
					Zones[-1].Neighbour.east= nZones

				if(k2 == 1):
					Zones[-1].Neighbour.west= ZonesToSplit[k]
				else:
					Zones[-1].Neighbour.west= nZones-2

				if(southToNorth):
					if(k == 0):
						Zones[-1].Neighbour.south= Zones[ZonesToSplit[0]].Neighbour.south
					else:
						Zones[-1].Neighbour.south= nZones - nSplits

					if(k == ns-1):
						Zones[-1].Neighbour.north= Zones[ZonesToSplit[-1]].Neighbour.north
					else:
						Zones[-1].Neighbour.north= nZones + (nSplits-2)

				else:
					if(k == 0):
						Zones[-1].Neighbour.north= Zones[ZonesToSplit[0]].Neighbour.north
					else:
						Zones[-1].Neighbour.north= nZones - nSplits
						
					if(k == ns-1):
						Zones[-1].Neighbour.south= Zones[ZonesToSplit[-1]].Neighbour.south
					else:
						Zones[-1].Neighbour.south= nZones + (nSplits-2)

				Zones[-1].magz	= np.array([magz[0], magz[1], magz[2]+cur], dtype='i4')
				Zones[-1].Nx	= Nx
				Zones[-1].Nz	= Nznew[k2]
				curb			= curb + Nznew[k2]
				cur				= cur  + Nznew[k2]

		for ko in range(nZonesOld):
			if(Zones[ko].Neighbour.west == ZonesToSplit[k]):
				Zones[ko].Neighbour.west = nZones-1
				
	
