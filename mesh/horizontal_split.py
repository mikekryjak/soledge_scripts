import types
import numpy			as np
from math 				import floor
from mesh.P_routines 	import copyP

def horizontal_split(Zones, ZonesToSplit, nSplits):
	
	nz		  =  ZonesToSplit
	ns		  =  len(ZonesToSplit)
	nZones	  =  len(Zones)
	nZonesOld =  nZones
	if((Zones[ZonesToSplit[0]].Neighbour.west<0) or (Zones[ZonesToSplit[0]].Neighbour.west ==  ZonesToSplit[-1])):
		westToEast =  True
	else:
		westToEast =  False

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
		iXnew	= np.linspace(0,Nx,nSplits+1).astype(np.int)
		for k2 in range(nSplits):
			Nxnew[k2] 		 = iXnew[k2+1]-iXnew[k2]
			Nxnew[nSplits-1] = Nx-np.sum(Nxnew[:nSplits-1])
			if(k2== 0):
				Zones[nz[k]].zb		= zb[:Nxnew[k2]+1, :]
				Zones[nz[k]].xb		= xb[:Nxnew[k2]+1, :]
				Zones[nz[k]].z		= z[:Nxnew[k2], :]
				Zones[nz[k]].x		= x[:Nxnew[k2], :]
				Zones[nz[k]].gridRc	= R[:Nxnew[k2], :]
				Zones[nz[k]].gridZc	= Z[:Nxnew[k2], :]
				Zones[nz[k]].gridR	= Rr[:Nxnew[k2]+1, :]
				Zones[nz[k]].gridZ	= Zr[:Nxnew[k2]+1, :]
				Zones[nz[k]].Br		= Br[:Nxnew[k2], :]
				Zones[nz[k]].Bz		= Bz[:Nxnew[k2], :]
				Zones[nz[k]].Bphi	= Bphi[:Nxnew[k2], :]
				Zones[nz[k]].Chi	= Chi[:Nxnew[k2], :]
				
				Zones[nz[k]].WestP.Br	= WestP.Br[:Nxnew[k2]]
				Zones[nz[k]].WestP.Bz	= WestP.Bz[:Nxnew[k2]]
				Zones[nz[k]].WestP.Bphi	= WestP.Bphi[:Nxnew[k2]]
				Zones[nz[k]].WestP.R	= WestP.R[:Nxnew[k2]]
				Zones[nz[k]].WestP.Z	= WestP.Z[:Nxnew[k2]]
				
				Zones[nz[k]].EastP.Br	= EastP.Br[:Nxnew[k2]]
				Zones[nz[k]].EastP.Bz	= EastP.Bz[:Nxnew[k2]]
				Zones[nz[k]].EastP.Bphi	= EastP.Bphi[:Nxnew[k2]]
				Zones[nz[k]].EastP.R	= EastP.R[:Nxnew[k2]]
				Zones[nz[k]].EastP.Z	= EastP.Z[:Nxnew[k2]]
				
				Zones[nz[k]].SouthP		= copyP(SouthP)
				
				if(k2== nSplits-1):
					Zones[nz[k]].NorthP =  copyP(NorthP)
					Zones[nz[k]].Neighbour.north= North
				else:
					Zones[nz[k]].NorthP			= types.SimpleNamespace()
					Zones[nz[k]].NorthP.Br		= Br[Nxnew[k2], :]
					Zones[nz[k]].NorthP.Bz		= Bz[Nxnew[k2], :]
					Zones[nz[k]].NorthP.Bphi	= Bphi[Nxnew[k2], :]
					Zones[nz[k]].NorthP.R		= R[Nxnew[k2], :]
					Zones[nz[k]].NorthP.Z		= Z[Nxnew[k2], :]
					Zones[nz[k]].Neighbour.north = nZones

				Zones[nz[k]].magz	= np.copy(magz)
				Zones[nz[k]].Nx		= Nxnew[k2]
				Zones[nz[k]].Nz		= Nz
				curb = Nxnew[k2]+1
				cur  = Nxnew[k2]
			else:
				nZones= nZones+1
				Zones.append(types.SimpleNamespace())
				Zones[-1].zb		= zb[curb-1:curb+Nxnew[k2],:]
				Zones[-1].xb		= xb[curb-1:curb+Nxnew[k2],:]
				Zones[-1].z			= z[cur:cur+Nxnew[k2],:]
				Zones[-1].x			= x[cur:cur+Nxnew[k2],:]
				Zones[-1].gridRc	= R[cur:cur+Nxnew[k2],:]
				Zones[-1].gridZc	= Z[cur:cur+Nxnew[k2],:]
				Zones[-1].gridR		= Rr[cur:cur+Nxnew[k2]+1,:]
				Zones[-1].gridZ		= Zr[cur:cur+Nxnew[k2]+1,:]
				Zones[-1].Br		= Br[cur:cur+Nxnew[k2],:]
				Zones[-1].Bz		= Bz[cur:cur+Nxnew[k2],:]
				Zones[-1].Bphi		= Bphi[cur:cur+Nxnew[k2],:]
				Zones[-1].Chi		= Chi[cur:cur+Nxnew[k2],:]
				
				Zones[-1].WestP			= types.SimpleNamespace()
				Zones[-1].WestP.Br		= WestP.Br[cur:cur+Nxnew[k2]]
				Zones[-1].WestP.Bz		= WestP.Bz[cur:cur+Nxnew[k2]]
				Zones[-1].WestP.Bphi	= WestP.Bphi[cur:cur+Nxnew[k2]]
				Zones[-1].WestP.R		= WestP.R[cur:cur+Nxnew[k2]]
				Zones[-1].WestP.Z		= WestP.Z[cur:cur+Nxnew[k2]]
				
				Zones[-1].EastP			= types.SimpleNamespace()
				Zones[-1].EastP.Br		= EastP.Br[cur:cur+Nxnew[k2]]
				Zones[-1].EastP.Bz		= EastP.Bz[cur:cur+Nxnew[k2]]
				Zones[-1].EastP.Bphi	= EastP.Bphi[cur:cur+Nxnew[k2]]
				Zones[-1].EastP.R		= EastP.R[cur:cur+Nxnew[k2]]
				Zones[-1].EastP.Z		= EastP.Z[cur:cur+Nxnew[k2]]
				
				Zones[-1].SouthP 		= types.SimpleNamespace()
				Zones[-1].SouthP.Br		= Br[cur-1,:]
				Zones[-1].SouthP.Bz		= Bz[cur-1,:]
				Zones[-1].SouthP.Bphi	= Bphi[cur-1,:]
				Zones[-1].SouthP.R		= R[cur-1,:]
				Zones[-1].SouthP.Z		= Z[cur-1,:]
				
				Zones[-1].Neighbour		= types.SimpleNamespace()
				if(k2 ==  nSplits-1):
					Zones[-1].NorthP =  copyP(NorthP)
					Zones[-1].Neighbour.north= North
				else:
					Zones[-1].NorthP 		= types.SimpleNamespace()
					Zones[-1].NorthP.Br		= Br[cur+Nxnew[k2],:]
					Zones[-1].NorthP.Bz		= Bz[cur+Nxnew[k2],:]
					Zones[-1].NorthP.Bphi	= Bphi[cur+Nxnew[k2],:]
					Zones[-1].NorthP.R		= R[cur+Nxnew[k2],:]
					Zones[-1].NorthP.Z		= Z[cur+Nxnew[k2],:]
					Zones[-1].Neighbour.north= nZones

				if(k2 == 1):
					Zones[-1].Neighbour.south= ZonesToSplit[k]
				else:
					Zones[-1].Neighbour.south= nZones-2

				if(westToEast):
					if(k == 0):
						if(Zones[ZonesToSplit[0]].Neighbour.west < -1):
							Zones[-1].Neighbour.west= Zones[ZonesToSplit[0]].Neighbour.west
						else:														#periodic
							Zones[-1].Neighbour.west= nZones + (ns-1)*(nSplits-1) - 1

					else:
						Zones[-1].Neighbour.west= nZones - nSplits

					if(k == ns-1):
						if(Zones[ZonesToSplit[-1]].Neighbour.east < -1):
							Zones[-1].Neighbour.east= Zones[ZonesToSplit[-1]].Neighbour.east
						else:
							Zones[-1].Neighbour.east= nZones - (ns-1)*(nSplits-1) - 1

					else:
						Zones[-1].Neighbour.east= nZones + (nSplits-2)

				else:
					if(k == 0):
						if(Zones[ZonesToSplit[0]].Neighbour.east < -1):
							Zones[-1].Neighbour.east = Zones[ZonesToSplit[0]].Neighbour.east
						else:
							Zones[-1].Neighbour.east = nZones + (ns-1)*(nSplits-1) - 1

					else:
						Zones[-1].Neighbour.east = nZones - nSplits

					if(k == ns-1):
						if(Zones[ZonesToSplit[-1]].Neighbour.west < -1):
							Zones[-1].Neighbour.west = Zones[ZonesToSplit[-1]].Neighbour.west
						else:
							Zones[-1].Neighbour.west = nZones - (ns-1)*(nSplits-1) - 1

					else:
						Zones[-1].Neighbour.west = nZones + (nSplits-2)

				Zones[-1].magz	= np.array([magz[0], magz[1]+cur, magz[2]], dtype='i4')
				Zones[-1].Nx	= Nxnew[k2]
				Zones[-1].Nz	= Nz
				curb	= curb+Nxnew[k2]
				cur		= cur+Nxnew[k2]

		for ko in range(nZonesOld):
			if(Zones[ko].Neighbour.south  == ZonesToSplit[k]):
				Zones[ko].Neighbour.south = nZones - 1
