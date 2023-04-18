import numpy as np

def mag_to_zones_pos(Config, MagCells):

	MagZones = Config.MagZones
	Zones	 = Config.Zones

	if(len(MagCells.shape) == 2):
		Cells = np.zeros_like(MagCells)
		iMagZones = np.unique(MagCells[:,0])
		for iMagZone in iMagZones:										#loop on magnetic zones
			kk = np.where(MagCells[:,0] == iMagZone)[0]
			ii = MagCells[kk,1]
			jj = MagCells[kk,2]
			for iZone in MagZones[iMagZone].list:						#loop on splitted zones					
				ll = np.where((ii >= Zones[iZone].magz[1]) & (ii < Zones[iZone].magz[1] + Zones[iZone].Nx) &
							  (jj >= Zones[iZone].magz[2]) & (jj < Zones[iZone].magz[2] + Zones[iZone].Nz))[0]
				if(len(ll) > 0):
					Cells[kk[ll],0] = iZone								#Splitted zone index
					Cells[kk[ll],1] = ii[ll] - Zones[iZone].magz[1]		#splitted coordinate	
					Cells[kk[ll],2] = jj[ll] - Zones[iZone].magz[2]		#splitted coordinate	
	else:
		kMag = MagCells[0]
		iMag = MagCells[1]
		jMag = MagCells[2]
		for iZone in MagZones[kMag].list:					#loop on splitted zones
			if( (iMag >= Zones[iZone].magz[1]) and (iMag < Zones[iZone].magz[1] + Zones[iZone].Nx) and
				(jMag >= Zones[iZone].magz[2]) and (jMag < Zones[iZone].magz[2] + Zones[iZone].Nz)):
				Cells = np.array([iZone, iMag - Zones[iZone].magz[1], jMag - Zones[iZone].magz[2]])
				break

	return Cells
