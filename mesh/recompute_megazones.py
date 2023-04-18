import types
import numpy	as np

def recompute_megazones(Config):
	
	Zones		= Config.Zones
	MagZones	= Config.MagZones
	Megazones	= []
	PMegazones	= []

#	Find Megazones

	for k in range(len(MagZones)): MagZones[k].list = np.empty((0),dtype='i4')

	for k in range(len(Zones)):
		MagZones[Zones[k].magz[0]].list = np.append(MagZones[Zones[k].magz[0]].list, k)
		Zones[k].mz  = -1
		Zones[k].pmz = -1
	
	allChecked  = False
	while(not allChecked):
		
		mz = -1
		for k in range(len(Zones)):
			if(Zones[k].mz < 0): 								#new Megazone
				Megazones.append(types.SimpleNamespace())
				mz						   += 1 
				start						= k
				Zones[start].mz				= mz
				
				Megazones[mz].list			= np.array([start])
				Megazones[mz].isperiodic	= False
				
				east						= Zones[start].Neighbour.east
				while((east > -1) and (east != start)):
					Megazones[mz].list		= np.append(Megazones[mz].list,east)
					Zones[east].mz			= mz
					east					= Zones[east].Neighbour.east

				if(east == start): Megazones[-1].isperiodic = True
				else:
					west = Zones[start].Neighbour.west
					while(west > -1):
						Megazones[mz].list	= np.append(west,Megazones[mz].list)
						Zones[west].mz		= mz
						west				= Zones[west].Neighbour.west
		
		allChecked  = True
		for k in range(len(Zones)):
			if(Zones[k].mz < 0):
				allChecked = False
				break


#	Find PMegazones

	pmz = -1
	for k in range(len(Zones)):
		if(Zones[k].Neighbour.south < -1): 								#new PMegazone
			PMegazones.append(types.SimpleNamespace())
			pmz						+= 1 
			PMegazones[pmz].list	= np.array([k])
			Zones[k].pmz			= pmz
			num						= k
			while(Zones[num].Neighbour.north > -1):
				num 			= Zones[num].Neighbour.north
				Zones[num].pmz	= pmz
				PMegazones[pmz].list		= np.append(PMegazones[pmz].list,num)
		
	
	for k in range(len(Zones)):
		Zones[k].MagNeighbour	= types.SimpleNamespace()

		if(Zones[k].Neighbour.north < -1):	Zones[k].MagNeighbour.north = 1
		else:								Zones[k].MagNeighbour.north = 0

		if(Zones[k].Neighbour.south < -1):	Zones[k].MagNeighbour.south = 1
		else:								Zones[k].MagNeighbour.south = 0

		if(Zones[k].Neighbour.west < -1):	Zones[k].MagNeighbour.west = 1
		else:								Zones[k].MagNeighbour.west = 0

		if(Zones[k].Neighbour.east < -1):	Zones[k].MagNeighbour.east = 1
		else:								Zones[k].MagNeighbour.east = 0
		
	Config.Megazones  = Megazones
	Config.PMegazones = PMegazones


#	Set MagZones data
#=================================================================
#	Build Br, Bz and Bphi on MagZone from values on splitted zones

	for k in range(len(MagZones)):

		Br		= np.empty((MagZones[k].gridZ.shape[0]-1, MagZones[k].gridZ.shape[1]-1), dtype = 'f8')
		Bz		= np.empty_like(Br)
		Bphi	= np.empty_like(Br)
		x		= np.empty_like(Br)
		z		= np.empty_like(Br)

		for sz in range(len(MagZones[k].list)):
			sz_k	= MagZones[k].list[sz]
			sz_si	= Zones[sz_k].magz[1]
			sz_sj	= Zones[sz_k].magz[2]
			sz_ei	= sz_si + Zones[sz_k].Br.shape[0]
			sz_ej	= sz_sj + Zones[sz_k].Br.shape[1]

			Br[sz_si:sz_ei, sz_sj:sz_ej] = Zones[sz_k].Br
			Bz[sz_si:sz_ei, sz_sj:sz_ej] = Zones[sz_k].Bz
			Bphi[sz_si:sz_ei, sz_sj:sz_ej] = Zones[sz_k].Bphi

			x[sz_si:sz_ei, sz_sj:sz_ej] = Zones[sz_k].x
			z[sz_si:sz_ei, sz_sj:sz_ej] = Zones[sz_k].z

		MagZones[k].gridRc = 0.25*(MagZones[k].gridR[:-1,:-1] + MagZones[k].gridR[:-1,1:] + MagZones[k].gridR[1:,1:] + MagZones[k].gridR[1:,:-1])
		MagZones[k].gridZc = 0.25*(MagZones[k].gridZ[:-1,:-1] + MagZones[k].gridZ[:-1,1:] + MagZones[k].gridZ[1:,1:] + MagZones[k].gridZ[1:,:-1])

		MagZones[k].Br	 = Br
		MagZones[k].Bz	 = Bz
		MagZones[k].Bphi = Bphi
		MagZones[k].x 	  = x
		MagZones[k].z 	  = z

		MagZones[k].Nx		= MagZones[k].gridRc.shape[0]
		MagZones[k].Nz		= MagZones[k].gridRc.shape[1]

