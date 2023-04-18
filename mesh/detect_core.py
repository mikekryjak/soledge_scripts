import types
from routines.globals	import CORE_NEIGHBOUR, EW_BORDER_NEIGHBOUR, NS_BORDER_NEIGHBOUR

def detect_core(Config):

	Zones	  = Config.Zones
	Megazones = Config.Megazones
	
	nZones = len(Zones)
	for k in range(nZones):
		Zones[k].MagNeighbour		= types.SimpleNamespace()
		Zones[k].MagNeighbour.north = 0
		Zones[k].MagNeighbour.south = 0
		Zones[k].MagNeighbour.east	= 0
		Zones[k].MagNeighbour.west	= 0

	for k in range(nZones):
		if(Zones[k].Neighbour.north < 0):
			Zones[k].Neighbour.north	= NS_BORDER_NEIGHBOUR
			Zones[k].MagNeighbour.north =  1

		if(Zones[k].Neighbour.south < 0):
			Zones[k].Neighbour.south	= NS_BORDER_NEIGHBOUR
			Zones[k].MagNeighbour.south =  1

		if(Zones[k].Neighbour.east < 0):
			Zones[k].Neighbour.east		= EW_BORDER_NEIGHBOUR
			Zones[k].MagNeighbour.east	=  1

		if(Zones[k].Neighbour.west < 0):
			Zones[k].Neighbour.west		= EW_BORDER_NEIGHBOUR
			Zones[k].MagNeighbour.west	=  1

	for k in range(len(Megazones)):
		nze = Megazones[k].list[-1]
		nz1 = Megazones[k].list[0]

		if(Zones[nze].Neighbour.east == nz1):						#periodic
			Megazones[k].isperiodic = True
			for k1 in range(len(Megazones[k].list)):
				nz = Megazones[k].list[k1]
				if(Zones[nz].Neighbour.south < 0):
					Zones[nz].Neighbour.south	 = CORE_NEIGHBOUR
					Zones[nz].MagNeighbour.south =  1
		else:
			Megazones[k].isperiodic = False



