import numpy				as np
from routines.globals		import CORE_NEIGHBOUR, DEBUG

def get_core_sep_megazones(Config):

	if(DEBUG > 1): print("get_core_sep_megazones")

	Zones	  = Config.Zones
	Megazones = Config.Megazones

#	Search Core zone

	k = 0
	while (k < len(Zones)):
		South = Zones[k].Neighbour.south
		if(South == CORE_NEIGHBOUR):
			CoreZone = k
			break
		k += 1

	CoreMegazone = Zones[CoreZone].mz

#	Search separatrix Megazone

	SepMegazone  = CoreMegazone

	north = Zones[CoreZone].Neighbour.north
	while ((north >= 0) and Megazones[Zones[north].mz].isperiodic):
		SepMegazone  = Zones[north].mz
		north = Zones[north].Neighbour.north

	if(DEBUG > 1): print("get_core_sep_megazones: completed")

	return CoreMegazone, SepMegazone