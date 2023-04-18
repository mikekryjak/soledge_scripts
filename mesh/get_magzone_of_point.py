import types
import numpy				as np
from mesh.get_magzone_path	import get_magzone_path, get_zone_path


# Find MagZone of point (R, Z)
# ####################################################################

def get_magzone_of_point(Config, Rpoint, Zpoint):

	MagZones = Config.MagZones
	kMagZone = -1

	for k in range(len(MagZones)):
		EdgePath = get_magzone_path(MagZones[k])
		if(EdgePath.contains_point([Rpoint, Zpoint])):
			kMagZone = k
			break

	return kMagZone


# Find Zone of point (R, Z)
# ####################################################################

def get_zone_of_point(Config, Rpoint, Zpoint):

	Zones = Config.Zones
	kZone = -1

	for k in range(len(Zones)):
		EdgePath = get_zone_path(Zones[k])
		if(EdgePath.contains_point([Rpoint, Zpoint])):
			kZone = k
			break

	return kZone
