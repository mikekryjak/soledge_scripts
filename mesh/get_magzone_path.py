import numpy			 as np
from matplotlib.path	 import Path
from math				 import sqrt

# Define the border path for a magzone 
# ####################################################################

def get_magzone_path(MagZone):

	MagZoneEdges_R = np.copy(MagZone.north.R)
	MagZoneEdges_Z = np.copy(MagZone.north.Z)
	if (sqrt((MagZoneEdges_R[-1]-MagZone.west.R[0])**2+(MagZoneEdges_Z[-1]-MagZone.west.Z[0])**2) < 1.e-10):
		MagZoneEdges_R = np.append(MagZoneEdges_R, MagZone.west.R)
		MagZoneEdges_Z = np.append(MagZoneEdges_Z, MagZone.west.Z)
	elif (sqrt((MagZoneEdges_R[-1]-MagZone.west.R[-1])**2+(MagZoneEdges_Z[-1]-MagZone.west.Z[-1])**2) < 1.e-10):
		MagZoneEdges_R = np.append(MagZoneEdges_R, MagZone.west.R[::-1])
		MagZoneEdges_Z = np.append(MagZoneEdges_Z, MagZone.west.Z[::-1])
	elif (sqrt((MagZoneEdges_R[0]-MagZone.west.R[0])**2+(MagZoneEdges_Z[0]-MagZone.west.Z[0])**2) < 1.e-10):
		MagZoneEdges_R = np.append(MagZoneEdges_R[::-1], MagZone.west.R[::-1])
		MagZoneEdges_Z = np.append(MagZoneEdges_Z[::-1], MagZone.west.Z[::-1])
	else:
		MagZoneEdges_R = np.append(MagZoneEdges_R[::-1], MagZone.west.R[::-1])
		MagZoneEdges_Z = np.append(MagZoneEdges_Z[::-1], MagZone.west.Z[::-1])

	if (sqrt((MagZoneEdges_R[-1]-MagZone.south.R[0])**2+(MagZoneEdges_Z[-1]-MagZone.south.Z[0])**2) < 1.e-10):
		MagZoneEdges_R = np.append(MagZoneEdges_R, MagZone.south.R)
		MagZoneEdges_Z = np.append(MagZoneEdges_Z, MagZone.south.Z)
	else:
		MagZoneEdges_R = np.append(MagZoneEdges_R, MagZone.south.R[::-1])
		MagZoneEdges_Z = np.append(MagZoneEdges_Z, MagZone.south.Z[::-1])

	if (sqrt((MagZoneEdges_R[-1]-MagZone.east.R[0])**2+(MagZoneEdges_Z[-1]-MagZone.east.Z[0])**2) < 1.e-10):
		MagZoneEdges_R = np.append(MagZoneEdges_R, MagZone.east.R)
		MagZoneEdges_Z = np.append(MagZoneEdges_Z, MagZone.east.Z)
	else:
		MagZoneEdges_R = np.append(MagZoneEdges_R, MagZone.east.R[::-1])
		MagZoneEdges_Z = np.append(MagZoneEdges_Z, MagZone.east.Z[::-1])

#	print("get_magzone_path: MagZoneEdges_R=",MagZoneEdges_R)
#	print("get_magzone_path: MagZoneEdges_Z=",MagZoneEdges_Z)
	MagZonePath = Path(np.array([MagZoneEdges_R,MagZoneEdges_Z]).T, closed=True)
	
	return MagZonePath


# Define the border path for a zone 
# ####################################################################

def get_zone_path(Zone):

	ZoneEdges_R = np.concatenate((Zone.gridR[:,0], Zone.gridR[-1,1:], Zone.gridR[-2::-1,-1], Zone.gridR[0,-2::-1]))
	ZoneEdges_Z = np.concatenate((Zone.gridZ[:,0], Zone.gridZ[-1,1:], Zone.gridZ[-2::-1,-1], Zone.gridZ[0,-2::-1]))

	ZonePath = Path(np.array([ZoneEdges_R,ZoneEdges_Z]).T, closed=True)
	
	return ZonePath

