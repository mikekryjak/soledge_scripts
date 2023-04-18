import types
import numpy 			as np
from mesh.P_routines	import copyP

def copy_zones(Zones):
	NewZones = []
	for k in range(len(Zones)):
		NewZones.append(types.SimpleNamespace())
		
	return NewZones

def copy_zone(Zones):

	NewZone = types.SimpleNamespace()

	NewZone.Nx		=  Zone.Nx
	NewZone.Nz		=  Zone.Nz

	NewZone.gridRc	=  np.copy(Zone.gridRc)
	NewZone.gridZc	=  np.copy(Zone.gridZc)
	NewZone.gridR	=  np.copy(Zone.gridR)
	NewZone.gridZ	=  np.copy(Zone.gridZ)
	NewZone.Br		=  np.copy(Zone.Br)
	NewZone.Bz		=  np.copy(Zone.Bz)
	NewZone.Bphi	=  np.copy(Zone.Bphi)
		
	NewZone.zb		=  np.copy(Zone.zb)
	NewZone.xb		=  np.copy(Zone.xb)
	NewZone.z		=  np.copy(Zone.z)
	NewZone.x		=  np.copy(Zone.x)
	NewZone.Chi		=  np.copy(Zone.Chi)
		
	NewZone.SouthP	=  copyP(Zone.SouthP)
	NewZone.NorthP	=  copyP(Zone.NorthP)
	NewZone.EastP	=  copyP(Zone.EastP)
	NewZone.WestP	=  copyP(Zone.WestP)

	NewZone.Neighbour		= types.SimpleNamespace()
	NewZone.Neighbour.West	=  Zone.Neighbour.west
	NewZone.Neighbour.East	=  Zone.Neighbour.east
	NewZone.Neighbour.North	=  Zone.Neighbour.north
	NewZone.Neighbour.South	=  Zone.Neighbour.south

	NewZone.mz				= Zones.mz
	NewZone.pmz				= Zones.pmz
	NewZone.magz			= np.copy(Zone.magz)

	return NewZone


def copy_magzones(MagZones):
	NewMagZones = []
	for k in range(len(MagZones)):
		NewMagZones.append(copy_magzone(MagZones[k]))

	return NewMagZones

def copy_magzone(MagZone):
	NewMagZone = types.SimpleNamespace()

	NewMagZone.gridR	=  np.copy(MagZone.gridR)
	NewMagZone.gridZ	=  np.copy(MagZone.gridZ)

	NewMagZone.Neighbour		= types.SimpleNamespace()
	NewMagZone.Neighbour.north	= MagZone.Neighbour.north
	NewMagZone.Neighbour.south	= MagZone.Neighbour.south
	NewMagZone.Neighbour.east	= MagZone.Neighbour.east
	NewMagZone.Neighbour.west	= MagZone.Neighbour.west
	NewMagZone.mz				= MagZone.mz
	NewMagZone.pmz				= MagZone.pmz
	if(hasattr(MagZone,'list')): NewMagZone.list = np.copy(MagZone.list)

	return NewMagZone
		