from tkinter 					import messagebox
from tkinter.messagebox 		import Message

import types
from tkinter					import messagebox
from tkinter.messagebox 		import Message
import numpy					as np

from mesh.delete_magzone		import delete_magzone
from mesh.get_magzone_path		import get_magzone_path
from routines.intersect_contour	import intersect_contour
from routines.part_contour		import part_contour
from routines.globals			import DEBUG, EXTERNAL_PLASMA_WALL

# find_megazones_outside_plasma Detect megazones and Pmegazones that are fully outside the plasma

def find_megazones_outside_plasma(self, Config): 

	if(DEBUG > 0): print("\nfind_megazones_outside_plasma\n")
	
	MagZones		= Config.MagZones
	MagMegazones	= Config.MagMegazones
	MagPMegazones	= Config.MagPMegazones

	MagZoneOut = np.ones(len(MagZones), dtype='i4')

#	Loop on walls

	for iPWall in range(len(Config.iPwalls)):
		iWall = Config.iPwalls[iPWall]
		Wall  = Config.Walls[iWall]
		if(Wall.Type == EXTERNAL_PLASMA_WALL):					# plasma inside this wall
			Rwall = Wall.Rwall
			Zwall = Wall.Zwall
			for iMagZone in range(len(MagZones)):				# Loop on zones
				MagZonePath = get_magzone_path(MagZones[iMagZone])

#				A zone is in the plasma if its edge is partly inside the wall or the wall is entirely inside the zone

				zoneedges_inWall = np.sum(np.where(Wall.WallPath.contains_points(MagZonePath.vertices),   1,0))
				wall_inZoneEdges = np.sum(np.where(MagZonePath.contains_points(np.array([Rwall,Zwall]).T),1,0))

				if((zoneedges_inWall > 0) or (wall_inZoneEdges > 0)): MagZoneOut[iMagZone] = 0

#	Now check for megazones
	MagMegazoneOutList = []
	for iMz in range(len(MagMegazones)):
		if(np.sum(MagZoneOut[MagMegazones[iMz].list]) == len(MagMegazones[iMz].list)): MagMegazoneOutList.append(iMz)

#	Now check for Pmegazones
	MagPMegazoneOutList = []
	for iPmz in range(len(MagPMegazones)):
		if(np.sum(MagZoneOut[MagPMegazones[iPmz].list]) == len(MagPMegazones[iPmz].list)): MagPMegazoneOutList.append(iPmz)

#	Finally produce the list of concerned zones
	MagZoneOutList = np.empty(0, dtype='i4')
	if(len(MagMegazoneOutList) + len(MagPMegazoneOutList) > 0):
		for iMz in MagMegazoneOutList:
			MagZoneOutList = np.append(MagZoneOutList, MagMegazones[iMz].list)

		for iPmz in MagPMegazoneOutList:
			MagZoneOutList = np.append(MagZoneOutList, MagPMegazones[iPmz].list)

		MagZoneOutList = np.unique(MagZoneOutList)			#unique and sorted

#		Ask for deleteing
	
		outMessage = "Zones ("
		for k in MagZoneOutList: outMessage +="{:d},".format(k)
		outMessage = outMessage[:-1]+") are fully outside the plasma domain.\nDo you want automatically delete them?"

		ret_status = messagebox.askyesno("Magnetic zones outside plasma", outMessage, default=messagebox.NO)
		if(not ret_status): return True

#		Delete zones outside plasma

		for k in range(len(MagZoneOutList)):
			kZone = MagZoneOutList[k]
			delete_magzone(Config, kZone)
			if(k < len(MagZoneOutList) - 1): MagZoneOutList[k+1:] -= 1			# Renumber the other zones to delete
			
	if(DEBUG > 0): print("\nfind_megazones_outside_plasma: completed\n")

	return True
