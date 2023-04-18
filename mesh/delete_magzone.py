import types
import numpy			as np

from mesh.copy_zones	import copy_magzone
from routines.globals	import DEBUG, CORE_NEIGHBOUR, EW_BORDER_NEIGHBOUR, NS_BORDER_NEIGHBOUR

# Delete a given zone and update all the related structures to make them consistent with the new zone list and numbering
# Replace deleteZone.m but it is very simplified in updating MagMegazones e  MagPMegazones 
# ####################################################################

def delete_magzone(Config, kZoneDel):

	if(DEBUG > 0): print("delete_magzone")

	MagZones		= Config.MagZones
	MagMegazones	= Config.MagMegazones
	MagPMegazones	= Config.MagPMegazones

#	first update to MagMegazones e  MagPMegazones references (also delete MagMegazone e  MagPMegazone if it becomes empy)

	iMz  = MagZones[kZoneDel].mz
	iPMz = MagZones[kZoneDel].pmz
	if(DEBUG > 0): print("\tZoneDel={:d}, Mz={:d}, PMz={:d}".format(kZoneDel+1, iMz+1,iPMz+1))

	if(len(MagMegazones[iMz].list) > 1):
		MagMegazones[iMz].list = np.setdiff1d(MagMegazones[iMz].list, kZoneDel, assume_unique=True)		#assume_unique=True to avoid output sorting
		iMzDel = -1
	else:
		del MagMegazones[iMz]
		iMzDel = iMz

	if(len(MagPMegazones[iPMz].list) > 1):
		MagPMegazones[iPMz].list = np.setdiff1d(MagPMegazones[iPMz].list, kZoneDel, assume_unique=True)	#assume_unique=True to avoid output sorting
		MagPMegazones[iPMz].isperiodic = False
		iPMzDel = -1
	else:
		del MagPMegazones[iPMz]
		iPMzDel = iPMz

#	Delete MagZones we want to suppress

	del MagZones[kZoneDel]

#	Change neighbours array and if needed MagMegazones e  MagPMegazones references  in MagZones

	if((DEBUG > 0) and ((iMzDel > -1) or (iPMzDel > -1))): print("\tMzDel={:d}, PMzDel={:d}".format(iMzDel+1, iPMzDel+1))
	for k in range(len(MagZones)):
		if  (MagZones[k].Neighbour.north == kZoneDel): MagZones[k].Neighbour.north  = -1
		elif(MagZones[k].Neighbour.north  > kZoneDel): MagZones[k].Neighbour.north -= 1

		if  (MagZones[k].Neighbour.south == kZoneDel): MagZones[k].Neighbour.south  = -1
		elif(MagZones[k].Neighbour.south  > kZoneDel): MagZones[k].Neighbour.south -= 1

		if  (MagZones[k].Neighbour.east  == kZoneDel): MagZones[k].Neighbour.east  = -1
		elif(MagZones[k].Neighbour.east  >  kZoneDel): MagZones[k].Neighbour.east -= 1

		if  (MagZones[k].Neighbour.west  == kZoneDel): MagZones[k].Neighbour.west  = -1
		elif(MagZones[k].Neighbour.west   > kZoneDel): MagZones[k].Neighbour.west -= 1

		if((iMzDel > -1)  and (MagZones[k].mz  > iMzDel)): MagZones[k].mz -= 1
		if((iPMzDel > -1) and (MagZones[k].pmz > iPMzDel)): MagZones[k].pmz -= 1

#	update to MagMegazones e  MagPMegazones references to MagZones

	for MagMegazone in MagMegazones:
		iFix = np.where(MagMegazone.list > kZoneDel)[0]
		if(len(iFix) > 0): MagMegazone.list[iFix] = MagMegazone.list[iFix] - 1

	for MagPMegazone in MagPMegazones:
		iFix = np.where(MagPMegazone.list > kZoneDel)[0]
		if(len(iFix) > 0): MagPMegazone.list[iFix] = MagPMegazone.list[iFix] - 1

	if(DEBUG > 0): print("delete_magzone: completed\n")

	
	return
