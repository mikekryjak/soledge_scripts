#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from routines.cli_routines import cli_present, cli_get_value, cli_get_value
	from old_mesh_to_new import old_mesh_to_new
	 
	if(cli_present("-h",sys.argv)):
		print("\nThis script converts old versions of mesh.h5 file to the last one\n")
		print("old_mesh_to_new options")
		print("\t-old_mesh <value>  Old mesh file name [d='mesh.h5']")
		print("\t-new_mesh <value>  New mesh file name [d='new_mesh.h5']")
		print()
		exit()

	old_mesh	 	= cli_get_value("-old_mesh",	sys.argv,"mesh.h5")
	new_mesh	 	= cli_get_value("-new_mesh",	sys.argv,"new_mesh.h5")
	old_mesh_to_new(old_mesh=old_mesh, new_mesh=new_mesh)
	exit()

#=======================================

# Function definition is here

import types
import numpy as np
import h5py
from routines.h5_routines				import h5_read
from files.load_soledge_mesh_file_V2	import load_soledge_mesh_file_V2
from files.load_soledge_mesh_file_V0	import load_soledge_mesh_file_V0
from files.save_soledge_mesh_file		import save_soledge_mesh_file
from routines.globals					import DEBUG

def old_mesh_to_new(old_mesh="mesh.h", new_mesh="new_mesh.h5"):

#	Load old mesh and old plasma

	if(DEBUG > 0): print("\told_mesh_to_new: Converting mesh file ",old_mesh)

	try:	
		if_mesh = h5py.File(old_mesh, "r")
	except:
		if(DEBUG > 0): print("\told_mesh_to_new: Unable to open ",old_mesh)
		return

	try:
		FileVersion	= h5_read(if_mesh, 'Version',	keep_array=False)
	except:
		FileVersion = 0
	if_mesh.close()

	if(FileVersion < 1):
		Config				= load_soledge_mesh_file_V0(old_mesh)

#		Update to new data

		Config.MagZones			= Config.Zones
		Config.MagMegazones		= Config.Megazones
		Config.MagPMegazones	= Config.PMegazones

		if(Config.Zones_OK):
			Config.MagZones_OK = Config.Zones_OK
			for k in range(len(Config.Zones)):
				Config.MagZones[k].list = np.array([k])
				Config.Zones[k].magz	= np.array([k,0,0])

		if(Config.transp_prof_OK == True):
			for k in range(len(Config.TranspCuts)):
				Config.TranspCuts[k].Flux.MagZones = Config.TranspCuts[k].Flux.Zones

	elif(FileVersion < 3):
		Config				= load_soledge_mesh_file_V2(old_mesh)

	save_soledge_mesh_file(new_mesh, Config)
	
	return
	
