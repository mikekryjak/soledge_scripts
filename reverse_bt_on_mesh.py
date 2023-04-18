#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from routines.cli_routines import cli_present, cli_get_value, cli_get_value
	from reverse_bt_on_mesh import reverse_bt_on_mesh
	 
	if(cli_present("-h",sys.argv)):
		print("\nThis script reverse Bt sign on mesh file\n")
		print("reverse_bt_on_mesh options")
		print("\t-old_mesh <value>  Old mesh file name [d='mesh.h5']")
		print("\t-new_mesh <value>  New mesh file name [d='mesh.h5']")
		print()
		exit()

	old_mesh	 	= cli_get_value("-old_mesh",	sys.argv,"mesh.h5")
	new_mesh	 	= cli_get_value("-new_mesh",	sys.argv,"mesh.h5")
	reverse_bt_on_mesh(old_mesh=old_mesh, new_mesh=new_mesh)
	exit()

#=======================================

# Function definition is here

import types
import numpy as np
import h5py
from routines.h5_routines				import h5_read
from files.load_soledge_mesh_file		import load_soledge_mesh_file
from files.save_soledge_mesh_file		import save_soledge_mesh_file
from routines.globals					import DEBUG

def reverse_bt_on_mesh(old_mesh="mesh.h5", new_mesh="mesh.h5"):

#	Load old mesh and old plasma

	if(DEBUG > 0): print("\treverse_bt_on_mesh: Converting mesh file ",old_mesh)

	Config	= load_soledge_mesh_file(old_mesh)

	Config.Bphi2	= -Config.Bphi2
	Config.Bphi2D	= -Config.Bphi2D

	Zones  = Config.Zones
	for k in range (len(Zones)):
		Zones[k].Bphi		 = -Zones[k].Bphi
		Zones[k].SouthP.Bphi = -Zones[k].SouthP.Bphi
		Zones[k].NorthP.Bphi = -Zones[k].NorthP.Bphi
		Zones[k].WestP.Bphi	 = -Zones[k].WestP.Bphi
		Zones[k].EastP.Bphi	 = -Zones[k].EastP.Bphi

	save_soledge_mesh_file(new_mesh, Config)
	
	return
	
