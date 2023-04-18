#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from routines.cli_routines import cli_present, cli_get_value, cli_get_value
	from b_from_g_to_t_on_mesh import b_from_g_to_t_on_mesh
	 
	if(cli_present("-h",sys.argv)):
		print("\nThis script transform B from Gauss to Telsla on mesh file\n")
		print("b_from_g_to_t_on_mesh options")
		print("\t-old_mesh <value>  Old mesh file name [d='mesh.h5']")
		print("\t-new_mesh <value>  New mesh file name [d='mesh.h5']")
		print()
		exit()

	old_mesh	 	= cli_get_value("-old_mesh",	sys.argv,"mesh.h5")
	new_mesh	 	= cli_get_value("-new_mesh",	sys.argv,"mesh.h5")
	b_from_g_to_t_on_mesh(old_mesh=old_mesh, new_mesh=new_mesh)
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

def b_from_g_to_t_on_mesh(old_mesh="mesh.h5", new_mesh="mesh.h5"):

#	Load old mesh and old plasma

	if(DEBUG > 0): print("\tb_from_g_to_t_on_mesh: Converting mesh file ",old_mesh)

	Config	= load_soledge_mesh_file(old_mesh)

	Config.Br2		= Config.Br2*1e-4
	Config.Br2D		= Config.Br2D*1e-4

	Config.Bz2		= Config.Bz2*1e-4
	Config.Bz2D		= Config.Bz2D*1e-4

	Config.Bphi2	= Config.Bphi2*1e-4
	Config.Bphi2D	= Config.Bphi2D*1e-4

	Zones  = Config.Zones
	for k in range (len(Zones)):
		Zones[k].Br		 	 = Zones[k].Br*1e-4
		Zones[k].SouthP.Br   = Zones[k].SouthP.Br*1e-4
		Zones[k].NorthP.Br 	 = Zones[k].NorthP.Br*1e-4
		Zones[k].WestP.Br	 = Zones[k].WestP.Br*1e-4
		Zones[k].EastP.Br	 = Zones[k].EastP.Br*1e-4

		Zones[k].Bz		 	 = Zones[k].Bz*1e-4
		Zones[k].SouthP.Bz   = Zones[k].SouthP.Bz*1e-4
		Zones[k].NorthP.Bz 	 = Zones[k].NorthP.Bz*1e-4
		Zones[k].WestP.Bz	 = Zones[k].WestP.Bz*1e-4
		Zones[k].EastP.Bz	 = Zones[k].EastP.Bz*1e-4

		Zones[k].Bphi		 = Zones[k].Bphi*1e-4
		Zones[k].SouthP.Bphi = Zones[k].SouthP.Bphi*1e-4
		Zones[k].NorthP.Bphi = Zones[k].NorthP.Bphi*1e-4
		Zones[k].WestP.Bphi	 = Zones[k].WestP.Bphi*1e-4
		Zones[k].EastP.Bphi	 = Zones[k].EastP.Bphi*1e-4

	save_soledge_mesh_file(new_mesh, Config)
	
	return
	
