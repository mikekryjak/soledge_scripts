#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from fix_t_imp	import fix_t_imp
	from routines.cli_routines		import cli_get_value, cli_present
	 
	if(cli_present("-h",sys.argv)):
		print("\nThis script make T_imp=Ti\n")
		print("fix_ti options")
		print("\t-path <value>     Path of mesh [d='.']")
		print()
		exit()

	path = cli_get_value("-path",	sys.argv,  ".")

	fix_t_imp(path=path)
	exit()

#=======================================

# Function definition is here

import os
import types
import numpy as np
from files.load_soledge_mesh_file		import load_soledge_mesh_file
from files.load_plasma_files			import load_plasma_files
from files.save_plasma_files			import save_plasma_files

def fix_t_imp(path="."):

#	Load number of zones and existing plasma files

	if((len(path) > 0) and (path[-1] != "/")): path = path + "/"
	Plasmas	= load_plasma_files(path, DeNorm=False)

	nZones	= len(Plasmas[0])

#	rescale

	for i in range(2,len(Plasmas)):
		for k in range(nZones):
			Plasmas[i][k].Temp = np.copy(Plasmas[1][k].Temp)
			Plasmas[i][k].pi_parallel = np.zeros_like(Plasmas[i][k].pi_parallel)
			Plasmas[i][k].Gamma	= np.zeros_like(Plasmas[i][k].Gamma)

	save_plasma_files(path+"Results/", Plasmas)

	return