#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from fix_gamma	import fix_gamma
	from routines.cli_routines		import cli_get_value, cli_present
	 
	if(cli_present("-h",sys.argv)):
		print("\nThis script set Gamma to zero\n")
		print("fix_ti options")
		print("\t-path <value>     Path of mesh [d='.']")
		print("\t-d_only           Set only Ti nont impurities T [d=false]")
		print()
		exit()

	path		= cli_get_value("-path",	sys.argv,  ".")
	d_only		= cli_present("-d_only",				sys.argv)

	fix_gamma(path=path, d_only=d_only)
	exit()

#=======================================

# Function definition is here

import os
import types
import numpy as np
from files.load_soledge_mesh_file		import load_soledge_mesh_file
from files.load_plasma_files			import load_plasma_files
from files.load_refpar_file				import load_refpar_file
from files.save_plasma_files			import save_plasma_files

def fix_gamma(path=".", d_only=0):

#	Load number of zones and existing plasma files

	if((len(path) > 0) and (path[-1] != "/")): path = path + "/"
	Plasmas	= load_plasma_files(path, DeNorm=False)
	
	if(d_only == 0): nPlasmas = len(Plasmas)
	else:			 nPlasmas = 2

	nZones	= len(Plasmas[0])

	for k in range(nZones):
		for i in range(nPlasmas):
			Plasmas[i][k].Gamma	= np.zeros_like(Plasmas[i][k].Gamma)

	save_plasma_files(path+"Results/", Plasmas[0:nPlasmas])

	return