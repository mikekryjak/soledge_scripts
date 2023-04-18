#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from fix_ti	import fix_ti
	from routines.cli_routines		import cli_get_value, cli_present
	 
	if(cli_present("-h",sys.argv)):
		print("\nThis script reduces Ti a factor of Te\n")
		print("fix_ti options")
		print("\t-path <value>     Path of mesh [d='.']")
		print("\t-ti_fact <value>  Max allowed ratio Ti/Te [d=2.]")
		print("\t-no_mask          no masked out of wall [d=false]")
		print("\t-d_only           Set only Ti nont impurities T [d=false]")
		print()
		exit()

	path			= cli_get_value("-path",	sys.argv,  ".")
	ti_fact		= cli_get_value("-ti_fact",	sys.argv,  2.0)
	no_mask	= cli_present("-no_mask",	sys.argv)
	d_only		= cli_present("-d_only",				sys.argv)

	fix_ti(path=path, ti_fact=ti_fact, no_mask=0, d_only=d_only)
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

def fix_ti(path=".", ti_fact=2., no_mask=0, d_only=0):

#	Load number of zones and existing plasma files

	if((len(path) > 0) and (path[-1] != "/")): path = path + "/"
	Plasmas	= load_plasma_files(path, DeNorm=False)
	if(no_mask != 0):
		Config	= load_soledge_mesh_file(path+"mesh.h5")
		Zones  	= Config.Zones
		RefPar		= load_refpar_file(path+"Results/")
	
	if(d_only == 0): nPlasmas = len(Plasmas)
	else:				   nPlasmas = 2

	nZones	= len(Plasmas[0])

#	rescale

	for i in range(1,nPlasmas):
		if(ti_fact != 0.):  print("Set T{:d}=min(T{:d},{:f}*Te)".format(i,i,ti_fact))
		else:				    print("Set T{:d}=Te".format(i))

		for k in range(nZones):
			if(ti_fact != 0.): Plasmas[i][k].Temp = np.where(Plasmas[i][k].Temp < ti_fact*Plasmas[0][k].Temp, Plasmas[i][k].Temp, ti_fact*Plasmas[0][k].Temp)
			else:			   Plasmas[i][k].Temp = np.copy(Plasmas[0][k].Temp)

	for k in range(nZones):
		for i in range(nPlasmas):
			Plasmas[i][k].pi_parallel = np.zeros_like(Plasmas[i][k].pi_parallel)
			Plasmas[i][k].Gamma	= np.zeros_like(Plasmas[i][k].Gamma)
			Plasmas[i][k].Alpham	= np.zeros_like(Plasmas[i][k].Alpham)
			Plasmas[i][k].Alphap	= np.zeros_like(Plasmas[i][k].Alphap)

	if(no_mask != 0):
		for k in range(nZones):
			Chi2 = extend_mat(Zones[k].Chi)
			for i in range(nPlasmas):
				Plasmas[i][k].Dens   		= np.where(Chi2  == 0, Plasmas[i][k].Dens, 1.e15/RefPar.n0)
				Plasmas[i][k].Temp   		= np.where(Chi2  == 0, Plasmas[i][k].Temp, 0.1/RefPar.T0eV)

	save_plasma_files(path+"Results/", Plasmas[0:nPlasmas])

	return