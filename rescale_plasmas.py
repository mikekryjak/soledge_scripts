#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from rescale_plasmas	import rescale_plasmas
	from routines.cli_routines		import cli_get_value, cli_present
	 
	if(cli_present("-h",sys.argv)):
		print("\nThis script rescale plasma parameters (n & T) in plasma file\n")
		print("rescale_plasmas options")
		print("\t-path <value>     Path of mesh [d='.']")
		print("\t-n_fact <value>   Requested factor for density [d=1.]")
		print("\t-t_fact <value>   Requested factor for temperature [d=1.]")
		print("\t-plasmas <value>  Plasmas to scale [d=[0,1]]")
		print("\t-keep_time        Keep original time (or set to 0) [d=false]")
		print()
		exit()

	path	 	= cli_get_value("-path",		sys.argv,  ".")
	n_fact	 	= cli_get_value("-n_fact",	sys.argv,  1.0)
	t_fact	 	= cli_get_value("-t_fact",	sys.argv,  1.0)
	plasmas	 	= cli_get_value("-plasmas",	sys.argv,  [0,1])
	keep_time	= cli_present("-keep_time",	sys.argv)

	if((n_fact == 1.) and (t_fact == 1.)):
		print("rescale_plasmas: Nothing to rescale!")
		exit()

	rescale_plasmas(path=path, n_fact=n_fact, t_fact=t_fact, plasmas=plasmas, keep_time=keep_time)
	exit()

#=======================================

# Function definition is here

import os
import types
import numpy as np
from files.load_soledge_mesh_file		import load_soledge_mesh_file
from files.load_plasma_files			import load_plasma_files
from files.save_plasma_files			import save_plasma_files

def rescale_plasmas(path=".", n_fact=1., t_fact=1., plasmas=[0,1], keep_time=0):

#	Load number of zones and existing plasma files

	if((len(path) > 0) and (path[-1] != "/")): path = path + "/"
	PlasmasData	= load_plasma_files(path, DeNorm=False)
	nZones	= len(PlasmasData[0])

#	rescale

	for iPlasma in plasmas:
		print("\tRescaling plasma #",iPlasma)
		Charge = PlasmasData[iPlasma][0].charge
		for k in range(nZones):

			if((n_fact != 1.) and (iPlasma != 0)):
				PlasmasData[0][k].Dens += Charge*(n_fact - 1.)*PlasmasData[iPlasma][k].Dens
				PlasmasData[iPlasma][k].Dens *= n_fact

			if(t_fact != 1.): PlasmasData[iPlasma][k].Temp *= t_fact

			PlasmasData[iPlasma][k].pi_parallel = np.zeros_like(PlasmasData[iPlasma][k].pi_parallel)		

	if(keep_time == 0):
		nPlasmasData = len(PlasmasData)
		for iPlasma in range(nPlasmasData): PlasmasData[iPlasma][0].tempus = 0.

	save_plasma_files(path+"Results/", PlasmasData)

	if(keep_time == 0):
		os.remove(path+"residuals")
		for iPlasma in range(nPlasmasData): os.remove(path+"balances_{:d}".format(iPlasma))

	return