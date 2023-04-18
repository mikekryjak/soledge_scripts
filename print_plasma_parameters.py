#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from print_plasma_parameters 	import print_plasma_parameters
	from routines.cli_routines	import *

	if(cli_present("-h",sys.argv)):
		print("\nThis script print parameters available in plasma file\n")
		print("print_plasma_parameters options")
		print("\t-path        Directory with simulation or list of directories[d='./']")
		print()
		exit()

	path	 	= cli_get_value("-path",				sys.argv,   "")
	print_plasma_parameters(path=path)
	exit()

#=======================================

# Function definition is here

import types
import os
import h5py
import numpy						as np

from files.load_plasma_files		import load_plasma_files
from files.load_ions_list			import load_ions_list

#==============================================================================
# This routine prints useful informations
#==============================================================================

def print_plasma_parameters(path=""):

	print("print_plasma_parameters\n")

	if((len(path) > 0) and (path[-1] != "/")): path = path + "/"

#	Read mesh

	Plasmas = load_plasma_files(path)
	ions	= load_ions_list(path)

	print("\nSoledge mesh parameters:")
	for iPlasma in range(len(Plasmas)):
		print("\tPlasma = ",ions[iPlasma])
		for iPar in range(len(Plasmas[iPlasma][0].VNames)):
			print("\t\t",Plasmas[iPlasma][0].VNames[iPar])

	print("\nSoledge triangles parameters:")
	for iPlasma in range(len(Plasmas)):
		print("\tPlasma = ",ions[iPlasma])
		for iPar in range(len(Plasmas[iPlasma][0].Triangles.VNames)):
			print("\t\t",Plasmas[iPlasma][0].Triangles.VNames[iPar])

	print("\nSoledge wall parameters:")
	for iPlasma in range(len(Plasmas)):
		print("\tPlasma = ",ions[iPlasma])
		for iPar in range(len(Plasmas[iPlasma][0].Wall.VNames)):
			print("\t\t",Plasmas[iPlasma][0].Wall.VNames[iPar])

	return

