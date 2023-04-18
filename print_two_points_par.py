#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from print_two_points_par 	import print_two_points_par
	from routines.cli_routines	import *

	if(cli_present("-h",sys.argv)):
		print("\nThis script parameters useful for two point model computation\n")
		print("print_two_points_par options")
		print("\t-path        Directory with simulation or list of directories[d='./']")
		print()
		exit()

	path	 	= cli_get_value("-path",				sys.argv,   "")
	print_two_points_par(path=path)
	exit()

#=======================================

# Function definition is here

import types
import os
import h5py
from math								import sqrt
import numpy							as np

from routines.h5_routines				import h5_read
from routines.globals					import DEBUG, KB, BALLOONING_NAMES

from mesh.get_rz_core_sep				import get_rz_core_sep
from files.load_soledge_mesh_file		import load_soledge_mesh_file

#==============================================================================
# This routine prints useful informations
#==============================================================================

def print_two_points_par(path=""):

	
	print("print_two_points_par")

	if((len(path) > 0) and (path[-1] != "/")): path = path + "/"

#	Read mesh

	Config = load_soledge_mesh_file(path+"mesh.h5")

#	Find separatrix

	RZcore, RZsep, CoreMegazone, SepMegazone, jSep = get_rz_core_sep(Config, core_and_sep = True)

	dl		= np.sqrt((RZsep[1:,0]-RZsep[:-1,0])**2 + (RZsep[1:,1]-RZsep[:-1,1])**2)
	Surface	= np.pi*np.sum(dl*(RZsep[1:,0]+RZsep[:-1,0]))
	Area	= 0.5*np.sum((RZsep[1:,1]+RZsep[:-1,1])*(RZsep[1:,0]-RZsep[:-1,0]))

	R0 		= 0.5*(np.max(RZsep[:,0])+np.min(RZsep[:,0]))
	a		= 0.5*(np.max(RZsep[:,0])-np.min(RZsep[:,0]))
	k		= (Surface/(4.*np.pi**2*a*R0))**2
	print("R0              = ",R0)
	print("a               = ",a)
	print("k               = ",k)
	print("Separatrix Surf = ",Surface)
	print("Separatrix Area = ",Area)

	print("print_two_points_par: Completed")

	return

