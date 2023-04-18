# Function definition is here

import types
import h5py
from routines.h5_routines import h5_read

#=========================================================
# This routine to write data to the H5DF mesh file
#=========================================================

def load_refpar_file(Path): 

	if((len(Path) > 0) and (Path[-1] != "/")): Path = Path + "/"

	RefPar = types.SimpleNamespace()
	if_refpar   = h5py.File(Path+"reference_parameters", "r")
	RefPar.T0eV	= h5_read(if_refpar, "T0eV",	keep_array=False)
	RefPar.n0	= h5_read(if_refpar, "n0",		keep_array=False)
	try:
		RefPar.c0	  = h5_read(if_refpar,"c0",	keep_array=False)
		RefPar.tau0	  = h5_read(if_refpar,"tau0",	keep_array=False)
		RefPar.R0	  = h5_read(if_refpar,"R0",	keep_array=False)
		RefPar.rs0	  = h5_read(if_refpar,"rs0",	keep_array=False)
		RefPar.T0	  = h5_read(if_refpar, "T0",		keep_array=False)
	except:
		RefPar.T0	  = RefPar.T0eV
		RefPar.c0	  = 1.
		RefPar.tau0	  = 1.
		RefPar.R0	  = 2
		RefPar.rs0	  = 0.1


	if_refpar.close()

	try:
		if_glob	    = h5py.File(Path+"globals", "r")
		RefPar.time	= h5_read(if_glob, "tempus",	keep_array=False)
		if_glob.close()
	except:
		RefPar.time	= 0.
	
	return RefPar

