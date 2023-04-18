# Function definition is here

import types
import h5py
from routines.h5_routines import h5_write

#=================================================================
# This routine to write reference parameters to the H5DF mesh file
#=================================================================

def save_refpar_file(PlasmaFilePath, RefPar): 

	if_refpar   = h5py.File(PlasmaFilePath+"reference_parameters", "w")
	h5_write(if_refpar, "/T0eV",	RefPar.T0eV)
	h5_write(if_refpar, "/n0",		RefPar.n0)
	if_refpar.close()


	if_global = h5py.File(PlasmaFilePath+"globals", "w")
	h5_write(if_global, "/tempus", RefPar.time)
	if_global.close()
	