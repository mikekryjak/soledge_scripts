# Function definition is here

import numpy as np
import os
import h5py
from routines.h5_routines		import h5_write
from routines.utils_routines	import extend_mat
from routines.globals			import DEBUG

#=========================================================
# This write plasma values H5DF file
#=========================================================

def save_plasma_files(Path, Plasmas, skip=0): 


	if(DEBUG > 0): print("save_plasma_files: Saving to ",Path)

	try:
		os.mkdir(Path)
	except OSError:
		pass

#	plasma files

	if((len(Path) > 0) and (Path[-1] != "/")): Path = Path + "/"
	nPlasmas = len(Plasmas)
	for i in range(skip, nPlasmas):
		if_plasma = h5py.File(Path+"plasma_{:d}".format(i), "w")
	
		for k in range(len(Plasmas[0])):
			zone = "zone{:d}".format(k+1)
			if_plasma.create_group(zone)
	
			Nx  = Plasmas[i][k].Alpham.shape[0]
			Nz  = Plasmas[i][k].Alpham.shape[1]
	
			h5_write(if_plasma, zone+ '/Nx', Nx)
			h5_write(if_plasma, zone+ '/Nz', Nz)

			if(Plasmas[i][k].Dens.shape[0] == Nx+2):
				h5_write(if_plasma, zone+ '/density', 	  	Plasmas[i][k].Dens, order = 'F')
				h5_write(if_plasma, zone+ '/temperature',		Plasmas[i][k].Temp, order = 'F')
				
				try:
					h5_write(if_plasma, zone+ '/pi_parallel', 	Plasmas[i][k].pi_parallel,	order = 'F')
				except:
					h5_write(if_plasma, zone+ '/pi_parallel', 	np.zeros_like(Plasmas[i][k].Temp),	order = 'F')

				h5_write(if_plasma, zone+ '/Gamma', 	  		Plasmas[i][k].Gamma,		order = 'F')
			else:
				h5_write(if_plasma, zone+ '/density', 	  	extend_mat(Plasmas[i][k].Dens), order = 'F')
				h5_write(if_plasma, zone+ '/temperature',		extend_mat(Plasmas[i][k].Temp), order = 'F')
				
				try:
					h5_write(if_plasma, zone+ '/pi_parallel', 	extend_mat(Plasmas[i][k].pi_parallel),	order = 'F')
				except:
					h5_write(if_plasma, zone+ '/pi_parallel', 	extend_mat(np.zeros_like(Plasmas[i][k].Temp)),	order = 'F')

				h5_write(if_plasma, zone+ '/Gamma', 	  		extend_mat(Plasmas[i][k].Gamma),		order = 'F')
	
			h5_write(if_plasma, zone+ '/alpham', 	  Plasmas[i][k].Alpham, order = 'F')
			h5_write(if_plasma, zone+ '/alphap', 	  Plasmas[i][k].Alphap, order = 'F')
		

		h5_write(if_plasma, '/charge',		Plasmas[i][0].charge)
		h5_write(if_plasma, '/mass',			Plasmas[i][0].mass)
		h5_write(if_plasma, '/tempus',		Plasmas[i][0].tempus)

		if_plasma.close()

	if_global = h5py.File(Path+"globals", "w")
	tempus = Plasmas[skip][0].tempus
	h5_write(if_global, "/tempus", tempus)
	if_global.close()
	
	return
