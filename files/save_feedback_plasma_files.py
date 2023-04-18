#!/usr/bin/python

# Function definition is here

import numpy as np
import types
import os
import h5py
from files.save_refpar_file		import save_refpar_file
from routines.h5_routines		import h5_write
from routines.utils_routines	import extend_mat
from routines.globals			import DEBUG

#=========================================================
# This routine to write data to the H5DF mesh file
#=========================================================

def save_feedback_plasma_files(PlasmaFilePath, Config): 

	if(DEBUG > 0): print("save_plasma_files: Saving to ",PlasmaFilePath)

	try:
		os.mkdir(PlasmaFilePath)
	except OSError:
		pass

	Zones			= Config.Zones
	FeedbackTransp	= Config.FeedbackTransp

#	plasma files
	ions = ["e","D+"]
	masses  = [0,2]
	charges = [-1,1]

	for k in range(len(Zones)):
		Zones[k].FeedbackValues[0] = np.where(Zones[k].Chi  == 0, Zones[k].FeedbackValues[0], 1.e15/FeedbackTransp.Data.NRef)
		Zones[k].FeedbackValues[1] = np.where(Zones[k].Chi  == 0, Zones[k].FeedbackValues[1], 1.e15/FeedbackTransp.Data.NRef)
		Zones[k].FeedbackValues[2] = np.where(Zones[k].Chi  == 0, Zones[k].FeedbackValues[2], 0.1/FeedbackTransp.Data.TRef)
		Zones[k].FeedbackValues[3] = np.where(Zones[k].Chi  == 0, Zones[k].FeedbackValues[3], 0.1/FeedbackTransp.Data.TRef)

	for Ion in range(2):
			
		if_plasma = h5py.File(PlasmaFilePath+"plasma_{:d}".format(Ion), "w")
	
		for k in range(len(Zones)):
			zone = "zone{:d}".format(k+1)
			if_plasma.create_group(zone)
	
			Nx  = Zones[k].FeedbackValues[0].shape[0]
			Nz  = Zones[k].FeedbackValues[0].shape[1]
	
			h5_write(if_plasma, zone+ '/Nx', Nx)
			h5_write(if_plasma, zone+ '/Nz', Nz)
	
			if(Ion == 0):
				h5_write(if_plasma, zone+ '/density', 	  	extend_mat(Zones[k].FeedbackValues[0]), order = 'F')
				h5_write(if_plasma, zone+ '/temperature', 	extend_mat(Zones[k].FeedbackValues[2]), order = 'F')
			else:
				h5_write(if_plasma, zone+ '/density', 	  	extend_mat(Zones[k].FeedbackValues[1]), order = 'F')
				h5_write(if_plasma, zone+ '/temperature', 	extend_mat(Zones[k].FeedbackValues[3]), order = 'F')
				
			ZeroMat=np.zeros((Nx+2,Nz+2), dtype='f8')
			h5_write(if_plasma, zone+ '/pi_parallel',	ZeroMat, order = 'F')
			h5_write(if_plasma, zone+ '/Gamma', 	  	ZeroMat, order = 'F')
	
			ZeroMat=np.zeros((Nx,Nz), dtype='f8')
			h5_write(if_plasma, zone+ '/alpham', 	  ZeroMat, order = 'F')
			h5_write(if_plasma, zone+ '/alphap', 	  ZeroMat, order = 'F')
			
		h5_write(if_plasma, '/charge',		charges[Ion])
		h5_write(if_plasma, '/mass',			masses[Ion])
		h5_write(if_plasma, '/tempus',		0.)

		if_plasma.close()

	RefPar		= types.SimpleNamespace()
	RefPar.n0	= FeedbackTransp.Data.NRef
	RefPar.T0eV = FeedbackTransp.Data.TRef
	RefPar.time = 0.

	save_refpar_file(PlasmaFilePath, RefPar)

	return