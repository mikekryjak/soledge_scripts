# Function definition is here

import numpy as np
import types
import os
import h5py
from files.load_refpar_file		import load_refpar_file
from routines.h5_routines		import h5_read
from routines.globals			import DEBUG

#=========================================================
# This routine read transport coefficients
#=========================================================

def load_transports_coefficients(Path="", nZones=0): 

	if((len(Path) > 0) and (Path[-1] != "/")): Path = Path + "/"

	if(nZones == 0):
		if_mesh = h5py.File(Path+"../mesh.h5", "r")
		nZones = h5_read(if_mesh,"NZones", keep_array=False)
		if_mesh.close()

	RefPar		= load_refpar_file(Path+"../Results/")
	D0=(RefPar.rs0**2.)/RefPar.tau0

#	Transports Coefficients files

	if(DEBUG > 0): print("load_transports_coefficients: Loading: ",Path+"transports_coefficients")
	if_transp = h5py.File(Path+"transports_coefficients", "r")

	TranspCoeffNames =  ["Chie", "D", "Nu", "Chii", "Vpinch"] 					#Balloning names in transport coefficient
	Zones = []
	for k in range(nZones):

		Zones.append(types.SimpleNamespace())
		zone = "zone{:d}".format(k+1)
		Zones[k].Ballooning = []
		for i in range(len(TranspCoeffNames)):
			if((TranspCoeffNames[i] != "Nu") or (TranspCoeffNames[i] != "Vpinch")):
				try:
					Zones[k].Ballooning.append(h5_read(if_transp,zone+"/"+ TranspCoeffNames[i], order='F')[1:-1,1:-1]*D0)
				except:
					Zones[k].Ballooning.append(np.copy(Zones[k].Ballooning[1]))

	TranspCoeff  = types.SimpleNamespace()
	TranspCoeff.Zones = Zones

	return TranspCoeff
