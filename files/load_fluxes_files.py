# Function definition is here

import numpy as np
import types
import os
import h5py
from files.load_refpar_file		import load_refpar_file
from routines.h5_routines		import h5_read
from routines.globals			import DEBUG, KB

#=========================================================
# This routine to write data to the H5DF mesh file
#=========================================================

def load_fluxes_files(Path="", nZones=0, iFluxes = [], DeNorm=True): 

	if((len(Path) > 0) and (Path[-1] != "/")): Path = Path + "/"

	if(DeNorm):
		RefPar		= load_refpar_file(Path+"../Results/")
		FluxEFact	= KB*RefPar.n0*RefPar.c0*RefPar.T0*RefPar.rs0/(2*np.pi*RefPar.R0)
		FluxGFact	= 1.
		FluxNFact	= RefPar.n0*RefPar.c0*RefPar.rs0/(2*np.pi*RefPar.R0)
	else:
		FluxEFact	= 1.
		FluxGFact	= 1.
		FluxNFact	= 1.


	if(nZones == 0):
		if_mesh = h5py.File(Path+"../mesh.h5", "r")
		nZones = h5_read(if_mesh,"NZones", keep_array=False)
		if_mesh.close()

	if(len(iFluxes) == 0):
		try:
			ions = load_ions_list(Path + "..")
			iFluxes = [k for k in range(len(ions))]
		except:
			iFluxes = [0,1]
		
#	Fluxes files

	Fluxes = []
	for i in range(len(iFluxes)):
		Fluxes.append([])

		if(DEBUG > 0): print("load_fluxes_files: Loading from ",Path+"fluxes_{:d}".format(iFluxes[i]))
		if_fluxes = h5py.File(Path+"fluxes_{:d}".format(iFluxes[i]), "r")
	
		for k in range(nZones):
			
			Fluxes[i].append(types.SimpleNamespace())
			zone = "zone{:d}".format(k+1)
			
			Fluxes[i][k].FluxE	  		= h5_read(if_fluxes, zone+ '/fluxE', order = 'F')*FluxEFact				#[Nx,Nz,4] (North, South,East,West)
			Fluxes[i][k].FluxG	  		= h5_read(if_fluxes, zone+ '/fluxG', order = 'F')*FluxGFact				#[Nx,Nz,4]
			Fluxes[i][k].Fluxn	  		= h5_read(if_fluxes, zone+ '/fluxn', order = 'F')*FluxNFact				#[Nx,Nz,4]

			Fluxes[i][k].Nx				= Fluxes[i][k].FluxE.shape[0]
			Fluxes[i][k].Nz				= Fluxes[i][k].FluxE.shape[1]

			Fluxes[i][k].Values = []
			Fluxes[i][k].Values.append(Fluxes[i][k].FluxE)			#0= FluxE
			Fluxes[i][k].Values.append(Fluxes[i][k].FluxG)			#1= FluxGE
			Fluxes[i][k].Values.append(Fluxes[i][k].Fluxn)			#2= Fluxn
	
		if_fluxes.close()

	return Fluxes
