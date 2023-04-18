#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from routines.cli_routines import cli_present, cli_get_value, cli_get_value
	from remesh_plasma import remesh_plasma
	 
	if(cli_present("-h",sys.argv)):
		print("\nThis script converts plasma files from one mesh to another\n")
		print("remesh_plasma options")
		print("\t-old_path <value>  Path of old mesh [d='old/']")
		print("\t-new_path <value>  Path of old mesh [d='new/']")
		print("\t-no_mask           no masked out of wall [d=false]")
		print("\t-d_only            Remesh only e- and D+ [d=false]")
		print("\t-keep_time         Keep original time (or set to 0) [d=false]")
		print()
		exit()

	old_path	= cli_get_value("-old_path",	sys.argv,"old/")
	new_path	= cli_get_value("-new_path",	sys.argv,"new/")
	d_only		= cli_present("-d_only",				sys.argv)
	no_mask		= cli_present("-no_mask",	sys.argv)
	keep_time	= cli_present("-keep_time",	sys.argv)
	remesh_plasma(old_path=old_path, new_path=new_path, d_only=d_only, no_mask=no_mask, keep_time=keep_time)
	exit()

#=======================================

# Function definition is here

import os
import glob
import types
import numpy as np
from scipy.interpolate					import griddata
from files.load_soledge_mesh_file		import load_soledge_mesh_file
from files.load_plasma_files			import load_plasma_files
from files.save_plasma_files			import save_plasma_files
from files.load_ions_list				import load_ions_list
from files.save_ions_list				import save_ions_list
from files.load_refpar_file				import load_refpar_file
from files.save_refpar_file				import save_refpar_file
from routines.utils_routines			import extend_mat

def remesh_plasma(old_path="old/", new_path="new/", d_only=0, no_mask=0, keep_time=0):

#	Load old mesh and old plasma
	
	if((len(old_path) > 0) and (old_path[-1] != "/")): old_path = old_path + "/"
	if((len(new_path) > 0) and (new_path[-1] != "/")):	new_path = new_path + "/"
	OldMeshFile	= old_path + "mesh.h5"
	OldConfig	= load_soledge_mesh_file(OldMeshFile)
	OldZones  	= OldConfig.Zones
	OldPlasmas	= load_plasma_files(old_path, len(OldZones), DeNorm=False)
	RefPar		= load_refpar_file(old_path+"Results/")
	Ions		= load_ions_list(old_path)

	nIons	  = len(OldPlasmas)
	if(d_only != 0): nIons = 2
	Ions = Ions[:nIons]

#load new mesh

	NewMeshFile = new_path + "mesh.h5"
	NewConfig	= load_soledge_mesh_file(NewMeshFile)
	NewZones 	= NewConfig.Zones

#	Compute number of odl points

	OldRZ, OldInOut, OldValues = prepare_arrays(OldZones)	
	NewRZ, NewInOut, NewValues = prepare_arrays(NewZones)	

	OldPlasmaValues = set_plasma_values(OldPlasmas)
	
	NewPlasmas = []
	for i in range(nIons):
		NewPlasmaValues = []
		for l in range(len(OldPlasmaValues[i])):
			NewPlasmaValues.append([])
			
#			prepare old values
			nTot  = 0
			for k in range(len(OldZones)):
				n  = OldZones[k].gridRc.size
				if(OldZones[0].gridRc.shape[0] == OldPlasmaValues[i][l][0].shape[0]):
					OldValues[nTot:nTot+n] = OldPlasmaValues[i][l][k].reshape(n)
				elif(OldZones[0].gridRc.shape[0]+2 == OldPlasmaValues[i][l][0].shape[0]):
					Nx = OldZones[k].gridRc.shape[0]
					Nz = OldZones[k].gridRc.shape[1]
					OldValues[nTot:nTot+n] = OldPlasmaValues[i][l][k][1:Nx+1, 1:Nz+1].reshape(n)
				else:
					print("\tError!!! mismach between Zone size and parameter")
					exit()

				nTot += n
				
			for k in range(2):
				NewValues[NewInOut[k]] = griddata(OldRZ[OldInOut[k],:], OldValues[OldInOut[k]], NewRZ[NewInOut[k],:], method='linear')
				iNaN=np.where(np.isnan(NewValues[NewInOut[k]]))[0]
				if(len(iNaN) > 0):
					NewValues[NewInOut[k][iNaN]] = griddata(OldRZ[OldInOut[k],:], OldValues[OldInOut[k]], NewRZ[NewInOut[k][iNaN],:], method='nearest')
				
			
			nTot  = 0
			for k in range(len(NewZones)):
				n = NewZones[k].gridRc.size
				if(OldZones[0].gridRc.shape[0] == OldPlasmaValues[i][l][0].shape[0]):
					NewPlasmaValues[l].append(np.copy(NewValues[nTot:nTot+n].reshape(NewZones[k].gridRc.shape)))
				elif(OldZones[0].gridRc.shape[0]+2 == OldPlasmaValues[i][l][0].shape[0]):
					NewPlasmaValues[l].append(extend_mat(np.copy(NewValues[nTot:nTot+n].reshape(NewZones[k].gridRc.shape))))
				else:
					print("\tError!!! mismach between Zone size and parameter")
					exit()			
				nTot += n
	
		NewPlasmas.append([])
		for k in range(len(NewZones)):
			NewPlasmas[i].append(types.SimpleNamespace())

			if(no_mask != 0):
				Chi2 = extend_mat(NewZones[k].Chi)
				NewPlasmas[i][k].Dens   		= np.where(Chi2  == 0, NewPlasmaValues[0][k], 1.e15/RefPar.n0)
				NewPlasmas[i][k].Temp   		= np.where(Chi2  == 0, NewPlasmaValues[1][k], 0.1/RefPar.T0eV)
				NewPlasmas[i][k].pi_parallel 	= np.where(Chi2  == 0, NewPlasmaValues[2][k], 0.)
				NewPlasmas[i][k].Gamma  		= np.where(Chi2  == 0, NewPlasmaValues[3][k], 0.)
				NewPlasmas[i][k].Alpham			= np.where(NewZones[k].Chi  == 0, NewPlasmaValues[4][k], 0.)
				NewPlasmas[i][k].Alphap			= np.where(NewZones[k].Chi  == 0, NewPlasmaValues[5][k], 0.)
				
			else:
				NewPlasmas[i][k].Dens   		= np.copy(NewPlasmaValues[0][k])
				NewPlasmas[i][k].Temp   		= np.copy(NewPlasmaValues[1][k])
				NewPlasmas[i][k].pi_parallel 	= np.copy(NewPlasmaValues[2][k])
				NewPlasmas[i][k].Gamma  		= np.copy(NewPlasmaValues[3][k])
				NewPlasmas[i][k].Alpham			= np.copy(NewPlasmaValues[4][k])
				NewPlasmas[i][k].Alphap			= np.copy(NewPlasmaValues[5][k])

		NewPlasmas[i][0].charge   = OldPlasmas[i][0].charge
		NewPlasmas[i][0].mass     = OldPlasmas[i][0].mass 
		if(keep_time == 0): 		NewPlasmas[i][0].tempus   = 0. 
		else:								NewPlasmas[i][0].tempus   = OldPlasmas[i][0].tempus

	save_plasma_files(new_path+"Results/", NewPlasmas)
	save_refpar_file(new_path+"Results/", RefPar)
	save_ions_list(new_path, ions=Ions)

	if(keep_time == 0):
		try:
			os.remove(new_path+"residuals")
		except:
			pass
		balances = glob.glob(new_path+"balances_*")
		for balance in balances: os.remove(balance)

	return


#============================================================
			
def set_plasma_values(Plasmas):
	nIons = len(Plasmas)
	
	PlasmaValues = []
	for i in range(nIons):
		PlasmaValues.append([])
		for l in range(6): PlasmaValues[i].append([])
		for k in range(len(Plasmas[0])):
			PlasmaValues[i][0].append(Plasmas[i][k].Dens)
			PlasmaValues[i][1].append(Plasmas[i][k].Temp)
			PlasmaValues[i][2].append(Plasmas[i][k].pi_parallel)
			PlasmaValues[i][3].append(Plasmas[i][k].Gamma)
			PlasmaValues[i][4].append(Plasmas[i][k].Alpham)
			PlasmaValues[i][5].append(Plasmas[i][k].Alphap)
			
	return PlasmaValues

#============================================================

def prepare_arrays(Zones):
	nTot  = 0
	for k in range(len(Zones)): nTot += Zones[k].gridRc.size
	RZ		= np.empty((nTot,2), dtype='f8')
	Chi		= np.empty(nTot, dtype='i4')
	Values	= np.empty(nTot, dtype='f8')

	nTot  = 0
	for k in range(len(Zones)):
		n = Zones[k].gridRc.size
		RZ[nTot:nTot+n,0]  = Zones[k].gridRc.reshape(n)
		RZ[nTot:nTot+n,1]  = Zones[k].gridZc.reshape(n)
		Chi[nTot:nTot+n] = Zones[k].Chi.reshape(n)
		nTot += n

	In  = np.where(Chi != 1)[0]
	Out = np.where(Chi == 1)[0]

	return RZ, (In,Out), Values