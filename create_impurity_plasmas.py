#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from create_impurity_plasmas	import create_impurity_plasmas
	from routines.cli_routines		import cli_get_value, cli_present
	 
	if(cli_present("-h",sys.argv)):
		print("\nThis script create and impurity and corresponding plasma files\n")
		print("create_impurity_plasmas options")
		print("\t-path <value>     Path of mesh [d='.']")
		print("\t-specie <value>   Impurity specie (He, Li, Be, B, C, N, O, Ne, Ar, Xe) [d='C']")
		print("\t-zeff <value>     Requested Zeff [d=1.]")
		print("\t-ne <value>       Requested electron density from impurity [d=5e18]")
		print("\t-factne <value>   Requested fraction of electron from impurity [d=0.05]")
		print("\t-ni <value>       Requested impurity density (sum over all ionization levels) [d=1e18]")
		print("\t-factni <value>   Requested fraction of D/H impurity density (sum over all ionization levels)  [d=0.01]")
		print("\t-clean_old        Remove previous impurities [d=false]")
		print("\t-keep_time        Keep original time (or set to 0) [d=false]")
		print()
		exit()

	path	 	= cli_get_value("-path",		sys.argv,  ".")
	specie	 	= cli_get_value("-specie",	sys.argv,  "C")
	zeff	 	= cli_get_value("-zeff",		sys.argv,  2.0)
	ne	 		= cli_get_value("-ne",			sys.argv,  5e18)
	factne	 	= cli_get_value("-factne",	sys.argv,  0.05)
	ni	 		= cli_get_value("-ni",			sys.argv,  1e18)
	factni	 	= cli_get_value("-factni",	sys.argv,  0.01)
	ip_ne		= cli_present("-ne",			sys.argv)
	ip_factne	= cli_present("-factne",		sys.argv)
	ip_ni		= cli_present("-ni",			sys.argv)
	ip_factni	= cli_present("-factni",		sys.argv)
	clean_old	= cli_present("-clean_old",		sys.argv)
	keep_time	= cli_present("-keep_time",		sys.argv)

	mode = 0
	if(ip_ne != 0): 	mode = 1
	if(ip_factne != 0): mode = 2
	if(ip_ni != 0): 	mode = 3
	if(ip_factni != 0): mode = 4

	create_impurity_plasmas(path=path, specie=specie, inp_zeff=zeff, inp_ne=ne, inp_factne=factne, inp_ni=ni, inp_factni=factni, mode=mode, clean_old=clean_old, keep_time=keep_time)
	exit()

#=======================================

# Function definition is here

import os
import types
import h5py
import numpy							as np
from files.load_soledge_mesh_file		import load_soledge_mesh_file
from files.load_refpar_file				import load_refpar_file
from files.load_plasma_files			import load_plasma_files
from files.save_plasma_files			import save_plasma_files
from files.load_ions_list				import load_ions_list
from files.save_ions_list				import save_ions_list
from mesh.get_core_sep_megazones		import get_core_sep_megazones
from routines.h5_routines				import h5_write
from routines.globals					import DEBUG

IONS		= ("He", "Li", "Be",  "B",  "C",   "N",    "O",   "Ne",  "Ar", "Xe")
ION_MASSES	= (4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007, 15.999, 20.180, 39.948, 131.29) 

def create_impurity_plasmas(path="./", specie="C", inp_zeff=2.0, inp_ne=5.e18, inp_factne=0.05, inp_ni=1.e18, inp_factni=0.01, mode=0, clean_old=0, keep_time=0):

#	Load number of zones and existing plasma files

	if((len(path) > 0) and (path[-1] != "/")): path = path + "/"
	RefPar	= load_refpar_file(path+"Results/")
	ions	= load_ions_list(path)

	Plasmas	= load_plasma_files(path,  DeNorm=False)
	nTotOldPlasmas = len(Plasmas)

	nZones		= len(Plasmas[0])
	if(clean_old != 0): 
		ions=ions[:2]
		Plasmas = Plasmas[:2]
		for k in range(nZones):
			Plasmas[0][k].Dens = Plasmas[1][k].Dens*Plasmas[1][0].charge

	nOldPlasmas = len(Plasmas)

	if(keep_time == 0):
		for k in range(nOldPlasmas): Plasmas[k][0].tempus	= 0.

#	Read corona distribution

	FilePath = os.path.dirname(os.path.abspath(__file__))
	CorData = np.loadtxt(FilePath + "/data/corona_"+specie+".csv", delimiter=",", skiprows=1)
	CorTe	= CorData[:,0]
	CorFra	= CorData[:,1:]
	nStages = CorFra.shape[1]

	if(mode == 0):													#For Zeff comppute average Ni fraction at separatrix

		Config = load_soledge_mesh_file(path+"mesh.h5")			#Read mesh
		CoreMegazone, SepMegazone = get_core_sep_megazones(Config)	#Find separatrix zone
		Megazones = Config.Megazones

		nSepZones 	= len(Megazones[SepMegazone].list)
		nzZones	  	= np.empty(nSepZones, dtype = 'i4')
		FactNiZones	= np.empty(nSepZones, dtype = 'f8')
		nNnorms		= np.empty(nStages-1, dtype = 'f8')
		
		for iZone in range(nSepZones):								#Compute average Ni fraction at seapratrix
			k = Megazones[SepMegazone].list[iZone]
			Temp = Plasmas[0][k].Temp*RefPar.T0eV

			Frac = np.empty((Temp.shape[1],nStages-1), dtype='f8')
			NZ   = np.zeros(Temp.shape[1], dtype='f8')
			NZ2  = np.zeros(Temp.shape[1], dtype='f8')
			for i in range(nStages-1):
				Frac[:,i] = np.interp(Temp[-2,:], CorTe, CorFra[:,i])
				NZ  += Frac[:,i]*(i+1)
				NZ2 += Frac[:,i]*(i+1)**2

			for i in range(nStages-1):
				nNnorms[i] = np.mean((inp_zeff - 1.)/(NZ2 - NZ*inp_zeff)*Frac[:,i])

			FactNiZones[iZone] = np.sum(nNnorms)
			nzZones[iZone]	   = Temp.shape[1]

		inp_factni = np.sum(FactNiZones*nzZones)/np.sum(nzZones)
		mode = 4													#switch to mode fraction of ion density
		if(DEBUG > 1):
			print("\tComputed inp_factni = {:0.4f}".format(inp_factni))

	for i in range(1,nStages):
		ions.append(specie+"{:d}+".format(i))

		Plasmas.append([])
		for k in range(nZones):
			Plasmas[-1].append(types.SimpleNamespace())

		Plasmas[-1][0].charge	= i
		Plasmas[-1][0].mass		= ION_MASSES[IONS.index(specie)]
		Plasmas[-1][0].tempus	= Plasmas[0][0].tempus

	for k in range(nZones):
		
		Temp = Plasmas[0][k].Temp*RefPar.T0eV
		Dens = Plasmas[1][k].Dens

		Frac = np.empty((Temp.shape[0],Temp.shape[1],nStages), dtype='f8')
		NZ= np.zeros_like(Temp)

		Frac[:,:,0] = np.interp(Temp, CorTe, CorFra[:,0])
		
		for i in range(1,nStages):
			Frac[:,:,i] = np.interp(Temp, CorTe, CorFra[:,i])
			NZ  += Frac[:,:,i]*i

		if(mode == 0):									#mode Zeff
			NZ2= np.zeros_like(Temp)
			for i in range(1,nStages):	NZ2 += Frac[:,:,i]*i**2
			Nnorm = Dens*(inp_zeff - 1.)/(NZ2 - NZ*inp_zeff)

		elif(mode == 1):								#electron density from impurity
			Nnorm = inp_ne/NZ/RefPar.n0

		elif(mode == 2):								#fraction of eletron density from impurity	
			Nnorm = inp_factne*Dens/NZ

		elif(mode == 3):								#total ion density from impurity	
			NFrac= np.zeros_like(Temp)
			for i in range(1,nStages):
				NFrac  += Frac[:,:,i]

			Nnorm = inp_ni/NFrac/RefPar.n0

		elif(mode == 4):								#fraction of ion density from impurity	
			NFrac= np.zeros_like(Temp)
			for i in range(1,nStages):
				NFrac  += Frac[:,:,i]

			Nnorm = inp_factni*Dens/NFrac
#			print("zone=",k," Dens.max=",Dens.max()*RefPar.n0*1e-19," Nnorm.max=",Nnorm.max())

		for i in range(nStages-1):

			Plasmas[nOldPlasmas+i][k].Dens	  = Nnorm*Frac[:,:,i+1]									#skip neutral stage
#			print("\tIplasma=",nOldPlasmas+i," Dens.max=",Plasmas[nOldPlasmas+i][k].Dens.max()*RefPar.n0*1e-19)
			Plasmas[nOldPlasmas+i][k].Temp	  = Plasmas[1][k].Temp
					
			Plasmas[nOldPlasmas+i][k].Gamma  = np.zeros_like(Plasmas[1][k].Gamma)
		
			Plasmas[nOldPlasmas+i][k].Alpham = Plasmas[1][k].Alpham
			Plasmas[nOldPlasmas+i][k].Alphap = Plasmas[1][k].Alphap

		Plasmas[0][k].Dens	  += Nnorm*NZ				#add electrons from ions to electron density

	
	save_plasma_files(path+"Results/", Plasmas)
	save_ions_list(path, ions)

	if(keep_time == 0):
#		set to zero tempus time

		tempus 	  = 0.										
		if_global = h5py.File(path+"Results/globals", "w")
		h5_write(if_global, "/tempus", tempus)
		if_global.close()

#	Remove residuals and balances files because they are incompatible with the new number of ions 

	try:												
		os.remove(path+"/residuals")												
		for k in range(nTotOldPlasmas): os.remove(path+"/balances_{:d}".format(k))
	except:
		pass

	return