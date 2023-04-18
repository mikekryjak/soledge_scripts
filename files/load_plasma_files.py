# Function definition is here

import numpy as np
import types
import os
import h5py
from math					 import sqrt
from files.load_refpar_file	 import load_refpar_file
from files.load_ions_list	 import load_ions_list
from routines.h5_routines	 import h5_read
from routines.globals		 import DEBUG, EV, KB
from files.load_soledge_mesh_file		import load_soledge_mesh_file

#=========================================================
# This routine to write data to the H5DF mesh file
#=========================================================

def load_plasma_files(Path="/", nZones=0, iPlasmas = [], Evolution=0, DeNorm=True, ToKnodes = 0): 

	if((len(Path) > 0) and (Path[-1] != "/")): Path = Path + "/"

	if(Evolution == 0): 	Results_dir = Path+"Results/"
	elif (Evolution == -1):	Results_dir = Path+"Progress/"
	elif (Evolution == -2):	Results_dir = Path+"Results_last/"
	else:					Results_dir = Path+"Evolution/"

	if(nZones == 0):
		if_mesh = h5py.File(Path+"mesh.h5", "r")
		nZones = h5_read(if_mesh,"NZones", keep_array=False)
		if_mesh.close()

	try:
#	if(True):

#		Load dvol to convert from triangles to quadrangles
		if_metric = h5py.File(Path+"Results/metric", "r")

		SolVolMet = []
		for k in range(nZones): 
			SolVolMet.append(h5_read(if_metric, "zone{:d}/dvol".format(k+1), order = 'F'))	#[Nx,Nz]
		if_metric.close()

#		Load TriKnots to convert from triangles to quadrangles

		if_triangles = h5py.File(Path+"triangles.h5", "r")
		RKnots = h5_read(if_triangles,"knots/R")*0.01												#cm to m
		ZKnots = h5_read(if_triangles,"knots/Z")*0.01												#cm to m

		TriKnots = h5_read(if_triangles,"triangles/tri_knots") - 1								#Knots index of each triangle (Matlab/Fortan to python indexes)
		TriToSol   = h5_read(if_triangles,"triangles/back_interp") - 1						#Soledge index of quadrangles contain triangles (Matlab/Fortan to python indexes)
		if_triangles.close()

		nTriangles	= TriToSol.shape[0]
		SolToTri	= []
		nSolToTri	= []
		for k in range(nZones):
			SolToTri.append(-np.ones((SolVolMet[k].shape[0],SolVolMet[k].shape[1],6), dtype='i4'))		#Index of triangles inside quadrangles
			nSolToTri.append(np.zeros((SolVolMet[k].shape[0],SolVolMet[k].shape[1]), dtype='i4'))
			ij = np.where(TriToSol[:,0] == k)[0]
			for iTri in range(len(ij)):
				ti = TriToSol[ij[iTri],1]
				tj = TriToSol[ij[iTri],2]
				SolToTri[-1][ti,tj,nSolToTri[-1][ti,tj]] = ij[iTri]
				nSolToTri[-1][ti,tj] += 1
			nSolMax = np.max(nSolToTri[-1])
			if(nSolMax > 0): SolToTri[-1] = SolToTri[-1][:,:,:nSolMax]+1												#+1 to skip zero size TriVolOff
			else:			 SolToTri[-1] = np.empty((0,0), dtype = 'f8') 
			
		IndTri = [0, 1, 2, 0]
		TriVol = np.zeros((TriKnots.shape[0]), dtype='f8')
		for k in range(3):
			TriVol +=  (RKnots[TriKnots[:,IndTri[k+1]]] - RKnots[TriKnots[:,IndTri[k]]])*(ZKnots[TriKnots[:,IndTri[k+1]]] + ZKnots[TriKnots[:,IndTri[k]]])
		TriVol *= np.pi*np.sum(RKnots[TriKnots[:,:]],axis=1)/3.
		TriVolOff = np.append(0.,np.abs(TriVol))
		TriVol = 0

		SolVol = []
		for k in range(nZones): 
			if(SolToTri[k].size>0):
				SolVol.append(np.sum(TriVolOff[SolToTri[k]], axis=2))
				SolVol[k] = np.where(SolVol[k] > 0., SolVol[k], SolVolMet[k])
			else:
				SolVol.append(np.copy(SolVolMet[k]))

		del SolVolMet
		"""
#		Check Volumes

		Config = load_soledge_mesh_file(Path+"mesh.h5")
		Zones	= Config.Zones

		for k in range(nZones): 
			if(len(SolToTri[k])>0):
				ii,jj=np.where(SolVol[k]/SolVolMet[k]>2)
				if(len(ii) > 0):
					for n in range(len(ii)):
						i = ii[n]
						j = jj[n]
						print("k,i,j           =",k,i,j)
						print("gridRc          =",Zones[k].gridRc[i,j])
						print("gridzc          =",Zones[k].gridZc[i,j])
						print("gridR           =",Zones[k].gridR[i:i+2,j:j+2])
						print("gridZ           =",Zones[k].gridZ[i:i+2,j:j+2])
						print("SolVolMet       =",SolVolMet[k][i,j])
						print("SolVol          =",SolVol[k][i,j])
						print("SolVol/SolVolMet=",SolVol[k][i,j]/SolVolMet[k][i,j])
						for nTri in range(nSolToTri[k][i,j]):
							iTri = SolToTri[k][i,j,nTri]-1
							print("\tRKnots =",RKnots[TriKnots[iTri,:]])
							print("\tZKnots =",ZKnots[TriKnots[iTri,:]])
							print("\tTriVol =",TriVolOff[iTri+1])
			"""
	except:
		pass

	try:
		if(ToKnodes > 0):													#Convert to Knots
			KnotsTri = TriKnots.reshape(-1)
			unique, KnotsCounts = np.unique(KnotsTri, return_counts=True)
			nKnots	 = unique.max()+1
			unique = 0
		else: 																#Converto to triangles
			KnotsTri	= TriKnots
			KnotsCounts =  0
			nKnots = -(TriKnots.max()+1)
	except:
		nKnots = 0

	try:
		ions = load_ions_list(Path)
	except:
		ions = ["e","D+"]
		masses  = [0,2]
		charges = [-1,1]

	if(len(iPlasmas) == 0):	iPlasmas = [k for k in range(len(ions))]

	if(DeNorm):
		RefPar		= load_refpar_file(Path+"Results")
		DenFact		= RefPar.n0
		TempFact	= RefPar.T0eV
		VFact		= RefPar.c0
		RadFact		= EV*RefPar.T0eV*RefPar.n0/RefPar.tau0
		PnFact		= EV
		SnFact		= RefPar.n0/RefPar.tau0
		SEFact		= KB*RefPar.T0*RefPar.n0/RefPar.tau0
		MFact		= sqrt(RefPar.T0eV)/RefPar.c0
	else:
		DenFact		= 1.
		TempFact	= 1.
		VFact		= 1.
		RadFact		= 1.
		PnFact		= 1.
		SnFact		= 1.
		SEFact		= 1.
		MFact		= 1.

#	Plasmas files

	if((Evolution >= -2) and (Evolution <= 0)):
		try:
			if_neutrals = h5py.File(Results_dir+"eirene_neutrals", "r")		
			NeutralsFile = True
			iNeutral 	 = 0
		except:
			NeutralsFile = False

	Plasmas = []
	for i in range(len(iPlasmas)):
		Plasmas.append([])

		if((Evolution >= -2) and (Evolution <= 0)):
			if(DEBUG > 0): print("\tload_plasma_files: Loading: ",Results_dir+"plasma_{:d}".format(iPlasmas[i]))
			if_plasma = h5py.File(Results_dir+"plasma_{:d}".format(iPlasmas[i]), "r")

		else:
			if(DEBUG > 0): print("\tload_plasma_files: Loading: ",Results_dir+"{:d}_plasma_{:d}".format(Evolution,iPlasmas[i]))
			if_plasma = h5py.File(Results_dir+"{:d}_plasma_{:d}".format(Evolution,iPlasmas[i]), "r")

		charge	= h5_read(if_plasma, '/charge',	keep_array=False)
		for k in range(nZones):

			Plasmas[i].append(types.SimpleNamespace())
			zone = "zone{:d}".format(k+1)
			
			Plasmas[i][k].Dens	  		= h5_read(if_plasma, zone+ '/density', order = 'F')*DenFact					#[Nx+2,Nz+2]
			Plasmas[i][k].Temp	  		= h5_read(if_plasma, zone+ '/temperature', order = 'F')*TempFact		#[Nx+2,Nz+2]
			try:
				Plasmas[i][k].pi_parallel	= h5_read(if_plasma, zone+ '/pi_parallel',   order = 'F')						#[Nx+2,Nz+2]
			except:
				Plasmas[i][k].pi_parallel	= np.zeros_like(Plasmas[i][k].Temp)

			Plasmas[i][k].Gamma 		= h5_read(if_plasma, zone+ '/Gamma', 	  		order = 'F')								#[Nx+2,Nz+2] (Gamma=velocity*density)

			Plasmas[i][k].Alpham		= h5_read(if_plasma, zone+ '/alpham',   order = 'F')									#[Nx,Nz]
			Plasmas[i][k].Alphap		= h5_read(if_plasma, zone+ '/alphap',   order = 'F')										#[Nx,Nz]

			try:
				Plasmas[i][k].Mach  	= h5_read(if_plasma, zone+ '/mach', order = 'F')												#[Nx+2,Nz+2]
			except:
				Plasmas[i][k].Mach  	= np.zeros_like(Plasmas[i][k].Temp)
			

			try:
				Plasmas[i][k].IRad 		= h5_read(if_plasma, zone+ '/rad', order = 'F')*RadFact									#[Nx,Nz]
			except:
				Plasmas[i][k].IRad 		= np.zeros_like(Plasmas[i][k].Alpham)

			try:
				Plasmas[i][k].Sn  		= h5_read(if_plasma, zone+ '/sources/Sn', order = 'F')*SnFact					#[Nx,Nz]
				Plasmas[i][k].SG  		= h5_read(if_plasma, zone+ '/sources/SG', order = 'F')								#[Nx,Nz]
				Plasmas[i][k].SE  		= h5_read(if_plasma, zone+ '/sources/SE', order = 'F')*SEFact					#[Nx,Nz]
			except:
				Plasmas[i][k].Sn  		= np.zeros_like(Plasmas[i][k].Alpham)
				Plasmas[i][k].SG  		= np.zeros_like(Plasmas[i][k].Alpham)
				Plasmas[i][k].SE  		= np.zeros_like(Plasmas[i][k].Alpham)																			#[Nx,Nz]

			Plasmas[i][k].Nx			= Plasmas[i][k].Alphap.shape[0]
			Plasmas[i][k].Nz			= Plasmas[i][k].Alphap.shape[1]

			Plasmas[i][k].Values = [Plasmas[i][k].Dens,   Plasmas[i][k].Temp,  Plasmas[i][k].Gamma,
									Plasmas[i][k].Mach,   Plasmas[i][k].IRad,  Plasmas[i][k].Dens*Plasmas[i][k].Temp**1.5, Plasmas[i][k].pi_parallel,
									Plasmas[i][k].Alpham, Plasmas[i][k].Alphap]		#0 = Dens;  #1 = Temp;  #2 = Gamma;  #3 = Mach
																													#4 = Rad;    #5 = pi_parallel;  #6 = Alpham;  #7 = Alphap
			if(k == 0):	Plasmas[i][0].VNames =  ["Dens","Temp","velocity","M","IRad","(Dens*Temp^1.5)","pi_parallel","Alpham","Alphap"]
			Plasmas[i][k].Values.extend([Plasmas[i][k].Sn, Plasmas[i][k].SG, Plasmas[i][k].SE])		#8= Sn; #9 = SG; #10 = SE
			if(k == 0):	Plasmas[i][0].VNames.extend(["Sn", "SG", "SE"])

			try:
				Plasmas[i][k].Zhdanov_min_n  		= h5_read(if_plasma, zone+ '/counters/Zhdanov_min_density', order = 'F')						#[Nx+2,Nz+2]
				Plasmas[i][k].clean_min_n  		= h5_read(if_plasma, zone+ '/counters/clean_min_density', order = 'F')						#[Nx+2,Nz+2]
				Plasmas[i][k].clean_min_t  		= h5_read(if_plasma, zone+ '/counters/clean_min_temperature', order = 'F')						#[Nx+2,Nz+2]
				Plasmas[i][k].coupling_min_t  		= h5_read(if_plasma, zone+ '/counters/coupling_min_density', order = 'F')						#[Nx+2,Nz+2]

				Plasmas[i][k].Values.extend([Plasmas[i][k].Zhdanov_min_n, Plasmas[i][k].clean_min_n, Plasmas[i][k].clean_min_t, Plasmas[i][k].coupling_min_t])
				if(k == 0):	Plasmas[i][0].VNames.extend(["count_Zhdanov_min_n", "count_clean_min_n", "count_clean_min_t", "count_coupling_min_t"])
				if(i > 1): 
					Plasmas[i][k].clean_t_excurs  		= h5_read(if_plasma, zone+ '/counters/clean_temperature_excursion', order = 'F')						#[Nx,Nz]
					Plasmas[i][k].Values.extend([Plasmas[i][k].clean_t_excurs])
					if(k == 0):	Plasmas[i][0].VNames.extend(["count_clean_t_excurs"])
			except:
				pass

		try:
			Plasmas[i][0].charge	= charge
			Plasmas[i][0].mass		= h5_read(if_plasma, '/mass',		keep_array=False)
			Plasmas[i][0].tempus	= h5_read(if_plasma, '/tempus',	keep_array=False)
			Plasmas[i][0].ion   	= ions[iPlasmas[i]]
		except:
			Plasmas[i][0].charge	= charges[iPlasmas[i]]
			Plasmas[i][0].mass		= masses[iPlasmas[i]]
			Plasmas[i][0].ion   	= ions[iPlasmas[i]]
			Plasmas[i][0].tempus	= 0.



#		Triangles and wall 
#		===========

#		Triangles data
		try:

#			Triangles data

			Plasmas[i][0].Triangles			= types.SimpleNamespace()

			Plasmas[i][0].Triangles.Dens	 = h5_read(if_plasma, 'triangles/density')*DenFact											#on knots
			Plasmas[i][0].Triangles.Temp	 = h5_read(if_plasma, 'triangles/temperature')*TempFact								#on knots
			Plasmas[i][0].Triangles.velocity = h5_read(if_plasma, 'triangles/velocity')*VFact											#on knots
			Plasmas[i][0].Triangles.IRad 	 = h5_read(if_plasma, 'triangles/radiation')*RadFact										#on knots
			Plasmas[i][0].Triangles.IRad  	 = np.where(Plasmas[i][0].Triangles.IRad < 1.,1., Plasmas[i][0].Triangles.IRad)

			Plasmas[i][0].Triangles.Values = [Plasmas[i][0].Triangles.Dens, Plasmas[i][0].Triangles.Temp, Plasmas[i][0].Triangles.velocity, Plasmas[i][0].Triangles.IRad]
			Plasmas[i][0].Triangles.VNames = ["Dens", "Temp", "velocity", "IRad"]

#			For electron there is Zeff

			if(Plasmas[i][0].charge == -1):
				Plasmas[i][0].Triangles.Zeff  = h5_read(if_plasma, 'triangles/Zeff')																#on Knots
				Plasmas[i][0].Triangles.SE    = h5_read(if_plasma, 'triangles/SE')*SEFact														#on triangles
				try:
					Plasmas[i][0].Triangles.Sn = h5_read(if_plasma, 'triangles/Sn')*SnFact														#on triangles
				except:
					Plasmas[i][0].Triangles.Sn = np.zeros_like(Plasmas[i][0].Triangles.SE)

				Plasmas[i][0].Triangles.Values.extend([Plasmas[i][0].Triangles.Zeff, Plasmas[i][0].Triangles.Sn, Plasmas[i][0].Triangles.SE])
				Plasmas[i][0].Triangles.VNames.extend(["Zeff ", "Sn", "SE"])
				
			if(Plasmas[i][0].charge == 1):
				Plasmas[i][0].Triangles.NRad  = -h5_read(if_plasma, 'triangles/NRad')															#on triangles
				Plasmas[i][0].Triangles.Sn	  = h5_read(if_plasma, 'triangles/Sn')*SnFact															#on triangles
				Plasmas[i][0].Triangles.SE	  = h5_read(if_plasma, 'triangles/SE')*SEFact															#on triangles
				try:
					Plasmas[i][0].Triangles.SG	  = h5_read(if_plasma, 'triangles/SG')																	#on triangles
				except:
					Plasmas[i][0].Triangles.SG = np.zeros_like(Plasmas[i][0].Triangles.SE)

				Plasmas[i][0].Triangles.Values.extend([Plasmas[i][0].Triangles.NRad, Plasmas[i][0].Triangles.Sn, Plasmas[i][0].Triangles.SE, Plasmas[i][0].Triangles.SG])
				Plasmas[i][0].Triangles.VNames.extend(["NRad",   "Sn",   "SE",   "SG"])

				if(NeutralsFile):
					iNeutral += 1
					try:
						Plasmas[i][0].Triangles.Nn	  = h5_read(if_neutrals, 'atomic_species/dens_{:d}'.format(iNeutral))*DenFact
						Plasmas[i][0].Triangles.Tn	  = h5_read(if_neutrals, 'atomic_species/T_{:d}'.format(iNeutral))

						Plasmas[i][0].Triangles.vxn	= h5_read(if_neutrals, 'atomic_species/Vx_{:d}'.format(iNeutral))
						Plasmas[i][0].Triangles.vyn	= h5_read(if_neutrals, 'atomic_species/Vy_{:d}'.format(iNeutral))
						Plasmas[i][0].Triangles.vzn	= h5_read(if_neutrals, 'atomic_species/Vz_{:d}'.format(iNeutral))

						Plasmas[i][0].Triangles.Values.extend([ Plasmas[i][0].Triangles.Nn,   Plasmas[i][0].Triangles.Tn,
																Plasmas[i][0].Triangles.vxn,  Plasmas[i][0].Triangles.vyn, Plasmas[i][0].Triangles.vzn])
						Plasmas[i][0].Triangles.VNames.extend(["Nn", "Tn", "vxn", "vyn", "vzn"])
						Plasmas[i][0].Triangles.Pn	  = Plasmas[i][0].Triangles.Nn*Plasmas[i][0].Triangles.Tn*PnFact
					except:
						Plasmas[i][0].Triangles.Pn	  = np.zeros_like(Plasmas[i][0].Triangles.NRad)

					try:
						Plasmas[i][0].Triangles.Nm	  = h5_read(if_neutrals, 'molecular_species/dens_{:d}'.format(iNeutral))*DenFact
						Plasmas[i][0].Triangles.Tm	  = h5_read(if_neutrals, 'molecular_species/T_{:d}'.format(iNeutral))

						Plasmas[i][0].Triangles.vxm	= h5_read(if_neutrals, 'molecular_species/Vx_{:d}'.format(iNeutral))
						Plasmas[i][0].Triangles.vym	= h5_read(if_neutrals, 'molecular_species/Vy_{:d}'.format(iNeutral))
						Plasmas[i][0].Triangles.vzm	= h5_read(if_neutrals, 'molecular_species/Vz_{:d}'.format(iNeutral))

						Plasmas[i][0].Triangles.Pn	  = (Plasmas[i][0].Triangles.Nn*Plasmas[i][0].Triangles.Tn+Plasmas[i][0].Triangles.Nm*Plasmas[i][0].Triangles.Tm)*PnFact

						Plasmas[i][0].Triangles.Values.extend([	Plasmas[i][0].Triangles.Nm,   Plasmas[i][0].Triangles.Tm,   Plasmas[i][0].Triangles.Pn,
																Plasmas[i][0].Triangles.vxm,  Plasmas[i][0].Triangles.vym, Plasmas[i][0].Triangles.vzm])
						Plasmas[i][0].Triangles.VNames.extend(["Nm", "Tm",  "Pn",  "vxm", "vym", "vzm"])


						Plasmas[i][0].Triangles.Nti	  = h5_read(if_neutrals, 'test_ion_species/dens_{:d}'.format(iNeutral))*DenFact

						Plasmas[i][0].Triangles.Values.append(Plasmas[i][0].Triangles.Nti)
						Plasmas[i][0].Triangles.VNames.append("Nti")
					except:
						pass
					Plasmas[i][0].Triangles.Values.append(Plasmas[i][0].Triangles.Pn)
					Plasmas[i][0].Triangles.VNames.append("Pn")

				else:
					try:
						Plasmas[i][0].Triangles.Nn	  = h5_read(if_plasma, 'triangles/Nn')*DenFact
						Plasmas[i][0].Triangles.Nti	  = h5_read(if_plasma, 'triangles/Nti')
						Plasmas[i][0].Triangles.Nm	  = h5_read(if_plasma, 'triangles/Nm')*DenFact
						Plasmas[i][0].Triangles.Tn	  = h5_read(if_plasma, 'triangles/Tn')
						Plasmas[i][0].Triangles.Tm	  = h5_read(if_plasma, 'triangles/Tm')
						Plasmas[i][0].Triangles.Pn	  = (Plasmas[i][0].Triangles.Nn*Plasmas[i][0].Triangles.Tn+Plasmas[i][0].Triangles.Nm*Plasmas[i][0].Triangles.Tm)*PnFact

						Plasmas[i][0].Triangles.vxn	= h5_read(if_plasma, 'triangles/vxn')
						Plasmas[i][0].Triangles.vyn	= h5_read(if_plasma, 'triangles/vyn')
						Plasmas[i][0].Triangles.vzn	= h5_read(if_plasma, 'triangles/vzn')

						Plasmas[i][0].Triangles.vxm	= h5_read(if_plasma, 'triangles/vxm')
						Plasmas[i][0].Triangles.vym	= h5_read(if_plasma, 'triangles/vym')
						Plasmas[i][0].Triangles.vzm	= h5_read(if_plasma, 'triangles/vzm')

						Plasmas[i][0].Triangles.Values.extend([ 
															   Plasmas[i][0].Triangles.Nn,   Plasmas[i][0].Triangles.Nm, Plasmas[i][0].Triangles.Nti,
															   Plasmas[i][0].Triangles.Tn,   Plasmas[i][0].Triangles.Tm,  Plasmas[i][0].Triangles.Pn,
															   Plasmas[i][0].Triangles.vxn,  Plasmas[i][0].Triangles.vyn, Plasmas[i][0].Triangles.vzn,
															   Plasmas[i][0].Triangles.vxm,  Plasmas[i][0].Triangles.vym, Plasmas[i][0].Triangles.vzm])
						Plasmas[i][0].Triangles.VNames.extend(["Nn", "Nm",   "Nti", "Tn",  "Tm", "Pn",
															   "vxn", "vyn", "vzn", "vxm", "vym", "vzm"])
					except:
						Plasmas[i][0].Triangles.Pn	  = np.zeros_like(Plasmas[i][0].Triangles.NRad)
						Plasmas[i][0].Triangles.Values.append(Plasmas[i][0].Triangles.Pn)
						Plasmas[i][0].Triangles.VNames.append("Pn")

			Triangles_OK = True
		except:
			Triangles_OK = False
			print("\tWARNING: Triangles data not available")

#			Triangles to Soledge mesh
#			================

		if(Triangles_OK):
			for k in range(nZones):

				Plasmas[i][k].IRad1	= tri_mesh_conv(Plasmas[i][0].Triangles.IRad,	SolToTri[k], SolVol[k], TriVolOff, TriKnots)
				Plasmas[i][k].Values.append(Plasmas[i][k].IRad1)
				if(k == 0):	Plasmas[i][0].VNames.append("IRad1")

#			For first stage of ions there is Neutral radiation 

				if(Plasmas[i][0].charge == -1):
					Plasmas[i][k].Zeff  = tri_mesh_conv(Plasmas[i][0].Triangles.Zeff,	SolToTri[k], SolVol[k], TriVolOff, TriKnots,Name="Zeff")
					Plasmas[i][k].SE1   = tri_mesh_conv(Plasmas[i][0].Triangles.SE,	    SolToTri[k], SolVol[k], TriVolOff, TriKnots)

					Plasmas[i][k].Values.extend([Plasmas[i][k].Zeff, Plasmas[i][k].SE1])
					if(k == 0):	Plasmas[i][0].VNames.extend(["Zeff ", "SE1"])

				if(Plasmas[i][0].charge == 1):
					Plasmas[i][k].NRad  = tri_mesh_conv(Plasmas[i][0].Triangles.NRad,SolToTri[k], SolVol[k], TriVolOff, TriKnots)
					Plasmas[i][k].Values.append(Plasmas[i][k].NRad)
					if(k == 0):	Plasmas[i][0].VNames.append("NRad")
					try:
#					if(True):
						Plasmas[i][k].Nn	= tri_mesh_conv(Plasmas[i][0].Triangles.Nn,  SolToTri[k], SolVol[k], TriVolOff, TriKnots)
						Plasmas[i][k].Tn	= tri_mesh_conv(Plasmas[i][0].Triangles.Tn,  SolToTri[k], SolVol[k], TriVolOff, TriKnots)
						Plasmas[i][k].Pn	= tri_mesh_conv(Plasmas[i][0].Triangles.Pn,  SolToTri[k], SolVol[k], TriVolOff, TriKnots)
						Plasmas[i][k].vxn	= tri_mesh_conv(Plasmas[i][0].Triangles.vxn, SolToTri[k], SolVol[k], TriVolOff, TriKnots)
						Plasmas[i][k].vyn	= tri_mesh_conv(Plasmas[i][0].Triangles.vyn, SolToTri[k], SolVol[k], TriVolOff, TriKnots)
						Plasmas[i][k].vzn	= tri_mesh_conv(Plasmas[i][0].Triangles.vzn, SolToTri[k], SolVol[k], TriVolOff, TriKnots)
						Plasmas[i][k].Values.extend([
													 Plasmas[i][k].Nn, Plasmas[i][k].Tn, Plasmas[i][k].Pn, 
													 Plasmas[i][k].vxn, Plasmas[i][k].vyn, Plasmas[i][k].vzn])
						if(k == 0):	Plasmas[i][0].VNames.extend([ "Nn",  "Tn",  "Pn", "vxn", "vyn", "vzn"])

						Plasmas[i][k].Nm	= tri_mesh_conv(Plasmas[i][0].Triangles.Nm,  SolToTri[k], SolVol[k], TriVolOff, TriKnots)
						Plasmas[i][k].Tm	= tri_mesh_conv(Plasmas[i][0].Triangles.Tm,  SolToTri[k], SolVol[k], TriVolOff, TriKnots)
						Plasmas[i][k].vxm	= tri_mesh_conv(Plasmas[i][0].Triangles.vxm, SolToTri[k], SolVol[k], TriVolOff, TriKnots)
						Plasmas[i][k].vym	= tri_mesh_conv(Plasmas[i][0].Triangles.vym, SolToTri[k], SolVol[k], TriVolOff, TriKnots)
						Plasmas[i][k].vzm	= tri_mesh_conv(Plasmas[i][0].Triangles.vzm, SolToTri[k], SolVol[k], TriVolOff, TriKnots)

						Plasmas[i][k].Values.extend([
													 Plasmas[i][k].Nm, Plasmas[i][k].Tm, 
													 Plasmas[i][k].vxm, Plasmas[i][k].vym, Plasmas[i][k].vzm])
						if(k == 0):	Plasmas[i][0].VNames.extend([ "Nm",  "Tm",  "vxm", "vym", "vzm"])

						Plasmas[i][k].Nti	= tri_mesh_conv(Plasmas[i][0].Triangles.Nti, SolToTri[k], SolVol[k], TriVolOff, TriKnots)
						Plasmas[i][k].Values.append( Plasmas[i][k].Nti)
						if(k == 0):	Plasmas[i][0].VNames.append("Nti")

					except:
						pass

#			If required convert from triangles to knot or viceversa
#			===============================

			for k in range(len(Plasmas[i][0].Triangles.VNames)):
				Plasmas[i][0].Triangles.Values[k] = tri_knots_conv(Plasmas[i][0].Triangles.Values[k], KnotsTri, KnotsCounts, nKnots)
				Plasmas[i][0].Triangles.__dict__[Plasmas[i][0].Triangles.VNames[k]] = Plasmas[i][0].Triangles.Values[k] 

#		Wall data
		try:
			Plasmas[i][0].Wall					= types.SimpleNamespace()
			Plasmas[i][0].Wall.ntri	  		= h5_read(if_plasma, 'wall/ntri') - 1 #python index
			Plasmas[i][0].Wall.side	  		= h5_read(if_plasma, 'wall/side') - 1 #python index
			Plasmas[i][0].Wall.swall	  		= h5_read(if_plasma, 'wall/swall')
			Plasmas[i][0].Wall.warea	  	= h5_read(if_plasma, 'wall/warea')
			Plasmas[i][0].Wall.Dens	  		= h5_read(if_plasma, 'wall/density')
			Plasmas[i][0].Wall.Temp	  		= h5_read(if_plasma, 'wall/temperature')
			Plasmas[i][0].Wall.velocity 	= h5_read(if_plasma, 'wall/velocity')
			Plasmas[i][0].Wall.fluxE 			= h5_read(if_plasma, 'wall/fluxE')

			Plasmas[i][0].Wall.Values = [Plasmas[i][0].Wall.Dens, Plasmas[i][0].Wall.Temp, 
										 Plasmas[i][0].Wall.velocity, Plasmas[i][0].Wall.fluxE]		#0 N; #1 T; #2 V; #3 FluxE
			Plasmas[i][0].Wall.VNames =  ["Dens","Temp","V","FluxE"]

			if(iPlasmas[i]>0):
				Plasmas[i][0].Wall.fluxn		= h5_read(if_plasma, 'wall/fluxn')
				Plasmas[i][0].Wall.Values.append(Plasmas[i][0].Wall.fluxn)				#4 Fluxn
				Plasmas[i][0].Wall.VNames.append("Fluxn")

			if(iPlasmas[i]==1):
				Plasmas[i][0].Wall.albedo	= h5_read(if_plasma, 'wall/albedo')
				Plasmas[i][0].Wall.jsat		= h5_read(if_plasma, 'wall/jsat')
				Plasmas[i][0].Wall.jsat_par	= Plasmas[i][0].Wall.Dens*np.abs(Plasmas[i][0].Wall.velocity)*EV
				Plasmas[i][0].Wall.qtot		= h5_read(if_plasma, 'wall/qtot')

				Plasmas[i][0].Wall.Values.extend([Plasmas[i][0].Wall.albedo,   Plasmas[i][0].Wall.jsat,
												  Plasmas[i][0].Wall.jsat_par, Plasmas[i][0].Wall.qtot])	# 5 Albedo; #6 Jsat perpendicular; #7 Jsat parallel; #8 qtot
				Plasmas[i][0].Wall.VNames.extend(["Albedo ","Jsat ","Jsat_par ","FluxEtot "])
			
			Wall_OK = True

		except:
			Wall_OK = False
			print("\tWARNING: Wall data not available")

		if_plasma.close()

	
	if((iPlasmas[0] == 0) and (len(iPlasmas) > 1)):				#For main specie compute Mach and Plasma pressure
		for k in range(nZones):
			for i in range(1,len(iPlasmas)):
				if(Plasmas[i][0].charge == 1):
					iFirst = i
					Plasmas[iFirst][k].TDens	= np.copy(Plasmas[iFirst][k].Dens)
					Plasmas[iFirst][k].TDense	= np.copy(Plasmas[iFirst][k].Dens)
					if(i == 1):
#					Plasmas[1][k].Mach  = Plasmas[1][k].velocity/np.sqrt((Plasmas[0][k].Temp+Plasmas[1][k].Temp)/Plasmas[1][0].mass)*MFact		
						Plasmas[1][k].Pp = EV*Plasmas[1][k].Dens*(Plasmas[0][k].Temp+Plasmas[1][k].Temp)*(1.+Plasmas[1][k].Mach**2)		
						Plasmas[1][k].Ep = EV*(Plasmas[0][k].Dens*Plasmas[0][k].Temp+Plasmas[1][k].Dens*Plasmas[1][k].Temp)		
						Plasmas[1][k].Values.extend([Plasmas[1][k].Pp, Plasmas[1][k].Ep])
						if(k == 0): Plasmas[i][0].VNames.extend([ "Pp", "Ep"])
						iLast = iFirst
				else:
					Plasmas[iFirst][k].TDens +=  Plasmas[i][k].Dens
					Plasmas[iFirst][k].TDense += Plasmas[i][k].Dens*Plasmas[i][0].charge

				if(((i > 1) and (Plasmas[i][0].charge == 1)) or (i+1 == len(iPlasmas))):
					Plasmas[iLast][k].FracDens		= np.where(Plasmas[0][k].Dens > 0., Plasmas[iLast][k].TDens/Plasmas[0][k].Dens, 0.)
					Plasmas[iLast][k].Cimp		= np.where(Plasmas[1][k].Dens > 0., Plasmas[iLast][k].TDens/Plasmas[1][k].Dens, 0.)
					Plasmas[iLast][k].FraceDens		= np.where(Plasmas[0][k].Dens > 0., Plasmas[iLast][k].TDense/Plasmas[0][k].Dens, 0.)
					Plasmas[iLast][k].Zave		= np.where( Plasmas[iLast][k].TDens > 0., Plasmas[iLast][k].TDense/Plasmas[iLast][k].TDens, 0.)
					Plasmas[iLast][k].Values.extend([Plasmas[iLast][k].TDens,Plasmas[iLast][k].TDense, Plasmas[iLast][k].FracDens, Plasmas[iLast][k].Cimp, Plasmas[iLast][k].FraceDens, Plasmas[iLast][k].Zave])
					if(k == 0): Plasmas[iLast][0].VNames.extend(["TDens","TDense","FracDens","Cimp","FraceDens", "Zave"])
					iLast = iFirst
	

#	Compute total dens & radiation for all ions and for each specie
	try:
#	if(True):
		if((iPlasmas[0] == 0) and (len(iPlasmas) > 1)):				#if reading at least e and first ion compute total rad
			for i in range(len(iPlasmas)):
				if(i == 0):
					Plasmas[0][0].Triangles.TRad  = np.zeros_like(Plasmas[0][0].Triangles.IRad)																	#on knots
					Plasmas[0][0].Triangles.TNRad = np.zeros_like(Plasmas[0][0].Triangles.IRad)																	#on knots
					Plasmas[0][0].Triangles.Values.extend([Plasmas[0][0].Triangles.TRad, Plasmas[0][0].Triangles.TNRad])						#on knots
					Plasmas[0][0].Triangles.VNames.extend(["TotRad ", "TotNRad "])
				elif(Plasmas[i][0].charge == 1):
					iFirst = i
					Plasmas[iFirst][0].Triangles.TDens	= np.copy(Plasmas[iFirst][0].Triangles.Dens)																#on knots
					Plasmas[iFirst][0].Triangles.TDense	= np.copy(Plasmas[iFirst][0].Triangles.Dens)															#on knots
					Plasmas[iFirst][0].Triangles.TRad  			= Plasmas[iFirst][0].Triangles.NRad + Plasmas[iFirst][0].Triangles.IRad
					Plasmas[0][0].Triangles.TRad		+= Plasmas[iFirst][0].Triangles.TRad
					Plasmas[0][0].Triangles.TNRad		+= Plasmas[iFirst][0].Triangles.NRad
					if(i == 1):
						Plasmas[i][0].Triangles.Mach = Plasmas[1][0].Triangles.velocity/np.sqrt((Plasmas[0][0].Triangles.Temp+Plasmas[1][0].Triangles.Temp)/Plasmas[1][0].mass)*MFact		
						Plasmas[i][0].Triangles.Pp	 = EV*Plasmas[1][0].Triangles.Dens*(Plasmas[0][0].Triangles.Temp+Plasmas[1][0].Triangles.Temp)*(1.+Plasmas[1][0].Triangles.Mach**2)		
						Plasmas[i][0].Triangles.Ep	 = EV*(Plasmas[0][0].Triangles.Dens*Plasmas[0][0].Triangles.Temp+Plasmas[1][0].Triangles.Dens*Plasmas[1][0].Triangles.Temp)		

						Plasmas[i][0].Triangles.Values.extend([Plasmas[i][0].Triangles.TRad,Plasmas[i][0].Triangles.Mach, Plasmas[i][0].Triangles.Pp, Plasmas[i][0].Triangles.Ep])
						Plasmas[i][0].Triangles.VNames.extend(["TotRad","M", "Pp", "Ep"])
						iLast = iFirst
					else:
						Plasmas[i][0].Triangles.Values.extend([Plasmas[i][0].Triangles.TRad])
						Plasmas[i][0].Triangles.VNames.extend(["TotRad"])
				else:
					Plasmas[iFirst][0].Triangles.TDens += Plasmas[i][0].Triangles.Dens
					Plasmas[iFirst][0].Triangles.TDense += Plasmas[i][0].Triangles.Dens*Plasmas[i][0].charge
					Plasmas[iFirst][0].Triangles.TRad  += Plasmas[i][0].Triangles.IRad
					Plasmas[0][0].Triangles.TRad	   += Plasmas[i][0].Triangles.IRad

				if(((i > 1) and (Plasmas[i][0].charge == 1)) or (i+1 == len(iPlasmas))):
					Plasmas[iLast][0].Triangles.FracDens = np.where(Plasmas[0][0].Triangles.Dens > 0., Plasmas[iLast][0].Triangles.TDens/Plasmas[0][0].Triangles.Dens, 0.)
					Plasmas[iLast][0].Triangles.Cimp = np.where(Plasmas[1][0].Triangles.Dens > 0., Plasmas[iLast][0].Triangles.TDens/Plasmas[1][0].Triangles.Dens, 0.)																		#Fraction of impurity density over D. It's 1 for D
					Plasmas[iLast][0].Triangles.FracDense = np.where(Plasmas[0][0].Triangles.Dens > 0., Plasmas[iLast][0].Triangles.TDense/Plasmas[0][0].Triangles.Dens, 0.)
					Plasmas[iLast][0].Triangles.Zave= np.where(Plasmas[iLast][0].Triangles.TDens > 0., Plasmas[iLast][0].Triangles.TDense/Plasmas[iLast][0].Triangles.TDens, 0.)
					Plasmas[iLast][0].Triangles.Values.extend([Plasmas[iLast][0].Triangles.TDens,Plasmas[iLast][0].Triangles.TDense, Plasmas[iLast][0].Triangles.FracDens, Plasmas[iLast][0].Triangles.Cimp, Plasmas[iLast][0].Triangles.FracDense, Plasmas[iLast][0].Triangles.Zave])
					Plasmas[iLast][0].Triangles.VNames.extend(["TDens","TDense","FracDens","Cimp","FracDense","Zave"])
					iLast = iFirst

	except:
		print("\tWARNING: Error triangles parameters computation")
					
	try:
#	if(True):
		if((iPlasmas[0] == 0) and (len(iPlasmas) > 1)):				#if reading at least e and first ion compute total rad
			for k in range(nZones):
				for i in range(len(iPlasmas)):
					if(i == 0):
						Plasmas[i][k].TRad  = np.zeros_like(Plasmas[i][k].IRad)
						Plasmas[i][k].TNRad = np.zeros_like(Plasmas[i][k].IRad)
						Plasmas[i][k].Values.extend([Plasmas[i][k].TRad, Plasmas[i][k].TNRad])
						if(k == 0): Plasmas[0][k].VNames.extend(["TotRad ", "TotNRad "])
					elif(Plasmas[i][0].charge == 1):
						iFirst = i
						Plasmas[iFirst][k].TRad  = Plasmas[iFirst][k].NRad + Plasmas[iFirst][k].IRad
						Plasmas[0][k].TRad		+= Plasmas[iFirst][k].TRad
						Plasmas[0][k].TNRad		+= Plasmas[iFirst][k].NRad
						if(i == 1):
							Plasmas[1][k].Values.append(Plasmas[iFirst][k].TRad)
							if(k == 0): Plasmas[i][0].VNames.append("TotRad")
						else:
							Plasmas[i][k].Values.append(Plasmas[iFirst][k].TRad)
							if(k == 0): Plasmas[i][0].VNames.append("TotRad")
					else:
						Plasmas[iFirst][k].TRad  += Plasmas[i][k].IRad
						Plasmas[0][k].TRad		 += Plasmas[i][k].IRad

	except:
		print("\tWARNING: Error on mesh total radiation computation")

#	Set name extensions	

	for i in range(len(iPlasmas)):

#		name extensions

		if(iPlasmas[i] == 0):	NameExt = "e"
		elif(iPlasmas[i] == 1):	NameExt = "i"
		else:					NameExt = Plasmas[i][0].ion[:-1]

#		Plasma values name

		for k in range(len(Plasmas[i][0].VNames)):
			if(Plasmas[i][0].VNames[k][-1] != " "):	Plasmas[i][0].VNames[k] += NameExt
			else:										Plasmas[i][0].VNames[k]  = Plasmas[i][0].VNames[k][:-1]

		if(Triangles_OK):									#		Set triangles values name
			for k in range(len(Plasmas[i][0].Triangles.VNames)):
				if(Plasmas[i][0].Triangles.VNames[k][-1] != " "):	Plasmas[i][0].Triangles.VNames[k] += NameExt
				else:												Plasmas[i][0].Triangles.VNames[k]  = Plasmas[i][0].Triangles.VNames[k][:-1]

		if(Wall_OK):											#		Set Wall values name		
			for k in range(len(Plasmas[i][0].Wall.VNames)): 
				if(Plasmas[i][0].Wall.VNames[k][-1]!= " "): Plasmas[i][0].Wall.VNames[k] += NameExt
				else:										 Plasmas[i][0].Wall.VNames[k]  = Plasmas[i][0].Wall.VNames[k][:-1]
						
	return Plasmas

#	Convert values on knots to values on triangles
#	===========================

def tri_knots_conv(values, KnotsTri, KnotsCounts, nKnots):
	if(nKnots== 0): return values

	if(nKnots > 0):								#convert to knots position
		if(len(values) == nKnots):				#values already on knots
			return values
		else:
			out_values = np.zeros(nKnots, dtype = 'f8')
			tri_values = np.array([values, values, values]).T.reshape(-1)
			for k in range(len(KnotsTri)): out_values[KnotsTri[k]] += tri_values[k]
			return out_values/KnotsCounts
	else:										#convert to triangles positions
		if(len(values) == -nKnots):				#values on knots
			return (values[KnotsTri[:,0]] + values[KnotsTri[:,1]] + values[KnotsTri[:,2]])/3.
		else:									#values already on triangles
			return values


#	Convert values on knots or triangle to values on mesh
#	================================

def tri_mesh_conv(values, SolToTri, SolVol, TriVolOff, TriKnots,Name="None"):
	if(len(values) != TriKnots.shape[0]):				#values on knots ==> convert to triangles positions
		TriValuesOff = np.append(0,(values[TriKnots[:,0]] + values[TriKnots[:,1]] + values[TriKnots[:,2]])/3.)
	else:
		TriValuesOff = np.append(0.,values)

	if(SolToTri.size > 0):
		SolValues = np.sum(TriValuesOff[SolToTri]*TriVolOff[SolToTri], axis=2)/SolVol
	else:
		SolValues = np.zeros_like(SolVol)

	return SolValues