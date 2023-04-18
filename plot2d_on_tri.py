#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot2d_on_tri		import plot2d_on_tri
	from routines.cli_routines	import *

	if(cli_present("-h",sys.argv)):
		print("\nThis script plot plasma parameters on eirene mesh\n")
		print("plot2d_on_tri options")
		print("\t-path <value>  Directory with simulation [d='']")
		print("\t-evolution     If > 0 plot evolution file [d=0]")
		print("\t-path_label    Labels to itentify runs in plot [d='']")
		print("\t-d_only        plot only e- and D+ [d=false]")
		print("\t-log_scale     Use log scale for z colors [d=false]")
		print("\t-diff          plot difference between two evolutions[d=false]")
		print("\t-one_plot      One plot on each figure [d=false]")
		print("\t-no_plot       No plots only print [d=false]")
		print("\t-no_samexy     No same xy scale [d=false]")
		print("\t-extra_walls   Flag to show extra walls [d=false]")
		print("\t-save <value>  Save figures on files, values=none/stat/png/ps/eps/pdf [D='none']")
		print("\t-xmin <value>  To set xmin on some plot by (FigNum,PosNum,value) [D=[]]")
		print("\t-xmax <value>  To set xmax on some plot by (FigNum,PosNum,value) [D=[]]")
		print("\t-ymin <value>  To set ymin on some plot by (FigNum,PosNum,value) [D=[]]")
		print("\t-ymax <value>  To set ymax on some plot by (FigNum,PosNum,value) [D=[]]")
		print("\t-zmin <value>  To set zmin on some plot by (FigNum,PosNum,value) [D=[]]")
		print("\t-zmax <value>  To set zmax on some plot by (FigNum,PosNum,value) [D=[]]")
		print()
		exit()

	path	 	= cli_get_value("-path",			sys.argv,"")
	evolution	= cli_get_value("-evolution",	sys.argv,	[])
	save		= cli_get_value("-save",			sys.argv, "none")
	path_label 	= cli_get_value("-path_label",	sys.argv, "")
	log_scale	= cli_present("-log_scale",		sys.argv)
	d_only		= cli_present("-d_only",			sys.argv)
	diff		= cli_present("-diff",			sys.argv)
	one_plot	= cli_present("-one_plot",		sys.argv)
	no_samexy	= cli_present("-no_samexy",		sys.argv)
	no_plot		= cli_present("-no_plot",			sys.argv)
	extra_walls	= cli_present("-extra_walls",	sys.argv)
	xmin	= cli_get_value("-xmin",	sys.argv,	[])
	xmax	= cli_get_value("-xmax",	sys.argv,	[])
	ymin	= cli_get_value("-ymin",	sys.argv,	[])
	ymax	= cli_get_value("-ymax",	sys.argv,	[])
	zmin	= cli_get_value("-zmin",	sys.argv,	[])
	zmax	= cli_get_value("-zmax",	sys.argv,	[])
	MinMax = [xmin, xmax, ymin, ymax, zmin, zmax]
	plot2d_on_tri(path=path, evolution=evolution, path_label=path_label, log_scale=log_scale, d_only=d_only, diff=diff, one_plot=one_plot, no_plot=no_plot, no_samexy=no_samexy, extra_walls=extra_walls, save=save, MinMax=MinMax)
	exit()

#=======================================

# Function definition is here

import h5py
import os
import numpy as np
import matplotlib.pyplot				as pyp
import matplotlib.tri						as tri
from matplotlib.backends.backend_pdf	import PdfPages
from matplotlib.colors 					import LogNorm
from matplotlib.path 						import Path
from routines.utils_walls				import plot2d_walls
from routines.h5_routines				import h5_read
from mesh.get_rz_core_sep			import get_rz_core_sep
from files.load_refpar_file				import load_refpar_file
from files.load_plasma_files			import load_plasma_files
from files.load_soledge_mesh_file	import load_soledge_mesh_file
from files.load_ions_list				import load_ions_list
from files.save_stat					import save_stat

#==============================================================================
# This routine plots parameters on eirene mesh
#==============================================================================

def plot2d_on_tri(path="", evolution=[], path_label="", log_scale=0, d_only=0, plot_neutrals=1, diff=0, one_plot=0, no_plot=0, no_samexy=0, extra_walls=0, save="none", MinMax=[[],[],[],[],[],[]]):

	print("plot2d_on_tri")

	save_norm = False
	if(save == "stat_norm"):
		save_norm = True
		save = "stat"

	if(diff != 0):
		if((evolution == 0) or (len(evolution) != 2)):
			print("\tWith diff option two evolutions must be provided")
			exit()
		log_scale = 0

	i_plot_file = 0
	blanks		= "                          "

#	Read simulation time

	if(len(evolution) == 0): evolution = [0]
	if((len(path) > 0) and (path[-1] != "/")): path = path + "/"

	if(len(path_label) < 0): path_label = os.path.basename(os.path.abspath(path))

	RefPar = load_refpar_file(path+"Results/")				#	load reference parameters

	ions   = load_ions_list(path)								#	load ions list

	Config = load_soledge_mesh_file(path+"mesh.h5")			#	load mesh
	nZones  = len(Config.Zones)

#	load Eirene mesh

	if_tri	 = h5py.File(path+"triangles.h5", "r")

	TriKnots = h5_read(if_tri,"triangles/tri_knots")
	TriKnots = TriKnots - 1 										#Matlab/Fortan to python indexes
	R		 = h5_read(if_tri,"knots/R")*0.01
	Z		 = h5_read(if_tri,"knots/Z")*0.01
	if_tri.close()

	if(save == "export"):
		if(len(path) > 0):	export_path=path+"/Export"
		else:						export_path="Export"
		try:
			os.mkdir(export_path)
		except OSError:
			pass		
		export_path +="/"
		np.savetxt(export_path+"TriKnots.txt", TriKnots, delimiter=",", fmt="%15d", comments="")
		np.savetxt(export_path+"Knots_R.txt", R, delimiter=",", fmt="%15.7e", comments="")
		np.savetxt(export_path+"Knots_Z.txt", Z, delimiter=",", fmt="%15.7e", comments="")	#Geometry done

	v1R = R[TriKnots[:,1]] - R[TriKnots[:,0]]
	v1Z = Z[TriKnots[:,1]] - Z[TriKnots[:,0]]
	v2R = R[TriKnots[:,2]] - R[TriKnots[:,0]]
	v2Z = Z[TriKnots[:,2]] - Z[TriKnots[:,0]]
	TriVols = np.pi*np.abs(v1R*v2Z-v1Z*v2R)*np.mean(R[TriKnots],axis=1)
	AllVolume  = np.sum(TriVols)
	v1R = v1Z = v2R = v2Z = 0

	RZcore, RZsep, CoreMegazone, SepMegazone, jSep = get_rz_core_sep(Config, core_and_sep = True, use_mag_zones = False)
	SepPath = Path(RZsep, closed=True)

	Radius = 1.e-2
	if(SepPath.contains_point([RZsep[0,0],RZsep[0,1]], radius= -Radius)): Radius = -Radius				#reverse radius if needed

	RZcore  = 0
	RZsep   = 1

	RTriCenter = np.sum(R[TriKnots], axis=1)/3
	ZTriCenter = np.sum(Z[TriKnots], axis=1)/3
	InSepTri   = np.where(SepPath.contains_points(np.array([RTriCenter,ZTriCenter]).T))[0]
	OnSepTri   = InSepTri[np.where(np.logical_not(SepPath.contains_points(np.array([RTriCenter[InSepTri],ZTriCenter[InSepTri]]).T, radius = -Radius)))[0]]
	InSepVolume  = np.sum(TriVols[InSepTri])
	OutSepVolume = AllVolume - InSepVolume
	OnSepVolume  = np.sum(TriVols[OnSepTri])

	if(diff == 0):
		Plasmas = load_plasma_files(path, nZones=nZones, Evolution=evolution[0], ToKnodes = diff)

		if(d_only != 0): Plasmas = Plasmas[:2]
	else:
		Evolutions= []
		for  iv in evolution:
			Plasmas = load_plasma_files(path, nZones=nZones, Evolution=iv, ToKnodes = diff)

			if(d_only != 0): Plasmas = Plasmas[:2]
			Evolutions.append(Plasmas)

		nPlasmas = len(Plasmas)
		nZones = len(Plasmas[0])
		for iPlasma in range(nPlasmas):
			nValues =min(len(Evolutions[0][iPlasma][0].Values), len(Evolutions[1][iPlasma][0].Values))
			nTValues =min(len(Evolutions[0][iPlasma][0].Triangles.Values), len(Evolutions[1][iPlasma][0].Triangles.Values))
			for itv in range(nTValues): Plasmas[iPlasma][0].Triangles.Values[itv] -= Evolutions[0][iPlasma][0].Triangles.Values[itv]
			"""
			for k in range(nZones):
				for iv in range(nValues): Plasmas[iPlasma][k].Values[iv] -= Evolutions[0][iPlasma][k].Values[iv]
			"""
	
		Evolutions = 0
		
	tempus = Plasmas[0][0].tempus

#	prepare for plot

	TripTriang = tri.Triangulation(R, Z, triangles=TriKnots)

	fid_stat2d=0
	if(save == "pdf"):
		pdf = PdfPages("plot2d_on_tri_t={:.3f}.".format(tempus)+save)   #pdf in one file only
	elif(save == "stat"):
		StatValues = [path_label, tempus]
		StatHeader =  "             path,             time,"
		StatFormats = ['"{:>15}",','{:17.4e},']
	elif(save == "print"):
		fid_stat2d = open("stat_2d_print_t={:.4f}.txt".format(tempus),'w')
	elif(save == "print_append"):
		fid_stat2d = open("stat_2d_print.txt",'a')


	print_write(fid_stat2d, "\n####################################################################")
	print_write(fid_stat2d, "\nTime = {:.4f}.txt".format(tempus))
	Fig = []
	Ax  = []
	Im  = []
	nPlasmas = len(Plasmas)
	for iPlasma in range(nPlasmas):
		print_write(fid_stat2d, "\n####################################################################")
		print_write(fid_stat2d, Plasmas[iPlasma][0].ion)

		if(iPlasma == 0):
			Labels = [["$n\ (*10^{19}\ m^{-3})$", "$T\ (keV)$",  "$Z_{eff}$", 
						"$Rad_{tot} (kW/m^3)$", "$Rad_{neu} (kW/m^3)$", "$S_{ne}\ (*10^{19}\ m^{-3}s^{-1})$"]]
			Pars   = [["Dens", "Temp","Zeff ","TotRad ","TotNRad ","Sn"]]
			PosPars = [[True, True, True, True, True,True]]
			NrowNcol = [[3,2]]
			Facts		  	  = [[1e-19,1.e-3, 1., 1e-3,1.e-3,1.e-19]]
			FactsUnits		  = [["10^19","10^3", "1.","10^3","10^3","10^19"]]
			MaxRange	      = [[1.e5,0., 0., 1.e5, 1.e5,0.]]
			PrintAllVolSumAve = [[1, 0, 0, 1, 1, 1]]								#1 for sum, 2 for ave, 3 for sum and ave
			PrintInSepSumAve  = [[3, 0, 2, 1, 1, 0]]
			PrintOutSepSumAve = [[0, 0, 0, 1, 1, 1]]
			PrintOnSepSumAve  = [[2, 2, 2, 0, 0, 1]]

		elif(iPlasma == 1):
			if(plot_neutrals):
				Labels = [["$n\ (*10^{19}\ m^{-3})$", "$T\ (keV)$",  "$M$", "$Rad_{tot} (kW/m^3)$", "$Rad_{neu} (kW/m^3)$", "$Rad_{ion} (kW/m^3)$"], \
						  ["$n_{n}\ (*10^{19}\ m^{-3})$","$n_{m}\ (*10^{19}\ m^{-3}$","$P_n\ (Pa)$", "$T_n\ (K)$", "$T_m\ (K)$", "$n_H/n_e$"]]
				Pars   = [["Dens", "Temp", "M", "TotRad", "NRad", "IRad"], \
						  ["Nn",  "Nm", "Pn", "Tn", "Tm","FracDens"]]
				PosPars = [[True, True, False, True, True, True], \
						  [True,    True, True, True, True, False]]
				NrowNcol = [[3,2],[3,2]]
				Facts  	 		  = [[1e-19, 1.e-3, 1., 1e-3,1.e-3,1.e-3],[1.e-19, 1.e-19,1, 1., 1., 1.]]
				FactsUnits 		  = [["10^19", "10^3", "1.", "10^3","10^3","10^3"],["10^19", "10^19","1.","1.","1.", "1."]]
				MaxRange		  = [[1.e5, 0., 0., 1.e5, 1.e5, 1.e5],[0., 0.,1.e3, 0.,0.,0.]]
				PrintAllVolSumAve = [[1, 0, 0, 1, 1, 1], [1, 1, 0, 0, 1, 0]]
				PrintInSepSumAve  = [[1, 0, 0, 1, 1, 1], [0, 0, 0, 0, 0, 2]]
				PrintOutSepSumAve = [[0, 0, 0, 1, 1, 1], [0, 0, 0, 0, 0, 0]]
				PrintOnSepSumAve  = [[2, 2, 0, 0, 0, 0], [0, 0, 0, 0, 0, 2]]
			else:
				Labels = [["$n\ (*10^{19}\ m^{-3})$", "$T\ (keV)$",  "$M$", "$Rad_{tot} (kW/m^3)$", "$Rad_{neu} (kW/m^3)$", "$Rad_{ion} (kW/m^3)$"], \
						  ["$S_{ni}\ (*10^{19}\ m^{-3}s^{-1})$","$S_{Ei}\ (MW/m^{3})$","$P_n\ (Pa)$", "$P_p\ (P)$", "$E_p\ (MJ/m^3)$", "$n_H/n_e$"]]
				Pars   = [["Dens", "Temp", "M", "TotRad", "NRad", "IRad"], \
						  ["Sn",    "SE", "Pn", "Pp", "Ep","FracDens"]]
				PosPars = [[True, True, False, True, True, True], \
						  [True,    True, True, True, True, False]]
				NrowNcol = [[3,2],[3,2]]
				Facts  	 		  = [[1e-19, 1.e-3, 1., 1e-3,1.e-3,1.e-3],[1.e-19, 1.e-6,1., 1.,1.e-6, 1.]]
				FactsUnits 		  = [["10^19", "10^3", "1.", "10^3","10^3","10^3"],["10^19", "10^6","1.", "1.","10^6", "1."]]
				MaxRange		  = [[1.e5, 0., 0., 1.e5, 1.e5, 1.e5],[0., 0.,1.e3, 0.,0.,0.]]
				PrintAllVolSumAve = [[1, 0, 0, 1, 1, 1], [0, 0, 0, 0, 1, 0]]
				PrintInSepSumAve  = [[3, 0, 0, 1, 1, 1], [0, 0, 0, 0, 0, 2]]
				PrintOutSepSumAve = [[0, 0, 0, 1, 1, 1], [0, 0, 0, 0, 0, 0]]
				PrintOnSepSumAve  = [[2, 2, 0, 0, 0, 0], [0, 0, 0, 0, 0, 2]]

		elif(Plasmas[iPlasma][0].charge == 1):
			if(plot_neutrals):
				Labels = [["$n\ (*10^{19}\ m^{-3})$", "$n\ (*10^{19}\ m^{-3})$",  "$n_{Tot}\ (*10^{19}\ m^{-3})$", "$Rad_{tot} (kW/m^3)$", "$Rad_{neu} (kW/m^3)$", "$Rad_{ion} (kW/m^3)$","$n_{tot}/n_e$", "$n_{Itot}/n_D$"]]
				Pars   = [["Nn", "Dens", "TDens", "TotRad", "NRad", "IRad","FracDens","Cimp"]]
				PosPars = [[True, True, True, True, True, True, True, True]]
				NrowNcol = [[4,2]]
				Facts  	 = [[1e-19,1e-19, 1.e-19, 1e-3,1.e-3,1.e-3, 1., 1.]]
				FactsUnits = [["10^19","10^3","10^19","10^3","10^3","10^3","1.","1."]]
				MaxRange	      = [[1.e3,1.e3, 1.e3, 1.e5, 1.e5,1.e5, 0.,0.]]
				PrintAllVolSumAve = [[1, 0, 1, 1, 1, 1, 0, 0]]
				PrintInSepSumAve  = [[2, 0, 2, 1, 1, 1, 2, 2]]
				PrintOutSepSumAve = [[0, 0, 0, 1, 1, 1, 0,0]]
				PrintOnSepSumAve  = [[0, 0, 0, 0, 0, 0, 2, 2]]
			else:
				Labels = [["$n\ (*10^{19}\ m^{-3})$", "$T\ (keV)$",  "$Sn\ (*10^{19}\ m^{-3}s^{-1})$", "$Rad_{tot} (kW/m^3)$", "$Rad_{neu} (kW/m^3)$", "$Rad_{ion} (kW/m^3)$","$Tot-n_Z$", "$n_{Itot}/n_D$"]]
				Pars   = [["Dens", "Temp", "Sn", "TotRad", "NRad", "IRad","FracDens","Cimp"]]
				PosPars = [[True, True, True, True, True, True, True, True]]
				NrowNcol = [[4,2]]
				Facts  	 = [[1e-19,1.e-3, 1.e-19, 1e-3,1.e-3,1.e-3,1.,1.]]
				FactsUnits = [["10^19","10^3","10^19","10^3","10^3","10^3","1.","1."]]
				MaxRange	      = [[1.e5,0., 0., 1.e5, 1.e5,1.e5, 0.,0.]]
				PrintAllVolSumAve = [[1, 0, 0, 1, 1, 1, 0, 0]]
				PrintInSepSumAve  = [[2, 0, 0, 1, 1, 1, 2, 2]]
				PrintOutSepSumAve = [[0, 0, 0, 1, 1, 1, 0, 2]]
				PrintOnSepSumAve  = [[0, 0, 0, 0, 0, 0, 2, 2]]

		else:
			Labels   = [["$n\ (*10^{19}\ m^{-3})$", "$T\ (keV)$",  "$Rad_{ion} (kW/m^3)$"]]
			Pars     = [["Dens", "Temp", "IRad"]]
			PosPars = [[True, True, True]]
			NrowNcol = [[2,2]]
			Facts  	 = [[1e-19,1.e-3, 1e-3]]
			FactsUnits = [["10^19","10^3","10^3"]]
			MaxRange	      = [[1.e5,0., 0., 1.e5]]
			PrintAllVolSumAve = [[1, 0, 1]]
			PrintInSepSumAve  = [[0, 0, 1]]
			PrintOutSepSumAve = [[0, 0, 1]]
			PrintOnSepSumAve  = [[0, 0, 0]]

		if(iPlasma == 0):	NameExt = "e"
		elif(iPlasma == 1):	NameExt = "i"
		else:				NameExt = Plasmas[iPlasma][0].ion[:-1]
		for iF in range(len(Pars)):
			for i in range(len(Pars[iF])): 
				if(Pars[iF][i][-1]!= " "): Pars[iF][i] += NameExt
				else:						Pars[iF][i]  = Pars[iF][i][:-1]

		if(no_plot == 0):
			FigNum = []
			PosNum= []
			FigPl = []
			AxPl = []
			if(one_plot != 1):
				for iF in range(len(Labels)):
					AxPl.append([])
					FigNum.append([])
					FigPl.append([])
					PosNum.append([])
					Fig.append(pyp.figure())
					Fig[-1].patch.set_facecolor('white')
					for i in range(len(Labels[iF])): 
						if((len(Ax) == 0) or (no_samexy != 0)):	Ax.append(Fig[-1].add_subplot(NrowNcol[iF][1],NrowNcol[iF][0],i+1))
						else:									Ax.append(Fig[-1].add_subplot(NrowNcol[iF][1],NrowNcol[iF][0],i+1, sharex = Ax[0], sharey = Ax[0]))
						AxPl[-1].append(Ax[-1])
						FigNum[-1].append(len(Fig))
						PosNum[-1].append(i+1)
						FigPl[-1].append(Fig[-1])

					Fig[-1].tight_layout(pad=2., w_pad=3., h_pad=3.)
			else:
				for iF in range(len(Labels)):
					AxPl.append([])
					FigPl.append([])
					FigNum.append([])
					PosNum.append([])
					for i in range(len(Labels[iF])):
						Fig.append(pyp.figure())
						Fig[-1].patch.set_facecolor('white')
						if((len(Ax) == 0) or (no_samexy != 0)):	Ax.append(Fig[-1].add_subplot(111))
						else:									Ax.append(Fig[-1].add_subplot(111, sharex = Ax[0], sharey = Ax[0]))

						AxPl[-1].append(Ax[-1])
						FigNum[-1].append(len(Fig))
						PosNum[-1].append(1)
						FigPl[-1].append(Fig[-1])
		
			for iF in range(len(Labels)):
				for i in range(len(Labels[iF])):
					AxPl[iF][i].set_aspect(1.)
					AxPl[iF][i].autoscale(enable=True, axis='both', tight=True)
					AxPl[iF][i].set_xlabel("$R\ (m)$")
					AxPl[iF][i].set_ylabel("$Z\ (m)$")
	#				if(i == 0): AxPl[iF][i].set_title(path_label+" @ t={:.3f} s".format(tempus))
					if(i == 0): AxPl[iF][i].set_title(Plasmas[iPlasma][0].ion)
					if(i == 1): AxPl[iF][i].set_title(path_label)
					if(i == 2): 
						if(diff == 0): AxPl[iF][i].set_title("Evol.={:d}".format(evolution[0]))
						else:			AxPl[iF][i].set_title("Ev_{:d}-Ev_{:d}".format(evolution[1],evolution[0]))


	#	Plot parameters

		for iF in range(len(Labels)):
			for i in range(len(Labels[iF])):
#				try:
					iPar = Plasmas[iPlasma][0].Triangles.VNames.index(Pars[iF][i])	

#					=========== Print section =============

#					Print total and average on whole plasma

					if(diff == 1):	TValues =  (Plasmas[iPlasma][0].Triangles.Values[iPar][TriKnots[:,0]] + Plasmas[iPlasma][0].Triangles.Values[iPar][TriKnots[:,1]] + Plasmas[iPlasma][0].Triangles.Values[iPar][TriKnots[:,2]])/3.
					else:				TValues =  Plasmas[iPlasma][0].Triangles.Values[iPar]

					if(PrintAllVolSumAve[iF][i] != 0):
						SumPar = np.sum(TriVols*TValues)
						if((PrintAllVolSumAve[iF][i] == 1) or (PrintAllVolSumAve[iF][i] == 1)):					
							Value   = SumPar
							PreName = "\tTotal        "
								
						if(PrintAllVolSumAve[iF][i] > 1):
							Value   = SumPar/AllVolume
							PreName = "\tAverage      "

						print_write(fid_stat2d, PreName+Pars[iF][i]+blanks[:15-len(Pars[iF][i])]+"= {:10.3e}".format(Value))
						if(save == "stat"):
							StatHeader += PreName[1:]+Pars[iF][i]+","
							StatValues.append(Value)
							StatFormats.append('{:17.4e},')
							if(save_norm):
								StatValues[-1] *=Facts[iF][i]
								StatHeader = StatHeader[:-1] + " (*" + FactsUnits[iF][i]+ "),"

#					Print total and average on plasma inside separatrix

					if(PrintInSepSumAve[iF][i] != 0):
						SumPar = np.sum(TriVols[InSepTri]*TValues[InSepTri])
						if((PrintInSepSumAve[iF][i] == 1) or (PrintInSepSumAve[iF][i] == 3)):
							Value   = SumPar
							PreName = "\tIn Sep. Tot  "

						if(PrintInSepSumAve[iF][i] > 1):
							Value   = SumPar/InSepVolume
							PreName = "\tIn Sep. Ave  "

						print_write(fid_stat2d, PreName+Pars[iF][i]+blanks[:15-len(Pars[iF][i])]+"= {:10.3e}".format(Value))

						if(save == "stat"):
							StatHeader += PreName[1:]+Pars[iF][i]+","
							StatValues.append(Value)
							StatFormats.append('{:17.4e},')
							if(save_norm):
								StatValues[-1] *=Facts[iF][i]
								StatHeader = StatHeader[:-1] + " (*" + FactsUnits[iF][i]+ "),"

#					Print total and average on plasma outside separatrix

					if(PrintOutSepSumAve[iF][i] != 0):
						SumPar = np.sum(TriVols*TValues) - np.sum(TriVols[InSepTri]*TValues[InSepTri])
						if((PrintOutSepSumAve[iF][i] == 1) or (PrintOutSepSumAve[iF][i] == 3)):
							Value   = SumPar
							PreName = "\tOut Sep. Tot "
						if(PrintOutSepSumAve[iF][i] > 1):
							Value   = SumPar/OutSepVolume
							PreName = "\tOut Sep. Ave "
						
						print_write(fid_stat2d, PreName+Pars[iF][i]+blanks[:15-len(Pars[iF][i])]+"= {:10.3e}".format(Value))
						if(save == "stat"):
							StatHeader += PreName[1:]+Pars[iF][i]+","
							StatValues.append(Value)
							StatFormats.append('{:17.4e},')
							if(save_norm):
								StatValues[-1] *=Facts[iF][i]
								StatHeader = StatHeader[:-1] + " (*" + FactsUnits[iF][i]+ "),"

#					Print total and average on plasma close to separatrix (but inside)

					if(PrintOnSepSumAve[iF][i] != 0):
						SumPar = np.sum(TriVols[OnSepTri]*TValues[OnSepTri])
						if((PrintOnSepSumAve[iF][i] == 1) or (PrintOnSepSumAve[iF][i] == 3)):
							Value   = SumPar
							PreName = "\tOn Sep. Tot  "
						if(PrintOnSepSumAve[iF][i] > 1):
							Value   = SumPar/OnSepVolume
							PreName = "\tOn Sep. Ave  "

						print_write(fid_stat2d, PreName+Pars[iF][i]+blanks[:15-len(Pars[iF][i])]+"= {:10.3e}".format(SumPar/OnSepVolume))
						if(save == "stat"):
							StatHeader += PreName[1:]+Pars[iF][i]+","
							StatValues.append(Value)
							StatFormats.append('{:17.4e},')
							if(save_norm):
								StatValues[-1] *=Facts[iF][i]
								StatHeader = StatHeader[:-1] + " (*" + FactsUnits[iF][i]+ "),"

#					=========== Plot section =============

					if(no_plot == 0):
						Vmin = Plasmas[iPlasma][0].Triangles.Values[iPar].min()
						Vmax = Plasmas[iPlasma][0].Triangles.Values[iPar].max()
						if((MaxRange[iF][i] != 0.) and (diff == 0)):
							MinValue = Vmax/MaxRange[iF][i]
							Plasmas[iPlasma][0].Triangles.Values[iPar] = np.where(Plasmas[iPlasma][0].Triangles.Values[iPar] > MinValue, Plasmas[iPlasma][0].Triangles.Values[iPar], MinValue)
						else:
							MinValue = Vmin

						MinSet,MinValue = search_min_max(FigNum[iF][i], PosNum[iF][i], MinMax[4])
						MaxSet,MaxValue = search_min_max( FigNum[iF][i], PosNum[iF][i], MinMax[5])
						if(not MinSet): MinValue = MinValue
						if(not MaxSet): MaxValue = Vmax

						if(diff == 1):
							Im.append(AxPl[iF][i].tricontourf(R, Z, TriKnots, Plasmas[iPlasma][0].Triangles.Values[iPar]*Facts[iF][i], cmap="jet"))
							if(Vmax > Vmin): AxPl[iF][i].tricontour(TripTriang, Plasmas[iPlasma][0].Triangles.Values[iPar]*Facts[iF][i], [0.], vmin=MinValue*Facts[iF][i], vmax=MaxValue*Facts[iF][i], colors='k')
						else:
							if((log_scale == 0) or (not PosPars[iF][i])):
								Im.append(AxPl[iF][i].tripcolor(TripTriang, Plasmas[iPlasma][0].Triangles.Values[iPar]*Facts[iF][i], vmin=MinValue*Facts[iF][i], vmax=MaxValue*Facts[iF][i]))
							else:
								if(MinValue > 0.):
									Im.append(AxPl[iF][i].tripcolor(TripTriang, Plasmas[iPlasma][0].Triangles.Values[iPar]*Facts[iF][i],norm=LogNorm(vmin=MinValue*Facts[iF][i], vmax=MaxValue*Facts[iF][i])))
								else:
									iPos = np.where(Plasmas[iPlasma][0].Triangles.Values[iPar] > 0.)[0]
									if(len(iPos) > 0):
										MinValue = np.min(Plasmas[iPlasma][0].Triangles.Values[iPar][iPos])
										VPos = np.where(Plasmas[iPlasma][0].Triangles.Values[iPar] > 0., Plasmas[iPlasma][0].Triangles.Values[iPar], MinValue)

										Im.append(AxPl[iF][i].tripcolor(TripTriang, VPos*Facts[iF][i],norm=LogNorm(vmin=MinValue*Facts[iF][i], vmax=MaxValue*Facts[iF][i])))
									else:
										Im.append(AxPl[iF][i].tripcolor(TripTriang, Plasmas[iPlasma][0].Triangles.Values[iPar]*Facts[iF][i]))

						plot2d_walls(AxPl[iF][i], Config.Walls, extra_wall=extra_walls)


						Xlims = AxPl[iF][i].get_xlim()
						Ylims = AxPl[iF][i].get_ylim()
						xMinSet,xMinValue = search_min_max(FigNum[iF][i], PosNum[iF][i], MinMax[0], Xlims[0])
						xMaxSet,xMaxValue = search_min_max( FigNum[iF][i], PosNum[iF][i], MinMax[1], Xlims[1])
						yMinSet,yMinValue = search_min_max(FigNum[iF][i], PosNum[iF][i], MinMax[2], Ylims[0])
						yMaxSet,yMaxValue = search_min_max( FigNum[iF][i], PosNum[iF][i], MinMax[3], Ylims[1])
						if(xMinSet or xMaxSet or yMinSet or yMaxSet):
							AxPl[iF][i].autoscale(enable=False, tight=True)
							AxPl[iF][i].set_xlim(xMinValue, xMaxValue)
							AxPl[iF][i].set_ylim(yMinValue, yMaxValue)

						cb = FigPl[iF][i].colorbar(Im[-1], ax=AxPl[iF][i])
						cb.set_label(Labels[iF][i])

#		=========== Export section =============

		if(save == "export"):
			
			ParExt="_"+NameExt
			if(iPlasma == 0):
				ExpPars   = ["Dens", "Temp"]
				ExpNames   = ["Dens", "Temp"]
			elif(iPlasma == 1):
				ExpPars   = ["Dens", "Temp","Nn", "Tn","vxn", "vyn", "vzn","Nm", "Tm",  "Pn",  "vxm", "vym", "vzm"]
				ExpNames   = ["Dens", "Temp","Nn", "Tn","vxn", "vyn", "vzn","Nm", "Tm",  "Pn",  "vxm", "vym", "vzm"]
			elif((iPlasma > 1) and (Plasmas[iPlasma][0].charge == 1)):
				ExpPars   = ["Nn", "Tn","vxn", "vyn", "vzn"]
				ExpNames   = ["Nn", "Tn","vxn", "vyn", "vzn"]
				ParExt=ParExt[:-1]
			else:
				ExpPars=[]

			for kPar in range(len(ExpPars)):			
				try:
						ePar = Plasmas[iPlasma][0].Triangles.VNames.index(ExpPars[kPar]+NameExt)
				except:
					print("WARNING: Not found ",ExpPars[kPar]+NameExt," for iPlasma = ", iPlasma)
					print("\tAvailable names=", Plasmas[iPlasma][0].Triangles.VNames)
					exit()
			
				np.savetxt(export_path+"/"+ExpNames[kPar]+ParExt+".txt", Plasmas[iPlasma][0].Triangles.Values[ePar] , delimiter=",", fmt="%15.7e", comments="")

	if(save != "none"):
		if(save == "stat"):
			save_stat("stat2d.csv", StatHeader, StatValues, StatFormats)
		elif((save == "print") or (save == "print_append")):
			fid_stat2d.close()
		elif((save == "export") or (no_plot != 0)):
			pass
		else:
			for i in range(len(Fig)):
				i_plot_file += 1
				if(one_plot != 1): Fig[i].set_size_inches(10.05,7.44)
				if(save == "pdf"):
					print("save pdf page=",i)
					pdf.savefig(Fig[i])
				else:
					Fig[i].savefig("plot2d_on_tri_t={:.3f}_{:d}.".format(tempus,i_plot_file)+save)
#		pyp.show(block=False)
#		pyp.close(0)
	elif(no_plot == 0):
		pyp.show()

	if(save == "pdf"):	pdf.close()

	print("plot2d_on_tri: Completed")

	return

def print_write(fid, string):
	print(string)
	if(fid != 0): fid.write(string+"\n")
	return
#
# This routine set to the minimum positive value all zero or negative values
#

def set_min_positive(Values):
	index_n = np.where(Values <= 0.)
	if(len(index_n[0]) > 0):
		index_p = np.where(Values > 0.)
		min_p   = np.min(Values[index_p])
		Values[index_n] = min_p

	return Values

def  search_min_max(FigNum, PosNum, MinMax, DefVal=None):
	nMinMax = int(len(MinMax)/3)
	if((nMinMax == 0) or (nMinMax*3 != len(MinMax))): return False, 0.

	for i in range(nMinMax):
		if((MinMax[3*i] == FigNum) and (MinMax[3*i+1] == PosNum) ): 
			
			return True, MinMax[3*i+2]

	if(DefVal != None): return True, DefVal
	else:				return False, 0.