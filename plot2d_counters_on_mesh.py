#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot2d_counters_on_mesh import plot2d_counters_on_mesh
	from routines.cli_routines import *

	if(cli_present("-h",sys.argv)):
		print("plot2d_counters_on_mesh options")
		print("-path        Directory with simulation [d='']")
		print("\t-evolution If > 0 plot evolution file [d=0]")
		print("-save        Save figures on files, values=none/png/ps/eps/pdf [D='none']")
		print("-no_shade    no shade 2d plot [d=false]")
		print("\t-no_mask   no masked out of wall [d=false]")
#		print("\t-log_scale Use log scale for z colors [d=false]")
		print("\t-d_only    plot only e- and D+ [d=false]")
		print("-one_plot    All figure on one plot [d=false]")
		exit()

	path		= cli_get_value("-path",		sys.argv,"")
	evolution	= cli_get_value("-evolution",	sys.argv,	[])
	save		= cli_get_value("-save",		sys.argv, "none")
	no_shade	= cli_present("-no_shade",	sys.argv)
	no_mask	= cli_present("-no_mask",	sys.argv)
	log_scale	= cli_present("-log_scale",		sys.argv)
	d_only		= cli_present("-d_only",			sys.argv)
	one_plot	= cli_present("-one_plot",	sys.argv)
	plot2d_counters_on_mesh(path=path, evolution=evolution, log_scale=log_scale, d_only=d_only, one_plot=one_plot, no_mask=no_mask, no_shade=no_shade, save=save)
	exit()

#=======================================

# Function definition is here

import h5py
import os
import numpy							as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot 				as pyp
from matplotlib.backends.backend_pdf	import PdfPages
from files.load_soledge_mesh_file		import load_soledge_mesh_file
from files.load_plasma_files			import load_plasma_files
from files.load_refpar_file				import load_refpar_file
from files.load_ions_list				import load_ions_list
from files.load_input_file				import load_input_file
from routines.utils_routines			import extend_mat1
from routines.utils_walls				import plot2d_walls
from routines.globals					import DEBUG, KB

#=========================================================
# This routine to plot M, ne/ni and Te/Ti results
#=========================================================

def plot2d_counters_on_mesh(path="", evolution=[], path_label="", log_scale=0, d_only=0, one_plot=0, no_shade=0, no_mask=0, diff=0, save="none"):

	print("plot2d_counters_on_mesh")

	shading = 'gouraud'
	if(no_shade == 1): shading = 'flat'

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
	input_params   = load_input_file(path)
	N_iterations = input_params.global_parameters.N_iterations
	N_save = input_params.global_parameters.N_save

	ions   = load_ions_list(path)								#	load ions list
	if(d_only != 0): ions = ions[0:2]
	iPlasmas = [i for i in range(len(ions))]

	Config = load_soledge_mesh_file(path+"mesh.h5")			#	load mesh
	Zones   = Config.Zones
	nZones  = len(Zones)

	if(diff == 0):
		Plasmas = load_plasma_files(path, nZones=nZones, Evolution=evolution[0], ToKnodes = 0, iPlasmas=iPlasmas)
	else:
		Evolutions= []
		for  iv in evolution:
			Plasmas = load_plasma_files(path, nZones=nZones, Evolution=iv, ToKnodes = 0, iPlasmas=iPlasmas)
			Evolutions.append(Plasmas)

		nPlasmas = len(Plasmas)
		nZones = len(Plasmas[0])
		for iPlasma in range(nPlasmas):
			nValues =min(len(Evolutions[0][iPlasma][0].Values), len(Evolutions[1][iPlasma][0].Values))
			nTValues =min(len(Evolutions[0][iPlasma][0].Triangles.Values), len(Evolutions[1][iPlasma][0].Triangles.Values))
			for itv in range(nTValues): Plasmas[iPlasma][0].Triangles.Values[itv] -= Evolutions[0][iPlasma][0].Triangles.Values[itv]

			for k in range(nZones):
				for iv in range(nValues): Plasmas[iPlasma][k].Values[iv] -= Evolutions[0][iPlasma][k].Values[iv]
				
		Evolutions = 0
		
	tempus = Plasmas[0][0].tempus

#	prepare for plot

	if(save == "pdf"):	pdf = PdfPages("plot2d_counters_on_mesh.py_t={:.3f}.".format(tempus)+save)   #pdf in one file only
	if(save == "stat"):
		StatValues = [path_label, tempus]
		StatHeader =  "             path,             time,"
		StatFormats = ['"{:>15}",','{:17.4e},']

	Fig = []
	Ax  = []
	Im  = []
	nPlasmas = len(Plasmas)
	for iPlasma in range(nPlasmas):
		print("\n####################################################################")
		print(Plasmas[iPlasma][0].ion)

		if(iPlasma < 2):
			Labels = [["Zhdanov_min_n", "clean_min_n",  "clean_min_t", "coupling_min_t"]]
			Pars   = [["count_Zhdanov_min_n",  "count_clean_min_n", "count_clean_min_t", "count_coupling_min_t"]]
			Labels = [["$c_{Z_min_n}/n_{iter}$", "$c_{min_n}/n_{iter}$",  "$c_{min_T}/n_{iter}$", "$c_{c_min_n}/n_{iter}$"]]
			PosPars = [[True, True, True, True]]
			NrowNcol = [[2,2]]
			Facts		  	  = [[1.,1., 1.,1.]]
			FactsUnits		  = [["1.","1.", "1."]]
			MaxRange	      = [[0., 0., 0., 0.]]
			PrintAllVolSumAve = [[0, 0, 0, 0]]								#1 for sum, 2 for ave, 3 for sum and ave
			PrintInSepSumAve  = [[0, 0, 0, 0]]
			PrintOutSepSumAve = [[0, 0, 0, 0]]
			PrintOnSepSumAve  = [[0, 0, 0, 0]]
			
		else:
			Labels = [["Zhdanov_min_n", "clean_min_n",  "clean_min_t", "coupling_min_t", "clean_t_excurs"]]
			Pars   = [["count_Zhdanov_min_n",  "count_clean_min_n", "count_clean_min_t", "count_coupling_min_t","count_clean_t_excurs"]]
			Labels = [["$c_{Z_min_n}/n_{iter}$", "$c_{min_n}/n_{iter}$",  "$c_{min_T}/n_{iter}$", "$c_{c_min_n}/n_{iter}$", "$c_{T_exurs}/n_{iter}$"]]
			PosPars = [[True, True, False, True, True, True], \
					  [True,    True, True, True, True, False]]
			NrowNcol = [[3,2]]
			Facts  	 		  = [[1.,1., 1.,1.,1.]]
			FactsUnits 		  =  [["1.","1.", "1.","1."]]
			MaxRange		  = [[0., 0., 0., 0.,0.]]
			PrintAllVolSumAve = [[0, 0, 0, 0,0]]
			PrintInSepSumAve  = [[0, 0, 0, 0,0]]
			PrintOutSepSumAve = [[0, 0, 0, 0,0]]
			PrintOnSepSumAve  = [[0, 0, 0, 0,0]]

		if(iPlasma == 0):	NameExt = "e"
		elif(iPlasma == 1):	NameExt = "i"
		else:				NameExt = Plasmas[iPlasma][0].ion[:-1]
		for iF in range(len(Pars)):
			for i in range(len(Pars[iF])): 
				if(Pars[iF][i][-1]!= " "): Pars[iF][i] += NameExt
				else:						Pars[iF][i]  = Pars[iF][i][:-1]

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
					Ax.append(Fig[-1].add_subplot(NrowNcol[iF][1],NrowNcol[iF][0],i+1))
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
					Ax.append(Fig[-1].add_subplot(111))

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
				if(i == 0): AxPl[iF][i].set_title(path_label)
				if(diff == 0):
					if(evolution[0] == 0):
						if(i == 1): AxPl[iF][i].set_title(Plasmas[iPlasma][0].ion)
					else:
						if(i == 1): AxPl[iF][i].set_title("Evol.={:d}".format(evolution[0]))
						if(i == 2): AxPl[iF][i].set_title(Plasmas[iPlasma][0].ion)
				else:
					if(i == 1): AxPl[iF][i].set_title("Ev_{:d}-Ev_{:d}".format(evolution[1],evolution[0]))
					if(i == 2): AxPl[iF][i].set_title(Plasmas[iPlasma][0].ion)


	#	Plot parameters

		nZones  = len(Plasmas[iPlasma])
		AVmax  = np.zeros((nZones), dtype='f8')
		AVmin  = np.zeros((nZones), dtype='f8')
		for iF in range(len(Labels)):
			for i in range(len(Labels[iF])):
#				try:
					iPar = Plasmas[iPlasma][0].VNames.index(Pars[iF][i])
					for k in range(nZones): Plasmas[iPlasma][k].Values[iPar] = Plasmas[iPlasma][k].Values[iPar]/N_iterations			#Normalize

#					=========== Print section =============

#					Print total and average on whole plasma

					if(PrintAllVolSumAve[iF][i] != 0):
						SumPar = np.sum(TriVols*Plasmas[iPlasma][0].Values[iPar])
						if((PrintAllVolSumAve[iF][i] == 1) or (PrintAllVolSumAve[iF][i] == 1)):					
							Value   = SumPar
							PreName = "\tTotal        "
								
						if(PrintAllVolSumAve[iF][i] > 1):
							Value   = SumPar/AllVolume
							PreName = "\tAverage      "

						print(PreName, Pars[iF][i],blanks[:15-len(Pars[iF][i])],"= {:.3e}".format(Value))
						if(save == "stat"):
							StatHeader += PreName[1:]+Pars[iF][i]+","
							StatValues.append(Value)
							StatFormats.append('{:17.4e},')
							if(save_norm):
								StatValues[-1] *=Facts[iF][i]
								StatHeader = StatHeader[:-1] + " (*" + FactsUnits[iF][i]+ "),"

#					Print total and average on plasma inside separatrix

					if(PrintInSepSumAve[iF][i] != 0):
						SumPar = np.sum(TriVols[InSepTri]*Plasmas[iPlasma][0].Triangles.Values[iPar][InSepTri])
						if((PrintInSepSumAve[iF][i] == 1) or (PrintInSepSumAve[iF][i] == 3)):
							Value   = SumPar
							PreName = "\tIn Sep. Tot  "

						if(PrintInSepSumAve[iF][i] > 1):
							Value   = SumPar/InSepVolume
							PreName = "\tIn Sep. Ave  "

						print(PreName,Pars[iF][i],blanks[:15-len(Pars[iF][i])],"= {:.3e}".format(Value))

						if(save == "stat"):
							StatHeader += PreName[1:]+Pars[iF][i]+","
							StatValues.append(Value)
							StatFormats.append('{:17.4e},')
							if(save_norm):
								StatValues[-1] *=Facts[iF][i]
								StatHeader = StatHeader[:-1] + " (*" + FactsUnits[iF][i]+ "),"

#					Print total and average on plasma outside separatrix

					if(PrintOutSepSumAve[iF][i] != 0):
						SumPar = np.sum(TriVols*Plasmas[iPlasma][0].Triangles.Values[iPar]) - np.sum(TriVols[InSepTri]*Plasmas[iPlasma][0].Triangles.Values[iPar][InSepTri])
						if((PrintOutSepSumAve[iF][i] == 1) or (PrintOutSepSumAve[iF][i] == 3)):
							Value   = SumPar
							PreName = "\tOut Sep. Tot "
						if(PrintOutSepSumAve[iF][i] > 1):
							Value   = SumPar/OutSepVolume
							PreName = "\tOut Sep. Ave "
						
						print(PreName,Pars[iF][i],blanks[:15-len(Pars[iF][i])],"= {:.3e}".format(Value))
						if(save == "stat"):
							StatHeader += PreName[1:]+Pars[iF][i]+","
							StatValues.append(Value)
							StatFormats.append('{:17.4e},')
							if(save_norm):
								StatValues[-1] *=Facts[iF][i]
								StatHeader = StatHeader[:-1] + " (*" + FactsUnits[iF][i]+ "),"

#					Print total and average on plasma close to separatrix (but inside)

					if(PrintOnSepSumAve[iF][i] != 0):
						SumPar = np.sum(TriVols[OnSepTri]*Plasmas[iPlasma][0].Triangles.Values[iPar][OnSepTri])
						if((PrintOnSepSumAve[iF][i] == 1) or (PrintOnSepSumAve[iF][i] == 3)):
							Value   = SumPar
							PreName = "\tOn Sep. Tot  "
						if(PrintOnSepSumAve[iF][i] > 1):
							Value   = SumPar/OnSepVolume
							PreName = "\tOn Sep. Ave  "

						print(PreName,Pars[iF][i],blanks[:15-len(Pars[iF][i])],"= {:.3e}".format(SumPar/OnSepVolume))
						if(save == "stat"):
							StatHeader += PreName[1:]+Pars[iF][i]+","
							StatValues.append(Value)
							StatFormats.append('{:17.4e},')
							if(save_norm):
								StatValues[-1] *=Facts[iF][i]
								StatHeader = StatHeader[:-1] + " (*" + FactsUnits[iF][i]+ "),"

#					=========== Plot section =============

					
					for k in range(nZones): AVmax[k] = Plasmas[iPlasma][k].Values[iPar].max()
					Vmax = AVmax.max()
					if((MaxRange[iF][i] != 0.) and (diff == 0)):
						Vmin = Vmax/MaxRange[iF][i]
						for k in range(nZones):  Plasmas[iPlasma][k].Values[iPar] = np.where(Plasmas[iPlasma][k].Values[iPar] > Vmin, Plasmas[iPlasma][k].Values[iPar], Vmin)
					else:
						for k in range(nZones): AVmin[k] = Plasmas[iPlasma][k].Values[iPar].min()
						Vmin = AVmin.min()

					for k in range(nZones):
						if(no_shade == 0):
							if(no_mask == 0):	MaskedValues = np.ma.masked_where(extend_mat1(Zones[k].Chi) == 1., Plasmas[iPlasma][k].Values[iPar][1:,1:]*Facts[iF][i])
							else:				MaskedValues = Plasmas[iPlasma][k].Values[iPar][1:,1:]*Facts[iF][i]
						else:
							if(no_mask == 0):	MaskedValues = np.ma.masked_where(extend_mat1(Zones[k].Chi) == 1., Plasmas[iPlasma][k].Values[iPar][1:,1:]*Facts[iF][i])
							else:				MaskedValues = Plasmas[iPlasma][k].Values[iPar][1:,1:]*Facts[iF][i]

						if(k == 0):	Im.append(AxPl[iF][i].pcolormesh(Zones[k].gridR, Zones[k].gridZ, MaskedValues, cmap=pyp.get_cmap('Reds'), shading=shading,  vmin = Vmin,  vmax = Vmax))
						else:			Im[-1] = AxPl[iF][i].pcolormesh(Zones[k].gridR, Zones[k].gridZ, MaskedValues, cmap=pyp.get_cmap('Reds'), shading=shading,  vmin = Vmin,  vmax = Vmax)

					plot2d_walls(AxPl[iF][i], Config.Walls)

					cb = FigPl[iF][i].colorbar(Im[-1], ax=AxPl[iF][i])
					cb.set_label(Labels[iF][i])

#				except:
#					print("WARNING: Not found ",Pars[iF][i]," for iPlasma = ", iPlasma)
#					print("\tIon = ",Plasmas[iPlasma][0].ion, " Charge=",Plasmas[iPlasma][0].charge)

	if(save != "none"):
		if(save == "stat"):
			save_stat("stat_2d.csv", StatHeader, StatValues, StatFormats)
		else:
			for i in range(len(Fig)):
				i_plot_file += 1
				if(one_plot != 1): Fig[i].set_size_inches(10.05,7.44)
				if(save == "pdf"):
					print("save pdf page=",i)
					pdf.savefig(Fig[i])
				else:
					Fig[i].savefig("plot2d_counters_on_mesh.py_t={:.3f}_{:d}.".format(tempus,i_plot_file)+save)
			pyp.show(block=False)
			pyp.close(0)
	else:
		pyp.show()

	if(save == "pdf"):	pdf.close()

	print("plot2d_counters_on_mesh.py: Completed")

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

def  search_min_max(FigNum, PosNum, MinMax):
	nMinMax = int(len(MinMax)/3)
	if((nMinMax == 0) or (nMinMax*3 != len(MinMax))): return False, 0.

	for i in range(nMinMax):
		if((MinMax[3*i] == FigNum) and (MinMax[3*i+1] == PosNum) ): 
			
			return True, MinMax[3*i+2]

	return False, 0.