#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot1d_custom import plot1d_custom
	from routines.cli_routines import *

	if(cli_present("-h",sys.argv)):
		print("\nThis script plot custom plot temporal files\n")
		print("plot_balances options")
		print("\t-name        Name of custom plot [d='./']")
		print("\t-path        Directory with simulation [d='./']")
		print("\t-last        If <= 0 ignored, <= 10 plot last time interval, > 10 plot last points [D=0]")
		print("\t-use_time    Plot versus time [d=false]")
		print("\t-d_only      plot only e- and D+ [d=false]")
		print("\t-log_scale   Use log scale for y axis [d=false]")
		print("\t-one_plot    One plot on each figure [d=false]")
		print("\t-no_samex    No same x scale [d=false]")
		print("\t-no_zero     No zero on y scale [d=false]")
		print("\t-use_pts     Use points instead than lines [d=false]")
		print("\t-save        Save figures on files, values=none/png/ps/eps/pdf [D='none']")
		print()
		exit()

	custom_name	 	= cli_get_value("-name",			sys.argv,"")
	path	 	= cli_get_value("-path",			sys.argv,"")
	last		= cli_get_value("-last",			sys.argv,	0.)
	use_time	= cli_present("-use_time",		sys.argv)
	d_only			= cli_present("-d_only",				sys.argv)
	log_scale	= cli_present("-log_scale",		sys.argv)
	one_plot	= cli_present("-one_plot",		sys.argv)
	no_samex	= cli_present("-no_samex",		sys.argv)
	no_zero	= cli_present("-no_zero",		sys.argv)
	use_pts    	= cli_present("-use_pts",		sys.argv)
	save		= cli_get_value("-save",			sys.argv, "none")
	plot1d_custom(custom_name=custom_name, path=path, last=last, use_time=use_time, use_pts=use_pts ,log_scale=log_scale, d_only=d_only, no_samex=no_samex, no_zero=no_zero, one_plot=one_plot, save=save)
	exit()

#=======================================

# Function definition is here

import h5py
import os
import numpy								as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot					as pyp
from matplotlib.backends.backend_pdf		import PdfPages
from files.load_ions_list				import load_ions_list
from files.load_custom_plot_files				import load_custom_plot_temporal_files

#==============================================================================
# This routine plots ne/ni and Te/Ti ionization and gas pressure on eirene mesh
#==============================================================================

def plot1d_custom(custom_name, path="", last=0., use_time=0, use_pts=0, log_scale=0, d_only=0, one_plot=0, no_samex=0, no_zero=0, save="none"):

	print("plot1d_custom")

	if((len(path) > 0) and (path[-1] != "/")): path = path + "/"

	ions	= load_ions_list(path)										#read ions list
	if(d_only == 0):
		nPlasmas = len(ions)
		Atoms  = [ions[1][:ions[1].find("+")-1]]
		lAtoms = [len(Atoms[-1])]
		for i in range(1,nPlasmas):
			if(Atoms[-1] != ions[i][:lAtoms[-1]]):
				Atoms.append(ions[i][:ions[i].find("+")-1])
				lAtoms.append(len(Atoms[-1]))
				
		nAtoms = len(Atoms)
	else:
		nPlasmas=2

	if(use_pts == 0): V_lines	= ["-b",  "-g",   "-r",  "-c", "-m",  "-y", "-k","--m", "--k"]
	else:			  V_lines	= [".b",  ".g",   ".r",  ".c", ".m",  ".y", "k.","om",  "ok"]

	Fig = []
	Ax  = []
	Im  = []
	for iPlasma in range(nPlasmas):
		if(iPlasma == 0):
			Labels = [["$n\ (*10^{19}\ m^{-3})$", "$T\ (keV)$",  "$Rad_{tot} (kW/m^3)$"]]
			Pars   = [["density", "temperature","rad"]]
			PosPars = [[True, True, True]]
			NrowNcol = [[3,1]]
			Facts		  	  = [[1e-19,1.e-3, 1e-3]]
			FactsUnits		  = [["10^19","10^3", "10^3"]]
			MaxRange	      = [[1.e5,0., 0.]]
			BottomZero	= [[True,True,True]]
			LogScales	= [["log", "log", "log"]]

		elif(iPlasma == 1):
			Labels = [["$n\ (*10^{19}\ m^{-3})$", "$T\ (keV)$",  "$M$", "$Rad (kW/m^3)$"]]
			Pars   = [["density", "temperature", "mach","rad"]]
			PosPars = [[True, True, True, True]]
			NrowNcol = [[4,1]]
			Facts  	 		  = [[1e-19, 1.e-3, 1., 1e-3]]
			FactsUnits 		  = [["10^19", "10^3", "1.", "10^3"]]
			MaxRange		  = [[1.e5, 0., 0., 0.]]
			BottomZero	= [[True,True,False,True]]
			LogScales	= [["log", "log", "linear","log"]]

		else:
			Labels   = [["$n\ (*10^{19}\ m^{-3})$", "$T\ (keV)$",  "$Rad_{ion} (kW/m^3)$"]]
			Pars     = [["density", "temperature", "rad"]]
			PosPars = [[True, True, True]]
			NrowNcol = [[3,1]]
			Facts  	 = [[1e-19,1.e-3, 1e-3]]
			FactsUnits = [["10^19","10^3","10^3"]]
			MaxRange	      = [[1.e5,0., 0., 1.e5]]
			BottomZero	= [[True,True,True]]
			LogScales	= [["log", "log", "log"]]

		"""
		if(iPlasma == 0):	NameExt = "e"
		elif(iPlasma == 1):	NameExt = "i"
		else:				NameExt = Plasmas[iPlasma][0].ion[:-1]
		for iF in range(len(Pars)):
			for i in range(len(Pars[iF])): 
				if(Pars[iF][i][-1]!= " "): Pars[iF][i] += NameExt
				else:						Pars[iF][i]  = Pars[iF][i][:-1]
		"""

		custom_data=load_custom_plot_temporal_files(custom_name, Path=path, iCustoms = [iPlasma])
		if(last > custom_data[0].Values.shape[0]): n_last = 0
		else: n_last = last

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
					if(((iPlasma==0) and (i == 0)) or (no_samex != 0)):	Ax.append(Fig[-1].add_subplot(NrowNcol[iF][0],NrowNcol[iF][1],i+1))
					else:																					Ax.append(Fig[-1].add_subplot(NrowNcol[iF][0],NrowNcol[iF][1],i+1,sharex = Ax[0]))
					AxPl[-1].append(Ax[-1])
					FigNum[-1].append(len(Fig))
					PosNum[-1].append(i+1)
					FigPl[-1].append(Fig[-1])
				Fig[-1].tight_layout(pad=2., w_pad=2., h_pad=2.)
		else:
			for iF in range(len(Labels)):
				AxPl.append([])
				FigPl.append([])
				FigNum.append([])
				PosNum.append([])
				for i in range(len(Labels[iF])):
					Fig.append(pyp.figure())
					Fig[-1].patch.set_facecolor('white')
					if(((iPlasma==0) and (i == 0)) or (no_samex != 0)):	Ax.append(Fig[-1].add_subplot(111))
					else:																					Ax.append(Fig[-1].add_subplot(111,sharex = Ax[0]))

					AxPl[-1].append(Ax[-1])
					FigNum[-1].append(len(Fig))
					PosNum[-1].append(1)
					FigPl[-1].append(Fig[-1])
	
		for iF in range(len(Labels)):
			for i in range(len(Labels[iF])):
				AxPl[iF][i].set_xlabel("$t\ (s)$")
				AxPl[iF][i].set_ylabel(Labels[iF][i])
				if(i == 0): AxPl[iF][i].set_title(os.path.basename(os.path.abspath(path)) + " " +ions[iPlasma])
				if(log_scale == 0):
					AxPl[iF][i].set_yscale('linear')
				else:
					AxPl[iF][i].set_yscale(LogScales[iF][i])

	#	Plot parameters

		for iF in range(len(Labels)):
			for i in range(len(Labels[iF])):
				try:
					iPar = custom_data[0].Names.index(Pars[iF][i])
				except:
					print("WARNING: Not found ",Pars[iF][i]," for iPlasma = ", iPlasma)
					print("\tAvailable names=",custom_data[0].Names)
					exit()
				if(n_last > 0): AxPl[iF][i].plot(custom_data[0].Values[-n_last:,0],custom_data[0].Values[-n_last:,iPar]*Facts[iF][i], V_lines[0])
				else:				AxPl[iF][i].plot(custom_data[0].Values[:,0],           custom_data[0].Values[:,iPar]*Facts[iF][i],            V_lines[0])
				AxPl[iF][i].autoscale(enable=True, axis='both', tight=True)
				if((log_scale == 0)  and (no_zero == 0) and BottomZero[iF][i]):  AxPl[iF][i].set_ylim(bottom=0.)

	if(save != "none"):
		for i in range(len(Fig)):
			i_plot_file += 1
			if(one_plot != 1): Fig[i].set_size_inches(10.05,7.44)
			if(save == "pdf"):
				print("save pdf page=",i)
				pdf.savefig(Fig[i])
			else:
				Fig[i].savefig("plot1d_custom_"+custom_name+"_{:d}.".format(iPlasma)+save)
		pyp.show(block=False)
		pyp.close()
	else:
		pyp.show()

	if(save == "pdf"):	pdf.close()

	print("plot1d_custom: Completed")

	return
