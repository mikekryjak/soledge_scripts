#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot2d_custom import plot2d_custom
	from routines.cli_routines import *

	if(cli_present("-h",sys.argv)):
		print("\nThis script plot custom plot temporal files\n")
		print("plot_balances options")
		print("\t-name        Name of custom plot [d='./']")
		print("\t-last        If <= 0 ignored, <= 10 plot last time interval, > 10 plot last points [D=0]")
		print("\t-path        Directory with simulation [d='./']")
		print("\t-d_only      plot only e- and D+ [d=false]")
		print("\t-one_plot    One plot on each figure [d=false]")
		print("\t-no_samexy   No same x&y scale [d=false]")
		print("\t-no_zero     No zero on y scale [d=false]")
		print("\t-save        Save figures on files, values=none/png/ps/eps/pdf [D='none']")
		print()
		exit()

	custom_name	 = cli_get_value("-name",	sys.argv,"")
	path	 	= cli_get_value("-path",	sys.argv,"")
	last		= cli_get_value("-last",	sys.argv,	0)
	d_only		= cli_present("-d_only",	sys.argv)
	log_scale	= cli_present("-log_scale",	sys.argv)
	one_plot	= cli_present("-one_plot",	sys.argv)
	no_samexy	= cli_present("-no_samexy",	sys.argv)
	no_zero		= cli_present("-no_zero",	sys.argv)
	save		= cli_get_value("-save",	sys.argv, "none")
	plot2d_custom(custom_name=custom_name, path=path, last=last, log_scale=log_scale, d_only=d_only, no_samexy=no_samexy, no_zero=no_zero, one_plot=one_plot, save=save)
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

from routines.utils_walls				import plot2d_walls
from matplotlib.colors 					import LogNorm
from mesh.get_rz_core_sep				import get_rz_core_sep
from files.load_ions_list				import load_ions_list
from files.load_custom_plot_files			import load_custom_plot_line_files
from files.load_soledge_mesh_file		import load_soledge_mesh_file

#==============================================================================
# This routine plots ne/ni and Te/Ti ionization and gas pressure on eirene mesh
#==============================================================================

def plot2d_custom(custom_name, path="", last=0, log_scale=0, d_only=0, one_plot=0, no_samexy=0, no_zero=0, save="none"):

	print("plot2d_custom")

	if((len(path) > 0) and (path[-1] != "/")): path = path + "/"

	Config = load_soledge_mesh_file(path+"mesh.h5")
	RZcore, RZsep, CoreMegazone, SepMegazone, jSep = get_rz_core_sep(Config, core_and_sep = True)
	del Config, CoreMegazone, SepMegazone, jSep

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

	Fig = []
	Ax  = []
	Im  = []
	for iPlasma in range(nPlasmas):
		if(iPlasma == 0):
			Labels = [["$n\ (*10^{19}\ m^{-3})$", "$T\ (keV)$",  "$Rad_{tot} (kW/m^3)$"]]
			Pars   = [["density", "temperature","rad"]]
			PosPars = [[True, True, True]]
			NrowNcol = [[2,2]]
			VPosPlot	= [[2,3,4]]																	#Position of plasma plot
			Facts		  	  = [[1e-19,1.e-3, 1e-3]]
			FactsUnits		  = [["10^19","10^3", "10^3"]]
			MaxRange	      = [[1.e5,0., 0.]]
			BottomZero	= [[True,True,True]]
			LogScales	= [["log", "log", "log"]]

		elif(iPlasma == 1):
			Labels = [["$n\ (*10^{19}\ m^{-3})$", "$T\ (keV)$",  "$M$", "$Rad (kW/m^3)$"]]
			Pars   = [["density", "temperature", "mach","rad"]]
			PosPars = [[True, True, False, True]]
			NrowNcol = [[2,2]]
			VPosPlot	= [[1,2,3,4]]																	#Position of plasma plot
			Facts  	 		  = [[1e-19, 1.e-3, 1., 1e-3]]
			FactsUnits 		  = [["10^19", "10^3", "1.", "10^3"]]
			MaxRange		  = [[1.e5, 0., 0., 0.]]
			BottomZero	= [[True,True,False,True]]
			LogScales	= [["log", "log", "linear","log"]]

		else:
			Labels   = [["$n\ (*10^{19}\ m^{-3})$", "$T\ (keV)$",  "$Rad_{ion} (kW/m^3)$"]]
			Pars     = [["density", "temperature", "rad"]]
			PosPars = [[True, True, True]]
			NrowNcol = [[2,2]]
			VPosPlot	= [[1,2,3]]																	#Position of plasma plot
			Facts  	 = [[1e-19,1.e-3, 1e-3]]
			FactsUnits = [["10^19","10^3","10^3"]]
			MaxRange	      = [[1.e5,0., 0., 1.e5]]
			BottomZero	= [[True,True,True]]
			LogScales	= [["log", "log", "log"]]

		custom_data=load_custom_plot_line_files(custom_name, Path=path, iCustoms = [iPlasma])
		if(last > custom_data[0].Values.shape[0]): n_last = 0
		else: n_last = last

		FigNum = []
		FigPl = []
		AxPl = []
		Im	= []
		if(one_plot != 1):
			for iF in range(len(Labels)):
				AxPl.append([])
				Im.append([])
				FigNum.append([])
				FigPl.append([])
				Fig.append(pyp.figure())
				Fig[-1].patch.set_facecolor('white')

				if((iF == 0) and (iPlasma == 0)): Ax.append(Fig[-1].add_subplot(NrowNcol[iF][0],NrowNcol[iF][1],1))
				for i in range(len(Labels[iF])):
					if(((iPlasma==0) and (i == 0)) or (no_samexy != 0)):	Ax.append(Fig[-1].add_subplot(NrowNcol[iF][0],NrowNcol[iF][1],VPosPlot[iF][i]))
					else:																					Ax.append(Fig[-1].add_subplot(NrowNcol[iF][0],NrowNcol[iF][1],VPosPlot[iF][i],sharex = Ax[1], sharey = Ax[1]))
					AxPl[-1].append(Ax[-1])
					Im[-1].append(0)
					FigNum[-1].append(len(Fig))
					FigPl[-1].append(Fig[-1])
				Fig[-1].tight_layout(pad=2., w_pad=2., h_pad=2.)
		else:
			if(iPlasma == 0): 
				Fig.append(pyp.figure())
				Fig[-1].patch.set_facecolor('white')
				Ax.append(Fig[-1].add_subplot(111))

			for iF in range(len(Labels)):
				AxPl.append([])
				FigPl.append([])
				Im.append([])
				FigNum.append([])
				for i in range(len(Labels[iF])):
					Fig.append(pyp.figure())
					Fig[-1].patch.set_facecolor('white')
					Ax.append(Fig[-1].add_subplot(111))
					AxPl[-1].append(Ax[-1])
					Im[-1].append(0)
					FigNum[-1].append(len(Fig))
					FigPl[-1].append(Fig[-1])
	
		for iF in range(len(Labels)):
			for i in range(len(Labels[iF])):
				AxPl[iF][i].set_xlabel("$l\ (m)$")
				AxPl[iF][i].set_ylabel("$t\ (s)$")
				if(i == 0): AxPl[iF][i].set_title(ions[iPlasma])

		if(iPlasma == 0): 
			Ax[0].set_title(os.path.basename(os.path.abspath(path)))
			Ax[0].set_aspect(1.)
			Ax[0].set_xlabel("$R\ (m)$")
			Ax[0].set_ylabel("$Z\ (m)$")

			plot2d_walls(Ax[0], Config.Walls)
			Ax[0].plot(RZcore[:,0],  RZcore[:,1],  'b-')
			Ax[0].plot(RZsep[:,0],   RZsep[:,1],  'g-')
			Ax[0].plot(custom_data[0].Rgeom,  custom_data[0].Zgeom,  'r-')

			dist = np.sqrt((custom_data[0].Rgeom[1:]-custom_data[0].Rgeom[:-1])**2+(custom_data[0].Zgeom[1:]-custom_data[0].Zgeom[:-1])**2)
			Lenght =np.cumsum(dist)-dist/2
			ExtLenght = np.append(np.append(-Lenght[0],Lenght),Lenght[-1]+dist[-1]/2)

			if(n_last > 0): Times = custom_data[0].Times[-n_last:]
			else:				Times = custom_data[0].Times[-n_last:]
			ExtTimes = 0.5*(Times[:-1]+Times[1:])
			ExtTimes = np.append(np.append(ExtTimes[0]+0.5*(Times[1]-Times[0])/2,ExtTimes),ExtTimes[-1]+0.5*(Times[-1]-Times[2])/2)

	#	Plot parameters

		for iF in range(len(Labels)):
			for i in range(len(Labels[iF])):
				try:
					iPar = custom_data[0].Names.index(Pars[iF][i])
				except:
					print("WARNING: Not found ",Pars[iF][i]," for iPlasma = ", iPlasma)
					print("\tAvailable names=",custom_data[0].Names)
					exit()

				if(n_last > 0): Values = custom_data[0].Values[:,:,iPar]
				else:				Values = custom_data[0].Values[:,-n_last:,iPar]
				MinValue = Values[:,:].min()
				MaxValue = Values[:,:].max()
				if((log_scale != 0) and PosPars[iF][i] and (MinValue > 0.)):
					Im[iF][i] = AxPl[iF][i].pcolormesh(ExtLenght, ExtTimes, Values[:,:]*Facts[iF][i], norm=LogNorm(vmin=MinValue*Facts[iF][i], vmax=MaxValue*Facts[iF][i]))
				else:
					Im[iF][i] = AxPl[iF][i].pcolormesh(ExtLenght, ExtTimes, Values[:,:]*Facts[iF][i])
				AxPl[iF][i].autoscale(enable=True, axis='both', tight=True)

				cb = FigPl[iF][i].colorbar(Im[iF][i], ax=AxPl[iF][i])
				cb.set_label(Labels[iF][i])

	if(save != "none"):
		for i in range(len(Fig)):
			i_plot_file += 1
			if(one_plot != 1): Fig[i].set_size_inches(10.05,7.44)
			if(save == "pdf"):
				print("save pdf page=",i)
				pdf.savefig(Fig[i])
			else:
				Fig[i].savefig("plot2d_custom_"+custom_name+"_{:d}.".format(iPlasma)+save)
#		pyp.show(block=False)
#		pyp.close()
	else:
		pyp.show()

	if(save == "pdf"):	pdf.close()

	print("plot2d_custom: Completed")

	return
