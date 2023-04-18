#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot1d_diags 			import plot1d_diags
	from routines.cli_routines	import *

	if(cli_present("-h",sys.argv) and cli_present("-exp_file",sys.argv) ):
		from routines.exp_data_routines			import get_exp_data_los_help, get_exp_data_rz_help, get_exp_data_spectra_help
		from files.load_exp_descr				import load_exp_descr_help
		load_exp_descr_help()
		get_exp_data_rz_help()
		get_exp_data_los_help()
		exit()

	if(cli_present("-h",sys.argv)):
		print("\nThis script compute derived data on diagnostic lines of sight\n")
		print("plot1d_diags options")
		print("\t-path           Directory with simulation or list of directories[d='']")
		print("\t-evolution      If > 0 plot evolution file (or evolutions files with list [3,4,5])  [d=0]")
		print("\t-exp_file       Name or list exp file with description of diagnostic data [d=['']]")
		print("\t-shot           Shot to compare [d=0]")
		print("\t-tstart         Start time experimental data to compare [d=0.]")
		print("\t-tend           End time experimental data to compare [d=0.]")
		print("\t-diag_file      Name or list of diag files for no experimental data[D=['']")
		print("\t-all_data       Plot also RZ data points inside grid inner boundary [d=false]")
		print("\t-path_label     Labels to itentify runs in plot [d='']")
		print("\t-label_data     Label RZ points and lines on 2D plot [d=false]")
		print("\t-log_scale      Use log scale for y axis [d=false]")
		print("\t-one_plot       One plot on each figure [d=false]")
		print("\t-save           Save figures/data on files, values=none/png/ps/eps/pdf/csv [D='none']")
		print()
		exit()

	path	 	= cli_get_value("-path",				sys.argv,   [""])
	evolution	= cli_get_value("-evolution",		sys.argv,	[])
	diag_file	= cli_get_value("-diag_file",		sys.argv,	[""])
	exp_files	= cli_get_value("-exp_file",			sys.argv,[""])
	shot		= cli_get_value("-shot",				sys.argv,	0)
	tstart		= cli_get_value("-tstart",			sys.argv,	0.)
	tend		= cli_get_value("-tend",				sys.argv,	0.)
	all_data	= cli_present("-all_data",			sys.argv)
	path_label 	= cli_get_value("-path_label",		sys.argv,   [""])
	label_data	= cli_present("-label_data",		sys.argv)
	log_scale	= cli_present("-log_scale",			sys.argv)
	one_plot	= cli_present("-one_plot",			sys.argv)
	save		= cli_get_value("-save",		sys.argv, "none")
	plot1d_diags(path=path, evolution=evolution, exp_files=exp_files, shot=shot, tstart=tstart, tend=tend, diag_file=diag_file, path_label=path_label, label_data=label_data, all_data=all_data, log_scale=log_scale, one_plot=one_plot, save=save)
	exit()

#=======================================

# Function definition is here

import types
import os
import h5py
from math								import sqrt, exp, log, log10
import numpy							as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot				as pyp
import matplotlib.tri					as tri
from matplotlib.colors 					import LogNorm
from matplotlib.path					import Path
from matplotlib.backends.backend_pdf	import PdfPages

from mesh.get_mesh_path					import get_chords_mesh_path
from eirene.get_triangles_path			import get_chords_triangle_path
from eirene.get_ext_plasma_triseq		import get_ext_plasma_triseq
from mesh.get_rz_core_sep				import get_rz_core_sep
from routines.exp_data_routines			import get_exp_data_los, get_exp_data_rz, get_exp_data_spectra
from routines.utils_routines			import extend_mat, extend_mat1
from routines.utils_walls				import plot2d_walls
from routines.globals					import DEBUG

from files.load_soledge_mesh_file		import load_soledge_mesh_file
from files.load_plasma_files			import load_plasma_files
from files.load_eirene_triangles		import load_eirene_triangles
from files.load_diag_file				import load_diag_file
from files.load_ions_list				import load_ions_list
from files.load_refpar_file				import load_refpar_file
from files.load_adas_15					import load_adas_15
from files.load_exp_descr				import load_exp_descr
from scipy								import interpolate


#==============================================================================
# This 
#==============================================================================

def plot1d_diags(path=[], evolution=[], exp_files=[""], shot=0, tstart=0., tend=0., diag_file=[], path_label=[], label_data=0, all_data=0, log_scale=0, one_plot=0, save="none"):

	print("plot1d_diags")

#	prepare for experimental or modeled data

	if((len(exp_files) > 0) and (len(exp_files[0]) > 0)):
		ExpDescrs = []
		Diags	  = []
		for exp_file in exp_files:
			ExpDescrs.append(load_exp_descr(exp_file))					#load diags description
			if((tstart > 0) and (tstart <= tend)):
				ExpDescrs[-1].tstart = tstart
				ExpDescrs[-1].tend   = tend
			if(shot != 0): ExpDescrs[-1].shot = shot

		for exp_descr in ExpDescrs:
			Diag = get_exp_data_rz(exp_descr)							#load diags data for rz diagnostics
			if(len(Diag)): Diags.append(Diag)
		nExpRZ = len(Diags)

		for exp_descr in ExpDescrs:
			Diag = get_exp_data_los(exp_descr)							#load diags data for los diagnostics
			if(len(Diag)): Diags.append(Diag)
		nExpLine = len(Diags) - nExpRZ

		for exp_descr in ExpDescrs:
			Diag = get_exp_data_spectra(exp_descr)						#load diags data for los spectrum
			if(len(Diag)): Diags.append(Diag)
		nExpSpectra = len(Diags) - nExpRZ - nExpLine
	
	else:
		Diags = [[]]
		nParameters = len(diag_file)
		for iDiag in range(nParameters):
			Diags[0].append(load_diag_file(diag_file[iDiag]))

#	Put together diags with same parameter

	FirstDiag = True
	Par_names = []
	for iExp in range(len(Diags)):
		for iDiag in range(len(Diags[iExp])):
			if(FirstDiag):
				Par_names = [Diags[0][0].par_name]
				Par_diags = [[[0],[0],0]]
				FirstDiag = False

			iData = indexes_upper(Par_names, Diags[iExp][iDiag].par_name)
			if((iExp != 0) or (iDiag != 0)):
				if(len(iData) == 0):											#New parameter
					Par_names.append(Diags[iExp][iDiag].par_name)				#Append parameter name
					Par_diags.append([[iExp],[iDiag],0])						#New experiment and diag number
				elif(len(iData) == 1):											#Parameter exist
					Par_diags[iData[0]][0].append(iExp)							#Append in parameter experiment
					Par_diags[iData[0]][1].append(iDiag)						#Append Diag
				else:
					print("\tError in scanning parameters")
					exit()

	nParameters = len(Par_names)

	if(nParameters == 0):
		print("\n\nNO DATA AVAILABLE!!!\n\n")
		exit()
	

#	prepare paths

	if(len(evolution) == 0): evolution = [0]
	if(len(path) == 0):		 path = [""]
	elif(len(path) > 1):	 evolution = [evolution[0]]

	if(len(path_label) < len(path)):
		path_label = []
		for in_path  in path: path_label.append(os.path.basename(os.path.abspath(in_path)))

	if((evolution[0] == -1) or (evolution[0] == 0)):	UseTri = True
	else:					 													UseTri = False
	nRuns		= max(len(path),len(evolution))
	i_plot_file = 0

#	Read reference parameters
	
	path0  = path[0]
	if((len(path0) > 0) and (path0[-1] != "/")): path0 = path0 + "/"

	RefPar = load_refpar_file(path0+"Results/")

	ions = load_ions_list(path0)

#	Read mesh

	Config = load_soledge_mesh_file(path0+"mesh.h5")

	if(UseTri):
		Eirene = load_eirene_triangles(path0+"triangles.h5")
		get_ext_plasma_triseq(Config, Eirene.Wall)

		Eirene.TriKnots		= np.array([Eirene.Triangles.p1, 			Eirene.Triangles.p2,			Eirene.Triangles.p3]).T
		Eirene.TriNeigh		= np.array([Eirene.Triangles.neigh1,		Eirene.Triangles.neigh2,		Eirene.Triangles.neigh3]).T
		Eirene.TriTypeneigh	= np.array([Eirene.Triangles.typeneigh1,	Eirene.Triangles.typeneigh2,	Eirene.Triangles.typeneigh3]).T
		Eirene.TriSeq		=  Eirene.WallTriangles.ntri[Eirene.Wall.TriSequences[Eirene.Wall.iTriSeqExtPlasma]]


	
	if(nRuns > 1):								#If multiple simulation disable plot2D
		Plot2D = False
	else:
		Plot2D = True
		RZcore = None
		if(UseTri): 
			TripTriang = tri.Triangulation(Eirene.RKnots, Eirene.ZKnots, triangles=Eirene.TriKnots)
		else:
			Zones  = Config.Zones
			nZones = len(Zones)
			for k in range(nZones): Zones[k].Chi2 = extend_mat(Zones[k].Chi)

	DataCenterMesh = types.SimpleNamespace()
	if(nExpRZ > 0): 							#prepare for rz_points 
		if(UseTri):
			DataCenterMesh.RCent	= (Eirene.RKnots[Eirene.Triangles.p1] + Eirene.RKnots[Eirene.Triangles.p2] + Eirene.RKnots[Eirene.Triangles.p3])/3.
			DataCenterMesh.ZCent	= (Eirene.ZKnots[Eirene.Triangles.p1] + Eirene.ZKnots[Eirene.Triangles.p2] + Eirene.ZKnots[Eirene.Triangles.p3])/3.
		else:
			Zones  = Config.Zones
			nZones = len(Zones)
			RCent  = np.empty(0, dtype='f8')
			ZCent  = np.empty(0, dtype='f8')
			kMesh  = np.empty(0, dtype='i4')
			iMesh  = np.empty(0, dtype='i4')
			jMesh  = np.empty(0, dtype='i4')
			for k in range(nZones):
				RCent = np.append(RCent, Zones[k].gridRc.reshape(-1))
				ZCent = np.append(ZCent, Zones[k].gridZc.reshape(-1))
				kMesh = np.append(kMesh, k*np.ones(Zones[k].gridZc.shape, dtype='i4'))

				jMat,iMat = np.meshgrid(np.arange(Zones[k].gridRc.shape[1],dtype='i4'),np.arange(Zones[k].gridRc.shape[0],dtype='i4'))
				iMesh = np.append(iMesh, iMat.reshape(-1))
				jMesh = np.append(jMesh, jMat.reshape(-1))
			del jMat,iMat
			DataCenterMesh.RCent = RCent
			DataCenterMesh.ZCent = ZCent
			DataCenterMesh.kMesh = kMesh
			DataCenterMesh.iMesh = iMesh
			DataCenterMesh.jMesh = jMesh

		DataCenterMesh.all_data =  all_data
		if(all_data == 0):
			Rcore, Zcore, CoreMegazone = get_rz_core_sep(Config, core_and_sep = False)
			DataCenterMesh.CorePath = Path(np.array([Rcore,Zcore]).T, closed=True)
			del Rcore, Zcore, CoreMegazone 

	if(nExpLine+nExpSpectra > 0): 											#Find mesh along line
		for iExp in range(nExpRZ,nExpRZ+nExpLine+nExpSpectra):
			for iDiag in range(len(Diags[iExp])):
				for iSub in range(len(Diags[iExp][iDiag].SubDiags)):
					if(UseTri): get_chords_triangle_path(Eirene, Diags[iExp][iDiag].SubDiags[iSub])
					else:		get_chords_mesh_path(Config, Diags[iExp][iDiag].SubDiags[iSub])
		

	for iExp in range(len(Diags)):												#Prepara Value array
		for iDiag in range(len(Diags[iExp])):
			for iSub in range(len(Diags[iExp][iDiag].SubDiags)):
				Diags[iExp][iDiag].SubDiags[iSub].Values   = []
				Diags[iExp][iDiag].Values2D = []



#	read plasma data and compute parameters on lines

	nParData = 0
	for iPh in range(len(path)):
		for iEv in range(len(evolution)):
			path0  = path[iPh]
			if((len(path0) > 0) and (path0[-1] != "/")): path0 = path0 + "/"
			Plasmas = load_plasma_files(path0, nZones=len(Config.Zones), Evolution=evolution[iEv])

			for iExp in range(len(Diags)):
				for iDiag in range(len(Diags[iExp])):
					get_plasma_parameter(Plasmas, Diags[iExp][iDiag], DataCenterMesh, UseTri, Plot2D)

#	Count plot plots to manage spectra which produce one plot for each channel

	nPlots = 0
	for iParam in range(nParameters):
		if(Diags[Par_diags[iParam][0][0]][Par_diags[iParam][1][0]].type == "SPECTRA"): 	
			for iData in range(len(Par_diags[iParam][0])):
				for SubDiag in Diags[Par_diags[iParam][0][iData]][Par_diags[iParam][1][iData]].SubDiags: Par_diags[iParam][2] += len(SubDiag.path.Leni)
		else:
			Par_diags[iParam][2] = 1
		nPlots += Par_diags[iParam][2]

#	Get core and seperatrix in case of not plotting 2D countour

	if((not Plot2D) or (nPlots > nParameters)):	RZcore, RZsep, CoreMegazone, SepMegazone, jSep = get_rz_core_sep(Config, core_and_sep = True, use_mag_zones = False)

#	Prepare for plotting and saving data

	if(save == "pdf"):	pdf = PdfPages("plot1d_diags_t={:.3f}.".format(RefPar.time)+save)   #pdf in one file only
	if(save == "csv"):
		if((len(evolution) > 1) or (len(path) > 1)):
			save = "none"
			print("plot1d_diags: cannot save to csv more than one profile!!!")
		else:
			csv = []

	shading = "flat"
	colors  = ['b','g','r','c','m','y','b','g','r','c','m','y']

	if(nPlots > 1):
		nRows	   = 2
		nCols	   = 3
		PlotPerFig = 2
	else:			
		nRows	   = 1
		nCols	   = 3
		PlotPerFig = 1
		
	nFigs = int(nPlots/PlotPerFig)
	if(nFigs*PlotPerFig < nPlots): nFigs += 1

	Fig = []
	Ax  = []
	Im  = []
	if(one_plot != 1):
		for i in range(nFigs):	
			Fig.append(pyp.figure())
			iRow = 0
			for k in range(min(PlotPerFig,nPlots-i*PlotPerFig)):
				Ax.append(pyp.subplot2grid((nRows,nCols),(iRow,0)))
				Ax[-1].locator_params(axis='x',nbins=4)
				Ax.append(pyp.subplot2grid((nRows,nCols),(iRow,1), colspan=2))
				Ax[-1].locator_params(axis='x',nbins=4)
				iRow += 1

			Fig[-1].tight_layout(pad=2., w_pad=3., h_pad=3.)
	else:
		for i in range(2*nPlots):
			Fig.append(pyp.figure())
			Ax.append(Fig[i].add_subplot(111))

	for figure in Fig:  figure.patch.set_facecolor('white')

	nPlot = -1
	for iParam in range(nParameters):
		for iPlotChan in range(Par_diags[iParam][2]):				#If necessary loop on channels (only for SPECTRUM diagnostics
			nPlot += 1

			ip = 2*nPlot

			if(len(path) == 1):
				if((len(evolution) != 1) or (evolution[0] == 0)):
					Ax[ip].set_title(path_label[0]+" @ t={:.3f} s".format(RefPar.time))
				else:
					Ax[ip].set_title(path_label[0]+":{:d} @ t={:.3f} s".format(evolution[0],RefPar.time))


#====================
#		2D plot
#====================UseTri = False

			Ax[ip].set_aspect(1.)
			Ax[ip].autoscale(enable=True, axis='both', tight=True)
			Ax[ip].set_xlabel("$R (m)$")
			Ax[ip].set_ylabel("$Z (m)$")

#		Plot SOLEDGE2D results
#		==============

			Diag = Diags[Par_diags[iParam][0][0]][Par_diags[iParam][1][0]]
			if(Diag.Plot2D):
				if(UseTri):																					#Triangles
					if(Diag.v_min != None): Vmin = Diag.v_min												#Set min
					else:					Vmin = Diag.Values2D.min()
					if(Diag.v_max != None): Vmax = Diag.v_max												#Set max
					else:					Vmax = Diag.Values2D.max()

					v_log = False
					if(Diag.v_log != None): v_log = Diag.v_log
					elif(log_scale != 0):	v_log = True

					if((not v_log) or (Vmin < 0.)):
						Im.append(Ax[ip].tripcolor(TripTriang, Diag.Values2D))
					else:
						Im.append(Ax[ip].tripcolor(TripTriang, np.where(Diag.Values2D > Vmin, Diag.Values2D, Vmin), norm=LogNorm(vmin=Vmin, vmax=Vmax)))
				else:																					#Quadrangles
					if(Diag.v_min != None): Vmin = Diag.v_min											#Set min
					else:
						Vmin = +1e40
						for iZone in range(nZones):
							Vmin = min(Vmin, np.min(Diag.Values2D[iZone]))

					if(Diag.v_max != None): Vmax = Diag.v_max											#Set Max
					else:
						Vmax= -1e40
						for iZone in range(nZones):
							Vmax = max(Vmax, np.max(Diag.Values2D[iZone]))

					for iZone in range(nZones):
						if(iZone == 0): Im.append(0)
						if(Diag.Values2D[iZone].size == Zones[iZone].Chi.size):		Values2D  = np.ma.masked_where(Zones[iZone].Chi== 1., Diag.Values2D[iZone])
						else:																							Values2D  = np.ma.masked_where(Zones[iZone].Chi== 1., Diag.Values2D[iZone][1:-1,1:-1])
						Im[-1] = Ax[ip].pcolormesh(Zones[iZone].gridR, Zones[iZone].gridZ, Values2D, shading=shading,  vmin = Vmin,  vmax = Vmax)

				cb = Ax[ip].get_figure().colorbar(Im[-1], ax=Ax[ip])
				cb.set_label(Diag.v_label)

			else:
				Ax[ip].plot(RZcore[:,0], RZcore[:,1],"g-",linewidth=1)
				Ax[ip].plot(RZsep[:,0],  RZsep[:,1], "g-",linewidth=1)


#		Plot diagnostics informations
#		=================

			kPlotChan = 0
			for iData in range(len(Par_diags[iParam][0])):
				iExp  = Par_diags[iParam][0][iData]
				iDiag = Par_diags[iParam][1][iData]
				Diag  = Diags[iExp][iDiag]

#			Plot RZ points positions (and chan name)

				if(Diag.type == "RZ_POINTS"):
					for iSub in range(len(Diag.SubDiags)):
						SubDiag = Diag.SubDiags[iSub]
						Ax[ip].plot(SubDiag.r_values, SubDiag.z_values,"r.")
						if(label_data == 1):
							for i in range(len(SubDiag.n_chans)):
								Ax[ip].text(SubDiag.r_values[i], SubDiag.z_values[i],"{:d}".format(SubDiag.n_chans[i]),backgroundcolor='w')
				elif(Diag.type == "LINES"):
#			Plot lines (and chan name)

					for iSub in range(len(Diag.SubDiags)):
						SubDiag = Diag.SubDiags[iSub]
						for iChan in range(len(SubDiag.path.Leni)):
							if(UseTri):	
								Xi = [SubDiag.path.Ri[iChan][0][0], SubDiag.path.Ri[iChan][-1][-1]]
								Zi = [SubDiag.path.Zi[iChan][0][0], SubDiag.path.Zi[iChan][-1][-1]]
								for k in range(len(SubDiag.path.Leni[iChan])):
#								Ax[ip].triplot(Eirene.RKnots, Eirene.ZKnots, Eirene.TriKnots, 'b-')
									Ax[ip].plot(SubDiag.path.Ri[iChan][k], SubDiag.path.Zi[iChan][k],"g.-",linewidth=2)
							else:
								Xi = [SubDiag.path.Xi[iChan][0,0], SubDiag.path.Xi[iChan][-1,1]]
								Zi = [SubDiag.path.Zi[iChan][0,0], SubDiag.path.Zi[iChan][-1,1]]
								Ax[ip].plot(Xi, Zi,"g.-",linewidth=2)

							if(label_data != 0):
								Ax[ip].text(0.5*(Xi[0]+Xi[1]),0.5*(Zi[0]+Zi[1]),"{:d}".format(SubDiag.n_chans[iChan]),backgroundcolor='w')
				elif(Diag.type == "SPECTRA"):
					for iSub in range(len(Diag.SubDiags)):
						SubDiag = Diag.SubDiags[iSub]
						kPlotChan += len(SubDiag.path.Leni)
						if((kPlotChan > iPlotChan) and (kPlotChan - len(SubDiag.path.Leni)  <= iPlotChan)):
							iChan = iPlotChan - (kPlotChan - len(SubDiag.path.Leni))
							if(UseTri):	
								Xi = [SubDiag.path.Ri[iChan][0][0], SubDiag.path.Ri[iChan][-1][-1]]
								Zi = [SubDiag.path.Zi[iChan][0][0], SubDiag.path.Zi[iChan][-1][-1]]
								for k in range(len(SubDiag.path.Leni[iChan])):
#								Ax[ip].triplot(Eirene.RKnots, Eirene.ZKnots, Eirene.TriKnots, 'b-')
									Ax[ip].plot(SubDiag.path.Ri[iChan][k], SubDiag.path.Zi[iChan][k],"g.-",linewidth=2)
							else:
								Xi = [SubDiag.path.Xi[iChan][0,0], SubDiag.path.Xi[iChan][-1,1]]
								Zi = [SubDiag.path.Zi[iChan][0,0], SubDiag.path.Zi[iChan][-1,1]]
								Ax[ip].plot(Xi, Zi,"g.-",linewidth=2)

							if(label_data != 0):
								Ax[ip].text(0.5*(Xi[0]+Xi[1]),0.5*(Zi[0]+Zi[1]),"{:d}".format(SubDiag.n_chans[iChan]),backgroundcolor='w')

			plot2d_walls(Ax[ip], Config.Walls)			

			if((Diag.r_min == None) or (Diag.r_max == None)):
				if(Diag.r_min != None): Ax[ip].set_xlim(left=Diag.r_min)
				if(Diag.r_max != None): Ax[ip].set_xlim(right=Diag.r_min)
			else: 						Ax[ip].set_xlim(left=Diag.r_min, right=Diag.r_max)

			if((Diag.z_min == None) or (Diag.z_min == None)):
				if(Diag.z_min != None): Ax[ip].set_ylim(bottom=Diag.z_min)
				if(Diag.z_max != None): Ax[ip].set_ylim(top=Diag.z_max)
			else: 						Ax[ip].set_ylim(bottom=Diag.z_min, top=Diag.z_max)

#====================
#		1D plot
#====================

			ip += 1
			if(Diag.title != None):	Ax[ip].set_title(Diag.title)
			else:					Ax[ip].set_title(Par_names[iParam])
			Ax[ip].autoscale(enable=True, axis='both', tight=True)

			Diag  = Diags[Par_diags[iParam][0][0]][Par_diags[iParam][1][0]]
			if(Diag.y_label != None): Ax[ip].set_ylabel(Diag.y_label)
			if(Diag.x_label != None): Ax[ip].set_xlabel(Diag.x_label)
			if(Diag.y_log != None):
				if(Diag.y_log):		Ax[ip].set_yscale('log')
				else:				Ax[ip].set_yscale('linear')
			else:
				if(log_scale == 0): Ax[ip].set_yscale('linear')
				else:				Ax[ip].set_yscale('log')
			if(Diag.x_log != None):
				if(Diag.x_log):		Ax[ip].set_xscale('log')
				else:				Ax[ip].set_xscale('linear')
				
	#		Plot SOLEDGE2D results
	#		==============

			if(nRuns == 1):									#Single simulation (one path & one evolution)
				x_values = np.empty(0, dtype = 'f8')
				y_values = np.empty(0, dtype = 'f8')
				kPlotChan = 0
				for iData in range(len(Par_diags[iParam][0])):
					iExp  = Par_diags[iParam][0][iData]
					iDiag = Par_diags[iParam][1][iData]
					Diag  = Diags[iExp][iDiag]
					for iSub in range(len(Diag.SubDiags)):
						SubDiag = Diag.SubDiags[iSub]
						plot_style = "bo"
						if(Diag.type == "SPECTRA"):
							iChan	   = iPlotChan - kPlotChan
							kPlotChan += len(SubDiag.path.Leni)
							if((kPlotChan > iPlotChan) and (iChan >= 0)):
								x_values = SubDiag.SpectrumWaves
								y_values = SubDiag.Values[0][:,iChan]
								SpectrumLines = SubDiag.SpectrumLines
								if(Diag.fwhm != None):
									x_values_label = np.copy(x_values)
									x_values, y_values, y_values_label = gaussfold(x_values, y_values, Diag.fwhm)
									plot_style = "b-"
								else:
									x_values_label = x_values
									y_values_label = y_values
								if(not Diag.use_wavelength):
									x_values 		= 1./x_values
									x_values_label	= 1./x_values_label
						else:
							x_values = np.append(x_values, SubDiag.x_values)
							y_values = np.append(y_values, SubDiag.Values[0])

				Ax[ip].plot(x_values, y_values,  plot_style, label="soledge2d")
				if((Diag.type == "SPECTRA") and (Diag.label_lines != 0)):
					for kLine in range(len(x_values_label)):
						if(y_values_label[kLine] > Diag.label_threshold):
							Ax[ip].text(x_values_label[kLine],y_values_label[kLine],"  " + SpectrumLines[kLine], clip_on=True, rotation="vertical",horizontalalignment="center", verticalalignment="bottom")

			
				if(save == "csv"):
					len_csv = len(csv)
					csv.append(types.SimpleNamespace())
					csv[-1].Name   = "ch_"+Par_names[iData]
					csv[-1].Values = np.copy(x_values)
					csv.append(types.SimpleNamespace())
					csv[-1].Name   = Par_names[iData]
					csv[-1].Values = np.copy(y_values/Diag.y_factor)
			else:												#Single simulation (one path & one evolution)
				for iRun in range(nRuns):
					x_values = np.empty(0, dtype = 'f8')
					y_values = np.empty(0, dtype = 'f8')
					for iData in range(len(Par_diags[iParam][0])):
						iExp  = Par_diags[iParam][0][iData]
						iDiag = Par_diags[iParam][1][iData]
						Diag  = Diags[iExp][iDiag]
						for iSub in range(len(Diag.SubDiags)):
							SubDiag = Diag.SubDiags[iSub]
							plot_style = "o"
							if(Diag.type == "SPECTRA"):
								iChan	   = iPlotChan - kPlotChan
								kPlotChan += len(SubDiag.path.Leni)
								if((kPlotChan > iPlotChan) and (iChan >= 0)):
									x_values = SubDiag.SpectrumWaves
									y_values = SubDiag.Values[iRun][:,iChan]
									SpectrumLines = SubDiag.SpectrumLines
									if(Diag.fwhm != None):
										x_values_label = np.copy(x_values)
										x_values, y_values, y_values_label = gaussfold(x_values, y_values, Diag.fwhm)
										plot_style = "-"
									else:
										x_values_label = x_values
										y_values_label = y_values
									if(not Diag.use_wavelength):
										x_values 		= 1./x_values
										x_values_label	= 1./x_values_label
							else:
								x_values = np.append(x_values, SubDiag.x_values)
								y_values = np.append(y_values, SubDiag.Values[iRun])

					if(len(path) > 1):	Ax[ip].plot(x_values, y_values, plot_style, color = colors[iRun], label=path_label[iRun])
					else:				Ax[ip].plot(x_values, y_values, plot_style, color = colors[iRun], label="{:d}".format(evolution[iRun]))
					if((Diag.type == "SPECTRA") and (Diag.label_lines != 0)):
						for kLine in range(len(SpectrumLines)):
							if(y_values_label[kLine] > Diag.label_threshold):
								Ax[ip].text(x_values_label[kLine],y_values_label[kLine],SpectrumLines[kLine], clip_on=True, rotation="vertical",horizontalalignment="center", verticalalignment="bottom")

#		Plot experiments data
#		============

			for iData in range(len(Par_diags[iParam][0])):
				HasData		= False
				HasDataErr	= False
				x_values	= np.empty(0, dtype = 'f8')
				y_values	= np.empty(0, dtype = 'f8')
				y_err		= np.empty(0, dtype = 'f8')

				iExp  = Par_diags[iParam][0][iData]
				iDiag = Par_diags[iParam][1][iData]
				Diag  = Diags[iExp][iDiag]
				for iSub in range(len(Diag.SubDiags)):
					SubDiag = Diag.SubDiags[iSub]
					if(SubDiag.data):
						HasData = True
						x_values = np.append(x_values, SubDiag.x_data)
						y_values = np.append(y_values, SubDiag.y_data)
						if(SubDiag.err_data):
							HasDataErr	= True
							y_err = np.append(y_err, SubDiag.y_err)
						else:
							y_err = np.append(y_err, np.zeros_like(SubDiag.y_data))

				if(iExp < nExpRZ):	iDescr = iExp
				else:				iDescr = iExp-nExpRZ
				if(HasData):
					if(HasDataErr):
						Ax[ip].errorbar(x_values, y_values*Diag.y_factor, yerr=y_err*Diag.y_factor, fmt=Diag.marker, color=Diag.color, label=Diag.label+" "+ExpDescrs[iDescr].label)
					else:
						Ax[ip].plot(x_values, y_values*Diag.y_factor, Diag.marker, color=Diag.color, label=Diag.label+" "+ExpDescrs[iDescr].label)

					if(save == "csv"):
						len_csv = len(csv)
						csv.append(types.SimpleNamespace())
						csv[-1].Name   = "ch_exp_"+Diag.par_name
						csv[-1].Values = np.copy(x_values)
						csv.append(types.SimpleNamespace())
						csv[-1].Name   = "exp_"+Diag.par_name
						csv[-1].Values = np.copy(y_values)
						if(HasDataErr):
							csv.append(types.SimpleNamespace())
							csv[-1].Name   = "err_"+Diag.par_name
							csv[-1].Values = np.copy(y_err)


			if((Diag.x_min == None) or (Diag.x_max == None)):
				if(Diag.x_min != None): Ax[ip].set_xlim(left=Diag.x_min)
				if(Diag.x_max != None): Ax[ip].set_xlim(right=Diag.x_min)
			else: 						Ax[ip].set_xlim(left=Diag.x_min, right=Diag.x_max)

			if((Diag.y_min == None) or (Diag.y_min == None)):
				if(Diag.y_min != None): Ax[ip].set_ylim(bottom=Diag.y_min)
				if(Diag.y_max != None): Ax[ip].set_ylim(top=Diag.y_max)
			else: 						Ax[ip].set_ylim(bottom=Diag.y_min, top=Diag.y_max)

			if((nRuns > 1) or HasData): 
				Ax[2*iDiag+1].legend(fontsize='small', loc='upper left')

	if(save != "none"):
		if(save == "csv"):
			maxLen = 0
			for ii in range(len(csv)): maxLen = max(maxLen, len(csv[ii].Values))
			save_cvs= np.zeros((maxLen,len(csv)), dtype='f8')
			Header = ""
			for ii in range(len(csv)):
				Header = Header + csv[ii].Name + ","
				save_cvs[:len(csv[ii].Values), ii] = csv[ii].Values
				if(len(csv[ii].Values) < maxLen): save_cvs[len(csv[ii].Values):, ii]  = np.nan
				
			np.savetxt("plot1d_diags_t={:.3f}_{:d}.".format(RefPar.time,i_plot_file)+save, save_cvs, header=Header, delimiter=",", fmt="%15.7e", comments="")

		else:
			for i in range(len(Fig)):
				i_plot_file += 1
				if(one_plot != 1): Fig[i].set_size_inches(10.05,7.44)
				if(save == "pdf"):
					pdf.savefig(Fig[i])
				else:
					Fig[i].savefig("plot1d_diags_t={:.3f}_{:d}.".format(RefPar.time,i_plot_file)+save)

		pyp.show(block=False)
		pyp.close()
	else:
		pyp.show()

	if(save == "pdf"):	pdf.close()

	print("plot1d_diags: Completed")

	return

#===================================================================

def indexes_upper(list_strs, str):
	return  [i for i, j in enumerate(list_strs) if j.upper() == str.upper()]

#===================================================================		
#===================================================================
def get_plasma_parameter(Plasmas, Diag, DataCenterMesh, UseTri, Compute2D):
	if(Diag.type == "RZ_POINTS"):	get_plasma_parameter_on_rz(Plasmas, Diag, DataCenterMesh, UseTri, Compute2D)
	else:								get_plasma_parameter_on_path(Plasmas, Diag, UseTri, Compute2D)

def get_plasma_parameter_on_rz(Plasmas, Diag, DataCenterMesh, UseTri, Compute2D):
	Diag.Plot2D = False

	if(Diag.y_factor  == None):	Diag.y_factor	 = 1.
	if(len(Diag.y_label) == 0):	Diag.y_label = Diag.par_name
	if(len(Diag.v_label) == 0):	Diag.v_label = Diag.par_name

	for iSub in range(len(Diag.SubDiags)):
		if(iSub+1 == len(Diag.SubDiags)): Diag.Plot2D = Compute2D
		SubDiag = Diag.SubDiags[iSub]
		SubDiag = Diag.SubDiags[iSub]

		if(DataCenterMesh.all_data == 0):
			InGrid = np.where(~DataCenterMesh.CorePath.contains_points(np.array([SubDiag.r_values,SubDiag.z_values]).T))[0]					#find points internal to the grid
			SubDiag.n_chans	 = SubDiag.n_chans[InGrid]
			SubDiag.x_values = SubDiag.x_values[InGrid]
			SubDiag.r_values = SubDiag.r_values[InGrid]
			SubDiag.z_values = SubDiag.z_values[InGrid]

			InGrid = np.where(~DataCenterMesh.CorePath.contains_points(np.array([SubDiag.r_data,SubDiag.z_data]).T))[0]					#find points internal to the grid
			SubDiag.x_data   = SubDiag.x_data[InGrid]
			SubDiag.r_data   = SubDiag.r_data[InGrid]
			SubDiag.z_data   = SubDiag.z_data[InGrid]
			SubDiag.y_data   = SubDiag.y_data[InGrid]
			if(SubDiag.err_data): SubDiag.y_err = SubDiag.y_err[InGrid]

		if(UseTri):
			if(iSub == 0):
				iPar = -1
				for i in range(len(Plasmas)):
					try:
						iPar = Plasmas[i][0].Triangles.VNames.index(Diag.par_name)
						iPlasma = i
						break
					except:
						pass

				if(iPar == -1):
					print("\tERROR: Not found ",Diag.par_name," in plasma triangles data\n\t\tcheck availables with script: print_plasma_parameters")
					Values1D = np.array([])
					Values2D = []
					return

#			Get value on nearest triangle

			Values1D = np.empty(len(SubDiag.r_values), dtype='f8')
			for i in range(len(SubDiag.r_values)):
				iTri = np.argmin((DataCenterMesh.RCent - SubDiag.r_values[i])**2 + (DataCenterMesh.ZCent - SubDiag.z_values[i])**2)
				Values1D[i] = Plasmas[iPlasma][0].Triangles.Values[iPar][iTri]*Diag.y_factor
				
			if(Diag.Plot2D):	Values2D = Plasmas[iPlasma][0].Triangles.Values[iPar]*Diag.y_factor
			else:				Values2D = 0
		else:
			if(iSub == 0):
				iPar = -1
				for i in range(len(Plasmas)):
					try:
						iPar = Plasmas[i][0].VNames.index(Diag.par_name)
						iPlasma = i
						break
					except:
						pass

				if(iPar == -1):
					print("\tERROR: Not found ",Diag.par_name," in plasma mesh data\n\t\tcheck availables with script: print_plasma_parameters")
					Values1D = np.array([])
					Values2D = []
					return

#			Get value on nearest triangle

			Values1D = np.empty(len(SubDiag.r_values), dtype='f8')
			for i in range(len(SubDiag.r_values)):
				iMesh = np.argmin((DataCenterMesh.RCent - SubDiag.r_values[i])**2 + (DataCenterMesh.ZCent - SubDiag.z_values[i])**2)
				Values1D[i] = Plasmas[iPlasma][DataCenterMesh.kMesh[iMesh]].Values[iPar][DataCenterMesh.iMesh[iMesh],DataCenterMesh.jMesh[iMesh]]

			Values1D *= Diag.y_factor
			Values2D = []
			if(Diag.Plot2D):
				for k in range(len(Plasmas[iPlasma])):
					Values2D.append(Plasmas[iPlasma][k].Values[iPar]*Diag.y_factor)

		Diag.SubDiags[iSub].Values.append(Values1D)

	Diag.Values2D = Values2D

	return


def get_plasma_parameter_on_path(Plasmas, Diag, UseTri, Compute2D):

	Diag.Plot2D = False
	for iSub in range(len(Diag.SubDiags)):
		if(iSub+1 == len(Diag.SubDiags)): Diag.Plot2D = Compute2D
		path = Diag.SubDiags[iSub].path

		if(UseTri):
			nCells = 0
			for iChan in range(len(path.Leni)): 
				for k in range(len(path.Leni[iChan])): nCells += path.Leni[iChan][k].size

			Leni = np.empty(nCells, dtype = 'f8')
			Ch   = np.empty(nCells, dtype = 'i4')
			Ti   = np.empty(nCells, dtype = 'i4')

			iCell = 0
			for iChan in range(len(path.Leni)):
				for k in range(len(path.Leni[iChan])):
					nCells = path.Leni[iChan][k].size
					Ch[iCell: iCell+nCells]	  = iChan
					Leni[iCell: iCell+nCells] = path.Leni[iChan][k]
					Ti[iCell: iCell+nCells]	  = path.Ti[iChan][k][:-1]
					iCell					 += nCells

			Coord = [Ti, Leni, Ch]
		else:
	#		Build arrays with cell intersected by lines

			nCells = 0
			for iChan in range(len(path.Leni)): nCells += path.Leni[iChan].size

			Leni = np.empty(nCells, dtype = 'f8')
			Ch   = np.empty(nCells, dtype = 'i4')
			Ki	 = np.empty(nCells, dtype = 'i4')
			Ii	 = np.empty(nCells, dtype = 'i4')
			Ji	 = np.empty(nCells, dtype = 'i4')

			iCell = 0
			for iChan in range(len(path.Leni)):
				nCells = path.Leni[iChan].size

				Leni[iCell: iCell+nCells] = path.Leni[iChan]
				Ch[iCell: iCell+nCells]	  = iChan
				Ki[iCell: iCell+nCells]	  = path.Ki[iChan]
				Ii[iCell: iCell+nCells]	  = path.Ii[iChan]
				Ji[iCell: iCell+nCells]	  = path.Ji[iChan]
				iCell					 += nCells

			Coord = [Ki, Ii, Ji, Leni, Ch]


		if(Diag.type == "LINES"):
			if(Diag.par_name   == "IntNe"): 		Values1D, Values2D = get_diag_values_IntNe(Plasmas, Diag, iSub, Coord, UseTri, Diag.Plot2D)
			elif(Diag.par_name == "BOLO"):			Values1D, Values2D = get_diag_values_BOLO(Plasmas, Diag, iSub, Coord, UseTri, Diag.Plot2D)
			elif((len(Diag.par_name[:6]) > 5) and (Diag.par_name[:6].upper() == "SPECT_")): Values1D, Values2D = get_diag_values_LINE(Plasmas, Diag, iSub, Coord, UseTri, Diag.Plot2D)
			Diag.SubDiags[iSub].Values.append(np.zeros(len(path.Leni), dtype="f8"))

			for iChan in range(len(path.Leni)):
				iCh = np.where(Ch == iChan)[0]
				if(len(iCh) > 0): Diag.SubDiags[iSub].Values[-1][iChan] += np.sum(Leni[iCh]*Values1D[iCh])
		else:
			Diag.Plot2D = False
			Values2D = get_diag_values_SPECTRUM(Plasmas, Diag, iSub, Coord, UseTri, Diag.Plot2D)

	Diag.Values2D = Values2D

	return 

#==========================================
# diagnostic section
#==========================================

#--------------------------
#	Electron density
#--------------------------

def get_diag_values_IntNe(Plasmas, Diag, iSub, Coords, UseTri, Get2D):

	print("\t\tget_diag_values_IntNe")

	if(iSub == 0):
		if(Diag.y_factor  == None):	Diag.y_factor	 = 1e-20
		else:						Diag.y_factor	 = 10**(int(log10(Diag.y_factor)))
		if(len(Diag.y_label) == 0):	Diag.y_label = "$n_e (*10^{{{:d}}}\ m^{-2})$".format(int(-log10(Diag.y_factor)))
		if(len(Diag.v_label) == 0):	Diag.v_label = "$n_e (*10^{{{:d}}}\ m^{-3})$".format(int(-log10(Diag.y_factor)))

	if(UseTri):
		Ti = Coords[0]
		if(Get2D):	Values2D = Plasmas[0][0].Triangles.Dens*Diag.y_factor
		else:		Values2D = 0
		Values1D = Plasmas[0][0].Triangles.Dens[Ti]*Diag.y_factor
	else:
		Ki = Coords[0]
		Ii = Coords[1]
		Ji = Coords[2]
		Values2D = []
		Values1D = np.empty(len(Coords[0]), dtype = 'f8')
		for k in range(len(Plasmas[0])):
			iK = np.where(Ki == k)[0]
			if(Get2D):	   	 Values2D.append(Plasmas[0][k].Dens*Diag.y_factor)
			if(len(iK) > 0): Values1D[iK] = Plasmas[0][k].Dens[Ii[iK], Ji[iK]]*Diag.y_factor

	return Values1D, Values2D

#--------------------------
#	Bolometer
#--------------------------

def get_diag_values_BOLO(Plasmas, Diag, iSub, Coords, UseTri, Get2D):

	print("\t\tget_diag_values_BOLO")
	if(iSub == 0):
		if(Diag.y_factor  == None):	Diag.y_factor	 = 1e-6
		else:						Diag.y_factor	 = 10**(int(log10(Diag.y_factor)))
		if(len(Diag.y_label) == 0):	Diag.y_label = "$Rad (*10^{{{:d}}}\ W/m^2)$".format(int(-log10(Diag.y_factor)))
		if(len(Diag.v_label) == 0):	Diag.v_label = "$Rad (*10^{{{:d}}}\ W/m^3)$".format(int(-log10(Diag.y_factor)))

	if(UseTri):
		Ti = Coords[0]
		if(Get2D): Values2D = Plasmas[0][0].Triangles.TRad*Diag.y_factor
		else:	   Values2D = 0
		Values1D = Plasmas[0][0].Triangles.TRad[Ti]*Diag.y_factor
	else:
		Ki = Coords[0]
		Ii = Coords[1]
		Ji = Coords[2]
		Values1D = np.empty(len(Ki), dtype = 'f8')
		Values2D = []
		for k in range(len(Plasmas[0])):
			if(Get2D):	Values2D.append(Plasmas[0][k].TRad*Diag.y_factor)
			iK = np.where(Ki == k)[0]
			if(len(iK) > 0): Values1D[iK] = Plasmas[0][k].TRad[Ii[iK], Ji[iK]]*Diag.y_factor

	return Values1D, Values2D


#--------------------------
#	Spectroscopy
#--------------------------

"""
JET available

BE1A (457nm line intensity from Be1)
BE1B (441 nm...)
BE2A (527 nm line intensity from Be2)
BE2B (436 nm...)
BE2C (467 nm...)
BE2D (483 nm...)
BEMO (line intensity from molecular Be)
C2A (427 nm line from C)
C2B (515 nm...)
C2C (589 nm...)
C3A (465 nm...)
C3B (570 nm...)
C4A (581 nm...)
C4B (502 nm...)
C5A (494 nm...)
C6A (529 nm...)bachmann
N2B (500 nm line intensity from N)
N2C (567 nm...)
N2D (452 nm...)
"""

def get_diag_values_LINE(Plasmas, Diag, iSub, Coords, UseTri, Get2D):


	ADF15_Folder_1 = ["pec96#h", "pec96#he", "pec96#li", "pec96#be", "pec93#b", \
				      "pec96#c", "pec96#n",  "pec96#o",   "pec96#ne", "pec40#ar"]
	ADF15_Folder_2 = ["", "", "", "", "", \
				      "pec93#c", "",  "",   "", ""]


	IONS_NAME  = ["H","D","HE","LI","BE","B", "C","N","O","NE","AR"]
	IONS_POS   = [0,      0,    1,   2,    3,   4,   5,   6,   7,     8,    9]
	IONS_IONIZ = [        1,    2,   3,    4,   5,   6,   7,   8,    10,   18]

	IONIZATION_NAMES1 = []
	for i in range(1,31): IONIZATION_NAMES1.append("{:d}".format(i))
	IONIZATION_NAMES2 = ["I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV", \
						"XVI","XVII","XVIII","XIX","XX","XXI","XXII","XXIII","XXIV","XXV","XXVI","XXVII", \
						"XXVIII","XXIX","XXX"]

	line_name = Diag.par_name[6:].upper()+"                    "

#	Find ion name from line emission name

	iIon = -1
	for i in range(len(IONS_NAME)):
		if(line_name[:len(IONS_NAME[i])] == IONS_NAME[i]):
			iPos = i
			iIon = IONS_POS[i]
			break

	if(iIon == -1):
		print("\n\nWARNING: Not found ion file for line emission ",Diag.par_name[6:],"\n\n")
		return [],[]


#	Find plasma ion correspondin to requested line emission

	LineAtom = line_name[:len(IONS_NAME[iPos])]
	nPlasmas   = len(Plasmas)
	PlasmaIons = [Plasmas[i][0].ion for i in range(nPlasmas)]

	iPlasmaIon1 = -1
	for i in range(1,nPlasmas):
		PlasmaAtom = PlasmaIons[i][:len(IONS_NAME[iPos])].upper()
		if((LineAtom == PlasmaAtom) or ((LineAtom == "H") and (PlasmaAtom == "D"))):
			iPlasmaIon1 = i												#first ionization stage
			break

	if(iPlasmaIon1 == -1):
		print("\n\nWARNING: Not found plasma file for line emission ",Diag.par_name[6:],"\n\n")
		return [],[]

#	Find requested ionization stage

	ioniz_name = line_name[len(IONS_NAME[iPos]):]
	if(iIon > 0):
		ioniz = -1
		for i in range(len(IONS_NAME)-1,-1,-1):
			if((ioniz_name[:len(IONIZATION_NAMES1[i])] == IONIZATION_NAMES1[i]) or (ioniz_name[:len(IONIZATION_NAMES2[i])] == IONIZATION_NAMES2[i])):
				ioniz = i
				break
		if(ioniz == -1):
			print("\n\nWARNING: Not found ionization stage for line emission ",Diag.par_name[6:],"\n\n")
			exit()
	else:
		LineAtom = "H"
		ioniz = 0

#	get requested line wavelength

	iWave = ioniz_name.find("-")
	OneLine = True
	if((iIon == 0) and (ioniz_name[iWave+1:iWave+6].upper() == "ALPHA")):
		Wave	= 656.2
	else:
		Wave  = eval(ioniz_name[iWave+1:iWave+5])/10.						#Wavelenght from A to nm
		if((len(ioniz_name) > iWave+5) and (ioniz_name[iWave+5] == "-")):
			Wave_Min = Wave
			Wave_Max = eval(ioniz_name[iWave+6:iWave+10])/10.
			OneLine	 = False

#	define available adas files

	ADF15_Files = [ADF15_Folder_1[iIon]+"/"+ADF15_Folder_1[iIon]+"_pju#"+LineAtom.lower()+"{:d}.dat".format(ioniz),
				   ADF15_Folder_1[iIon]+"/"+ADF15_Folder_1[iIon]+"_vsu#"+LineAtom.lower()+"{:d}.dat".format(ioniz)]
	if(len(ADF15_Folder_2[iIon]) > 0): ADF15_Files.append(ADF15_Folder_2[iIon]+"/"+ADF15_Folder_2[iIon]+"_pju#"+LineAtom.lower()+"{:d}.dat".format(ioniz))

#	Find first adas file with requested wavelength

	iWaves = []
	for ADF15_File in ADF15_Files:
		Adas = load_adas_15(ADF15_File)
		if(Adas.has_data):
			if(OneLine):
				iWave  = np.argmin(np.abs(Adas.wavels - Wave*10))
				if(abs(Adas.wavels[iWave] - Wave*10) < 5): iWaves = np.where(Adas.wavels == Adas.wavels[iWave])[0]
			else:
				iWaves = np.where((Adas.wavels >= Wave_Min) & (Adas.wavels <= Wave_Max))[0]
			if(len(iWaves) != 0): break


	if(len(iWaves) == 0):
		if(OneLine):
			print("\n\nWARNING: Not found wavelength for line emission ",Diag.par_name[6:])
			print("         requested wavelength =",Wave*10,"\n\n")
			exit()
		else:
			print("\n\nWARNING: Not found wavelength for line emission ",Diag.par_name[6:])
			print("         min requested wavelength =",Wave_Min*10)
			print("         max requested wavelength =",Wave_Max*10,"\n\n")
			exit()

		if(Adas.has_data):
			print("\n\n        available wavelengths    =",np.unique(np.sort(Adas.wavels))*10,"\n\n")
		exit()
#		return [],[]

	if(iSub == 0):
		if(OneLine): line_label = "I_{"+"{:.1f}".format(Adas.wavels[iWave]/10.)+"}"
		else:		 line_label	= "I_{"+"{:.1f}-{:.1f}".format(np.min(Adas.wavels[iWaves])/10.,np.max(Adas.wavels[iWaves])/10.)+"}"
		if(Diag.y_factor  == None):	Diag.y_factor	 = 1.
		else:						Diag.y_factor	 = 10**(int(log10(Diag.y_factor)))
		if(len(Diag.y_label) == 0):	Diag.y_label = "$"+line_label+"\ (ph*s^{-1}*m^{-2}*sr^{-1})$"
		if(len(Diag.v_label) == 0):	Diag.v_label = "$"+line_label+"\ (ph*s^{-1}*m^{-3}*sr^{-1})$"

	ValuesFact			 = Diag.y_factor/(4.*np.pi)										#I am  not sure about to divide or not by 4*pi


#	Ready to compute emission

	nWaves		= len(iWaves)
	FitPecs		= []			#fit pec on log log scale
	TransTypes	= []			#Transaction types [excit,  recom, chexc, ioniz]
	MinNe		= np.empty(nWaves, dtype='f8')
	MaxNe		= np.empty(nWaves, dtype='f8')
	MinTe		= np.empty(nWaves, dtype='f8')
	MaxTe		= np.empty(nWaves, dtype='f8')
	for i in range(nWaves):
		iw = iWaves[i]
		TransTypes.append(Adas.blocks[iw].ctype)

		LogNe = np.log(Adas.blocks[iw].teda)
		LogTe = np.log(Adas.blocks[iw].teta)
		FitPecs.append(interpolate.RectBivariateSpline(LogNe,LogTe,Adas.blocks[iw].pec))

		MinNe[i] = np.min(Adas.blocks[iw].teda)
		MaxNe[i] = np.max(Adas.blocks[iw].teda)
		MinTe[i] = np.min(Adas.blocks[iw].teta)
		MaxTe[i] = np.max(Adas.blocks[iw].teta)

	iPlasmaIon = iPlasmaIon1 + ioniz - 1
	if(UseTri):
		Ti = Coords[0]
		eTriangles = Plasmas[0][0].Triangles
		if(Get2D):	
			Values2D = np.zeros_like(eTriangles.Dens)
			for i in range(nWaves):
				LogTe  = np.log(np.where(np.where(eTriangles.Temp > MinTe[i], eTriangles.Temp, MinTe[i]) < MaxTe[i], \
																				  eTriangles.Temp, MaxTe[i]))
				LogNe  = np.log(np.where(np.where(eTriangles.Dens > MinTe[i], eTriangles.Dens, MinNe[i]) < MaxNe[i], \
																				  eTriangles.Dens, MaxNe[i]))
				Pec = FitPecs[i].ev(LogNe,LogTe)
				if(TransTypes[i] == "excit"):
					if(ioniz == 0):
						Values2D += eTriangles.Dens*Plasmas[iPlasmaIon+1][0].Triangles.Nn*Pec
					else:
						Values2D += eTriangles.Dens*Plasmas[iPlasmaIon][0].Triangles.Dens*Pec
				elif(TransTypes[i] == "recom"):
					Values2D += eTriangles.Dens*Plasmas[iPlasmaIon+1][0].Triangles.Dens*Pec
				elif(TransTypes[i] == "chexc"):
					if(ioniz == 1):
						Values2D += eTriangles.Dens*Plasmas[iPlasmaIon][0].Triangles.Nn*Pec
					else:
						Values2D += eTriangles.Dens*Plasmas[iPlasmaIon-1][0].Triangles.Dens*Pec
				elif(TransTypes[i] == "ioniz"):
					Values2D += eTriangles.Dens*Plasmas[iPlasmaIon][0].Triangles.Dens*Pec

			Values2D  = np.where(Values2D > 0., Values2D*ValuesFact, 0.)

		else:
			Values2D = 0

		Values1D = np.zeros(Ti.shape, dtype='f8')
		for i in range(nWaves):
			LogTe  = np.log(np.where(np.where(eTriangles.Temp[Ti] > MinTe[i], eTriangles.Temp[Ti], MinTe[i]) < MaxTe[i], \
																		  eTriangles.Temp[Ti], MaxTe[i]))
			LogNe  = np.log(np.where(np.where(eTriangles.Dens[Ti] > MinTe[i], eTriangles.Dens[Ti], MinNe[i]) < MaxNe[i], \
																		  eTriangles.Dens[Ti], MaxNe[i]))

			Pec = FitPecs[i].ev(LogNe,LogTe)
			if(TransTypes[i] == "excit"):
				if(ioniz == 0):
					Values1D += eTriangles.Dens[Ti]*Plasmas[iPlasmaIon+1][0].Triangles.Nn[Ti]*Pec
				else:
					Values1D += eTriangles.Dens[Ti]*Plasmas[iPlasmaIon][0].Triangles.Dens[Ti]*Pec
			elif(TransTypes[i] == "recom"):
				Values1D += eTriangles.Dens[Ti]*Plasmas[iPlasmaIon+1][0].Triangles.Dens[Ti]*Pec
			elif(TransTypes[i] == "chexc"):
				if(ioniz == 1):
					Values1D += eTriangles.Dens[Ti]*Plasmas[iPlasmaIon+1][0].Triangles.Nn[Ti]*Pec
				else:
					Values1D += eTriangles.Dens[Ti]*Plasmas[iPlasmaIon-1][0].Triangles.Dens[Ti]*Pec
			elif(TransTypes[i] == "ioniz"):
				Values1D += eTriangles.Dens[Ti]*Plasmas[iPlasmaIon][0].Triangles.Dens[Ti]*Pec

		Values1D  = np.where(Values1D > 0., Values1D*ValuesFact, 0.)

	else:
		Ki = Coords[0]
		Ii = Coords[1]
		Ji = Coords[2]

		Values1D = np.zeros(len(Ki), dtype = 'f8')
		Values2D = []
		for k in range(len(Plasmas[0])):
			if(Get2D):
				Values2D.append(np.zeros_like(Plasmas[0][k].Temp))

				for i in range(nWaves):
					TeLow = np.where(Plasmas[0][k].Temp > MinTe[i], Plasmas[0][k].Temp, MinTe[i]+0.01*(MaxTe[i]-MinTe[i]))
					NeLow = np.where(Plasmas[0][k].Dens > MinTe[i], Plasmas[0][k].Dens, MinNe[i]+0.01*(MaxNe[i]-MinNe[i]))
					LogTe  = np.log(np.where(TeLow < MaxTe[i], TeLow, MaxTe[i]-0.01*(MaxTe[i]-MinTe[i])).reshape(-1))
					LogNe  = np.log(np.where(NeLow < MaxNe[i], NeLow, MaxNe[i]-0.01*(MaxNe[i]-MinNe[i])).reshape(-1))

					Pec = FitPecs[i].ev(LogNe,LogTe).reshape(Plasmas[0][k].Temp.shape)
					if(TransTypes[i] == "excit"):
						if(ioniz == 0):
							Values2D[-1] += Plasmas[0][k].Dens*Plasmas[iPlasmaIon+1][k].Dens0*Pec
						else:
							Values2D[-1] += Plasmas[0][k].Dens*Plasmas[iPlasmaIon][k].Dens*Pec
					elif(TransTypes[i] == "recom"):
						Values2D[-1] += Plasmas[0][k].Dens*Plasmas[iPlasmaIon+1][k].Dens*Pec
					elif(TransTypes[i] == "chexc"):
						Values2D[-1] += Plasmas[1][k].Dens0*Plasmas[iPlasmaIon+1][k].Dens*Pec
					elif(TransTypes[i] == "ioniz"):
						Values2D[-1] += Plasmas[0][k].Dens*Plasmas[iPlasmaIon][k].Dens*Pec

				Values2D[-1]  = np.where(Values2D[-1] > 0., Values2D[-1]*ValuesFact, 0.)

			iK = np.where(Ki == k)[0]
			if(len(iK) > 0):
				for i in range(nWaves):
					iw = iWaves[i]
					TeLow   = np.where(Plasmas[0][k].Temp[Ii[iK], Ji[iK]] > MinTe[i], Plasmas[0][k].Temp[Ii[iK], Ji[iK]] , MinTe[i]+0.01*(MaxTe[i]-MinTe[i]))
					NeLow = np.where(Plasmas[0][k].Dens[Ii[iK], Ji[iK]] > MinNe[i], Plasmas[0][k].Dens[Ii[iK], Ji[iK]] , MinNe[i]+0.01*(MaxNe[i]-MinNe[i]))
					LogTe  = np.log(np.where(TeLow < MaxTe[i], TeLow, MaxTe[i]-0.01*(MaxTe[i]-MinTe[i])))
					LogNe  = np.log(np.where(NeLow < MaxNe[i], NeLow, MaxNe[i]-0.01*(MaxNe[i]-MinNe[i])))

					Pec = FitPecs[i].ev(LogNe,LogTe)
					if(TransTypes[i] == "excit"):
						if(ioniz == 0):
							Values1D[iK] += Plasmas[0][k].Dens[Ii[iK], Ji[iK]]*Plasmas[iPlasmaIon+1][k].Dens0[Ii[iK], Ji[iK]]*Pec
						else:
							Values1D[iK] += Plasmas[0][k].Dens[Ii[iK], Ji[iK]]*Plasmas[iPlasmaIon][k].Dens[Ii[iK], Ji[iK]]*Pec
					elif(TransTypes[i] == "recom"):
						Values1D[iK] += Plasmas[0][k].Dens[Ii[iK], Ji[iK]]*Plasmas[iPlasmaIon+1][k].Dens[Ii[iK], Ji[iK]]*Pec
					elif(TransTypes[i] == "chexc"):
						Values1D[iK] += Plasmas[1][k].Dens0[Ii[iK], Ji[iK]]*Plasmas[iPlasmaIon+1][k].Dens[Ii[iK], Ji[iK]]*Pec
					elif(TransTypes[i] == "ioniz"):
						Values1D[iK] += Plasmas[0][k].Dens[Ii[iK], Ji[iK]]*Plasmas[iPlasmaIon][k].Dens[Ii[iK], Ji[iK]]*Pec

				Values1D  = np.where(Values1D > 0., Values1D*ValuesFact, 0.)

	return Values1D, Values2D

#	###################################
#	Compute full spectrum available from one or more adas files
#	####################################

def get_diag_values_SPECTRUM(Plasmas, Diag, iSub, Coords, UseTri, Get2D):


	IONS_NAME  = ["H","D","He","Li","Be","B", "C","N","O","Ne","Ar"]
	IONS_POS   = [0,      0,    1,   2,    3,   4,   5,   6,   7,     8,    9]
	IONS_IONIZ = [        1,    2,   3,    4,   5,   6,   7,   8,    10,   18]

	IONIZATION_NAMES1 = []
	for i in range(1,31): IONIZATION_NAMES1.append("{:d}".format(i))
	IONIZATION_NAMES2 = ["I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV", \
						"XVI","XVII","XVIII","XIX","XX","XXI","XXII","XXIII","XXIV","XXV","XXVI","XXVII", \
						"XXVIII","XXIX","XXX"]
	if(iSub == 0):
		if(Diag.y_factor  == None):	Diag.y_factor	 = 1.
		else:						Diag.y_factor	 = 10**(int(log10(Diag.y_factor)))
		if(len(Diag.y_label) == 0):	Diag.y_label = "$\ (ph*s^{-1}*m^{-2}*sr^{-1})$"
		if(len(Diag.v_label) == 0):	Diag.v_label = "$\ (ph*s^{-1}*m^{-3}*sr^{-1})$"

	ValuesFact    = Diag.y_factor/(4.*np.pi)										#I am  not sure about to divide or not by 4*pi
	path 	   = Diag.SubDiags[iSub].path

	nPlasmas   = len(Plasmas)
	PlasmaIons = [Plasmas[i][0].ion for i in range(nPlasmas)]

	SpectrumLines  = []
	SpectrumWaves  = np.empty((0,len(path.Leni)),dtype='f8')
	SpectrumValues = np.empty((0,len(path.Leni)),dtype='f8')
	for iAdas in range(len(Diag.adas_files)):												#Loop on selected adas folder
		adas_file	= Diag.adas_files[iAdas]												#files whitout extension

		ion_name	= adas_file[::-1]														#get ion name from adas file
		ion_name	= ion_name[ion_name.find("#")-1::-1]
		iStar		= ion_name.find("*")
		if(iStar > -1):																		#all ionizations are requested
			AllIoniz  = True
			ion_name = ion_name[:iStar].upper()
		else:																				#single ionization requested
			AllIoniz  = False
			LenName = 2
			if(ord(ion_name[1]) < 65): LenName = 1
			if((ion_name[LenName] == "0") and (len(ion_name) > LenName+1)):	adas_ioniz_values = [eval(ion_name[LenName+1:])]
			else:																adas_ioniz_values = [eval(ion_name[LenName:])]
			ion_name = ion_name[:LenName].upper()
			adas_file_ions = [adas_file+".dat"]

		iIon = -1
		for kIon in range(len(IONS_NAME)):
			if(ion_name == IONS_NAME[kIon].upper()):
				iIon = IONS_POS[kIon]
				if(AllIoniz):
					adas_file_ions = []
					adas_ioniz_values   = []
					for k in range(IONS_IONIZ[iIon]):
						adas_ioniz_values.append(k)
						adas_file_ions.append(adas_file[:-1]+"{:d}.dat".format(k))
				break

		if(iIon == -1):
			print("\n\nWARNING: Unknow ion ",ion_name,"\n\n")
			return [],[]

		iPlasmaIon1 = -1																#Find plasma files
		for i in range(1,nPlasmas):
			if(ord(PlasmaIons[i][1]) < 65): PlasmaAtom = PlasmaIons[i][:1].upper()
			else:							PlasmaAtom = PlasmaIons[i][:2].upper()
			if((ion_name == PlasmaAtom) or ((ion_name == "H") and (PlasmaAtom == "D"))):
				iPlasmaIon1 = i															#first ionization stage
				break
		if(iPlasmaIon1 == -1):
			print("\n\nWARNING: Not found plasma file for ion ",ion_name,"\n\n")
			return [],[]

		for iAdasIon in range(len(adas_file_ions)):										#Loop on ions (or single Ion)
			Adas  = load_adas_15(adas_file_ions[iAdasIon])
			if(not Adas.has_data):
				print("\tNo data for file = ",adas_file_ions[iAdasIon])
				continue

			ioniz	   = adas_ioniz_values[iAdasIon]
			iPlasmaIon = iPlasmaIon1 + ioniz - 1
			AdasLines  = []
			AdasWaves  = np.empty((0),dtype='f8')
			AdasValues = np.empty((0,len(path.Leni)),dtype='f8')
			ChanValues = np.zeros((1,len(path.Leni)),dtype='f8')
			for iw in range(len(Adas.wavels)):
				if(((Diag.w_min != None) and (Diag.w_min > Adas.wavels[iw])) or
				   ((Diag.w_max != None) and (Diag.w_max < Adas.wavels[iw]))): continue 

				TransTypes = Adas.blocks[iw].ctype														#Transaction types [excit,  recom, chexc, ioniz]

				LogNe = np.log(Adas.blocks[iw].teda)
				LogTe = np.log(Adas.blocks[iw].teta)
				FitPecs =	interpolate.RectBivariateSpline(LogNe,LogTe,Adas.blocks[iw].pec)			#fit pec on log log scale

				MinNe = np.min(Adas.blocks[iw].teda)
				MaxNe = np.max(Adas.blocks[iw].teda)
				MinTe = np.min(Adas.blocks[iw].teta)
				MaxTe = np.max(Adas.blocks[iw].teta)

				if(UseTri):
					Ti	 = Coords[0]
					Leni = Coords[1]
					Ch	 = Coords[2]
					eTriangles = Plasmas[0][0].Triangles

					LogTe  = np.log(np.where(np.where(eTriangles.Temp[Ti] > MinTe, eTriangles.Temp[Ti], MinTe) < MaxTe, \
																				   eTriangles.Temp[Ti], MaxTe))
					LogNe  = np.log(np.where(np.where(eTriangles.Dens[Ti] > MinTe, eTriangles.Dens[Ti], MinNe) < MaxNe, \
																				   eTriangles.Dens[Ti], MaxNe))
					Pec = FitPecs.ev(LogNe,LogTe)
					if(TransTypes == "excit"):
						if(ioniz == 0):
							Values1D = eTriangles.Dens[Ti]*Plasmas[iPlasmaIon+1][0].Triangles.Nn[Ti]*Pec
						else:
							Values1D = eTriangles.Dens[Ti]*Plasmas[iPlasmaIon][0].Triangles.Dens[Ti]*Pec
					elif(TransTypes == "recom"):
						Values1D = eTriangles.Dens[Ti]*Plasmas[iPlasmaIon+1][0].Triangles.Dens[Ti]*Pec
					elif(TransTypes == "chexc"):
						Values1D = Plasmas[1][0].Triangles.Nn[Ti]*Plasmas[iPlasmaIon+1][0].Triangles.Dens[Ti]*Pec
					elif(TransTypes == "ioniz"):
						Values1D = eTriangles.Dens[Ti]*Plasmas[iPlasmaIon][0].Triangles.Dens[Ti]*Pec

					Values1D  = np.where(Values1D > 0., Values1D*ValuesFact, 0.)

				else:
					Ki	 = Coords[0]
					Ii	 = Coords[1]
					Ji	 = Coords[2]
					Leni = Coords[3]
					Ch	 = Coords[4]

					Values1D = np.zeros(len(Ki), dtype = 'f8')
					for k in range(len(Plasmas[0])):
						iK = np.where(Ki == k)[0]
						if(len(iK) > 0):
							TeLow  = np.where(Plasmas[0][k].Temp[Ii[iK], Ji[iK]] > MinTe, Plasmas[0][k].Temp[Ii[iK], Ji[iK]] , MinTe+0.01*(MaxTe-MinTe))
							NeLow  = np.where(Plasmas[0][k].Dens[Ii[iK], Ji[iK]] > MinNe, Plasmas[0][k].Dens[Ii[iK], Ji[iK]] , MinNe+0.01*(MaxNe-MinNe))
							LogTe  = np.log(np.where(TeLow < MaxTe, TeLow, MaxTe-0.01*(MaxTe-MinTe)))
							LogNe  = np.log(np.where(NeLow < MaxNe, NeLow, MaxNe-0.01*(MaxNe-MinNe)))

							Pec = FitPecs.ev(LogNe,LogTe)
							if(TransTypes == "excit"):
								if(ioniz == 0):
									Values1D[iK] += Plasmas[0][k].Dens[Ii[iK], Ji[iK]]*Plasmas[iPlasmaIon+1][k].Dens0[Ii[iK], Ji[iK]]*Pec
								else:
									Values1D[iK] += Plasmas[0][k].Dens[Ii[iK], Ji[iK]]*Plasmas[iPlasmaIon][k].Dens[Ii[iK], Ji[iK]]*Pec
							elif(TransTypes == "recom"):
								Values1D[iK] += Plasmas[0][k].Dens[Ii[iK], Ji[iK]]*Plasmas[iPlasmaIon+1][k].Dens[Ii[iK], Ji[iK]]*Pec
							elif(TransTypes == "chexc"):
								Values1D[iK] += Plasmas[1][k].Dens0[Ii[iK], Ji[iK]]*Plasmas[iPlasmaIon+1][k].Dens[Ii[iK], Ji[iK]]*Pec
							elif(TransTypes == "ioniz"):
								Values1D[iK] += Plasmas[0][k].Dens[Ii[iK], Ji[iK]]*Plasmas[iPlasmaIon][k].Dens[Ii[iK], Ji[iK]]*Pec

						Values1D  = np.where(Values1D > 0., Values1D*ValuesFact, 0.)

				for iChan in range(len(path.Leni)):
					iCh = np.where(Ch == iChan)[0]
					if(len(iCh) > 0): ChanValues[0,iChan] = np.sum(Leni[iCh]*Values1D[iCh])

				if(iw != 0): iwPre = np.where(AdasWaves == Adas.wavels[iw])[0] 									#Wavelength is already present
				else:		 iwPre = []
				if(len(iwPre) == 0):
					if(Diag.label_lines == 0):		AdasLines.append("")
					elif(Diag.label_lines == 1):	AdasLines.append(IONS_NAME[kIon]+IONIZATION_NAMES2[adas_ioniz_values[iAdasIon]])
					else:							AdasLines.append(IONS_NAME[kIon]+IONIZATION_NAMES2[adas_ioniz_values[iAdasIon]]+" {:0.1f}".format(Adas.wavels[iw]))
					AdasWaves  = np.append(AdasWaves, Adas.wavels[iw])
					AdasValues = np.append(AdasValues, ChanValues, axis=0)
				else:
					iwPre = iwPre[0]
					AdasValues[iwPre,:] += ChanValues[0,:]
			if(Diag.threshold_max > 0.):
				i_above = np.where(np.max(AdasValues,axis=1) > Diag.threshold_max)[0]
				AdasWaves  = AdasWaves[i_above]
				AdasLines  = [AdasLines[iline] for iline in i_above]
				AdasValues = AdasValues[i_above,:]
			elif(Diag.threshold_min > 0.):
				i_above = np.where(np.min(AdasValues,axis=1) > Diag.threshold_min)[0]
				AdasWaves  = AdasWaves[i_above]
				AdasLines  = [AdasLines[iline] for iline in i_above]
				AdasValues = AdasValues[i_above,:]

			SpectrumWaves = np.append(SpectrumWaves,AdasWaves)
			SpectrumLines.extend(AdasLines)
			SpectrumValues = np.append(SpectrumValues, AdasValues, axis=0)

	Diag.SubDiags[iSub].SpectrumWaves = SpectrumWaves					#Same for all simulation (Amstrong)
	Diag.SubDiags[iSub].SpectrumLines = SpectrumLines					#
	Diag.SubDiags[iSub].Values.append(SpectrumValues)					#Append to manage more than one simulation
	Values2D = 0

	return Values2D

#===================================================================0
# NAME:		gaussfold
# PURPOSE:	Smoothes a plot by convolving with a Gaussian profile.
#			Main purpose is to convolve a spectrum (flux against wavelength) with a given instrument resolution.
#			Also applicable e.g. to smooth ligthcurves or in time-series analysis. 
#
# CALLING SEQUENCE:
#		smoothedLambda, smoothedFlux,  smoothedFluxOnLambda= gaussfold(Lambda, Flux, fwhm, lammin=lammin, lammin=lammax) 
#
# INPUTS:	Lambda	= In ascending order sorted array containing the values of the x-axis (e.g. the wavelength- or frequencygrid).
#						Irregularly spaced grids are acepted.
#			Llux		= Array ( same size as lambda containing the values of the y-axis (e.g. the flux).
# 			fwhm	= FWHM in units of the x-axis (e.g. Angstroem) of the Gaussian profile.
#
# KEYWORD PARAMETERS:
#			lammin &
#			lammax   = Defines the x-axis range. All y-data within this range will be smoothed.CAUTION: improtant in case of large arrays (memory!)
#			DEFAULTS:	lammin = MIN(lambda) - fwhm 
# 						lammax = MAX(lambda) + fwhm 
# OUTPUTS:
#			smoothedLambda			= New array of Lambda regularly spaced 
#			smoothedFlux				= New smoothed flux on regularly spaced smoothedLambda
#			smoothedFluxOnLambda	= Array same size as Lambda  containing the smoothed flux at Lambda values

from math import sqrt, log, exp

def gaussfold(lam, flux, fwhm, lammin=None, lammax=None):
	if(lammin == None): lammin = np.min(lam) - fwhm
	if(lammax == None): lammax = np.max(lam) + fwhm

	dlambda   = fwhm/17.
	smoothedLambda  = lammin + dlambda*np.arange(int((lammax-lammin)/dlambda+1))

#	set available values in thinner regular grid

	smoothedFlux = np.zeros_like(smoothedLambda)
	ismoothedLambda = ((lam-lammin)/dlambda+0.5).astype(int)
	iOk = np.where((ismoothedLambda > -1) & (ismoothedLambda < len(smoothedLambda)))[0]
	smoothedFlux[ismoothedLambda[iOk]] = flux[iOk]

	fwhm_pix = fwhm/dlambda
	window	 = int(17*fwhm_pix)

#	## get a 1D gaussian profile
	sigma   = fwhm_pix/(2.0*sqrt(2.0*log(2.0)))
	x_gauss = np.arange(window, dtype='f8') - (window-1)/2.
	gauss	= np.exp(-np.power(x_gauss, 2.)/(2*np.power(sigma, 2.)))
	gauss	= gauss/np.sum(gauss)

#	## convolve input spectrum with the gauss profile
	smoothedFlux = np.convolve(smoothedFlux, gauss)[int((window-1)/2.):len(smoothedLambda)+int((window-1)/2.)]

#	Return data also at the original wavelegths
	smoothedFluxOnLambda  = np.interp(lam, smoothedLambda, smoothedFlux)

	return smoothedLambda, smoothedFlux, smoothedFluxOnLambda
