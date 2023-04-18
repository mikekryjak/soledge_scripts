#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot1d_sp_on_wall 			import plot1d_sp_on_wall
	from routines.cli_routines	import *

	if(cli_present("-h",sys.argv) and cli_present("-exp_file",sys.argv) ):
		from routines.exp_data_routines			import get_exp_data_point_help
		from files.load_exp_descr						import load_exp_descr_help
		load_exp_descr_help()
		get_exp_data_point_help()
		exit()

	if(cli_present("-h",sys.argv)):
		print("\nThis script plot plasma parameters on strike points\n")
		print("plot1d_sp_on_wall options")
		print("\t-path        Directory with simulation or list of directories[d='./']")
		print("\t-evolution   If > 0 plot evolution file (or evolutions files with list [3,4,5])  [d=0]")
		print("\t-mod_file    Text (csv, or tsv) file with model data to compare [d='']")
		print("\t-exp_file    Exp file with description of diagnostic data [d='']")
		print("\t-shot        Shot to compare [d=0]")
		print("\t-tstart      Start time experimental data to compare [d=0.]")
		print("\t-tend        End time experimental data to compare [d=0.]")
		print("\t-deltas      Amplitude around strike points [d=0.05]")
		print("\t-d_from_sep  Max distance from separatrix, if >0 use d_from_sep instead then deltas [d=0.0]")
		print("\t-path_label  Labels to itentify runs in plot [d='']")
		print("\t-no_labels   Skip labels in plots [d=false]")
		print("\t-log_scale   Use log scale for y axis [d=false]")
		print("\t-diff        Plot difference between two evolutions[d=false]")
		print("\t-extra_walls Flag to show extra walls [d=false]")
		print("\t-one_plot    One plot on each figure [d=false]")
		print("\t-save        Save figures on files, values=none/png/ps/eps/pdf/cvs/stat/stat_norm [D='none']")
		print()
		exit()

	path	 	= cli_get_value("-path",				sys.argv,   [""])
	evolution	= cli_get_value("-evolution",		sys.argv,	[])
	mod_file	= cli_get_value("-mod_file",			sys.argv,"")
	exp_files	= cli_get_value("-exp_file",			sys.argv,[""])
	shot		= cli_get_value("-shot",				sys.argv,	0)
	tstart		= cli_get_value("-tstart",			sys.argv,	0.)
	tend		= cli_get_value("-tend",				sys.argv,	0.)
	d_from_sep	= cli_get_value("-d_from_sep",		sys.argv,   0.0)
	deltas		= cli_get_value("-deltas",			sys.argv,	0.1)
	path_label 	= cli_get_value("-path_label",		sys.argv,   [""])
	no_labels	= cli_present("-no_labels",			sys.argv)
	save		= cli_get_value("-save",				sys.argv, "none")
	log_scale	= cli_present("-log_scale",			sys.argv)
	diff		= cli_present("-diff",				sys.argv)
	extra_walls	= cli_present("-extra_walls",		sys.argv)
	one_plot	= cli_present("-one_plot",			sys.argv)
	plot1d_sp_on_wall(path=path, evolution=evolution, mod_file=mod_file, exp_files=exp_files, shot=shot, tstart=tstart, tend=tend, deltas=deltas, d_from_sep=d_from_sep, path_label=path_label, no_labels=no_labels, log_scale=log_scale, diff=diff, extra_walls=extra_walls, one_plot=one_plot, save=save)
	exit()

#=======================================

# Function definition is here

import types
import os
import h5py
from math													import sqrt
import numpy											as np
import scipy.interpolate							as interp
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot							as pyp
from matplotlib.backends.backend_pdf	import PdfPages

from routines.h5_routines							import h5_read
from routines.set_profile_from_filedata		import set_profile_from_filedata
from routines.exp_data_routines				import get_exp_data_point
from routines.utils_routines						import argsort_mat
from routines.utils_walls							import plot2d_walls
from routines.get_strike_points					import get_strike_points
from routines.intersect_contour				import intersect_2contours
from routines.globals									import DEBUG, KB, BALLOONING_NAMES

from mesh.get_rz_core_sep						import get_rz_core_sep
from mesh.get_rho_in_out_core_sep		import get_rho_in_out_core_sep
from mesh.find_zones_intersections			import find_zones_intersections
from mesh.compute_mesh_intersections	import compute_mesh_intersections
#from mesh.compute_target_from_mesh_intersections	import compute_target_from_mesh_intersections

from eirene.get_wall_triangles					import get_wall_triangles

from files.load_soledge_mesh_file				import load_soledge_mesh_file
from files.load_eirene_triangles					import load_eirene_triangles
from files.load_plasma_files						import load_plasma_files
from files.load_exp_data							import load_exp_data
from files.load_ions_list								import load_ions_list
from files.load_text_data							import load_text_data
from files.load_refpar_file							import load_refpar_file
from files.load_exp_descr							import load_exp_descr
from files.save_stat									import save_stat


#==============================================================================
# This routine plots strike points 
#==============================================================================

def plot1d_sp_on_wall(path=[], evolution=[], mod_file="", exp_files=[], shot=0, tstart=0., tend=0., deltas=0.1, d_from_sep=0., path_label=[], no_labels=0, log_scale=0, diff=0, extra_walls=0, one_plot=0, save="none"):

	print("plot1d_sp_on_wall")

	save_norm = False
	if(save == "stat_norm"):
		save_norm = True
		save = "stat"
	
	if(diff != 0):
		if((evolution == 0) or (len(evolution) != 2)):
			print("\tWith diff option two evolutions must be provided")
			exit()
		log_scale	 = 0
		mod_file	 = ""
		exp_files	 = ""

#	prepare for experimental data

	exp_data_ok = False
	if((len(exp_files) > 0) and (len(exp_files[0]) > 0)):
		exp_descr = []
		Diags	  = []
		for exp_file in exp_files:
			if(len(exp_file) > 0):
				exp_descr.append(load_exp_descr(exp_file))					#load diags description
				if((tstart > 0) and (tstart <= tend)):
					exp_descr[-1].tstart = tstart
					exp_descr[-1].tend   = tend
				if(shot != 0): exp_descr[-1].shot = shot
		exp_data_ok = True

#	preapre paths

	if(len(evolution) == 0): evolution = [0]
	if(len(path) == 0):		 path = [""]
	elif(len(path) > 1):	 evolution = [evolution[0]]

	if(((len(path) > 1) or (len(evolution) > 1)) and (save == "stat")):
		print("\tSave stat is possible with one evolution or one path only")
		exit()

	if((len(path_label) < len(path)) or (len(path_label[0]) == 0)):
		path_label = []
		for in_path  in path: path_label.append(os.path.basename(os.path.abspath(in_path)))
		HasPathLabel=False
	else:
		HasPathLabel=True

#	Read reference parameters
	
	path0  = path[0]
	if((len(path0) > 0) and (path0[-1] != "/")): path0 = path0 + "/"

	RefPar = load_refpar_file(path0+"Results/")

	ions = load_ions_list(path0)

#	Read soledge and eirene mesh

	Config = load_soledge_mesh_file(path0+"mesh.h5")

	Eirene = load_eirene_triangles(path0+"triangles.h5")

	RStrikes, ZStrikes, InXPoints = get_strike_points(Config)
	nInXPoints = RStrikes.shape[0]

	Plasmas = load_plasma_files(path0, nZones=len(Config.Zones), iPlasmas = [0], Evolution=evolution[0])

	if((len(path) < 2) and (len(evolution) < 2)): Tempus = Plasmas[0][0].tempus
	else:										  Tempus = RefPar.time

	RWalls, ZWalls, iKnots = get_wall_triangles(Eirene,  nTri=Plasmas[0][0].Wall.ntri, Side=Plasmas[0][0].Wall.side)

#	Sort strike points from internal mid-plane anticlock
#	print("k, RWalls, ZWalls")
#	for k in range(len(RWalls)): print(k, RWalls[k], ZWalls[k])


	Rcore, Zcore, CoreMegazone = get_rz_core_sep(Config, core_and_sep = False)
	R0 		= 0.5*(Rcore.min() + Rcore.max())
	iMax = np.argmax(Rcore)
	Z0 = Zcore[iMax]

#	define walls arrays starting from z at plasma center and going anti clock

	iWLeft = np.where(RWalls < R0)[0]
	iWmid  = iWLeft[np.argmin(np.abs(ZWalls[iWLeft] - Z0))]

	ImWalls = np.arange(len(ZWalls), dtype='i4')
	ImWalls = np.roll(ImWalls,-iWmid)
	RmWalls = np.roll(RWalls,-iWmid)
	ZmWalls = np.roll(ZWalls,-iWmid)
	ns8 	= int(len(ZmWalls)/8)
	WallDirection = 1.
	if(ZmWalls[ns8] > ZmWalls[-ns8]):
		WallDirection = -1.
		RmWalls = RmWalls[::-1]
		ZmWalls = ZmWalls[::-1]
		ImWalls = ImWalls[::-1]

#	Find strike points on new wall array

	iSegs = np.empty((nInXPoints, 2), dtype='i4')
	for k in range(nInXPoints):
		for i in range(2):
			iSegs[k,i]	= np.argmin(np.sqrt((RmWalls - RStrikes[k,i])**2 + (ZmWalls - ZStrikes[k,i])**2))


#		if(	((iSegs[k,0] > iSegs[k,1]) and (0.5*(ZStrikes[k,0] + ZStrikes[k,1]) < Z0)) or
#			((iSegs[k,0] < iSegs[k,1]) and (0.5*(ZStrikes[k,0] + ZStrikes[k,1]) > Z0))):			#swap strike points
		if(iSegs[k,0] > iSegs[k,1]):			#swap strike points
			iSegs[k,:]	  = iSegs[k,::-1]
			RStrikes[k,:] = RStrikes[k,::-1]
			ZStrikes[k,:] = ZStrikes[k,::-1]

	iXpoint, iStrike = argsort_mat(iSegs)										#get SP from inner bottom to inner top

#	Find triangle at strike points

	iTri = np.empty((nInXPoints, 2), dtype='i4')
	for k in range(nInXPoints):
		for i in range(2):
			iTri[k,i]	= np.argmin(np.sqrt((RWalls - RStrikes[k,i])**2 + (ZWalls - ZStrikes[k,i])**2))

	iTri	 = iTri[iXpoint, iStrike].reshape(-1)				#sorte by SPs

#	Define max wall interval to plot

	if(d_from_sep > 0.):
				
		RZLine		= np.array([[R0, Z0], [R0+100., Z0]])
		dist_ew, SpCell_ew, SpPos_ew = compute_target_from_mesh_intersections(Config, Eirene, RZLine, d_from_sep, east_west=1)
		iSeq_ew = np.where(Plasmas[0][0].Wall.ntri == Eirene.WallTriangles.ntri[SpCell_ew[-1,3]])[0][0]

		dist_we, SpCell_we, SpPos_we = compute_target_from_mesh_intersections(Config, Eirene, RZLine, d_from_sep, east_west=-1)
		iSeq_we = np.where(Plasmas[0][0].Wall.ntri == Eirene.WallTriangles.ntri[SpCell_we[-1,3]])[0][0]
		if(nInXPoints == 1):
			ConfigType = 0																			#SND/XD or SXD configuaration
			deltas_sp1 = Plasmas[0][0].Wall.swall[iTri[0]-1] -Plasmas[0][0].Wall.swall[iSeq_ew]
			deltas_sp2 = Plasmas[0][0].Wall.swall[iSeq_we] - Plasmas[0][0].Wall.swall[iTri[1]-1]
			deltas_sp = np.array([deltas_sp1, deltas_sp2])
		elif(nInXPoints == 2):
			if((X_points[InXPoints[0]].Z - Z0)*(X_points[InXPoints[1]].Z - Z0) < 0):
				ConfigType = 1																		#DND configuaration					

				deltas_sp3 = Plasmas[0][0].Wall.swall[iTri[2]-1] -Plasmas[0][0].Wall.swall[iSeq_ew]
				deltas_sp2 = Plasmas[0][0].Wall.swall[iSeq_we] - Plasmas[0][0].Wall.swall[iTri[1]-1]

				RZLine		= np.array([[R0, Z0], [0., Z0]])
				dist_ew, SpCell_ew, SpPos_ew = compute_target_from_mesh_intersections(Config, Eirene, RZLine, d_from_sep, east_west=1)
				iSeq_ew = np.where(Plasmas[0][0].Wall.ntri == Eirene.WallTriangles.ntri[SpCell_ew[-1,3]])[0][0]
				deltas_sp1 = Plasmas[0][0].Wall.swall[iTri[0]-1] -Plasmas[0][0].Wall.swall[iSeq_ew]

				dist_we, SpCell_we, SpPos_we = compute_target_from_mesh_intersections(Config, Eirene, RZLine, d_from_sep, east_west=-1)
				iSeq_we = np.where(Plasmas[0][0].Wall.ntri == Eirene.WallTriangles.ntri[SpCell_we[-1,3]])[0][0]
				deltas_sp4 = Plasmas[0][0].Wall.swall[iSeq_we] - Plasmas[0][0].Wall.swall[iTri[3]-1]

				deltas_sp = np.array([deltas_sp1, deltas_sp2, deltas_sp3, deltas_sp4])

			else:
				RZcore, RZsep, CoreMegazone, SepMegazone, jSep = get_rz_core_sep(Config, core_and_sep = True)
				dist_from_sep =np.empty(2,dtype='f8')
				for k in range(nInXPoints):
					dist_from_sep[k] = sqrt(np.min((RZsep[:,0] - X_points[InXPoints[k]].R)**2+(RZsep[:,1] - X_points[InXPoints[k]].Z)**2))
				iInXSep = np.argmin(dist_from_sep)
				if(X_points[InXPoints[iInXSep]].psi > X_points[InXPoints[1-iInXSep]].psi):
					ConfigType = 2																	#SFD+ configuaration					
				
					deltas_sp1 = Plasmas[0][0].Wall.swall[iTri[0]-1] -Plasmas[0][0].Wall.swall[iSeq_ew]
					deltas_sp4 = Plasmas[0][0].Wall.swall[iSeq_we] - Plasmas[0][0].Wall.swall[iTri[3]-1]
					deltas_sp = np.array([deltas_sp1, deltas_sp1, deltas_sp4, deltas_sp4])
				else:
					iSeq = np.where(Plasmas[0][0].Wall.ntri == Eirene.WallTriangles.ntri[SpCell_we[0,3]])[0][0]
					if(abs(iSeq - iTri[0] -1) < abs(iSeq - iTri[2] -1)):						
						ConfigType = 3																#SFD- configuaration low field side
						deltas_sp1 = Plasmas[0][0].Wall.swall[iTri[0]-1] -Plasmas[0][0].Wall.swall[iSeq_we]
						deltas_sp2 = 0.
						deltas_sp4 = 0. 
						for i_wall_tri in SpCell_ew[:,3]:
							iSeq = np.where(Plasmas[0][0].Wall.ntri == Eirene.WallTriangles.ntri[i_wall_tri])[0][0]
							if(abs(iSeq - iTri[1] -1) < abs(iSeq - iTri[3] -1)):
								deltas_sp2 = Plasmas[0][0].Wall.swall[iSeq] - Plasmas[0][0].Wall.swall[iTri[1]-1]
							else:
								deltas_sp4 = Plasmas[0][0].Wall.swall[iSeq] - Plasmas[0][0].Wall.swall[iTri[3]-1] 

						deltas_sp2 = deltas_sp2 + deltas_sp4
						deltas_sp = np.array([deltas_sp1, deltas_sp2, deltas_sp2, deltas_sp2])
					else:
						ConfigType = 4																#SFD- configuaration high field side
						deltas_sp4 = Plasmas[0][0].Wall.swall[iTri[3]-1] -Plasmas[0][0].Wall.swall[iSeq_ew]
						deltas_sp1 = 0.
						deltas_sp3 = 0. 
						for i_wall_tri in SpCell_we[:,3]:
							iSeq = np.where(Plasmas[0][0].Wall.ntri == Eirene.WallTriangles.ntri[i_wall_tri])[0][0]
							if(abs(iSeq - iTri[2] -1) < abs(iSeq - iTri[0] -1)):
								deltas_sp3 = Plasmas[0][0].Wall.swall[iTri[2]-1] -Plasmas[0][0].Wall.swall[iSeq] 
							else:
								deltas_sp1 = Plasmas[0][0].Wall.swall[iTri[0]-1] -Plasmas[0][0].Wall.swall[iSeq] 

						deltas_sp1 = deltas_sp1 + deltas_sp3
						deltas_sp = np.array([deltas_sp1, deltas_sp1, deltas_sp1, deltas_sp4])
	else:
		deltas_sp = np.ones(nInXPoints*2, dtype = 'f8')*deltas


	iTriLow  = np.empty_like(iTri)
	iTriHigh = np.empty_like(iTri)
	for k in range(len(iTri)):
		swall_sp = Plasmas[0][0].Wall.swall - Plasmas[0][0].Wall.swall[iTri[k]-1]
		iTriLow[k]  = np.argmax(swall_sp >= -deltas_sp[k])
		iTriHigh[k] = np.argmin(swall_sp <=  deltas_sp[k])

	for iX in range(len(iXpoint)):
		if(sqrt((RWalls[iTri[iX]]-RStrikes[iXpoint[iX],iStrike[iX]])**2 + (ZWalls[iTri[iX]]-ZStrikes[iXpoint[iX],iStrike[iX]])**2) > 1e-2):
			print("\tERROR in strike point position for xPoint #:",iX+1)
			print("\t\tRw=",RWalls[iTri[iX]]," Zw=",ZWalls[iTri[iX]])
			print("\t\tRs=",RStrikes[iXpoint[iX],iStrike[iX]]," Zs=",ZStrikes[iXpoint[iX],iStrike[iX]])
			print("\t\td =",sqrt((RWalls[iTri[iX]]-RStrikes[iXpoint[iX],iStrike[iX]])**2 + (ZWalls[iTri[iX]]-ZStrikes[iXpoint[iX],iStrike[iX]])**2) )
			exit()

#	Prepare for plotting and saving data

	i_plot_file = 0

	if(save == "pdf"):	pdf = PdfPages("plot1d_sp_on_wall_t={:.3f}.".format(Tempus)+save)   #pdf in one file only
	if(save == "csv"):
		if((len(evolution) > 1) or (len(path) > 1)):
			save = "none"
			print("plot1d_sp_on_wall: cannot save to csv more than one profile!!!")
		else:
			csv = []
#			csv.append(types.SimpleNamespace())
#			csv[-1].Name   = "Dist (m)"
#			csv[-1].Values = dist


	xLabels    = [["$s\ (m)$",	"$s\ (m)$",   "$s\ (m)$"], \
				  ["$s\ (m)$",	 "$s\ (m)$",   "$s\ (m)$","$s\ (m)$", "$s\ (m)$", "$s\ (m)$", "$s\ (m)$"]]
	yLabels    = [["$n\ (*10^{19}\ m^{-3})$", "$T\ (eV)$",  "$\Gamma_E\ (MW/m^{2})$"], \
				  ["$n\ (*10^{19}\ m^{-3})$", "$T\ (eV)$", "$\Gamma_E\ (MW/m^{2})$",	"$\Gamma_n\ (*10^{22}\ m^{-2})$","$\Gamma_{Etot}\ (MW/m^{2})$", "$J_{sat}\ (kA/m^{2})$", "$J_{sat-par}\ (kA/m^{2})$"]]
	yLogLabels = [["$Ln(n)\ (*10^{19}\ m^{-3})$", "$Ln(T)\ (eV)$",  "$Ln(\Gamma_E)\/(MW/m^{2})$"], \
				  ["$Ln(n)\ (*10^{19}\ m^{-3})$", "$Ln(T)\ (eV)$", "$Ln(\Gamma_E)\/(MW/m^{2})$",	"$Ln(\Gamma_n)\/(*10^{22}\ m^{-2})$", "$Ln(\Gamma_{Etot})\ (MW/m^{2})$", "$Ln(J_{sat})\/(kA/m^{2})$", "$Ln(J_{sat-par})\/(kA/m^{2})$"]]

	ValueFacts	= [[1e-19,1.,1e-6],[1e-19,1.,1e-6,1e-22,1e-6,1e-3,1e-3]]
	iwValues	= [[0,1,3],[0,1,3,4,8,6,7]]															#index in plasma file of wall parameter

	LogScales	= [["log", "log", "log"],["log", "log", "log", "log", "log", "log", "log"]]

	colors = ['b','g','r','c','m','y','b','g','r','c','m','y']
	TitlePlots = [0]

	if(save == "stat"):
		iSaveStat	  = [[1,1,0],[0,1,0,0,2,2,2]]															#0=no stat, 1=max, 2=max & integral
		nSaveStat	  = 9
		kSaveStatBase  = [2,4]
		if(save_norm):
			StatHeaderForm  = "ne_max_sp{:d} (*1e19 m-3),Te_max_sp{:d} (eV),Ti_max_sp{:d} (eV),Ptot_max_sp{:d} (MW/m^2),Ptot_int_sp{:d} (MW),Jsat_max_sp{:d} (kA/m^2),Jsat_int_sp{:d} (kA),Jsat_par_max_sp{:d} (kA/m^2),Jsat_par_int_sp{:d} (kA/m^3),"
			StatFact = [1.e-19, 1., 1., 1e-6, 1e-6, 1e-3, 1e-3, 1e-3,1e-3]
		else:
			StatHeaderForm  = "       ne_max_sp{:d},       Te_max_sp{:d},       Ti_max_sp{:d},     Ptot_max_sp{:d},     Ptot_int_sp{:d},     Jsat_max_sp{:d},     Jsat_int_sp{:d}, Jsat_par_max_sp{:d}, Jsat_par_int_sp{:d},"
			StatFact = [1, 1., 1., 1., 1., 1., 1., 1.,1.]

		StatHeader	   =  "             path,             time,"
		StatFormats		= ['"{:>15}",','{:17.4e},']
		StatValues		= ["", 0.]
		for iX in range(1,len(iXpoint)+1): 
			StatHeader += StatHeaderForm.format(iX, iX, iX, iX, iX, iX, iX, iX, iX)
			for i in range(nSaveStat):
				StatFormats.append('{:17.4e},')
				StatValues.append(StatFact[i])

		StatHeader  = StatHeader[:-1]
		StatFormats[-1] = StatFormats[-1][:-1]

		StatValues[0] = path_label[0]
		StatValues[1] = Tempus
		if(len(mod_file) > 0):
			StatValuesMod = np.zeros((nSaveStat*len(iXpoint)+2), dtype='f8')
			StatValuesMod[0] = path_label[0]
			StatValuesMod[1] = Tempus
			kSaveStat = 0
			for iX in range(1,len(iXpoint)+1): 
				for i in range(nSaveStat): 
					StatValuesMod[kSaveStat] = StatFact[i]
					kSaveStat += 1			

#	Prepare plotting structures

	nRows		= 2
	nCols		= 3
	PlotPerFig	= nRows*nCols
	nPLots		= len(TitlePlots) + (len(xLabels[0]) + len(xLabels[1]))*len(iXpoint)
	nFigs		= int(nPLots/PlotPerFig)
	if(nFigs*PlotPerFig < nPLots): nFigs += 1

	Fig = []
	Ax  = []
	if(one_plot != 1):
		for iF in range(nFigs):	
			Fig.append(pyp.figure())
			for i in range(min(PlotPerFig,nPLots-iF*PlotPerFig)): 
				Ax.append(Fig[-1].add_subplot(nRows,nCols,i+1))
				Ax[-1].locator_params(axis='x',nbins=4)

			Fig[-1].tight_layout(pad=2., w_pad=3., h_pad=3.)
	else:
		for i in range(nPLots):
			Fig.append(pyp.figure())
			Ax.append(Fig[i].add_subplot(111))

	for figure in Fig:  figure.patch.set_facecolor('white')

	for ip in TitlePlots: 
		Ax[ip].set_aspect(1.)
		Ax[ip].autoscale(enable=True, axis='both', tight=True)
		Ax[ip].set_xlabel("$R\ (m)$")
		Ax[ip].set_ylabel("$Z\ (m)$")
		if(len(path) == 1): Ax[ip].set_title(path_label[0]+" @ t={:.3f} s".format(Tempus))

	ip = 0
	for iPlasma in range(2):
		for iX in range(len(iXpoint)):
			for i in range(len(xLabels[iPlasma])):
				ip += 1
				if((i > 0) or (((len(evolution) != 1) or (evolution[0] == 0)) and (diff == 0))):
					Ax[ip].set_title(ions[iPlasma])
				elif(diff == 0):
					Ax[ip].set_title(ions[iPlasma]+" evol. {:d}".format(evolution[0]))
				else:
					Ax[ip].set_title(ions[iPlasma]+" diff {:d}-{:d}".format(evolution[1],evolution[0]))

				Ax[ip].autoscale(enable=True, axis='both', tight=True)
				Ax[ip].set_xlabel(xLabels[iPlasma][i]+" @ sp{:d}".format(iX+1))
				if(log_scale == 0):
					Ax[ip].set_ylabel(yLabels[iPlasma][i])
					Ax[ip].set_yscale('linear')
				else:
					Ax[ip].set_ylabel(yLogLabels[iPlasma][i])
					Ax[ip].set_yscale(LogScales[iPlasma][i])


#	Plot Xpoints

	for ip in TitlePlots:
		plot2d_walls(Ax[ip], Config.Walls, extra_wall=extra_walls)
		for k in range(nInXPoints):
			Ax[ip].contour(Config.r2D, Config.z2D, Config.flux2D, colors="r", linestyles='-', linewidths=2, levels=[Config.X_points[InXPoints[k]].psi])

		for iX in range(len(iXpoint)):
#			Ax[ip].text(RStrikes[iXpoint[iX],iStrike[iX]],ZStrikes[iXpoint[iX],iStrike[iX]], "{:d}".format(iX+1), horizontalalignment="center", size="large", color="g")
			Ax[ip].text(RStrikes[iXpoint[iX],iStrike[iX]],ZStrikes[iXpoint[iX],iStrike[iX]], "{:d}".format(iX+1), size="large", color="g")
			if(iStrike[iX] == 0): Ax[ip].plot(RWalls[iTriLow[iX]:iTriHigh[iX]+1],ZWalls[iTriLow[iX]:iTriHigh[iX]+1], 'r-', linewidth=2)
			else: 				  Ax[ip].plot(RWalls[iTriLow[iX]:iTriHigh[iX]+1],ZWalls[iTriLow[iX]:iTriHigh[iX]+1], 'g-', linewidth=2)
#		Load Plasma parameters

	ip = 0
	for iPlasma in range(2):
		Values	 = []
		for i in range(len(xLabels[iPlasma])):
			Values.append([])

		for iPh in range(len(path)):
			for iEv in range(len(evolution)):
				path0  = path[iPh]
				if((len(path0) > 0) and (path0[-1] != "/")): path0 = path0 + "/"
				Plasmas = load_plasma_files(path0, nZones=len(Config.Zones), iPlasmas = [iPlasma], Evolution=evolution[iEv])

				for i in range(len(xLabels[iPlasma])):
					Values[i].append(Plasmas[0][0].Wall.Values[iwValues[iPlasma][i]]*ValueFacts[iPlasma][i])

		if(diff != 0):
			for i in range(len(Values)): Values[i][0] = Values[i][1] - Values[i][0]

		mod_data_ok = False
		if(len(mod_file) > 0):
			print("reading=",mod_file)
			Headers, TextData = load_text_data(mod_file)
			xDataName  = "DIST"
			File_xValues = []
			File_yValues = []
			Exp_Data	 = []
			for i in range(len(iwValues[iPlasma])):
				File_xValues.append([])
				File_yValues.append([])
				for iX in range(len(iXpoint)):
					Exp_Data[-1].append([])
					Empty = np.empty(0, dtype='i4')
					xRange = np.array([-deltas,deltas])
					yDataName = Plasmas[0][0].Wall.VNames[iwValues[iPlasma][i]]
					xValue, yValue = set_profile_from_filedata(Headers, TextData, xDataName, yDataName, xRange, Empty, extrapolate = False, NameExt="sp{:d}".format(iX+1))
					if(len(yValue) > 0):
						File_xValues[-1].append(xValue)
						File_yValues[-1].append(yValue)
						mod_data_ok = True
					else:
						File_xValues[-1].append([])
						File_yValues[-1].append([])

#		Read exp_file. data

		if(exp_data_ok):
			Exp_Data  = []
			Par_xName = "DIST"
			Lengths	  = np.array([[-deltas],[deltas]])			#This works because DIST  is the first in LENGTH_TYPES and so the first column in Lengths
			for i in range(len(iwValues[iPlasma])):
				Exp_Data.append([])
				for iX in range(len(iXpoint)):
					Exp_Data[-1].append([])
					if(iwValues[iPlasma][i] > -1):
						Par_yName = Plasmas[0][0].Wall.VNames[iwValues[iPlasma][i]]
						for iExp in range(len(exp_descr)):
							diag_data = get_exp_data_point(Par_xName, Par_yName, exp_descr[iExp], Lengths, DiagType = "WALL", NameExt="sp{:d}".format(iX+1))
							Exp_Data[-1][-1].append(diag_data)
					else:
						Exp_Data[-1][-1].append([]) 
#		Plot parameters

		for iX in range(len(iXpoint)):
			swall_sp = Plasmas[0][0].Wall.swall[iTriLow[iX]:iTriHigh[iX]+1] - Plasmas[0][0].Wall.swall[iTri[iX]-1]
#			print("iX=",iX+1," sc=",Plasmas[0][0].Wall.swall[iTri[iX]-1])

			if(iStrike[iX] == 0): swall_sp = -WallDirection*swall_sp
			else:				  swall_sp =  WallDirection*swall_sp

			if((save == "csv") and (iPlasma == 0)):
				csv.append(types.SimpleNamespace())
				csv[-1].Name   ="dist_sp{:d}".format(iX+1)
				csv[-1].Values = np.copy(swall_sp)

			if(save == "stat"):
				kSaveStat    = kSaveStatBase[iPlasma] + nSaveStat*iX
				kSaveStatMod = kSaveStatBase[iPlasma] + nSaveStat*iX

			for i in range(len(xLabels[iPlasma])):
				ip += 1
				if((len(evolution) < 2) or (diff != 0)):
					if((len(path) < 2) or (diff != 0)):
						Ax[ip].plot(swall_sp,  Values[i][0][iTriLow[iX]:iTriHigh[iX]+1],  'b.-')

#						store extrema values

						if((save == "stat") and (iSaveStat[iPlasma][i] > 0)):
							StatValues[kSaveStat] = StatValues[kSaveStat]*np.max(Values[i][0][iTriLow[iX]:iTriHigh[iX]+1])/ValueFacts[iPlasma][i]
							kSaveStat += 1
						if((save == "stat") and (iSaveStat[iPlasma][i] > 1)):
							StatValues[kSaveStat] =  StatValues[kSaveStat]*abs(0.5*np.sum((swall_sp[1:]-swall_sp[:-1])* \
															 (Values[i][0][iTriLow[iX]:iTriHigh[iX]] + Values[i][0][iTriLow[iX]+1:iTriHigh[iX]+1]))/ValueFacts[iPlasma][i])
							kSaveStat += 1

						if(save == "csv"):
							csv.append(types.SimpleNamespace())
							i_slash= yLabels[iPlasma][i][1:-1].index("\ ")
							if(i_slash >= 0):
								if(iPlasma == 0):	csv[-1].Name   = yLabels[iPlasma][i][1:i_slash+1]+"_e_sp{:d}".format(iX+1) + yLabels[iPlasma][i][i_slash+2:-1].replace("\ ","")
								else:						csv[-1].Name   = yLabels[iPlasma][i][1:i_slash+1]+"_i_sp{:d}".format(iX+1) + yLabels[iPlasma][i][i_slash+2:-1].replace("\ ","")
							else:
								if(iPlasma == 0):	csv[-1].Name   = yLabels[iPlasma][i][1:-1]+"_e_sp{:d}".format(iX+1)
								else:						csv[-1].Name   = yLabels[iPlasma][i][1:-1]+"_i_sp{:d}".format(iX+1)

							csv[-1].Values = np.copy(Values[i][0][iTriLow[iX]:iTriHigh[iX]+1])
					else:
						for iPh in range(len(path)):
							if(no_labels == 0): Ax[ip].plot(swall_sp, Values[i][iPh][iTriLow[iX]:iTriHigh[iX]+1], '-', color = colors[iPh], label=path_label[iPh])
							else:						Ax[ip].plot(swall_sp, Values[i][iPh][iTriLow[iX]:iTriHigh[iX]+1], '-', color = colors[iPh])
							Ax[ip].plot(swall_sp, Values[i][iPh][iTriLow[iX]:iTriHigh[iX]+1], '-', color = colors[iPh])

				else:
					for iEv in range(len(evolution)):
						Ax[ip].plot(swall_sp, Values[i][iEv][iTriLow[iX]:iTriHigh[iX]+1],  '-', color = colors[iEv], label="{:d}".format(evolution[iEv]))
						Ax[ip].plot(swall_sp, Values[i][iEv][iTriLow[iX]:iTriHigh[iX]+1], '-', color = colors[iEv])

				if(mod_data_ok and (len(File_xValues[i][iX]) > 0)):
					Ax[ip].plot(File_xValues[i][iX], File_yValues[i][iX]*ValueFacts[iPlasma][i], "ro")

					if(save == "csv"):
						csv.append(types.SimpleNamespace())
						csv[-1].Name   = File_xNames[i][iX]+"_sp{:d}".format(iX+1)
						csv[-1].Values = np.copy(File_xValues[i][iX])
						csv.append(types.SimpleNamespace())
						csv[-1].Name   = File_yNames[i][iX]+"_sp{:d}".format(iX+1)
						csv[-1].Values = np.copy( File_yValues[i][iX])

#				store extrema values

				if(mod_data_ok):
					if((save == "stat") and (iSaveStat[iPlasma][i] > 0)):
						if(len(File_xValues[i][iX]) > 0): StatValuesMod[kSaveStatMod] = np.max(File_yValues[i][iX])
						kSaveStatMod += 1
					if((save == "stat") and (iSaveStat[iPlasma][i] > 1)):
						if(len(File_xValues[i][iX]) > 0): StatValuesMod[kSaveStatMod] = abs(0.5*np.sum((File_xValues[i][iX][1:]-File_xValues[i][iX][:-1])* \
																									 (File_yValues[i][iX][:-1] + File_yValues[i][iX][1:])))
						kSaveStatMod += 1

				if(exp_data_ok and (len(Exp_Data[i][iX]) > 0)):
					for iExp in range(len(Exp_Data[i][iX])):
						for id in range(len(Exp_Data[i][iX][iExp])):
							Exp_Data_Diag = Exp_Data[i][iX][iExp][id]
							if(iStrike[iX] == 0): xValues = -Exp_Data_Diag.xValues
							else:				  xValues =  Exp_Data_Diag.xValues

							if(Exp_Data_Diag.errors):
								Ax[ip].errorbar(xValues, Exp_Data_Diag.yValues*ValueFacts[iPlasma][i], yerr=Exp_Data_Diag.yErrors*ValueFacts[iPlasma][i], fmt = Exp_Data_Diag.marker, color=Exp_Data_Diag.color, label=Exp_Data_Diag.label)
							else:
								Ax[ip].plot(xValues, Exp_Data_Diag.yValues*ValueFacts[iPlasma][i], marker = Exp_Data_Diag.marker, color=Exp_Data_Diag.color, linestyle='none', label=Exp_Data_Diag.label)

							if(save == "csv"):
								csv.append(types.SimpleNamespace())
								csv[-1].Name   = "Exp_"+Exp_Data_Diag.xName+"_sp{:d}".format(iX+1)
								csv[-1].Values = np.copy(xValues)
								csv.append(types.SimpleNamespace())
								csv[-1].Name   = "Exp_"+Exp_Data_Diag.yName+"_sp{:d}".format(iX+1)
								csv[-1].Values = np.copy(Exp_Data_Diag.yValues)
								if(Exp_Data_Diag.errors):
									csv.append(types.SimpleNamespace())
									csv[-1].Name   = "Exp_"+Exp_Data_Diag.yName+"_err"+"_sp{:d}".format(iX+1)
									csv[-1].Values = np.copy(Exp_Data_Diag.yErrors)

	xv = 0.
	for i in range(1,len(Ax)):
#		pyp.setp(Ax[k].yaxis.get_ticklabels(), rotation=90.)
		if(len(np.where(TitlePlots == i)[0]) == 0):
			Ax[i].axvline(x=xv, color='k', linestyle='dashed')
			if(((len(path) > 1) or (len(evolution) > 1)) and (no_labels == 0)): Ax[i].legend(fontsize='small', loc='lower left')

	if(save != "none"):
		if(save == "stat"):
			save_stat("stat_sp_on_wall.csv", StatHeader, StatValues, StatFormats)

			if(mod_data_ok):
				np.savetxt("data_wall_sp_mod.csv", StatValuesMod.reshape(1,len(StatValuesMod)), header=StatHeader, delimiter=",", fmt="%13.4e", comments="")

		elif(save == "csv"):
			maxLen = 0
			for ii in range(len(csv)): maxLen = max(maxLen, len(csv[ii].Values))
			save_cvs= np.zeros((maxLen,len(csv)), dtype='f8')
			Header = ""
			for ii in range(len(csv)):
				Header = Header + csv[ii].Name + ","
				save_cvs[:len(csv[ii].Values), ii] = csv[ii].Values
				if(len(csv[ii].Values) < maxLen): save_cvs[len(csv[ii].Values):, ii]  = np.nan

				if(HasPathLabel): 
					np.savetxt("sp_on_wall_{:s}.{:s}".format(path_label[0], save), save_cvs, header=Header, delimiter=",", fmt="%15.7e", comments="")
				else:
					np.savetxt("plot1d_sp_on_wall_t={:.3f}_{:d}.".format(Tempus,i_plot_file)+save, save_cvs, header=Header, delimiter=",", fmt="%15.7e", comments="")
			
	
		elif(save != "none"):
			for i in range(len(Fig)):
				i_plot_file += 1
				if(one_plot != 1): Fig[i].set_size_inches(10.05,7.44)
				if(save == "pdf"):
					pdf.savefig(Fig[i])
				else:
					Fig[i].savefig("plot1d_sp_on_wall_t={:.3f}_{:d}.".format(Tempus,i_plot_file)+save)


		pyp.show(block=False)
		pyp.close()
	else:
		pyp.show()

	if(save == "pdf"):	pdf.close()

		
	print("plot1d_sp_on_wall: Completed")

	return

#=================================================================

from math								import sqrt, asin
import numpy							as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot				as pyp

from routines.utils_walls				import get_in_out_walls

from mesh.get_rho_in_out_core_sep		import get_rho_in_out_core_sep
from mesh.find_zones_intersections		import find_zones_intersections
from mesh.compute_mesh_intersections	import compute_mesh_intersections


def compute_target_from_mesh_intersections(Config, Eirene, RZLine, d_from_sep=0.01, east_west=0, l_pol=1):

	Zones	= Config.Zones

	WallTriangles = Eirene.WallTriangles

	Cut = find_zones_intersections(Config, RZLine)
	Lengths, IntRZ, IntCEll = compute_mesh_intersections(Config, Cut, also_pos=True, use_mag_zones=False)
	in_wall, out_wall = get_in_out_walls(Config, IntRZ[:,0], IntRZ[:,1])
	Lengths = Lengths[in_wall,:]
	IntRZ	= IntRZ[in_wall,:]
	IntCEll	= IntCEll[in_wall,:]

	dist = Lengths[:,0] - Lengths[0,0]
	Rho, In_Sep, Out_Sep, RZcore, RZsep = get_rho_in_out_core_sep(Config, IntRZ[:,0], IntRZ[:,1])
	
	if(len(In_Sep) > 0):
		if(Out_Sep[-1] < In_Sep[0]):
			Out_Sep = np.append(Out_Sep, In_Sep[0])
		else:
			Out_Sep = np.append(In_Sep[-1], Out_Sep)

	Ri, Zi, is1, is2  = intersect_2contours(RZsep[:,0], RZsep[:,1], IntRZ[:,0], IntRZ[:,1])
	if(len(Ri)==0):
		pyp.plot(RZsep[:,0], RZsep[:,1],'k-')
		pyp.plot(IntRZ[:,0], IntRZ[:,1],'r-')
		pyp.show()

	dsep = sqrt((IntRZ[0,0] - Ri[0])**2 + (IntRZ[0,1] - Zi[0])**2)
	dist -= dsep

	if(d_from_sep == 0): d_from_sep = 2*dist.max()
	OutSep	= np.where((dist > 0) & (dist < d_from_sep))[0]
	dist	= dist[OutSep]
	IntRZ	= IntRZ[OutSep,:]
	IntCEll = IntCEll[OutSep,:]

	SpCell = np.empty((IntCEll.shape[0],4), dtype = 'i4')
	SpPos  = np.empty((IntCEll.shape[0],3), dtype = 'f8')

	for iCell in range(IntCEll.shape[0]):											#find nearest S
		iZone 	= IntCEll[iCell,0]
		ix		= IntCEll[iCell,1]
		iTheta	= IntCEll[iCell,2]

#		Find zones along poloidal coordinate

		if(Zones[iZone].Chi[ix,-1] == 1):
			iThEast = np.array([np.min(np.where(Zones[iZone].Chi[ix,iTheta:] == 1)[0])+iTheta])
			East = -1
		else:
			iThEast = np.array([Zones[iZone].Chi.shape[1]])
			East = Zones[iZone].Neighbour.east
		
		nThetaPts = iThEast[0] - iTheta

		if(Zones[iZone].Chi[ix,0] == 1):
			iThWest = np.array([np.max(np.where(Zones[iZone].Chi[ix,:iTheta] == 1)[0])])
			West = -1
		else:
			iThWest = np.array([0])
			West = Zones[iZone].Neighbour.west

		iThetaOff  = iTheta - iThWest[0]
		nThetaPts  = iThEast[0] - iThWest[0]
		iZones	   = np.array([iZone])

#		Look East

		while (East > -1):
			iZones = np.append(iZones, East)
			iThWest = np.append(iThWest,0)
			if(Zones[East].Chi[ix,-1] == 1):
				iThEast = np.append(iThEast, np.min(np.where(Zones[East].Chi[ix,:] == 1)[0]))
				if(iThEast[-1] == 0):
					iZones  = iZones[:-1]
					iThWest = iThWest[:-1]
					iThEast = iThEast[:-1]
					East	= -2
				else: East = -1
			else:
				iThEast = np.append(iThEast, Zones[East].Chi.shape[1])
				East 	 = Zones[East].Neighbour.east

			if(East > -2):
				nThetaPts += iThEast[-1]

#		Look West

		while (West > -1):
			iZones = np.append(West, iZones)
			iThEast = np.append(Zones[West].Chi.shape[1], iThEast)
			if(Zones[West].Chi[ix,0] == 1):
				iThWest = np.append(np.max(np.where(Zones[West].Chi[ix,:] == 1)[0])+1, iThWest)
				if(iThWest[0] == Zones[West].Chi.shape[1]):
					iThWest = iThWest[1:]
					iZones  = iZones[1:]
					iThEast = iThEast[1:]
					West = -2
				else: West = -1
			else:
				iThWest = np.append(0, iThWest)
				West = Zones[West].Neighbour.west

			if(West > -2):
				nThetaPts += iThEast[0] - iThWest[0]
				iThetaOff += iThEast[0] - iThWest[0]

		Rpol  = np.empty((nThetaPts), dtype = 'f8')
		Zpol  = np.empty((nThetaPts), dtype = 'f8')
#		Bfact = np.empty((nThetaPts), dtype = 'f8')
		jOff = 0
		for k in range(len(iZones)):
			Rpol[jOff: jOff + iThEast[k] - iThWest[k]] = Zones[iZones[k]].gridRc[ix, iThWest[k]:iThEast[k]]
			Zpol[jOff: jOff + iThEast[k] - iThWest[k]] = Zones[iZones[k]].gridZc[ix, iThWest[k]:iThEast[k]]
#			Bfact[jOff: jOff + iThEast[k] - iThWest[k]] = np.sqrt(Zones[iZones[k]].Bphi[ix, iThWest[k]:iThEast[k]]**2 + Zones[iZones[k]].Br[ix, iThWest[k]:iThEast[k]]**2 + Zones[iZones[k]].Bz[ix, iThWest[k]:iThEast[k]]**2)/ \
#														  np.sqrt(Zones[iZones[k]].Br[ix, iThWest[k]:iThEast[k]]**2 + Zones[iZones[k]].Bz[ix, iThWest[k]:iThEast[k]]**2)

			jOff += iThEast[k] - iThWest[k]

		if(l_pol == 0):
#			Lpara1 = np.append(0., np.cumsum(np.sqrt((Rpol[1:]-Rpol[:-1])**2 + (Zpol[1:]-Zpol[:-1])**2)*(Bfact[1:]+Bfact[:-1])*0.5))

			dl = np.empty((nThetaPts), dtype = 'f8')
			jOff = 0
			for k in range(len(iZones)):
				dtheta = (Zones[iZones[k]].zb[ix, 1:] - Zones[iZones[k]].zb[ix, :-1])*2.*np.pi
				dlZone	   = 2.*dtheta/Gmet[iZones[k]][ix+1, 1:-1]
				
				dl[jOff: jOff + iThEast[k] - iThWest[k]] = dlZone[iThWest[k]:iThEast[k]]
				jOff += iThEast[k] - iThWest[k]
			dlZone = 0
			Lpara = np.cumsum(dl)
			dl	  = 0
		else:
			Lpara = np.append(0., np.cumsum(np.sqrt((Rpol[1:]-Rpol[:-1])**2 + (Zpol[1:]-Zpol[:-1])**2)))

		Lpara = Lpara - Lpara[iThetaOff]

		if(east_west == 0):
			if(Lpara[-1] >= -Lpara[0]): east_west = -1
			else:						east_west = 1

		if(east_west == -1):
			SpCell[iCell,0] = iZones[-1]						#SpCell[k_zone, i_grid, j_grid, ind_wall_tri]
			SpCell[iCell,1] = ix
			SpCell[iCell,2] = iThEast[-1]-1
			SpPos[iCell,2]	= -Lpara[-2]						#SpPos[R_target, Z_target,Line_lenght]
		else:
			SpCell[iCell,0] = iZones[0]
			SpCell[iCell,1] = ix
			SpCell[iCell,2] = iThWest[0]
			SpPos[iCell,2]	= Lpara[1]


		"""
		gridR = SpCell[iCell,0]].gridR
		gridZ = SpCell[iCell,0]].gridZ
		jsp = SpCell[iCell,2]
		R1	= gridR[ix, jsp]; R2 = gridR[ix, jsp+1]; R3 = gridR[ix+1, jsp+1]; R4 = gridR[ix+1, jsp]
		Z1	= gridZ[ix, jsp]; Z2 = gridZ[ix, jsp+1]; Z3 = gridZ[ix+1, jsp+1]; Z4 = gridZ[ix+1, jsp]
		DeltaRTheta = 0.5*((R2-R1) + (R3-R4))
		DeltaZTheta = 0.5*((Z2-Z1) + (Z3-Z4))
		DeltaLTheta = sqrt(DeltaRTheta**2 + DeltaRTheta**2)
		DeltaRx		= 0.5*((R4-R1) + (R3-R2))
		DeltaZx		= 0.5*((Z4-Z1) + (Z3-Z2))
		DeltaLPhi	= sqrt(DeltaRx**2 + DeltaRx**2)
		DeltaLPhi	= (DeltaRx*DeltaRTheta + DeltaZx*DeltaZTheta)/DeltaLTheta
		"""

		SpPos[iCell,0]	= Zones[SpCell[iCell,0]].gridRc[ix, SpCell[iCell,2]]
		SpPos[iCell,1]	= Zones[SpCell[iCell,0]].gridZc[ix, SpCell[iCell,2]]

		iZoneEnd	= SpCell[iCell,0]
		iThetaEnd	= SpCell[iCell,2]
		for iOff in range(5):
			if(iOff == 0):
				aa = np.where((WallTriangles.k == iZoneEnd) & (WallTriangles.i == ix) & (WallTriangles.j == iThetaEnd))[0]
			else:
				if(iThetaEnd+iOff < Zones[iZoneEnd].gridZc.shape[1]):
					aa = np.where((WallTriangles.k == iZoneEnd) & (WallTriangles.i == ix) & (WallTriangles.j == iThetaEnd+iOff))[0]
				else:
					aa = np.where((WallTriangles.k == Zones[iZoneEnd].Neighbour.east) & (WallTriangles.i == ix) & (WallTriangles.j == iThetaEnd+iOff-Zones[iZoneEnd].gridZc.shape[1]))[0]
				if(len(aa) > 0): break

				if(iThetaEnd-iOff >= 0):
					aa = np.where((WallTriangles.k == iZoneEnd) & (WallTriangles.i == ix) & (WallTriangles.j == iThetaEnd-iOff))[0]
				else:
					aa = np.where((WallTriangles.k == Zones[iZoneEnd].Neighbour.west) & (WallTriangles.i == ix) & (WallTriangles.j == Zones[Zones[iZoneEnd].Neighbour.west].gridZc.shape[1]+iThetaEnd-iOff))[0]
				
			if(len(aa) > 0): break

		if(len(aa) > 1):
			print("Attention more than one wall triangle in soledge quadrangle: triangles=",aa+1)

		SpCell[iCell,3] = aa[0]

	return dist, SpCell, SpPos
