#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot1d_on_pol_mesh_cut 		import plot1d_on_pol_mesh_cut
	from routines.cli_routines			import *

	if(cli_present("-h",sys.argv)):
		print("\nThis script plot plasma parameters on soledge mesh along a poloidal flux surface\n")
		print("plot1d_on_pol_mesh_cut options")
		print("\t-path        Directory with simulation or list of directories[d='./']")
		print("\t-evolution   If > 0 plot evolution file (or evolutions files with list [3,4,5])  [d=0]")
		print("\t-path_label  Labels to itentify runs in plot [d='']")
		print("\t-rz0_line    Core coordinate direction line (Es.: [3.,0.]) [D=[]")
		print("\t-theta_line  Angular (deg) direction line [D=[0.]")
		print("\t-d_from_sep  Distance from separatrix [d=0.01]")
		print("\t-l_pol       Use polidal length instead than parallel length [d=false]")
		print("\t-log_scale   Use log scale for y axis [d=false]")
		print("\t-d_only      plot only e- and D+ [d=false]")
		print("\t-diff        Plot difference between two evolutions[d=false]")
		print("\t-no_samex    No same x scale [d=false]")
		print("\t-path_label  Labels to itentify runs in plot [d='']")
		print("\t-no_label    Skip labels in plots [d=false]")
		print("\t-no_zero     No zero on y scale [d=false]")
		print("\t-extra_walls Flag to show extra walls [d=false]")
		print("\t-one_plot    One plot on each figure [d=false]")
		print("\t-save        Save figures on files, values=none/png/ps/eps/pdf [D='none']")
		print()
		exit()

	path	 	= cli_get_value("-path",			sys.argv,   [""])
	evolution	= cli_get_value("-evolution",		sys.argv,	[])
	path_label 	= cli_get_value("-path_label",		sys.argv,   [""])
	rz0_line	= cli_get_value("-rz0_line",		sys.argv,	[])
	theta_line	= cli_get_value("-theta_line",		sys.argv,	0.)
	d_from_sep	= cli_get_value("-d_from_sep",		sys.argv, 0.01)
	save		= cli_get_value("-save",			sys.argv, "none")
	l_pol		= cli_present("-l_pol",				sys.argv)
	log_scale	= cli_present("-log_scale",			sys.argv)
	d_only		= cli_present("-d_only",			sys.argv)
	diff		= cli_present("-diff",				sys.argv)
	no_samex	= cli_present("-no_samex",			sys.argv)
	no_zero		= cli_present("-no_zero",			sys.argv)
	path_label 	= cli_get_value("-path_label",		sys.argv,   [""])
	no_labels	= cli_present("-no_labels",			sys.argv)
	extra_walls		= cli_present("-extra_walls",		sys.argv)
	one_plot	= cli_present("-one_plot",			sys.argv)
	plot1d_on_pol_mesh_cut(path=path, evolution=evolution, rz0_line=rz0_line, theta_line=theta_line, d_from_sep=d_from_sep, l_pol=l_pol, log_scale=log_scale, d_only=d_only, no_samex=no_samex, no_zero=no_zero, path_label=path_label, no_labels=no_labels, diff=diff, extra_walls=extra_walls, one_plot=one_plot, save=save)
	exit()

#=======================================

# Function definition is here

import types
import os
import h5py
from math								import sqrt
import numpy							as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot				as pyp
from matplotlib.backends.backend_pdf	import PdfPages

from routines.h5_routines				import h5_read
from routines.intersect_contour			import intersect_2contours
from routines.utils_walls				import get_in_out_walls, plot2d_walls, get_dmax_points_walls
from routines.set_profile_from_filedata	import set_profile_from_filedata
from routines.globals					import DEBUG, KB, BALLOONING_NAMES

from mesh.get_rz_core_sep				import get_rz_core_sep
from mesh.get_rho_in_out_core_sep		import get_rho_in_out_core_sep
from mesh.find_zones_intersections		import find_zones_intersections
from mesh.compute_mesh_intersections	import compute_mesh_intersections

from files.load_soledge_mesh_file		import load_soledge_mesh_file
from files.load_plasma_files			import load_plasma_files
from files.load_fluxes_files			import load_fluxes_files
from files.load_exp_data				import load_exp_data
from files.load_ions_list				import load_ions_list
from files.load_text_data				import load_text_data
from files.load_refpar_file				import load_refpar_file

#==============================================================================
# This routine plots ne/ni and Te/Ti soledge mesh inpolidal direction
#==============================================================================

def plot1d_on_pol_mesh_cut(path=[], evolution=[], rz0_line = [2.,0.], theta_line=5., d_from_sep=0.01, l_pol=0, log_scale=0, d_only=0, no_samex=0, no_zero=0, path_label=[], no_labels=0, diff=0, extra_walls=0, one_plot=0, save="none"):

	matplotlib_ver = matplotlib.__version__
	
	print("plot1d_on_pol_mesh_cut")

	if(diff != 0):
		if((evolution == 0) or (len(evolution) != 2)):
			print("\tWith diff option two evolutions must be provided")
			exit()
		log_scale = 0

	if(len(evolution) == 0): evolution = [0]
	if(len(path) == 0):		 path = [""]
	elif(len(path) > 1):	 evolution = [evolution[0]]

	i_plot_file = 0

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
	if(d_only != 0): ions = ions[0:2]
	iPlasmas = [i for i in range(len(ions))]

#	Read mesh

	Config = load_soledge_mesh_file(path0+"mesh.h5")
	Zones	= Config.Zones

#	Read Metric

	if(l_pol == 0):
		if_metric = h5py.File(path0+"Results/metric", "r")
		Gmet = []
		for k in range(len(Zones)):
			zone = "zone{:d}".format(k+1)
			Gmet.append(h5_read(if_metric, zone+ '/G', order = 'F'))									#[Nx+2,Nz+2]

		if_metric.close()
	

#	Find mesh along line

	if(len(rz0_line) == 0):
		Rcore, Zcore, CoreMegazone = get_rz_core_sep(Config, core_and_sep = False)
		rz0_line = [0.5*(Rcore.min() + Rcore.max()), 0.]

	rMax		= 6*get_dmax_points_walls(Config, rz0_line[0],rz0_line[1])
	theta_line	= theta_line*np.pi/180.
	RZLine		= np.array([[rz0_line[0],						  rz0_line[1]], \
							[rz0_line[0]+rMax*np.cos(theta_line), rz0_line[1]+rMax*np.sin(theta_line)]])


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

#	Select cell

	iCell = np.argmin(np.abs(dist-d_from_sep))
	iZone 	= IntCEll[iCell,0]
	ix		= IntCEll[iCell,1]
	iTheta	= IntCEll[iCell,2]

#	Find zones along poloidal coordinate

	if(Zones[iZone].Chi[ix,-1] == 1):
		iThEast = np.array([np.min(np.where(Zones[iZone].Chi[ix,iTheta:] == 1)[0])+iTheta])
		East = -1
	else:
		iThEast = np.array([Zones[iZone].Chi.shape[1]])
		East = Zones[iZone].Neighbour.east
	
	if(Zones[iZone].Chi[ix,0] == 1):
		iThWest = np.array([np.max(np.where(Zones[iZone].Chi[ix,:iTheta] == 1)[0])])
		West = -1
	else:
		iThWest = np.array([0])
		West = Zones[iZone].Neighbour.west

	iThetaOff  = iTheta - iThWest[0]
	nThetaPts  = iThEast[0] - iThWest[0]
	iZones	   = np.array([iZone])

#	Look East

	while (East > -1):
		iZones = np.append(iZones, East)
		iThWest = np.append(iThWest,0)
		if(Zones[East].Chi[ix,-1] == 1):
			iThEast = np.append(iThEast, np.min(np.where(Zones[East].Chi[ix,:] == 1)[0]))
			East = -1
		else:
			iThEast = np.append(iThEast, Zones[East].Chi.shape[1])
			East 	 = Zones[East].Neighbour.east
		nThetaPts += iThEast[-1]

#	Look West

	while (West > -1):
		iZones = np.append(West, iZones)
		iThEast = np.append(Zones[West].Chi.shape[1], iThEast)
		if(Zones[West].Chi[ix,0] == 1):
			iThWest = np.append(np.max(np.where(Zones[West].Chi[ix,:] == 1)[0])+1, iThWest)
			West = -1
		else:
			iThWest = np.append(0, iThWest)
			West = Zones[West].Neighbour.west
		iThetaOff += iThEast[0] - iThWest[0]
		nThetaPts += iThEast[0] - iThWest[0]

	Rpol = np.empty((nThetaPts), dtype = 'f8')
	Zpol = np.empty((nThetaPts), dtype = 'f8')
	jOff = 0
	for k in range(len(iZones)):
		Rpol[jOff: jOff + iThEast[k] - iThWest[k]] = Zones[iZones[k]].gridRc[ix, iThWest[k]:iThEast[k]]
		Zpol[jOff: jOff + iThEast[k] - iThWest[k]] = Zones[iZones[k]].gridZc[ix, iThWest[k]:iThEast[k]]
		jOff += iThEast[k] - iThWest[k]

	if(l_pol == 0):
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
		Lpara = np.sqrt((Rpol[1:]-Rpol[:-1])**2 + (Zpol[1:]-Zpol[:-1])**2)
		Lpara = np.append(0., np.cumsum(np.sqrt((Rpol[1:]-Rpol[:-1])**2 + (Zpol[1:]-Zpol[:-1])**2)))

	Lpara = Lpara - Lpara[iThetaOff]

#	Read and plot parameters

	if(save == "pdf"):	pdf = PdfPages("plot1d_mesh_t={:.3f}.".format(RefPar.time)+save)   #pdf in one file only

	colors = ['b','g','r','c','m','y','b','g','r','c','m','y']

#	Load Plasma parameters

	Evolutions = []
	for iPh in range(len(path)):
		for iEv in range(len(evolution)):
			path0  = path[iPh]
			if((len(path0) > 0) and (path0[-1] != "/")): path0 = path0 + "/"
			Plasmas = load_plasma_files(path0, nZones=len(Config.Zones), Evolution=evolution[iEv], iPlasmas=iPlasmas)

			Evolutions.append(Plasmas)

#	Plot parameters
#	#######################

	if(l_pol == 0):
		xLabel =  "$L_{par}\ (m)$"
	else:
		xLabel	=  "$L_{pol}\ (m)$"

	nPlasmas = len(Evolutions[0])
	for iPlasma in range(nPlasmas):

		if(iPlasma == 0):
			Labels = [["$n\ (*10^{19}\ m^{-3})$", "$T\ (keV)$",  "$Z_{eff}$", 
						"$Rad_{tot} (kW/m^3)$", "$Rad_{neu} (kW/m^3)$"]]
			Pars   = [["Dens", "Temp","Zeff ","TotRad ","TotNRad "]]
			ParIsPositive = [[True, True, True, True, True]]
			NrowNcol = [[3,2]]
			Facts		  	  = [[1e-19,1.e-3, 1., 1e-3,1.e-3]]
			FactsUnits		  = [["10^19","10^3", "1.","10^3","10^3"]]
			MaxRange	      = [[1.e5,0., 0., 1.e5, 1.e5]]

		elif(iPlasma == 1):
			Labels = [["$n\ (*10^{19}\ m^{-3})$", "$T\ (keV)$",  "$M$", "$Rad_{tot} (kW/m^3)$", "$Rad_{neu} (kW/m^3)$", "$Rad_{ion} (kW/m^3)$"], \
					  ["$S_{ni}\ (*10^{19}\ m^{-3}s^{-1})$","$S_{Ei}\ (MW/m^{3})$","$P_n\ (Pa)$", "$P_p\ (P)$", "$PE_p\ (MJ/m^3)$", "$n_H/n_e$"]]
			Pars   = [["Dens", "Temp", "M", "TotRad", "NRad", "IRad"], \
					  ["Sn",    "SE", "Pn", "Pp", "Ep","FracDens"]]
			ParIsPositive = [[True, True, False, True, True, True], \
					  [True,    True, True, True, True, False]]
			NrowNcol = [[3,2],[3,2]]
			Facts  	 		  = [[1e-19, 1.e-3, 1., 1e-3,1.e-3,1.e-3],[1.e-19, 1.e-6,1., 1.,1.e-6, 1.]]
			FactsUnits 		  = [["10^19", "10^3", "1.", "10^3","10^3","10^3"],["10^19", "10^6","1.", "1.","10^6", "1."]]
			MaxRange		  = [[1.e5, 0., 0., 1.e5, 1.e5, 1.e5],[0., 0.,1.e3, 0.,0.,0.]]

		elif(Plasmas[iPlasma][0].charge == 1):
			Labels = [["$n\ (*10^{19}\ m^{-3})$", "$T\ (keV)$",  "$Sn\ (*10^{19}\ m^{-3}s^{-1})$", "$Rad_{tot} (kW/m^3)$", "$Rad_{neu} (kW/m^3)$", "$Rad_{ion} (kW/m^3)$","$Tot-n_Z\ (*10^{19}\ m^{-3})$", "$C_{imp}$"]]
			Pars   = [["Dens", "Temp", "Sn", "TotRad", "NRad", "IRad","FracDens","Cimp"]]
			ParIsPositive = [[True, True, True, True, True, True, True, True]]
			NrowNcol = [[4,2]]
			Facts  	 = [[1e-19,1.e-3, 1.e-19, 1e-3,1.e-3,1.e-3, 1.e-19, 1.]]
			FactsUnits = [["10^19","10^3","10^19","10^3","10^3","10^3","1.","1."]]
			MaxRange	      = [[1.e5,0., 0., 1.e5, 1.e5,1.e5, 0.,0.]]

		else:
			Labels   = [["$n\ (*10^{19}\ m^{-3})$", "$T\ (keV)$",  "$Rad_{ion} (kW/m^3)$"]]
			Pars     = [["Dens", "Temp", "IRad"]]
			ParIsPositive = [[True, True, True]]
			NrowNcol = [[2,2]]
			Facts  	 = [[1e-19,1.e-3, 1e-3]]
			FactsUnits = [["10^19","10^3","10^3"]]
			MaxRange	      = [[1.e5,0., 0., 1.e5]]

		if(iPlasma == 0):	NameExt = "e"
		elif(iPlasma == 1):	NameExt = "i"
		else:				NameExt = Plasmas[iPlasma][0].ion[:-1]
		for iF in range(len(Pars)):
			for i in range(len(Pars[iF])): 
				if(Pars[iF][i][-1]!= " "): Pars[iF][i] += NameExt
				else:						Pars[iF][i]  = Pars[iF][i][:-1]


		if(iPlasma == 0):

#			Plot poloidal
#			#######################

			Fig = []
			Ax  = []
			Fig.append(pyp.figure())
			if(one_plot != 1):
				Ax.append(Fig[-1].add_subplot(NrowNcol[0][1], NrowNcol[0][0], 1))
			else:
				Ax.append(Fig[i].add_subplot(111))

			ip = 0
			Ax[ip].set_aspect(1.)
			Ax[ip].autoscale(enable=True, axis='both', tight=True)
			Ax[ip].set_xlabel("$R\ (m)$")
			Ax[ip].set_ylabel("$Z\ (m)$")
			if(len(path) == 1): Ax[0].set_title(path_label[0])

			plot2d_walls(Ax[ip], Config.Walls, extra_wall=extra_walls)
			Ax[ip].plot(RZcore[:,0],  RZcore[:,1],  'b-')
			Ax[ip].plot(RZsep[:,0],   RZsep[:,1],  'g-')

			Ax[ip].plot(IntRZ[Out_Sep,0], IntRZ[Out_Sep,1], 'g.-')
			Ax[ip].plot(IntRZ[In_Sep,0],  IntRZ[In_Sep,1],  'b.-')

			Ax[ip].plot(Rpol, Zpol, 'r.-')

			Ax[ip].text(0.5*(Rcore.min() + Rcore.max()), 0., "d={:.3f}".format(d_from_sep), horizontalalignment="center",verticalalignment="center")

			Ax[ip].plot(Rpol[0],  Zpol[0],  "ro")
			Ax[ip].plot(Rpol[-1], Zpol[-1], "ro")

			Ax[ip].text(Rpol[0],  Zpol[0],  "{:d}".format(int(Lpara[0])),  horizontalalignment="center",verticalalignment="top")
			Ax[ip].text(Rpol[-1], Zpol[-1], "{:d}".format(int(Lpara[-1])), horizontalalignment="center",verticalalignment="top")

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
				if((iPlasma == 0) and (iF == 0)):
					skip_first=1
				else:
					skip_first=0
					Fig.append(pyp.figure())
				Fig[-1].patch.set_facecolor('white')
				for i in range(len(Labels[iF])): 
					if(((skip_first != 0) and (i == 0)) or (no_samex != 0)): Ax.append(Fig[-1].add_subplot(NrowNcol[iF][1], NrowNcol[iF][0], i+1+skip_first))
					else: 									  Ax.append(Fig[-1].add_subplot(NrowNcol[iF][1], NrowNcol[iF][0], i+1+skip_first, sharex = Ax[1]))
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
				AxPl[iF][i].autoscale(enable=True, axis='both', tight=True)
#				if(i == 0): AxPl[iF][i].set_title(path_label+" @ t={:.3f} s".format(tempus))
				if(diff == 0):
					if(evolution[0] == 0):
						if(i == 1): AxPl[iF][i].set_title(Plasmas[iPlasma][0].ion)
					else:
						if(i == 1): AxPl[iF][i].set_title("Evol.={:d}".format(evolution[0]))
						if(i == 2): AxPl[iF][i].set_title(Plasmas[iPlasma][0].ion)
				else:
					if(i == 1): AxPl[iF][i].set_title("Ev_{:d}-Ev_{:d}".format(evolution[1],evolution[0]))
					if(i == 2): AxPl[iF][i].set_title(Plasmas[iPlasma][0].ion)

				AxPl[iF][i].set_xlabel(xLabel)
				AxPl[iF][i].set_ylabel(Labels[iF][i])
				if((log_scale != 0) and  (diff == 0)  and ParIsPositive[iF][i]):
					AxPl[iF][i].set_yscale('log')
				else:
					AxPl[iF][i].set_yscale('linear')

	#	Plot parameters

		for iF in range(len(Labels)):
			for i in range(len(Labels[iF])):
				try:
					iPar = Evolutions[0][iPlasma][0].VNames.index(Pars[iF][i])	
				except:
					print("\tWARNING: not found parameter ",Pars[iF][i])
					print("\t\tVNames=",Evolutions[0][iPlasma][0].VNames)
					continue

				if((len(Evolutions) < 2) or (diff != 0)):
					Values = get_plasma_parameter_on_pol(Evolutions[0][iPlasma], iPar, ix, iZones, iThWest, iThEast, nThetaPts)*Facts[iF][i]
					if(diff != 0): Values  = get_plasma_parameter_on_pol(Evolutions[1][iPlasma],   iPar, ix, iZones, iThWest, iThEast, nThetaPts)*Facts[iF][i] - Values
				else:
					Values = []
					for iv in range(len(Evolutions)):
						Values.append(get_plasma_parameter_on_pol(Evolutions[iv][iPlasma], iPar, ix, iZones, iThWest, iThEast, nThetaPts)*Facts[iF][i])

				if((len(evolution) < 2) or (diff != 0)):
					if((len(path) < 2) or (diff != 0)):
						AxPl[iF][i].plot(Lpara,  Values,  'b-')
					else:
						for iPh in range(len(path)):
							if(no_labels == 0):	AxPl[iF][i].plot(Lpara,  Values[iPh],  '-', color = colors[iPh], label=path_label[iPh])
							else:				AxPl[iF][i].plot(Lpara,  Values[iPh],  '-', color = colors[iPh])
				else:
					for iEv in range(len(Evolutions)):
						AxPl[iF][i].plot(Lpara,   Values[iEv], '-', color = colors[iEv], label="{:d}".format(evolution[iEv]))

				if((log_scale == 0) and (diff == 0) and (no_zero == 0) and ParIsPositive[iF][i]):
					AxPl[iF][i].set_ylim(bottom=0.)

				AxPl[iF][i].axvline(x=0., color='k', linestyle='dashed')
				if((len(path) > 1) or (len(evolution) > 1)): 
					AxPl[iF][i].legend(fontsize='small', loc='lower left')

	if(save != "none"):
		for i in range(len(Fig)):
			i_plot_file += 1
			if(one_plot != 1): Fig[i].set_size_inches(10.05,7.44)
			if(save == "pdf"):
				pdf.savefig(Fig[i])
			else:
				Fig[i].savefig("plot1d_on_pol_mesh_cut_t={:.3f}_{:d}.".format(RefPar.time,i_plot_file)+save)

#		pyp.show(block=False)
#		pyp.close()
	else:
		pyp.show()

	if(save == "pdf"):	pdf.close()

	print("plot1d_on_pol_mesh_cut: Completed")

	return


#===================================================================

def get_plasma_parameter_on_pol(Plasma, iPar, ix, iZones, iThWest, iThEast, nThetaPts):

	Par = np.zeros(nThetaPts, dtype='f8')

	if(iPar > -1):
		jOff = 0
		for k in range(len(iZones)):
			if(Plasma[iZones[k]].Nz < Plasma[iZones[k]].Values[iPar].shape[1]):
				Par[jOff: jOff + iThEast[k] - iThWest[k]] = Plasma[iZones[k]].Values[iPar][ix, iThWest[k]+1:iThEast[k]+1]		#[Nx+2,Nz+2] parameter
			else:
				Par[jOff: jOff + iThEast[k] - iThWest[k]] = Plasma[iZones[k]].Values[iPar][ix, iThWest[k]:iThEast[k]]				#[Nx,Nz] parameter

			jOff += iThEast[k] - iThWest[k]

	return Par
