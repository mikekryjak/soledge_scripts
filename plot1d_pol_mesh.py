#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot1d_pol_mesh 		import plot1d_pol_mesh
	from routines.cli_routines	import *

	if(cli_present("-h",sys.argv)):
		print("\nThis script plot plasma parameters on soledge mesh along a poloidal flux surface\n")
		print("plot1d_pol_mesh options")
		print("\t-path        Directory with simulation or list of directories[d='./']")
		print("\t-evolution   If > 0 plot evolution file (or evolutions files with list [3,4,5])  [d=0]")
		print("\t-rz0_line    Core coordinate direction line (Es.: [3.,0.]) [D=[]")
		print("\t-theta_line  Angular (deg) direction line [D=[0.]")
		print("\t-d_from_sep  Distance from separatrix [d=0.01]")
		print("\t-l_pol       Use polidal length instead than parallel length [d=false]")
		print("\t-log_scale   Use log scale for y axis [d=false]")
		print("\t-one_plot    One plot on each figure [d=false]")
		print("\t-save        Save figures on files, values=none/png/ps/eps/pdf [D='none']")
		print()
		exit()

	path	 	= cli_get_value("-path",				sys.argv,   [""])
	evolution	= cli_get_value("-evolution",		sys.argv,	[])
	rz0_line	= cli_get_value("-rz0_line",			sys.argv,	[])
	theta_line	= cli_get_value("-theta_line",		sys.argv,	0.)
	d_from_sep	= cli_get_value("-d_from_sep",		sys.argv, 0.01)
	save		= cli_get_value("-save",				sys.argv, "none")
	l_pol		= cli_present("-l_pol",				sys.argv)
	log_scale	= cli_present("-log_scale",			sys.argv)
	one_plot	= cli_present("-one_plot",			sys.argv)
	plot1d_pol_mesh(path=path, evolution=evolution, rz0_line=rz0_line, theta_line=theta_line, d_from_sep=d_from_sep, l_pol=l_pol, log_scale=log_scale, one_plot=one_plot, save=save)
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
from routines.utils_walls				import get_in_out_walls, plot2d_walls
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

def plot1d_pol_mesh(path=[], evolution=[], rz0_line = [2.,0.], theta_line=5., d_from_sep=0.01, l_pol=0, log_scale=0, rho_scale=0, plot_fluxes=0, one_plot=0, save="none"):

	matplotlib_ver = matplotlib.__version__
	
	print("plot1d_pol_mesh")

	if(len(evolution) == 0): evolution = [0]
	if(len(path) == 0):		 path = [""]
	elif(len(path) > 1):	 evolution = [evolution[0]]

	print("path=",path)
	i_plot_file = 0

#	Read reference parameters

	
	path0  = path[0]
	if((len(path0) > 0) and (path0[-1] != "/")): path0 = path0 + "/"

	RefPar = load_refpar_file(path0+"Results/")

	ions = load_ions_list(path0)

#	Read mesh

	Config = load_soledge_mesh_file(path0+"mesh.h5")
	Zones	= Config.Zones

#	Read Metric

	if(l_pol == 0):
		if_metric = h5py.File(path0+"Results/metric", "r")
		Gmet = []
		for k in range(len(Zones)):
			zone = "zone{:d}".format(k+1)
			Gmet.append(h5_read(if_metric, zone+ '/G', order = 'F'))

		if_metric.close()
	

#	Find mesh along line

	if(len(rz0_line) == 0):
		Rcore, Zcore, CoreMegazone = get_rz_core_sep(Config, core_and_sep = False)
		rz0_line = [0.5*(Rcore.min() + Rcore.max()), 0.]

	rMax		= 6*get_dmax_points_walls(Config, rz0_line[0], z0_line[1], plasma_wall=True, eirene_wall=False, extra_wall=False)
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
			dlZone	   = -2.*dtheta/Gmet[iZones[k]][ix+1, 1:-1]
			dl[jOff: jOff + iThEast[k] - iThWest[k]] = dlZone[iThWest[k]:iThEast[k]]
			jOff += iThEast[k] - iThWest[k]
		dlZone = 0
		Lpara = np.cumsum(dl)
		dl	  = 0
	else:
		Lpara = np.append(0., np.cumsum(np.sqrt((Rpol[1:]-Rpol[:-1])**2 + (Zpol[1:]-Zpol[:-1])**2)))

	Lpara = Lpara - Lpara[iThetaOff]

#	Read and plot parameters

	if(save == "pdf"):	pdf = PdfPages("plot1d_mesh_t={:.3f}.".format(RefPar.time)+save)   #pdf in one file only


	if(l_pol == 0):
		xLabels		= [[ "$L_{par}\ (m)$",	 "$l_{par}\ (m)$", "$l_{par}\ (m)$"], \
					   [ "$l_{par}\ (m)$",	 "$l_{par}\ (m)$", "$l_{par}\ (m)$",  "$l_{par}\ (m)$"]]
	else:
		xLabels		= [[ "$L_{pol}\ (m)$",	 "$l_{pol}\ (m)$", "$l_{pol}\ (m)$"], \
					   [ "$l_{pol}\ (m)$",	 "$l_{pol}\ (m)$", "$l_{pol}\ (m)$",  "$l_{pol}\ (m)$"]]
		
	yLabels     = [[ "$n\ (*10^{19}\ m^{-3})$",	  "$T\ (eV)$",	"$Pp\ (kP)$" ], \
				   [ "$n\ (*10^{19}\ m^{-3})$",	  "$T\ (eV)$",	"$Mach\ number$",	 "$Rad\ (kW)$" ]]

	yLogLabels   = [[ "$Ln(n)\ (*10^{19}\ m^{-3})$",	  "$Ln(T)\ (eV)$",	"$Ln(Pp)\ (kP)$" ], \
					[ "$Ln(n)\ (*10^{19}\ m^{-3})$",	  "$Ln(T)\ (eV)$",	"$Mach\ number$",	 "$Rad\ (kW)$" ]]

	LogScales	= [["log", "log", "log"], ["log", "log", "linear", "log"]]

	iValues		= [[0,1,-1],[0,1,3,4]]

	Facts		= [[1e-19, 1., 1.e16,],[1e-19, 1., 1., 1.e-3]]
	BottomZero	= [[True,True,True],[True,True,False,True]]
	VPosPlot	= [[1,2,3],[4,5,6,7]]
	EPosPlot	= [0]

#	repeat for all ions like the first one

	if(len(ions) > 2):
		for iPlasma in range(2,len(ions)):
			xLabels.append(xLabels[1]) 
			yLabels.append(yLabels[1]) 
			yLogLabels.append(yLogLabels[1]) 
			LogScales.append(LogScales[1]) 
			iValues.append(iValues[1]) 
			Facts.append(Facts[1]) 
			BottomZero.append(BottomZero[1]) 
			VPosPlot.append(VPosPlot[1][:]) 
			for i in range(len(VPosPlot[-1])): VPosPlot[-1][i] += VPosPlot[-2][-1] + 1

	colors = ['b','g','r','c','m','y','b','g','r','c','m','y']

#	Load Plasma parameters

	Values = []
	for iPh in range(len(path)):
		for iEv in range(len(evolution)):
			Values.append([])
			for iPlasma in range(len(ions)):
				Values[-1].append([])
				path0  = path[iPh]
				if((len(path0) > 0) and (path0[-1] != "/")): path0 = path0 + "/"
				Plasmas = load_plasma_files(path0, nZones=len(Config.Zones), iPlasmas = [iPlasma], Evolution=evolution[iEv])

				for i in range(len(iValues[iPlasma])):
					Values[-1][iPlasma].append(get_plasma_parameter_on_pol(Plasmas[0], iValues[iPlasma][i], ix, iZones, iThWest, iThEast, nThetaPts)*Facts[iPlasma][i])

#	Compute some parameters

			Values[-1][0][2] = 1.6e-19*Values[-1][0][0]*(Values[-1][0][1] + Values[-1][1][1])*(1.+ Values[-1][1][2]**2)*Facts[0][2]			#Plasma_pressure Pp=1.6e-19*density*(Te+Ti)*(1+M**2)


#	Prepare for plotting

	nRows		= 2
	nCols		= 3
	PlotPerFig	= nRows*nCols
	nPLots		= len(EPosPlot) + len(VPosPlot[0]) + len(VPosPlot[1])*(len(ions)-1)
	nFigs		= int(nPLots/PlotPerFig)
	if(nFigs*PlotPerFig < nPLots): nFigs += 1

	Fig = []
	Ax  = []
	if(one_plot != 1):
		for i in range(nFigs):	
			Fig.append(pyp.figure())
			for i in range(min(PlotPerFig,nPLots-i*PlotPerFig)):
				Ax.append(Fig[-1].add_subplot(nRows,nCols,i+1))
				Ax[-1].locator_params(axis='x',nbins=4)

			Fig[-1].tight_layout(pad=2., w_pad=3., h_pad=3.)
	else:
		for i in range(nPLots):
			Fig.append(pyp.figure())
			Ax.append(Fig[i].add_subplot(111))

	for ip in EPosPlot: 
		Ax[ip].set_aspect(1.)
		Ax[ip].autoscale(enable=True, axis='both', tight=True)
		Ax[ip].set_xlabel("$R\ (m)$")
		Ax[ip].set_ylabel("$Z\ (m)$")
		if(len(path) == 1):
			Ax[ip].set_title(os.path.basename(os.path.abspath(path[0]))+" @ t={:.3f} s".format(RefPar.time))

	for iPlasma in range(len(ions)):
		for i in range(len(VPosPlot[iPlasma])):
			ip = VPosPlot[iPlasma][i]
			if((i > 0) or (len(evolution) != 1) or (evolution[0] == 0)):
				Ax[ip].set_title(ions[iPlasma])
			else:
				Ax[ip].set_title(ions[iPlasma]+" evol. {:d}".format(evolution[0]))

			Ax[ip].autoscale(enable=True, axis='both', tight=True)
			Ax[ip].set_xlabel(xLabels[iPlasma][i])
			if(log_scale == 0):
				Ax[ip].set_ylabel(yLabels[iPlasma][i])
				Ax[ip].set_yscale('linear')
			else:
				Ax[ip].set_ylabel(yLogLabels[iPlasma][i])
				Ax[ip].set_yscale(LogScales[iPlasma][i])

#	Plot poloidal

	for ip in EPosPlot: 
		plot2d_walls(Ax[ip], Config.Walls)

		Ax[ip].plot(RZcore[:,0],  RZcore[:,1],  'b-')
		Ax[ip].plot(RZsep[:,0],   RZsep[:,1],  'g-')

		Ax[ip].plot(IntRZ[Out_Sep,0], IntRZ[Out_Sep,1], 'g.-')
		Ax[ip].plot(IntRZ[In_Sep,0],  IntRZ[In_Sep,1],  'b.-')

		Ax[ip].plot(Rpol, Zpol, 'r.-')

		Ax[ip].text(0.5*(Rcore.min() + Rcore.max()), 0., "d={:.3f}".format(d_from_sep), horizontalalignment="center",verticalalignment="center")

#		Plot parameters

	for iPlasma in range(len(ions)):
		for i in range(len(VPosPlot[iPlasma])):
			if(len(evolution) < 2):
				if(len(path) < 2):
					Ax[VPosPlot[iPlasma][i]].plot(Lpara,  Values[0][iPlasma][i],  'b.-')
				else:
					for iPh in range(len(path)):
						Ax[VPosPlot[iPlasma][i]].plot(Lpara,  Values[iPh][iPlasma][i],  '-', color = colors[iPh], label=path[iPh])
			else:
				for iEv in range(len(evolution)):
					Ax[VPosPlot[iPlasma][i]].plot(Lpara,  Values[iEv][iPlasma][i],  '-', color = colors[iEv], label="{:d}".format(evolution[iEv]))

			if((log_scale == 0) and BottomZero[iPlasma][i]):
				Ax[VPosPlot[iPlasma][i]].set_ylim(bottom=0.)

		xv = 0.
		for i in range(len(VPosPlot[iPlasma])):
#			pyp.setp(Ax[k].yaxis.get_ticklabels(), rotation=90.)
			Ax[VPosPlot[iPlasma][i]].axvline(x=xv, color='k', linestyle='dashed')
			if((len(path) > 1) or (len(evolution) > 1)): 
				Ax[VPosPlot[iPlasma][i]].legend(fontsize='small', loc='lower left')

	if(save != "none"):
		for i in range(len(Fig)):
			i_plot_file += 1
			if(one_plot != 1): Fig[i].set_size_inches(10.05,7.44)
			if(save == "pdf"):
				pdf.savefig(Fig[i])
			else:
				Fig[i].savefig("plot1d_pol_mesh_t={:.3f}_{:d}.".format(RefPar.time,i_plot_file)+save)

		pyp.show(block=False)
		pyp.close()
	else:
		pyp.show()

	if(save == "pdf"):	pdf.close()

	print("plot1d_pol_mesh: Completed")

	return


#===================================================================

def get_plasma_parameter_on_pol(Plasma, iPar, ix, iZones, iThWest, iThEast, nThetaPts):

	Par = np.zeros(nThetaPts, dtype='f8')

	if(iPar > -1):
		jOff = 0
		for k in range(len(iZones)):
			if(Plasma[iZones[k]].Nz < Plasma[iZones[k]].Values[iPar].shape[1]):
				Par[jOff: jOff + iThEast[k] - iThWest[k]] = Plasma[iZones[k]].Values[iPar][ix, iThWest[k]+1:iThEast[k]+1]
			else:
				Par[jOff: jOff + iThEast[k] - iThWest[k]] = Plasma[iZones[k]].Values[iPar][ix, iThWest[k]:iThEast[k]]

			jOff += iThEast[k] - iThWest[k]

	return Par
