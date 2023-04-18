#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot1d_sp_on_mesh_cut 		import plot1d_sp_on_mesh_cut
	from routines.cli_routines	import *

	if(cli_present("-h",sys.argv)):
		print("\nThis script plot plasma parameters on soledge mesh along a poloidal flux surface\n")
		print("plot1d_sp_on_mesh_cut options")
		print("\t-path        Directory with simulation or list of directories[d='./']")
		print("\t-evolution   If > 0 plot evolution file (or evolutions files with list [3,4,5])  [d=0]")
		print("\t-rz0_line    Core coordinate direction line (Es.: [3.,0.]) [D=[]")
		print("\t-theta_line  Angular (deg) direction line [D=[0.]")
		print("\t-d_from_sep  Max distance from separatrix [d=0.01]")
		print("\t-east_west   Stryke points on east=-1, west=1 or shorted distance from cut [d=0]")
		print("\t-par_fluxes  plot parallel fluxes [d=false]")
		print("\t-path_label  Labels to itentify runs in plot [d='']")
		print("\t-no_labels   Skip labels in plots [d=false]")
		print("\t-d_only      plot only e- and D+ [d=false]")
		print("\t-log_scale   Use log scale for y axis [d=false]")
		print("\t-one_plot    One plot on each figure [d=false]")
		print("\t-save        Save figures on files, values=none/png/ps/eps/pdf/csv [D='none']")
		print()
		exit()

	path	 	= cli_get_value("-path",				sys.argv,   [""])
	evolution	= cli_get_value("-evolution",		sys.argv,	[])
	rz0_line	= cli_get_value("-rz0_line",			sys.argv,	[])
	theta_line	= cli_get_value("-theta_line",		sys.argv,	0.)
	d_from_sep	= cli_get_value("-d_from_sep",		sys.argv, 0.01)
	east_west	= cli_get_value("-east_west",		sys.argv, 0)
	save		= cli_get_value("-save",				sys.argv, "none")
	par_fluxes	= cli_present("-par_fluxes",		sys.argv)
	path_label 	= cli_get_value("-path_label",		sys.argv,   [""])
	no_labels	= cli_present("-no_labels",			sys.argv)
	d_only		= cli_present("-d_only",				sys.argv)
	log_scale	= cli_present("-log_scale",			sys.argv)
	one_plot	= cli_present("-one_plot",			sys.argv)
	plot1d_sp_on_mesh_cut(path=path, evolution=evolution, rz0_line=rz0_line, theta_line=theta_line, d_from_sep=d_from_sep, east_west=east_west, par_fluxes=par_fluxes, path_label=path_label, no_labels=no_labels, d_only=d_only, log_scale=log_scale, one_plot=one_plot, save=save)
	exit()

#=======================================

# Function definition is here

import types
import os
import h5py
from math								import sqrt, asin
import numpy							as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot				as pyp
from matplotlib.path 					import Path
from matplotlib.backends.backend_pdf	import PdfPages

from routines.h5_routines				import h5_read
from routines.intersect_contour			import intersect_2contours
from routines.utils_walls				import get_in_out_walls, get_min_max_walls, plot2d_walls, get_dmax_points_walls
from routines.set_profile_from_filedata	import set_profile_from_filedata
from routines.globals					import DEBUG, KB, BALLOONING_NAMES

from mesh.get_rz_core_sep				import get_rz_core_sep
from mesh.get_rho_in_out_core_sep		import get_rho_in_out_core_sep
from mesh.find_zones_intersections		import find_zones_intersections
from mesh.compute_mesh_intersections	import compute_mesh_intersections
from eirene.get_ext_plasma_triseq		import get_ext_plasma_triseq

from files.load_soledge_mesh_file		import load_soledge_mesh_file
from files.load_eirene_triangles		import load_eirene_triangles
from files.load_plasma_files			import load_plasma_files
from files.load_fluxes_files			import load_fluxes_files
from files.load_exp_data				import load_exp_data
from files.load_ions_list				import load_ions_list
from files.load_text_data				import load_text_data
from files.load_refpar_file				import load_refpar_file

#==============================================================================
# This routine plots ne/ni and Te/Ti soledge mesh inpolidal direction
#==============================================================================

def plot1d_sp_on_mesh_cut(path=[], evolution=[], rz0_line = [2.,0.], theta_line=5., d_from_sep=0.01, east_west=0, par_fluxes=0, path_label=[], no_labels=0, d_only=0, log_scale=0, one_plot=0, save="none"):

	
	print("plot1d_sp_on_mesh_cut")

#	Read reference parameters

	if(len(evolution) == 0): evolution = [0]
	if(len(path) == 0):		 path = [""]
	elif(len(path) > 1):	 evolution = [evolution[0]]

	if((len(path_label) < len(path)) or (len(path_label[0]) == 0)):
		path_label = []
		for in_path  in path: path_label.append(os.path.basename(os.path.abspath(in_path)))
		HasPathLabel=False
	else:
		HasPathLabel=True

	path0  = path[0]
	if((len(path0) > 0) and (path0[-1] != "/")): path0 = path0 + "/"

	RefPar = load_refpar_file(path0+"Results/")

	ions = load_ions_list(path0)
	if(d_only != 0): ions = ions[0:2]

#	Read mesh

	Config = load_soledge_mesh_file(path0+"mesh.h5")
	Zones	= Config.Zones

	Eirene = load_eirene_triangles(path0+"triangles.h5")
	get_ext_plasma_triseq(Config, Eirene.Wall)

	WallTriangles = Eirene.WallTriangles

#	Read Metric

	l_pol = 0
	if(l_pol == 0):
		if_metric = h5py.File(path0+"Results/metric", "r")
		Gmet = []
		for k in range(len(Zones)):
			zone = "zone{:d}".format(k+1)
			Gmet.append(h5_read(if_metric, zone+ '/G', order = 'F'))

		if_metric.close()

	X_points	= Config.X_points
	if(len(X_points) > 1):
		X_points_R = np.array([X_points[k].R for k in range(len(X_points))])
		X_points_Z = np.array([X_points[k].Z for k in range(len(X_points))])
		InXPoints, pts_out = get_in_out_walls(Config, X_points_R, X_points_Z)
	else:
		InXPoints = np.array([0])

	nInXPoints = len(InXPoints)
		
#	Find mesh along line

	if(len(rz0_line) == 0):
		Rcore, Zcore, CoreMegazone = get_rz_core_sep(Config, core_and_sep = False)
		rz0_line = [0.5*(Rcore.min() + Rcore.max()), Zcore[np.argmax(Rcore)]]

	rMax		= 6*get_dmax_points_walls(Config, rz0_line[0], rz0_line[1], plasma_wall=True, eirene_wall=False, extra_wall=False)
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

	if(d_from_sep == 0): d_from_sep = 2*dist.max()
	OutSep	= np.where((dist > 0) & (dist < d_from_sep))[0]
	dist	= dist[OutSep]
	IntRZ	= IntRZ[OutSep,:]
	IntCEll = IntCEll[OutSep,:]

	SpCell = np.empty((IntCEll.shape[0],4), dtype = 'i4')
	SpPos  = np.empty((IntCEll.shape[0],5), dtype = 'f8')

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
			SpCell[iCell,0] = iZones[-1]
			SpCell[iCell,1] = ix
			SpCell[iCell,2] = iThEast[-1]-1
			SpPos[iCell,2]	= -Lpara[-2]
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

#		Compute angle of incidence on target

		iWall = Eirene.Wall.iTriSeqExtPlasma
		iSeq = np.where(Eirene.Wall.TriSequences[iWall] == SpCell[iCell,3])[0][0]

		Dr   = Eirene.Wall.R12[iWall] [iSeq,1] - Eirene.Wall.R12[iWall] [iSeq,0]
		Dz   = Eirene.Wall.Z12[iWall] [iSeq,1] - Eirene.Wall.Z12[iWall] [iSeq,0]
		Br   = Zones[SpCell[iCell,0]].Br[ix, SpCell[iCell,2]]
		Bz   = Zones[SpCell[iCell,0]].Bz[ix, SpCell[iCell,2]]
		Bphi = Zones[SpCell[iCell,0]].Bphi[ix, SpCell[iCell,2]]
		
		SpPos[iCell,3] = abs(asin((-Dz*Br + Dr*Bz)/(sqrt(Dr**2+Dz**2)*sqrt(Br**2+Bz**2+Bphi**2))))

#		Flux line dist at cut

		DrF =   Zones[iZone].gridR[ix+1, iTheta] - Zones[iZone].gridR[ix, iTheta]
		DzF =   Zones[iZone].gridZ[ix+1, iTheta] - Zones[iZone].gridZ[ix, iTheta]

		DrL =   Zones[iZone].gridR[ix, iTheta+1] - Zones[iZone].gridR[ix, iTheta]
		DzL =   Zones[iZone].gridZ[ix, iTheta+1] - Zones[iZone].gridZ[ix, iTheta]

		FdistCut =  abs((DrL*DzF -DzL*DrF)/sqrt(DrL**2+DzL**2))

#		Flux dist at target

		DrF =   Zones[SpCell[iCell,0]].gridR[ix+1, SpCell[iCell,2]] - Zones[SpCell[iCell,0]].gridR[ix, SpCell[iCell,2]]
		DzF =   Zones[SpCell[iCell,0]].gridZ[ix+1, SpCell[iCell,2]] - Zones[SpCell[iCell,0]].gridZ[ix, SpCell[iCell,2]]

		DrL =   Zones[SpCell[iCell,0]].gridR[ix, SpCell[iCell,2]+1] - Zones[SpCell[iCell,0]].gridR[ix, SpCell[iCell,2]]
		DzL =   Zones[SpCell[iCell,0]].gridZ[ix, SpCell[iCell,2]+1] - Zones[SpCell[iCell,0]].gridZ[ix, SpCell[iCell,2]]

		FdistSp = abs((DrL*DzF -DzL*DrF)/sqrt(DrL**2+DzL**2))

		SpPos[iCell,4]  = FdistSp/FdistCut

#		print("Dr, Dz=",Dr, Dz)
#		print("Br, Bz, Bphi=",Br, Bz, Bphi)
#		print("SpPos[iCell,3]=",SpPos[iCell,3])

#	Read and plot parameters

	dist *= 1000.		#change dist in mm

	if(save == "pdf"):	pdf = PdfPages("plot1d_sp_on_mesh_cut_t={:.3f}.".format(RefPar.time)+save)   #pdf in one file only
	if(save == "csv"):
		if(len(path) > 1):
			save = "none"
			print("plot1d_sp_on_mesh_cut: cannot save to csv more than one profile!!!")
		else:
			csv = []
			csv.append(types.SimpleNamespace())
			csv[-1].Name   = "Dist (mm)"
			csv[-1].Values = dist


	xLabels	   		= [["$l\ (mm)$", "$l\ (mm)$", "$l\ (mm)$"],["$l\ (mm)$", "$l\ (mm)$", "$l\ (mm)$", "$l\ (mm)$"]]
	yLabels    		= [["$n\ (*10^{19}\ m^{-3})$", "$T\ (eV)$",	"$\Gamma_E\ (MW/m^{2})$"], \
					   ["$n\ (*10^{19}\ m^{-3})$", "$T\ (eV)$",	"$\Gamma_n\ (*10^{22}\ m^{-2})$", "$\Gamma_E\ (MW/m^{2})$"]]
	yLogLabels		= [["$Ln(n)\ (*10^{19}\ m^{-3})$", "$Ln(T)\ (eV)$",	"$Ln(\Gamma_E)\/(MW/m^{2})$"], \
					   ["$Ln(n)\ (*10^{19}\ m^{-3})$", "$Ln(T)\ (eV)$",	"$Ln(\Gamma_n)\/(*10^{22}\ m^{-2})$", "$Ln(\Gamma_E)\/(MW/m^{2})$"]]
	FLabels    		= [["n_{:s} (*10^19 m^-3)", "T_{:s} (eV)",	 "Flux_E_{:s} (MW/m^2)"], \
					   ["n_{:s} (*10^19 m^-3)", "T_{:s} (eV)",	"Flux_n_{:s} (*10^22 m^-2)", "Flux_E_{:s} (MW/m^2)"]]
	LogScales		= [["log", "log", "linear", "log"],["log", "log", "log"]]
	iPlasmaValues	= [[0,  1, 3],[0,  1, 4, 3]]
	Facts			= [[1e-19, 1., 1e-6],[1e-19, 1., 1.e-22,1e-6]]
	BottomZero		= [[True,True,True],[True,True,True,True]]

	VPosPlot	= [[1,2,9],[6,7,8,10]]											#Position of plasma plot
	EPosPlot    = [0,3,4,5]																#Position of extra plots


#	repeat for all ions like the first one

	if(len(ions) > 2):
		for iPlasma in range(2,len(ions)):
			xLabels.append(xLabels[1][:]) 
			yLabels.append(yLabels[1][:]) 
			yLogLabels.append(yLogLabels[1][:]) 
			LogScales.append(LogScales[1][:]) 
			iPlasmaValues.append(iPlasmaValues[1][:]) 
			Facts.append(Facts[1][:]) 
			BottomZero.append(BottomZero[1][:]) 
			VPosPlot.append(VPosPlot[1][:]) 
			for i in range(len(VPosPlot[-1])): VPosPlot[-1][i] = VPosPlot[-2][-1] + i+1

	xLabels[1].append("$l\ (mm)$")											#Add total output power
	yLabels[1].append("$\Gamma_{Etot}\ (MW/m^{2})$") 
	FLabels[1].append("Flux_Etot_{:s} (MW/m^2)") 
	yLogLabels[1].append("$Ln(\Gamma_{Etot})\/(MW/m^{2})$") 
	LogScales[1].append("log") 
	iPlasmaValues[1].append(8) 
	Facts[1].append(1e-6) 
	BottomZero[1].append(True) 
	VPosPlot[1].append(11) 

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

				if((iPh == 0) and (iEv == 0) and (iPlasma == 0)):
					for k in range(IntCEll.shape[0]):
						SpCell[k,3] = np.where(Plasmas[0][0].Wall.ntri == WallTriangles.ntri[SpCell[k,3]])[0][0]

				for i in range(len(iPlasmaValues[iPlasma])):
					ParValues = np.zeros(IntCEll.shape[0], dtype = 'f8')
					if(iPlasmaValues[iPlasma][i] > -1):
						for k in range(IntCEll.shape[0]): ParValues[k] = Plasmas[0][0].Wall.Values[iPlasmaValues[iPlasma][i]][SpCell[k,3]]
						if(((iPlasmaValues[iPlasma][i] == 3) or (iPlasmaValues[iPlasma][i] == 8)) and (par_fluxes == 1)): ParValues = ParValues/np.sin(SpPos[iCell,3])					#plot field parallel fluxes 
						Values[-1][-1].append(ParValues)


#	Prepare for plotting

	nRows		= 2
	nCols		= 3
	PlotPerFig	= nRows*nCols
	nPLots		= VPosPlot[-1][-1]+1
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

	for figure in Fig:  figure.patch.set_facecolor('white')

	for i in range(len(EPosPlot)):
		ip = EPosPlot[i]
		if(i == 0):															#plot section and cut
			if(len(path) == 1):
				Ax[ip].set_title(os.path.basename(os.path.abspath(path[0]))+" @ t={:.3f} s".format(RefPar.time))
			Ax[ip].set_aspect(1.)
			Ax[ip].set_xlabel("$R\ (m)$")
			Ax[ip].set_ylabel("$Z\ (m)$")

			plot2d_walls(Ax[ip], Config.Walls)

			Ax[ip].plot(RZcore[:,0],  RZcore[:,1],  'b-')
			Ax[ip].plot(RZsep[:,0],   RZsep[:,1],  'g-')

			Ax[ip].plot(IntRZ[:,0],  IntRZ[:,1],  'b.-')

			for k in range(nInXPoints):
				Ax[ip].contour(Config.r2D, Config.z2D, Config.flux2D, colors="r", linestyles='-', linewidths=2, levels=[Config.X_points[InXPoints[k]].psi])

			Ax[ip].plot(SpPos[:,0],  SpPos[:,1],  'go')					#Plot stryke points positions

			r_min, r_max, z_min, z_max = get_min_max_walls(Config)
			Ax[ip].set_xlim(r_min, r_max)
			Ax[ip].set_ylim(z_min, z_max)

		elif(i == 1):														#connections lenghs
			Ax[ip].autoscale(enable=True, axis='both', tight=True)
			Ax[ip].set_xlabel(xLabels[0][0])
			Ax[ip].set_ylabel("$L_{par}\ (m)$")
			SpPos[:,2] = np.abs(SpPos[:,2])
			Ax[ip].plot(dist,  SpPos[:,2],  'b.-')
			Ax[ip].set_ylim(bottom=0.)

			if(save == "csv"):
				csv.append(types.SimpleNamespace())
				csv[-1].Name   = "L_par (m)"
				csv[-1].Values = np.copy(SpPos[:,2])

		elif(i == 2):														#flux expansion
			Ax[ip].autoscale(enable=True, axis='both', tight=True)
			Ax[ip].set_xlabel(xLabels[0][0])
			Ax[ip].set_ylabel("$\\Gamma_{exp}$")
			Ax[ip].plot(dist, SpPos[:,4],  'b.-')
			Ax[ip].set_ylim(bottom=0., top=SpPos[:,4].max()*1.1)

			if(save == "csv"):
				csv.append(types.SimpleNamespace())
				csv[-1].Name   = "Flux_exp"
				csv[-1].Values = np.copy(SpPos[:,4])
		elif(i == 3):														#strike points angles
			Ax[ip].autoscale(enable=True, axis='both', tight=True)
			Ax[ip].set_xlabel(xLabels[0][0])
			Ax[ip].set_ylabel("$\\alpha_i\ (^{\\circ})$")
			Ax[ip].plot(dist, SpPos[:,3]*180./np.pi,  'b.-')
			Ax[ip].set_ylim(bottom=0., top=SpPos[:,3].max()*180./np.pi*1.1)

			if(save == "csv"):
				csv.append(types.SimpleNamespace())
				csv[-1].Name   = "Alpha_i (deg)"
				csv[-1].Values = np.copy(SpPos[:,3])
		else:
			print("ERROR: extraplot non defined")



	for iPlasma in range(len(ions)):
		for i in range(len(VPosPlot[iPlasma])):
			ip = VPosPlot[iPlasma][i]
			Ax[ip].autoscale(enable=True, axis='both', tight=True)
			Ax[ip].set_xlabel(xLabels[iPlasma][i])
			if(log_scale == 0):
				Ax[ip].set_ylabel(yLabels[iPlasma][i])
				Ax[ip].set_yscale('linear')
			else:
				Ax[ip].set_ylabel(yLogLabels[iPlasma][i])
				Ax[ip].set_yscale(LogScales[iPlasma][i])

#	Plot parameters

	for iPlasma in range(len(ions)):
		for i in range(len(VPosPlot[iPlasma])):
			ip = VPosPlot[iPlasma][i]

			if(iPlasmaValues[iPlasma][i] != -1):	Title = ions[iPlasma]
			Ax[ip].set_title(Title)
			if(len(evolution) < 2):
				if(len(path) < 2):
					Ax[ip].plot(dist,  Values[0][iPlasma][i]*Facts[iPlasma][i],  'b.-')
					if(save == "csv"):
						csv.append(types.SimpleNamespace())
						csv[-1].Name   = FLabels[iPlasma][i].format(Title[:-1])
						csv[-1].Values = np.copy(Values[0][iPlasma][i]*Facts[iPlasma][i])
				else:
					for iPh in range(len(path)):
						if(no_labels == 0): Ax[ip].plot(dist, Values[iPh][iPlasma][i]*Facts[iPlasma][i],  '-', color = colors[iPh], label=path[iPh])
						else:					    Ax[ip].plot(dist, Values[iPh][iPlasma][i]*Facts[iPlasma][i],  '-', color = colors[iPh])
			else:
				for iEv in range(len(evolution)):
					Ax[ip].plot(dist, Values[iEv][iPlasma][i]*Facts[iPlasma][i],  '-', color = colors[iEv], label="{:d}".format(evolution[iEv]))

			if((log_scale == 0) and (BottomZero[iPlasma][i])): Ax[ip].set_ylim(bottom=0.)

		xv = 0.
		for i in range(len(VPosPlot[iPlasma])):
#			pyp.setp(Ax[k].yaxis.get_ticklabels(), rotation=90.)
			Ax[VPosPlot[iPlasma][i]].axvline(x=xv, color='k', linestyle='dashed')
			if(((len(path) > 1) or (len(evolution) > 1)) and (no_labels == 0)): Ax[VPosPlot[iPlasma][i]].legend(fontsize='small', loc='lower left')

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
				
			if(HasPathLabel): 
				np.savetxt("sp_on_mesh_cut_{:s}.{:s}".format(path_label[0], save), save_cvs, header=Header[:-1], delimiter=",", fmt="%15.7e", comments="")
			else:
				np.savetxt("plot1d_sp_on_mesh_cut_t={:.3f}.".format(RefPar.time)+save, save_cvs, header=Header[:-1], delimiter=",", fmt="%15.7e", comments="")


		else:
			i_plot_file = 0
			for i in range(len(Fig)):
				if(one_plot != 1): Fig[i].set_size_inches(10.05,7.44)
				if(save == "pdf"):
					pdf.savefig(Fig[i])
				else:
					i_plot_file += 1
					Fig[i].savefig("plot1d_sp_on_mesh_cut_t={:.3f}_{:d}.".format(RefPar.time,i_plot_file)+save)

		pyp.show(block=False)
		pyp.close()
	else:
		pyp.show()

	if(save == "pdf"):	pdf.close()

	print("plot1d_sp_on_mesh_cut: Completed")

	return

