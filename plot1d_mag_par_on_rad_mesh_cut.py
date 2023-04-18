#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot1d_mag_par_on_rad_mesh_cut 		import plot1d_mag_par_on_rad_mesh_cut
	from routines.cli_routines	import *

	if(cli_present("-h",sys.argv)):
		print("\nThis script plot plasma parameters on soledge mesh along a poloidal flux surface\n")
		print("plot1d_mag_par_on_rad_mesh_cut options")
		print("\t-path        Directory with simulation or list of directories[d='./']")
		print("\t-rz0_line    Core coordinate direction line (Es.: [3.,0.]) [D=[]")
		print("\t-theta_line  Angular (deg) direction line [D=[0.]")
		print("\t-east_west   Stryke points on east=-1, west=1 or shorted distance from cut [d=0]")
		print("\t-rho_scale   Use rho_pol scale for x axis [d=false]")
		print("\t-psi_scale   Use psi_pol norm. scale for x axis [d=false]")
		print("\t-l_pol       Use polidal length instead than parallel length [d=false]")
		print("\t-log_scale   Use log scale for y axis [d=false]")
		print("\t-one_plot    One plot on each figure [d=false]")
		print("\t-save        Save figures on files, values=none/png/ps/eps/pdf [D='none']")
		print()
		exit()

	path	 	= cli_get_value("-path",				sys.argv,   [""])
	rz0_line	= cli_get_value("-rz0_line",			sys.argv,	[])
	theta_line	= cli_get_value("-theta_line",		sys.argv,	0.)
	east_west	= cli_get_value("-east_west",		sys.argv, 0)
	rho_scale	= cli_present("-rho_scale",			sys.argv)
	psi_scale	= cli_present("-psi_scale",			sys.argv)
	l_pol		= cli_present("-l_pol",				sys.argv)
	log_scale	= cli_present("-log_scale",			sys.argv)
	one_plot	= cli_present("-one_plot",			sys.argv)
	save		= cli_get_value("-save",				sys.argv, "none")
	plot1d_mag_par_on_rad_mesh_cut(path=path, rz0_line=rz0_line, theta_line=theta_line, east_west=east_west, rho_scale=rho_scale, psi_scale=psi_scale, l_pol=l_pol, log_scale=log_scale, one_plot=one_plot, save=save)
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
from matplotlib.backends.backend_pdf	import PdfPages

from routines.h5_routines				import h5_read
from routines.intersect_contour			import intersect_2contours
from routines.utils_walls				import get_in_out_walls
from routines.set_profile_from_filedata	import set_profile_from_filedata
from routines.globals					import DEBUG, KB, BALLOONING_NAMES

from mesh.get_rz_core_sep				import get_rz_core_sep
from mesh.get_rho_in_out_core_sep		import get_rho_in_out_core_sep
from mesh.find_zones_intersections		import find_zones_intersections
from mesh.compute_mesh_intersections	import compute_mesh_intersections

from files.load_soledge_mesh_file		import load_soledge_mesh_file
from files.load_eirene_triangles		import load_eirene_triangles

#==============================================================================
# This routine plots ne/ni and Te/Ti soledge mesh inpolidal direction
#==============================================================================

def plot1d_mag_par_on_rad_mesh_cut(path=[], rz0_line = [2.,0.], theta_line=5.,  east_west=0, rho_scale=0, psi_scale=0, l_pol=0, log_scale=0, one_plot=0, save="none"):


	matplotlib_ver = matplotlib.__version__
	
	print("plot1d_mag_par_on_rad_mesh_cut")


	i_plot_file = 0

#	Read and plot parameters


#	Load Mesh  parameters

	if(len(path) == 0):	path = [""]

	Dists = []
	Values = []
	for iPh in range(len(path)):
		Values.append([])

		path0  = path[iPh]
		if((len(path0) > 0) and (path0[-1] != "/")): path0 = path0 + "/"

#		Read mesh

		Config = load_soledge_mesh_file(path0+"mesh.h5")
		Zones	= Config.Zones

		Eirene = load_eirene_triangles(path0+"triangles.h5")
		WallTriangles = Eirene.WallTriangles

#		Read Metric

		if(l_pol == 0):
			if_metric = h5py.File(path0+"Results/metric", "r")
			Gmet = []
			for k in range(len(Zones)):
				zone = "zone{:d}".format(k+1)
				Gmet.append(h5_read(if_metric, zone+ '/G', order = 'F'))

			if_metric.close()
		
#		Find mesh along line

		if(len(rz0_line) == 0):
			Rcore, Zcore, CoreMegazone = get_rz_core_sep(Config, core_and_sep = False)
			rz0_line = [0.5*(Rcore.min() + Rcore.max()), Zcore[np.argmax(Rcore)]]

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
		if(rho_scale != 0): xName = "Rho_pol"
		else:  				xName = "Psi_pol"
		Rho, In_Sep, Out_Sep, RZcore, RZsep = get_rho_in_out_core_sep(Config, IntRZ[:,0], IntRZ[:,1], rho_type = xName)
		
		if(len(In_Sep) > 0):
			if(Out_Sep[-1] < In_Sep[0]):
				Out_Sep = np.append(Out_Sep, In_Sep[0])
			else:
				Out_Sep = np.append(In_Sep[-1], Out_Sep)

		if((rho_scale != 0) or (psi_scale != 0)):
			dist = Rho
			dsep = 0.
		else:
			xName = "Dist (mm)"
			Ri, Zi, is1, is2  = intersect_2contours(RZsep[:,0], RZsep[:,1], IntRZ[:,0], IntRZ[:,1])
			if(len(Ri)==0):
				pyp.plot(RZsep[:,0], RZsep[:,1],'k-')
				pyp.plot(IntRZ[:,0], IntRZ[:,1],'r-')
				pyp.show()

			dsep = sqrt((IntRZ[0,0] - Ri[0])**2 + (IntRZ[0,1] - Zi[0])**2)
			dist -= dsep

		OutCels		   = np.where(dist >= 0.)[0]

		Dists.append(dist[OutCels])

		LConLength	   = np.zeros((len(OutCels)), dtype='f8')
		TargetAngle	   = np.zeros((len(OutCels)), dtype='f8')
		FLExpans	   = np.zeros((len(OutCels)), dtype='f8')
		FPerpExpans	   = np.zeros((len(OutCels)), dtype='f8')

		DLStart		   = np.empty((len(OutCels)), dtype='f8')
		DLEnd		   = np.empty((len(OutCels)), dtype='f8')
		DPerpStart	   = np.empty((len(OutCels)), dtype='f8')
		DPerpEnd	   = np.empty((len(OutCels)), dtype='f8')

		NOutCels	   = 0
		for iCell in range(OutCels.shape[0]):

#			Select cell

			iZone 	= IntCEll[OutCels[iCell],0]
			ix		= IntCEll[OutCels[iCell],1]
			iTheta	= IntCEll[OutCels[iCell],2]

#			Find zones along poloidal coordinate

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

			DR1				   = Zones[iZone].gridR[ix+1, iTheta+1] - Zones[iZone].gridR[ix+1, iTheta]
			DZ1				   = Zones[iZone].gridZ[ix+1, iTheta+1] - Zones[iZone].gridZ[ix+1, iTheta]
			DR2				   = Zones[iZone].gridR[ix+1, iTheta+1] - Zones[iZone].gridR[ix, iTheta+1]
			DZ2			   	   = Zones[iZone].gridZ[ix+1, iTheta+1] - Zones[iZone].gridZ[ix, iTheta+1]
			DLStart[iCell]	   = sqrt(DR2**2+DZ2**2)
			DPerpStart[iCell]  = abs(DR1*DZ2 - DZ1*DR2)/sqrt(DR1**2+DZ1**2)

#			Look East

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

#			Look West

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

			if(east_west == 0):
				if(Lpara[-1] >= -Lpara[0]): east_west = -1
				else:						east_west = 1

			if(east_west == -1):
				iZoneEnd		  = iZones[-1]
				iThetaEnd		  = iThEast[-1]-1
				LConLength[iCell] = Lpara[-2]
			else:
				iZoneEnd		  = iZones[0]
				iThetaEnd		  = iThWest[0]
				LConLength[iCell] = -Lpara[1]
			
			DR1				= Zones[iZoneEnd].gridR[ix+1, iThetaEnd+1] - Zones[iZoneEnd].gridR[ix+1, iThetaEnd]
			DZ1				= Zones[iZoneEnd].gridZ[ix+1, iThetaEnd+1] - Zones[iZoneEnd].gridZ[ix+1, iThetaEnd]
			DR2				= Zones[iZoneEnd].gridR[ix+1, iThetaEnd+1] - Zones[iZoneEnd].gridR[ix, iThetaEnd+1]
			DZ2			   	= Zones[iZoneEnd].gridZ[ix+1, iThetaEnd+1] - Zones[iZoneEnd].gridZ[ix, iThetaEnd+1]
			DLEnd[iCell] = sqrt(DR2**2+DZ2**2)
			DPerpEnd[iCell] = abs(DR1*DZ2 - DZ1*DR2)/sqrt(DR1**2+DZ1**2)

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

			if(len(aa) == 0):
				print("Attention no wall triangle in soledge quadrangle: k,i,j=",iZoneEnd,ix,iThetaEnd)
				aa = np.where((WallTriangles.k == iZoneEnd) & (WallTriangles.i == ix))[0]
				print("Found in soledge quadrangle with j=",WallTriangles.j[aa+[0]])
				break
		
			if(len(aa) > 1):
				print("Attention more than one wall triangle in soledge quadrangle: triangles=",aa+1)

			NOutCels += 1

	#		Compute angle of incidence on target

			iSeq = np.where(Eirene.Wall.TriSequences[Eirene.Wall.iTriSeqExtPlasma] == aa[0])[0][0]

			Dr   = Eirene.Wall.R12[0][iSeq,1] - Eirene.Wall.R12[0][iSeq,0]
			Dz   = Eirene.Wall.Z12[0][iSeq,1] - Eirene.Wall.Z12[0][iSeq,0]
			Br   = Zones[iZoneEnd].Br[ix, iThetaEnd]
			Bz   = Zones[iZoneEnd].Bz[ix, iThetaEnd]
			Bphi = Zones[iZoneEnd].Bphi[ix, iThetaEnd]
			
			TargetAngle[iCell] = abs(asin((-Dz*Br + Dr*Bz)/(sqrt(Dr**2+Dz**2)*sqrt(Br**2+Bz**2+Bphi**2))))


		FLExpans[:NOutCels] 	= DLEnd[:NOutCels]/DLStart[:NOutCels]
		FPerpExpans [:NOutCels] = DPerpEnd[:NOutCels]/DPerpStart[:NOutCels]

		Values[-1].append(TargetAngle*180./np.pi)
		Values[-1].append(np.abs(LConLength))
		Values[-1].append(FLExpans)
		Values[-1].append(FPerpExpans)

#	Prepare for plotting
	
	if(save == "pdf"):	pdf = PdfPages("plot1d_mesh_t={:.3f}.".format(RefPar.time)+save)   #pdf in one file only

	if(save == "csv"):
		if(len(path) > 1):
			save = "none"
			print("plot1d_sp_on_mesh_cut: cannot save to csv more than one profile!!!")
		else:
			csv = []
			csv.append(types.SimpleNamespace())

			csv[-1].Name   = xName
			csv[-1].Values = Dists[0]


	if(rho_scale != 0):
		xLabels    = ["$\\rho_{pol}$", "$\\rho_{pol}$", "$\\rho_{pol}$", "$\\rho_{pol}$"]
	elif(psi_scale != 0):
		xLabels    = ["$\\tilde{\\psi}_{pol}$", "$\\tilde{\\psi}_{pol}$", "$\\tilde{\\psi}_{pol}$", "$\\tilde{\\psi}_{pol}$"]
	else:
		xLabels    = ["$l\ (m)$", "$l\ (m)$", "$l\ (m)$", "$l\ (m)$"]

	if(l_pol == 0):
		yLabels		= [ "$\\alpha_i\ (^{\\circ})$", "$L_{par}\ (m)$",	 "$DlExp$",	 "$DperpExp$"]
	else:
		yLabels		= [ "$\\alpha_i\ (^{\\circ})$", "$L_{pol}\ (m)$",	 "$DlExp$",	 "$DperpExp$"]

	LogScales		= ["log", "log", "log", "log"]

	FLabels    		= ["Alpha_i (deg)", "ConLength (m)",	"SideExp", "FluxExp"]


	BottomZero	= [True, True, True,True]
	VPosPlot	= [0,1,2,3]

	colors = ['b','g','r','c','m','y','b','g','r','c','m','y']

#	Prepare for plotting

	nRows		= 2
	nCols		= 2
	PlotPerFig	= nRows*nCols
	nPLots		= len(VPosPlot)
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

	for i in range(len(VPosPlot)):
		ip = VPosPlot[i]
		if((VPosPlot[i] == 0) and (len(path) == 1)):
			Ax[ip].set_title(os.path.basename(os.path.abspath(path[0])))

		Ax[ip].autoscale(enable=True, axis='both', tight=True)
		Ax[ip].set_xlabel(xLabels[i])
		if(log_scale == 0):
			Ax[ip].set_ylabel(yLabels[i])
			Ax[ip].set_yscale('linear')
		else:
			Ax[ip].set_ylabel(yLogLabels[i])
			Ax[ip].set_yscale(LogScales[i])


#	Plot parameters

	for i in range(len(VPosPlot)):
		if(len(path) < 2):
			Ax[VPosPlot[i]].plot(Dists[0],  Values[0][i],  'b.-')
			if(save == "csv"):
				csv.append(types.SimpleNamespace())
				csv[-1].Name   = FLabels[i]
				csv[-1].Values = np.copy(Values[0][i])

		else:
			for iPh in range(len(path)):
				Ax[VPosPlot[i]].plot(Dists[iPh],  Values[iPh][i],  '-', color = colors[iPh], label=path[iPh])

		if((log_scale == 0) and BottomZero[i]): Ax[VPosPlot[i]].set_ylim(bottom=0.)

	xv = 0.
	for i in range(len(VPosPlot)):
#		pyp.setp(Ax[k].yaxis.get_ticklabels(), rotation=90.)
		Ax[VPosPlot[i]].axvline(x=xv, color='k', linestyle='dashed')
		if(len(path) > 1): Ax[VPosPlot[i]].legend(fontsize='small', loc='lower left')

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
				
			if(east_west == -1):
				np.savetxt("plot1d_mag_par_on_rad_mesh_cut_east.txt", save_cvs, header=Header[:-1], delimiter=",", fmt="%15.7e", comments="")
			else:
				np.savetxt("plot1d_mag_par_on_rad_mesh_cut_west.txt", save_cvs, header=Header[:-1], delimiter=",", fmt="%15.7e", comments="")
		else:
			for i in range(len(Fig)):
				i_plot_file += 1
				if(one_plot != 1): Fig[i].set_size_inches(10.05,7.44)
				if(save == "pdf"):
					pdf.savefig(Fig[i])
				else:
					Fig[i].savefig("plot1d_mag_par_on_rad_mesh_cut_t={:.3f}_{:d}.".format(RefPar.time,i_plot_file)+save)

		pyp.show(block=False)
		pyp.close()
	else:
		pyp.show()

	if(save == "pdf"):	pdf.close()

	print("plot1d_mag_par_on_rad_mesh_cut: Completed")

	return
