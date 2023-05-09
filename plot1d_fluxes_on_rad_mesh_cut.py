#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot1d_fluxes_on_rad_mesh_cut 			import plot1d_fluxes_on_rad_mesh_cut
	from routines.cli_routines	import *

	if(cli_present("-h",sys.argv)):
		print("plot1d_fluxes_on_rad_mesh_cut options")
		print("\t-path        Directory with simulation or list of directories[d='./']")
		print("\t-rz0_line    Core coordinate direction line (Es.: [3.,0.]) [D=[]")
		print("\t-theta_line  Angular (deg) direction line [D=[0.]")
		print("\t-no_labels   Skip labels in plots [d=false]")
		print("\t-d_only      plot only e- and D+ [d=false]")
		print("\t-rho_scale   Use rho_pol scale for x axis [d=false]")
		print("\t-psi_scale   Use psi_pol norm. scale for x axis [d=false]")
		print("\t-log_scale   Use log scale for y axis [d=false]")
		print("\t-one_plot    One plot on each figure [d=false]")
		print("\t-save        Save figures on files, values=none/png/ps/eps/pdf [D='none']")
		print()
		exit()

	path	 	= cli_get_value("-path",				sys.argv,   [""])
	rz0_line	= cli_get_value("-rz0_line",			sys.argv,	[])
	theta_line	= cli_get_value("-theta_line",		sys.argv,	0.)
	save		= cli_get_value("-save",				sys.argv, "none")
	no_labels	= cli_present("-no_labels",			sys.argv)
	d_only		= cli_present("-d_only",				sys.argv)
	rho_scale	= cli_present("-rho_scale",			sys.argv)
	psi_scale	= cli_present("-psi_scale",			sys.argv)
	log_scale	= cli_present("-log_scale",			sys.argv)
	one_plot	= cli_present("-one_plot",			sys.argv)
	plot1d_fluxes_on_rad_mesh_cut(path=path, rz0_line=rz0_line, theta_line=theta_line, no_labels=no_labels, d_only=d_only, log_scale=log_scale, rho_scale=rho_scale, psi_scale=psi_scale,one_plot=one_plot, save=save)
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
from routines.globals					import DEBUG, KB, BALLOONING_NAMES

from mesh.get_rz_core_sep				import get_rz_core_sep
from mesh.get_rho_in_out_core_sep		import get_rho_in_out_core_sep
from mesh.find_zones_intersections		import find_zones_intersections
from mesh.compute_mesh_intersections	import compute_mesh_intersections

from files.load_soledge_mesh_file		import load_soledge_mesh_file
from files.load_fluxes_files			import load_fluxes_files
from files.load_ions_list				import load_ions_list
from files.load_refpar_file				import load_refpar_file
from files.load_input_file				import load_input_file
from files.load_transports_coefficients	import load_transports_coefficients

#==============================================================================
# This routine plots fluxes on mesh cut
#==============================================================================

def plot1d_fluxes_on_rad_mesh_cut(path=[], rz0_line = [2.,0.], theta_line=5., log_scale=0, rho_scale=0, psi_scale=0, no_labels=0, d_only=0, one_plot=0, save="none"):

	print("plot1d_fluxes_on_rad_mesh_cut")

	if(len(path) == 0):		 path = [""]

#	Read reference parameters

	
	path0  = path[0]
	if((len(path0) > 0) and (path0[-1] != "/")): path0 = path0 + "/"

	RefPar = load_refpar_file(path0+"Results/")

	ions = load_ions_list(path0)
	if(d_only != 0): ions = ions[0:2]

	input = load_input_file(path0)

#	Read mesh

	Config = load_soledge_mesh_file(path0+"mesh.h5")

#	Find mesh along line

	if(len(rz0_line) == 0):
		Rcore, Zcore, CoreMegazone = get_rz_core_sep(Config, core_and_sep = False)

		if(Zcore.min()*Zcore.max() < 0.):
			rz0_line = [0.5*(Rcore.min() + Rcore.max()), 0.]
		else:
			iMax = np.argmax(Rcore)
			rz0_line = [0.5*(Rcore.min() + Rcore.max()), Zcore[iMax]]

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
#	print("Lengths=",Lengths)

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
		xSep = 1.
	else:
		xName = "DIST"
		Ri, Zi, is1, is2  = intersect_2contours(RZsep[:,0], RZsep[:,1], IntRZ[:,0], IntRZ[:,1])
		if(len(Ri)==0):
			pyp.plot(RZsep[:,0], RZsep[:,1],'k-')
			pyp.plot(IntRZ[:,0], IntRZ[:,1],'r-')
			pyp.show()

		dsep = sqrt((IntRZ[0,0] - Ri[0])**2 + (IntRZ[0,1] - Zi[0])**2)
		dist -= dsep
		xSep = 0.


#	Read and plot parameters

	if(save == "pdf"):	pdf = PdfPages("plot1d_fluxes_on_rad_mesh_cut_t={:.3f}.".format(RefPar.time)+save)   #pdf in one file only

	FluxDir	= [["North", "South", "East", "West", "S-N+W-E"],["North", "South","S-N"]]
	if(rho_scale != 0):
		xLabels    = [["$\\rho_{pol}$",	 "$\\rho_{pol}$",  "$\\rho_{pol}$"],["$\\rho_{pol}$",	 "$\\rho_{pol}$",  "$\\rho_{pol}$","$\\rho_{pol}$",	 "$\\rho_{pol}$",  "$\\rho_{pol}$"]]
	elif(psi_scale != 0):
		xLabels    = [["$\\tilde{\\psi}_{pol}$",	 "$\\tilde{\\psi}_{pol}$",  "$\\tilde{\\psi}_{pol}$"],["$\\tilde{\\psi}_{pol}$",	 "$\\tilde{\\psi}_{pol}$",  "$\\tilde{\\psi}_{pol}$","$\\tilde{\\psi}_{pol}$",	 "$\\tilde{\\psi}_{pol}$",  "$\\tilde{\\psi}_{pol}$"]]
	else:
		xLabels    = [["$l\ (m)$",	 "$l\ (m)$",  "$l\ (m)$",  "$l\ (m)$"], ["$l\ (m)$",	 "$l\ (m)$",  "$l\ (m)$","$l\ (m)$",	 "$l\ (m)$",  "$l\ (m)$"]]

	yLabels 	= [["$CellFluxE\ (MW/m^2)$",	 "$FluxE\ (MW)$",	 "$CellFluxN\ (*10^{20}\ m^{-2})$", "$FluxN\ (*10^{20})$"], \
				   ["$CellFluxE\ (MW/m^2)$",	 "$FluxE\ (MW)$",	 "$CellFluxN\ (*10^{20}\ m^{-2})$", "$FluxN\ (*10^{20})$", "$CellFluxG$",	"$FluxG$"]]
	yLogLabels  = [["$Ln(CellFluxE)\ (MW/m^2)$", "$Ln(SurfFluxE)\ (MW)$", "$Ln(vFluxN)\ (*10^{20}\ m^{-2})", "$Ln(SurfFluxN\ (*10^{20})"], \
				   ["$Ln(CellFluxE)\ (MW/m^2)$", "$Ln(SurfFluxE)\ (MW)$", "$Ln(vFluxN)\ (*10^{20}\ m^{-2})", "$Ln(SurfFluxN\ (*10^{20})", "$Ln(CellFluxG)",  "$Ln(SurfFluxG)"]]

	LogScales	= [["log", "log", "log", "log"],["log", "log", "log", "log", "log", "log"]]

	iValues		= [[0,2],[0,2,1]]
	Facts		= [[1e-6, 1e-6, 1e-20, 1e-20],[1e-6, 1e-6, 1e-20, 1e-20, 1., 1.]]
	BottomZero	= [[False, False, False, False],[False, False, False, False, False, False]]
	VPosPlot	= [[1,4,2,5],[6,9,7,10,8,11]]											#Position of plasma plot
	EPosPlot    = [0,3]																	#Position of extra plots

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
			for i in range(len(VPosPlot[-1])): VPosPlot[-1][i] = VPosPlot[1][i] + (iPlasma - 1)*len(VPosPlot[1])

	colors = ['b','g','r','c','m','y','b','g','r','c','m','y']

	Transp = np.zeros((len(BALLOONING_NAMES), len(dist)), dtype='f8')
	for k in range(len(BALLOONING_NAMES)):
		if(k < len(BALLOONING_NAMES)-1):
			if((input.ballooning_parameters.ballooning_model == 2) and (Config.transp_values_OK)):
				Transp[k,:] = get_transp_parameter_on_mesh(Config, IntCEll, k)
			else:
				Transp[k,:] = 1.

			if(BALLOONING_NAMES[k] == "Chie"):
				Transp[k,:] *= input.transport_parameters.chie_p
			elif(BALLOONING_NAMES[k] == "D"):
				Transp[k,:] *= input.transport_parameters.Dn_p[0]
			elif(BALLOONING_NAMES[k] == "Nu"):
				Transp[k,:] *= input.transport_parameters.nu_p[0]
			elif(BALLOONING_NAMES[k] == "Chi"):
				Transp[k,:] *= input.transport_parameters.chii_p[0]
		else:
			if(BALLOONING_NAMES[k] == "Vpinch"):
				if((input.transport_parameters.pinch_model[0] == 3) and (Config.transp_values_OK)):
					Transp[k,:] = get_transp_parameter_on_mesh(Config, IntCEll, k)*input.transport_parameters.v_pinch[0]
				elif(input.transport_parameters.pinch_model[0] == 2):
					if(TranspCoeffRes == None): TranspCoeffRes = load_transports_coefficients(path0+"Results/", nZones=len(Config.Zones))
					Transp[iPh,k,:] = get_transp_parameter_on_mesh(TranspCoeff, IntCEll, k)
				elif(input.transport_parameters.pinch_model[0] == 1):
					Transp[k,:] = input.transport_parameters.v_pinch[0]
				elif(input.transport_parameters.pinch_model[0] == 0):
					Transp[iPh,k,:] = 0.
				else:
					print("\tUNKNOWN PINCH MODEL = ",input[-1].transport_parameters.pinch_model[0] )
					exit()

		TranspCoeffMesh = TranspCoeffRes =None
			
#	Load Flux parameters

	FluxValues = []
	for iPh in range(len(path)):
		path0  = path[iPh]
		if((len(path0) > 0) and (path0[-1] != "/")): path0 = path0 + "/"

		FluxValues.append([])
		for iPlasma in range(len(ions)):
			FluxValues[-1].append([])
			Fluxes = load_fluxes_files(path0+"Results/", nZones=len(Config.Zones), iFluxes = [iPlasma])

			for i in range(len(iValues[iPlasma])):
				if(iValues[iPlasma][i] > -1):
					FluxValues[-1][-1].append(get_flux_parameter_on_mesh(Fluxes[0], Config, iValues[iPlasma][i], IntCEll))
					FluxValues[-1][-1].append(get_flux_parameter_on_flux_surfaces(Fluxes[0], Config, iValues[iPlasma][i], IntCEll, RefPar))
				else:
					FluxValues[-1][-1].append([])
					FluxValues[-1][-1].append([])

#	Prepare for plotting and saving data

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
			for k in range(min(PlotPerFig,nPLots-i*PlotPerFig)):
				Ax.append(Fig[-1].add_subplot(nRows,nCols,k+1))
				Ax[-1].locator_params(axis='x',nbins=4)

			Fig[-1].tight_layout(pad=2., w_pad=3., h_pad=3.)
	else:
		for i in range(nPLots):
			Fig.append(pyp.figure())
			Ax.append(Fig[i].add_subplot(111))

	for iPlasma in range(len(ions)):
		for i in range(len(VPosPlot[iPlasma])):
			ip = VPosPlot[iPlasma][i]
			Ax[ip].set_title(ions[iPlasma])

			Ax[ip].autoscale(enable=True, axis='both', tight=True)
			Ax[ip].set_xlabel(xLabels[iPlasma][i])
			if(log_scale == 0):
				Ax[ip].set_ylabel(yLabels[iPlasma][i])
				Ax[ip].set_yscale('linear')
			else:
				Ax[ip].set_ylabel(yLogLabels[iPlasma][i])
				Ax[ip].set_yscale(LogScales[iPlasma][i])


#	Draw extra plots

	for i in range(len(EPosPlot)):
		ip = EPosPlot[i]
		Ax[ip].autoscale(enable=True, axis='both', tight=True)
		if(i == 0):															#plot section and cut
			if(len(path) == 1):
				Ax[ip].set_title(os.path.basename(os.path.abspath(path[0]))+" @ t={:.3f} s".format(RefPar.time))
			Ax[ip].set_aspect(1.)
			Ax[ip].set_xlabel("$R\ (m)$")
			Ax[ip].set_ylabel("$Z\ (m)$")

			plot2d_walls(Ax[ip], Config.Walls, plasma_wall=True, eirene_wall=True, extra_wall=True)

			Ax[ip].plot(RZcore[:,0],  RZcore[:,1],  'b-')
			Ax[ip].plot(RZsep[:,0],   RZsep[:,1],  'g-')

			Ax[ip].plot(IntRZ[Out_Sep,0], IntRZ[Out_Sep,1], 'g.-')
			Ax[ip].plot(IntRZ[In_Sep,0],  IntRZ[In_Sep,1],  'b.-')

		elif(i == 1):														#plot balloning
			Ax[ip].set_xlabel(xLabels[0][0])
			for k in range(len(BALLOONING_NAMES)):
				Ax[ip].plot(dist[In_Sep],  Transp[k,In_Sep],  'b.-', color = colors[k], label=BALLOONING_NAMES[k])
				Ax[ip].plot(dist[Out_Sep], Transp[k,Out_Sep], 'g.-', color = colors[k])
				Ax[ip].axvline(x=xSep, color='k', linestyle='dashed')

			Ax[ip].set_ylim(bottom=0.)
		else:
			print("ERROR: extraplot non defined")

	Ax[EPosPlot[1]].legend(fontsize='small', loc='upper left')

#	Plot parameters

	for iPh in range(len(path)):
		for iPlasma in range(len(ions)):
			iDir = 0
			for i in range(len(VPosPlot[iPlasma])):
				ip = VPosPlot[iPlasma][i]
				for iF in range(FluxValues[iPh][iPlasma][i].shape[1]):
					Ax[ip].plot(dist[In_Sep],  FluxValues[iPh][iPlasma][i][In_Sep,  iF]*Facts[iPlasma][i],  '.-', color = colors[iF], label=FluxDir[iDir][iF])
					Ax[ip].plot(dist[Out_Sep], FluxValues[iPh][iPlasma][i][Out_Sep, iF]*Facts[iPlasma][i], '.-', color = colors[iF])
				if(iDir == 0): iDir += 1
				else:		   iDir  = 0

		if((log_scale == 0) and BottomZero[iPlasma][i]): Ax[ip].set_ylim(bottom=0.)


		for i in range(len(VPosPlot[iPlasma])):
			ip = VPosPlot[iPlasma][i]
#			pyp.setp(Ax[k].yaxis.get_ticklabels(), rotation=90.)
			Ax[ip].axvline(x=xSep, color='k', linestyle='dashed')
			Ax[ip].legend(fontsize='small', loc='upper right')

	if(save != "none"):
		i_plot_file = 0
		for i in range(len(Fig)):
			i_plot_file += 1
			if(one_plot != 1): Fig[i].set_size_inches(10.05,7.44)
			if(save == "pdf"):
				pdf.savefig(Fig[i])
			else:
				Fig[i].savefig("plot1d_fluxes_on_rad_mesh_cut_t={:.3f}_{:d}.".format(RefPar.time,i_plot_file)+save)

		pyp.show(block=False)
		pyp.close()
	else:
		pyp.show()

	if(save == "pdf"):	pdf.close()

	print("plot1d_fluxes_on_rad_mesh_cut: Completed")

	return


#===================================================================

def get_transp_parameter_on_mesh(Config, IntCell, iBall):

	Zones = Config.Zones

	Transp = np.empty((IntCell.shape[0]), dtype='f8')

	iZones = np.unique(IntCell[:,0])
	for iZone in iZones:
		index = np.where(IntCell[:,0] == iZone)[0]
		Transp[index] = Zones[iZone].Ballooning[iBall][IntCell[index,1],IntCell[index,2]]

	return Transp


#===================================================================

def get_flux_parameter_on_mesh(Fluxes, Config, iPar, IntCell):

	Zones		= Config.Zones

	Par = np.empty((IntCell.shape[0], 5), dtype='f8')

	iZones = np.unique(IntCell[:,0])
	for iZone in iZones:
		index = np.where(IntCell[:,0] == iZone)[0]
		ii = IntCell[index,1]
		jj = IntCell[index,2]
		Par[index,:4] = Fluxes[iZone].Values[iPar][ii,jj,:4]

		ii1			  = ii + 1
		jj1			  = jj + 1
		SideLen		  = np.sqrt((Zones[iZone].gridR[ii1,jj1] - Zones[iZone].gridR[ii1,jj])**2 + (Zones[iZone].gridZ[ii1,jj1] - Zones[iZone].gridZ[ii1,jj])**2)			#North
		TotLen		  = np.copy(SideLen)
		Par[index,4]  = -Par[index,0]*SideLen
		SideLen		  = np.sqrt((Zones[iZone].gridR[ii,jj1]  - Zones[iZone].gridR[ii,jj])**2  + (Zones[iZone].gridZ[ii,jj1]  - Zones[iZone].gridZ[ii,jj])**2)			#South
		TotLen		 += SideLen
		Par[index,4] += Par[index,1]*SideLen
		SideLen		  = np.sqrt((Zones[iZone].gridR[ii1,jj1] - Zones[iZone].gridR[ii,jj1])**2 + (Zones[iZone].gridZ[ii1,jj1] - Zones[iZone].gridZ[ii,jj1])**2)			#East
		TotLen		 += SideLen
		Par[index,4] -= Par[index,2]*SideLen
		SideLen		  = np.sqrt((Zones[iZone].gridR[ii1,jj]  - Zones[iZone].gridR[ii,jj])**2  + (Zones[iZone].gridZ[ii1,jj]  - Zones[iZone].gridZ[ii,jj])**2)			#West
		TotLen		 += SideLen
		Par[index,4] += Par[index,3]*SideLen
		Par[index,4] /= TotLen

	return Par


#===================================================================

def get_flux_parameter_on_flux_surfaces(Fluxes, Config, iPar, IntCell, RefPar):

	Zones		= Config.Zones
	Megazones	= Config.Megazones

	Par = np.zeros((IntCell.shape[0], 3), dtype='f8')

	iZones = np.unique(IntCell[:,0])
	for iZone in iZones:

		index = np.where(IntCell[:,0] == iZone)[0]
		ii = IntCell[index,1]
		ii1 = ii + 1

		for kZone in Megazones[Zones[iZone].mz].list[:]:
			LenSouth2R	  = np.sqrt((Zones[kZone].gridR[ii ,1:] - Zones[kZone].gridR[ii, :-1])**2 + (Zones[kZone].gridZ[ii ,1:] - Zones[kZone].gridZ[ii, :-1])**2)* \
								    (Zones[kZone].gridR[ii ,1:] + Zones[kZone].gridR[ii, :-1])*(1.-Zones[kZone].Chi[ii,:])
			LenNorth2R	  = np.sqrt((Zones[kZone].gridR[ii1,1:] - Zones[kZone].gridR[ii1,:-1])**2 + (Zones[kZone].gridZ[ii1,1:] - Zones[kZone].gridZ[ii1,:-1])**2)* \
								    (Zones[kZone].gridR[ii1,1:] + Zones[kZone].gridR[ii1,:-1])*(1.-Zones[kZone].Chi[ii,:])

			Par[index,0] += np.sum(Fluxes[kZone].Values[iPar][ii,:,0]*LenNorth2R,axis=1)				#North
			Par[index,1] += np.sum(Fluxes[kZone].Values[iPar][ii,:,1]*LenSouth2R,axis=1)				#South

			InOut	 = Zones[kZone].Chi[ii,1:] - Zones[kZone].Chi[ii,:-1]
			InOutMax = np.max(InOut, axis=1)															#East jump up
			InOutMin = np.min(InOut, axis=1)															#West jump down
			
			if(InOutMax.max() > 0):
				iigt, jgt = np.where(InOut > 0)															#Index jumps plasma to wall at East
				igt		 = ii[iigt]
				igt1 	 = igt + 1
				jgt1	 = jgt + 1
				LenEast	 = np.sqrt((Zones[kZone].gridR[igt1,jgt1] - Zones[kZone].gridR[igt,jgt1])**2 + (Zones[kZone].gridZ[igt1,jgt1] - Zones[kZone].gridZ[igt,jgt1])**2)* \
							       (Zones[kZone].gridR[igt1,jgt1] + Zones[kZone].gridR[igt,jgt1])
				Par[index[iigt],2] -= Fluxes[kZone].Values[iPar][igt,jgt,2]*LenEast

			if(InOutMin.min() < 0):
				iilt, jlt = np.where(InOut < 0)															#Index jumps  dowm wall to plasma West
				ilt		  = ii[iilt]
				ilt1 	  = ilt + 1
				jlt		  = jlt + 1																			#move to first cell with chi = 0
				LenEast	  = np.sqrt((Zones[iZone].gridR[ilt1,jlt] - Zones[iZone].gridR[ilt,jlt])**2 + (Zones[iZone].gridZ[ilt1,jlt] - Zones[iZone].gridZ[ilt,jlt])**2)* \
								    (Zones[iZone].gridR[ilt1,jlt] + Zones[iZone].gridR[ilt,jlt])
				Par[index[iilt],2] += Fluxes[kZone].Values[iPar][ilt,jlt,3]*LenEast
		Par[index,2] += Par[index,1] - Par[index,0]


	Par *= np.pi		#2*pi*R (alredy multiply by 2*R)

	return Par
