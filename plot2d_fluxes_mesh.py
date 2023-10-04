#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot2d_fluxes_mesh import plot2d_fluxes_mesh
	from routines.cli_routines import *

	if(cli_present("-h",sys.argv)):
		print("plot2d_fluxes_mesh options")
		print("-path        Directory with simulation [d='']")
		print("-save        Save figures on files, values=none/png/ps/eps/pdf [D='none']")
		print("-no_shade    no shade 2d plot [d=false]")
		print("\t-no_mask   no masked out of wall [d=false]")
#		print("\t-log_scale Use log scale for z colors [d=false]")
		print("\t-d_only    plot only e- and D+ [d=false]")
		print("\t-no_samexy No same xy scale [d=false]")
		print("-one_plot    All figure on one plot [d=false]")
		exit()

	path		= cli_get_value("-path",			sys.argv,"")
	save		= cli_get_value("-save",			sys.argv, "none")
	no_shade	= cli_present("-no_shade",		sys.argv)
	no_mask		= cli_present("-no_mask",			sys.argv)
	log_scale	= cli_present("-log_scale",		sys.argv)
	d_only		= cli_present("-d_only",			sys.argv)
	no_samexy	= cli_present("-no_samexy",		sys.argv)
	one_plot	= cli_present("-one_plot",		sys.argv)
	plot2d_fluxes_mesh(path=path, log_scale=log_scale, d_only=d_only, one_plot=one_plot, no_mask=no_mask, no_shade=no_shade, no_samexy=no_samexy, save=save)
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
from files.load_fluxes_files			import load_fluxes_files
from files.load_refpar_file				import load_refpar_file
from files.load_ions_list				import load_ions_list
from routines.utils_routines			import extend_mat1
from routines.utils_walls				import plot2d_walls
from routines.globals					import DEBUG, KB

#=========================================================
# This routine to plot M, ne/ni and Te/Ti results
#=========================================================

def plot2d_fluxes_mesh(path="", log_scale=0, d_only=0, one_plot=0, no_shade=0, no_mask=0, no_samexy=0, save="none"):
	no_shade =1
	no_mask =1
	print("plot2d_fluxes_mesh")

#	Read references data

	if((len(path) > 0) and (path[-1] != "/")): path = path + "/"

	shading = 'gouraud'
	if(no_shade == 1): shading = 'flat'

#	Load plasma results

#	Read mesh

	Config = load_soledge_mesh_file(os.path.join(path,"mesh.h5"))

	Zones	= Config.Zones
	nZones	= len(Zones)

	RefPar = load_refpar_file(path+"Results/")

	ions = load_ions_list(path)

	if(d_only != 0): ions = ions[0:2]

	Titles 	= [["$\Gamma_E\ North\ (MW/m^2)$",				"$\Gamma_E\ South\ (MW/m^2)$",	 			"$\Gamma_E\ East\ (MW/m^2)$",  			"$\Gamma_E\ West\ (MW/m^2)$", 
				"$\Gamma_E\ S-N\ (MW/m^2)$", 				"$\Gamma_E\ W-E\ (MW/m^2)$", 				"$\Gamma_E\ S-N+W-E\ (MW/m^2)$", \
				"$\Gamma_n\ North\ *10^{20}\ m^{-2})$", 	"$\Gamma_n\ South\ *10^{20}\ m^{-2})$",	"$\Gamma_n\ East\ *10^{20}\ m^{-2})$",	"$\Gamma_n\ West\ *10^{20}\ m^{-2})$", 
				"$\Gamma_n\ S-N\ (*10^{20}\ m^{-2})$",		"$\Gamma_n\ W-E\ (*10^{20}\ m^{-2})$",		"$\Gamma_n\ S-N+W-E\ (*10^{20}\ m^{-2})$"]]

	EFluxFact 	= 1e-6*KB*RefPar.n0*RefPar.c0*RefPar.T0*RefPar.rs0/(2*np.pi*RefPar.R0)
	nFluxFact 	= 1e-20*RefPar.n0*RefPar.c0*RefPar.rs0/(2*np.pi*RefPar.R0)
	FluxFacts	= [[EFluxFact, EFluxFact, EFluxFact, EFluxFact, EFluxFact, EFluxFact, EFluxFact, nFluxFact, nFluxFact, nFluxFact, nFluxFact, nFluxFact, nFluxFact, nFluxFact]]

	iValues		= [[0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2]]
	iPars		= [[0, 1, 2, 3,-1,-2,-3, 0, 1, 2, 3,-1,-2,-3]]
	PosPlots	= [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13]]
	if(len(ions) > 1):
		for iPlasma in range(1,len(ions)):
			Titles.append(Titles[0]) 
			iValues.append(iValues[0]) 
			iPars.append(iPars[0]) 
			FluxFacts.append(FluxFacts[0]) 
			PosPlots.append(PosPlots[0][:]) 											#[:] slice to force real copy
			for i in range(len(PosPlots[-1])): PosPlots[-1][i] += PosPlots[-2][-1] + 1

#	Preparing for plot

	if(save == "pdf"):	pdf = PdfPages("plot2d_fluxes_mesh_t={:.3f}.".format(RefPar.time)+save)   #pdf in one file only
	else:				i_plot_file = 0

#	Prepare for plotting

	nRows		= 2
	nCols		= 3
	PlotPerFig	= nRows*nCols
	nPLots		= len(Titles[0])*len(ions)
	nFigs		= int(nPLots/PlotPerFig)
	if(nFigs*PlotPerFig < nPLots): nFigs += 1

	Fig		= []
	Ax		= []
	AxFig	= []
	Im		= []
	if(one_plot != 1):
		for i in range(nFigs):	
			Fig.append(pyp.figure())
			for i in range(min(PlotPerFig,nPLots-i*PlotPerFig)):
				if((len(Ax) == 0) or (no_samexy != 0)):	Ax.append(Fig[-1].add_subplot(nRows,nCols,i+1))
				else:									Ax.append(Fig[-1].add_subplot(nRows,nCols,i+1, sharex = Ax[0], sharey = Ax[0]))
				AxFig.append(Fig[-1])
				Im.append(0)
			Fig[-1].tight_layout()
	else:
		for i in range(nPLots):
			Fig.append(pyp.figure())
			if((len(Ax) == 0) or (no_samexy != 0)):	Ax.append(Fig[i].add_subplot(111))
			else:									Ax.append(Fig[i].add_subplot(111, sharex = Ax[0], sharey = Ax[0]))
			Im.append(0)
			AxFig.append(Fig[-1])

	for iPlasma in range(len(ions)):
		for i in range(len(PosPlots[iPlasma])):
			ip = PosPlots[iPlasma][i]
			if(i == 0): Ax[ip].set_title(os.path.basename(os.path.abspath(path))+" @ t={:.3f} s".format(RefPar.time))
			Ax[ip].set_aspect(1.)
			Ax[ip].set_xlabel(ions[iPlasma])

#	PLot all	

	for iPlasma in range(len(ions)):
		Fluxes = load_fluxes_files(path+"Results/", nZones=nZones, iFluxes = [iPlasma])
		Values = []
		for i in range(len(PosPlots[iPlasma])): 
			Values.append([])
			for k in range(nZones):
				if(iValues[iPlasma][i] > -1):
					if(iPars[iPlasma][i] > -1):
						Values[-1].append(Fluxes[0][k].Values[iValues[iPlasma][i]][:,:,iPars[iPlasma][i]])
					elif(iPars[iPlasma][i] == -1):
						Values[-1].append(Fluxes[0][k].Values[iValues[iPlasma][i]][:,:,1]-Fluxes[0][k].Values[iValues[iPlasma][i]][:,:,0])	#S-N
					elif(iPars[iPlasma][i] == -2):
						Values[-1].append(Fluxes[0][k].Values[iValues[iPlasma][i]][:,:,3]-Fluxes[0][k].Values[iValues[iPlasma][i]][:,:,2])	#W-E
					elif(iPars[iPlasma][i] == -3):
						Values[-1].append(Fluxes[0][k].Values[iValues[iPlasma][i]][:,:,1]-Fluxes[0][k].Values[iValues[iPlasma][i]][:,:,0]+Fluxes[0][k].Values[iValues[iPlasma][i]][:,:,3]-Fluxes[0][k].Values[iValues[iPlasma][i]][:,:,2]) #S-N+W-E

		Vmin	= np.empty((nZones, len(PosPlots[iPlasma])), dtype='f8')
		Vmax	= np.empty((nZones, len(PosPlots[iPlasma])), dtype='f8')
		nn		= 0
		for k in range(nZones):
			ii,jj  = np.where(Zones[k].Chi != 1.)
			if(len(ii) > 0):
				for i in range(len(PosPlots[iPlasma])):
					Vmin[nn,i]	= np.min(Values[i][k][ii,jj])*FluxFacts[iPlasma][i]
					Vmax[nn,i]	= np.max(Values[i][k][ii,jj])*FluxFacts[iPlasma][i]
				nn	+= 1

		Vmin = np.min(Vmin[:nn,:], axis=0)
		Vmax = np.max(Vmax[:nn,:], axis=0)
				
		for i in range(len(PosPlots[iPlasma])): 
			ip = PosPlots[iPlasma][i]
			for k in range(nZones):
				if(no_shade == 0):
					if(no_mask == 0):	MaskedValues = np.ma.masked_where(extend_mat1(Zones[k].Chi) == 1., extend_mat1(Values[i][k]*FluxFacts[iPlasma][i]))
					else:				MaskedValues = extend_mat1(Values[i][k]*FluxFacts[iPlasma][i])
				else:
					if(no_mask == 0):	MaskedValues = np.ma.masked_where(Zones[k].Chi == 1., Values[i][k]*FluxFacts[iPlasma][i])
					else:				MaskedValues = Values[i][k]*FluxFacts[iPlasma][i]

				Im[ip] = Ax[ip].pcolormesh(Zones[k].gridR, Zones[k].gridZ, MaskedValues, shading=shading,  vmin = Vmin[i],  vmax = Vmax[i])

			plot2d_walls(Ax[ip], Config.Walls)

			cb = AxFig[ip].colorbar(Im[ip], ax=Ax[ip])
			cb.set_label(Titles[iPlasma][i])

	if(save != "none"):
		for i in range(len(Fig)):
			if(one_plot != 1): Fig[i].set_size_inches(20.,15.)
			if(save == "pdf"):
				pdf.savefig(Fig[i])
			else:
				i_plot_file += 1
				Fig[i].savefig("plot2d_fluxes_mesh_{:.3f}_{:d}.".format(RefPar.time,i_plot_file)+save)
		pyp.show(block=False)
		pyp.close()
	else:
		pyp.show()

	if(save == "pdf"):	pdf.close()

	print("plot2d_fluxes_mesh: Completed")

	return

