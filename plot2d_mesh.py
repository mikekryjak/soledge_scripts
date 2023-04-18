#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot2d_mesh import plot2d_mesh
	from routines.cli_routines import *

	if(cli_present("-h",sys.argv)):
		print("\nThis script plot plasma parameters on soledge mesh\n")
		print("plot2d_mesh options")
		print("\t-path        Directory with simulation [d='']")
		print("\t-evolution   If > 0 plot evolution file [d=0]")
		print("\t-save        Save figures on files, values=none/png/ps/eps/pdf [D='none']")
		print("\t-plot_all    Plot Rad & M for impurities also [d=false]")
		print("\t-extra_walls Flag to show extra walls [d=false]")
		print("\t-no_shade    no shade 2d plot [d=false]")
		print("\t-no_mask     no masked out of wall [d=false]")
		print("\t-no_samexy   No same xy scale [d=false]")
		print("\t-one_plot    All figure on one plot [d=false]")
		print()
		exit()

	path		= cli_get_value("-path",			sys.argv,"")
	evolution	= cli_get_value("-evolution",	sys.argv,	0)
	save		= cli_get_value("-save",			sys.argv, "none")
	plot_all	= cli_present("-plot_all",		sys.argv)
	extra_walls	= cli_present("-extra_walls",	sys.argv)
	no_shade	= cli_present("-no_shade",		sys.argv)
	no_mask		= cli_present("-no_mask",			sys.argv)
	no_samexy	= cli_present("-no_samexy",		sys.argv)
	one_plot	= cli_present("-one_plot",		sys.argv)
	plot2d_mesh(path=path, evolution=evolution, plot_all=plot_all, one_plot=one_plot, extra_walls=extra_walls, no_shade=no_shade, no_mask=no_mask, no_samexy=no_samexy, save=save)
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
from files.load_plasma_files			import load_plasma_files
from files.load_soledge_mesh_file		import load_soledge_mesh_file
from files.load_refpar_file				import load_refpar_file
from files.load_ions_list				import load_ions_list
from routines.h5_routines 				import h5_read
from routines.utils_routines			import extend_mat1
from routines.utils_walls				import plot2d_walls

#=========================================================
# This routine to plot M, ne/ni and Te/Ti results
#=========================================================

def plot2d_mesh(path="", evolution=0, plot_all=0, one_plot=0, extra_walls=0, no_shade=0, no_mask=0, no_samexy=0, save="none"):

	print("plot2d_mesh")

	shading = 'gouraud'
	if(no_shade == 1): shading = 'flat'

#	Read references data

	if((len(path) > 0) and (path[-1] != "/")): path = path + "/"
	RefPar	= load_refpar_file(path+"Results")
	ions	= load_ions_list(path)

	Config = load_soledge_mesh_file(path+"mesh.h5")			#	load mesh
	Zones	= Config.Zones
	nZones  = len(Zones)


	Plasmas = load_plasma_files(path, nZones=nZones, Evolution=evolution)

#	Preparing for plot

	if(save == "pdf"):	pdf = PdfPages("plot2d_mesh_t={:.3f}.".format(RefPar.time)+save)   #pdf in one file only
	else:				i_plot_file = 0

#	Electrons and Main specie plasmas
#     ####################################

	Labels  = [["$n_e\ (10^{19}\ m^{-3})$", "$T_e\ (keV)$", "$n_i\ (10^{19})\ m^{-3}$", "$T_i\ (keV)$", "$R_{rad}\ (kW)$", "$M$"], 
			   ["$S_{ne}\ (*10^{19}\ m^{-3}s^{-1})$","$S_{Ge}$","$S_{Ee}\ (kW/m^{3})$",
				"$S_{ni}\ (*10^{19}\ m^{-3}s^{-1})$","$S_{Gi}$","$S_{Ei}\ (kW/m^{3})$"]]
	Pars   = [["Dense", "Tempe","Densi", "Tempi","TotRad","Mi"],["Sne","SGe","SEe","Sni","SGi","SEi"]]

	iPLot  = [[1,2,4,5,3,6],[1,2,3,4,5,6]]
	NrowNcol = [[3,2],[3,2]]

	Facts  = [[1e-19,1.e-3,1e-19,1.e-3,1.e-3,1.,1.],[1e-19,1.,1.e3,1e-19,1.,1.e3]]

	Im  = []

	FigNum = []
	PosNum= []
	FigPl = []
	AxPl = []
	Im	  = []
	if(one_plot != 1):
		for iF in range(len(Labels)):
			AxPl.append([])
			FigPl.append([])
			PosNum.append([])
			Im.append([])
			Fig= pyp.figure()
			for i in range(len(Labels[iF])): 
				if((len(AxPl[0]) == 0) or (no_samexy != 0)): AxPl[-1].append(Fig.add_subplot(NrowNcol[iF][1],NrowNcol[iF][0],iPLot[iF][i]))
				else:										 AxPl[-1].append(Fig.add_subplot(NrowNcol[iF][1],NrowNcol[iF][0],iPLot[iF][i], sharex = AxPl[0][0], sharey = AxPl[0][0]))
				FigPl[-1].append(Fig)
				Im[-1].append(0)

			Fig.patch.set_facecolor('white')
			Fig.tight_layout(pad=2., w_pad=3., h_pad=3.)
	else:
		for iF in range(len(Labels)):
			AxPl.append([])
			FigNum.append([])
			FigPl.append([])
			PosNum.append([])
			Im.append([])
			for i in range(len(Labels[iF])):
				Fig=pyp.figure()
				Fig.patch.set_facecolor('white')
				if((len(AxPl[0]) == 0) or (no_samexy != 0)):	AxPl[-1].append(Fig.add_subplot(111))
				else:										AxPl[-1].append(Fig.add_subplot(111, sharex = AxPl[0][0], sharey = AxPl[0][0]))
				FigPl[-1].append(Fig)
				Im[-1].append(0)

	for iF in range(len(Labels)):
		for i in range(len(Labels[iF])):
			AxPl[iF][i].set_aspect(1.)
			AxPl[iF][i].autoscale(enable=True, axis='both', tight=True)
			AxPl[iF][i].set_xlabel("R (m)")
			AxPl[iF][i].set_ylabel("Z (m)")
			if(i == 1): AxPl[iF][i].set_title(os.path.basename(os.path.abspath(path))+" @ t={:.3f} s".format(RefPar.time))

#	Plot parameters

	VminZ	= np.empty((nZones), dtype='f8')
	VmaxZ	= np.empty((nZones), dtype='f8')

	nPlasmas = len(Plasmas)
	for iF in range(len(Labels)):
		for i in range(len(Labels[iF])):
			ToFoundPar = True
			for iPlasma in range(nPlasmas):
				if(ToFoundPar):
					try:
						iPar = Plasmas[iPlasma][0].VNames.index(Pars[iF][i])
					except:
#						print("Not found parameter: ",Pars[iF][i]," in Plasma=",iPlasma)
						continue

					ToFoundPar = False
					if(evolution == 0):	AxPl[iF][i].set_title(ions[iPlasma])
					else:				AxPl[iF][i].set_title(ions[iPlasma]+" evol. {:d}".format(evolution))

#					Find min & max on all zones
					nn = 0
					Values = []
					for iZone in range(nZones):
						if(Plasmas[iPlasma][iZone].Values[iPar].shape[0] == Zones[iZone].Chi.shape[0]):
							Values.append(Plasmas[iPlasma][iZone].Values[iPar]*Facts[iF][i])
						elif(Plasmas[iPlasma][iZone].Values[iPar].shape[0] == Zones[iZone].Chi.shape[0] + 2):
							Values.append(Plasmas[iPlasma][iZone].Values[iPar][1:-1,1:-1]*Facts[iF][i])
						else:
							print("iZone,Zones[iZone].gridR.shape,Values.shape=",iZone,Zones[iZone].Chi.shape,Plasmas[iPlasma][iZone].Values[iPar].shape)

						ii,jj  = np.where(Zones[iZone].Chi != 1.)
						if(len(ii) > 0):
							VminZ[nn]	= np.min(Values[-1][ii,jj])
							VmaxZ[nn]	= np.max(Values[-1][ii,jj])
							if(Pars[iF][i]=="Tempi"):
								iijj		= np.argmax(Values[-1][ii,jj])
							nn		   += 1
					Vmin = np.min(VminZ[:nn], axis=0)
					Vmax = np.max(VmaxZ[:nn], axis=0)

					for iZone in range(nZones):
						if(no_shade == 0):
							if(no_mask == 0):	MaskedValues = np.ma.masked_where(extend_mat1(Zones[iZone].Chi) == 1., extend_mat1(Values[iZone]))
							else:				MaskedValues = extend_mat1(Values[iZone])		
						else:
							if(no_mask == 0):	MaskedValues = np.ma.masked_where(Zones[iZone].Chi == 1., Values[iZone])
							else:				MaskedValues = Values[iZone]
						Im[iF][i] = AxPl[iF][i].pcolormesh(Zones[iZone].gridR, Zones[iZone].gridZ, MaskedValues, shading=shading,  vmin = Vmin,  vmax = Vmax)

					plot2d_walls(AxPl[iF][i], Config.Walls, extra_wall=extra_walls)

					cb = FigPl[iF][i].colorbar(Im[iF][i], ax=AxPl[iF][i])
					cb.set_label(Labels[iF][i])

	if(save != "none"):
		for iF in range(len(Labels)):
			if(one_plot != 1): FigPl[iF][0].set_size_inches(10.05,7.44)
			if(save == "pdf"):
				pdf.savefig(Fig[i])
			else:
				i_plot_file += 1
				Fig[i].savefig("plot2d_mesh_{:.3f}_{:d}.".format(RefPar.time,i_plot_file)+save)
		pyp.show(block=False)
		pyp.close()
	else:
		pyp.show()


	if(save == "pdf"):	pdf.close()

	print("plot2d_mesh: Completed")

	return

