#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot1d_tri 			import plot1d_tri
	from routines.cli_routines	import *

	if(cli_present("-h",sys.argv)):
		print("\nThis script plot plasma parameters on eirene mesh along a line starting from plasma core\n")
		print("plot1d_tri options")
		print("\t-path        Directory with simulation [d='./']")
		print("\t-evolution   If > 0 plot evolution file (or evolutions files with list [3,4,5])  [d=0]")
		print("\t-rz0_line    Core coordinate direction line (Es.: [3.,0.]) [D=[]")
		print("\t-theta_line  Angular (deg) direction line [D=[0.]")
		print("\t-mod_file    Text (csv, or tsv) file with model data to compare [d='']")
		print("\t-exp         HDF5 file with experiment data to compare [d='']")
		print("\t-shot        Shot to compare [d=0]")
		print("\t-time        Time to compare [d=0.]")
		print("\t-plot_all    Plot Rad & M for impurities also [d=false]")
		print("\t-rho_scale   Use rho scale for x axis [d=false]")
		print("\t-log_scale   Use log scale for y axis [d=false]")
		print("\t-one_plot    One plot on each figure [d=false]")
		print("\t-save        Save figures on files, values=none/png/ps/eps/pdf [D='none']")
		print()
		exit()

	path	 	= cli_get_value("-path",				sys.argv,"./")
	evolution	= cli_get_value("-evolution",			sys.argv,	[])
	rz0_line	= cli_get_value("-rz0_line",			sys.argv,	[])
	theta_line	= cli_get_value("-theta_line",		sys.argv,	0.)
	mod_file	= cli_get_value("-mod_file",			sys.argv,"")
	exp			= cli_get_value("-exp",				sys.argv,"")
	shot		= cli_get_value("-shot",				sys.argv,	0)
	time		= cli_get_value("-time",				sys.argv,	0.)
	save		= cli_get_value("-save",				sys.argv, "none")
	plot_all	= cli_present("-plot_all",			sys.argv)
	rho_scale	= cli_present("-rho_scale",			sys.argv)
	log_scale	= cli_present("-log_scale",			sys.argv)
	one_plot	= cli_present("-one_plot",			sys.argv)
	plot1d_tri(path=path, evolution=evolution, rz0_line=rz0_line, theta_line=theta_line, mod_file=mod_file, exp=exp, shot=shot, time=time, log_scale=log_scale, rho_scale=rho_scale, plot_all=plot_all, one_plot=one_plot, save=save)
	exit()

#=======================================

# Function definition is here

import os
import h5py
from math								import sqrt
import numpy							as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot				as pyp
from matplotlib.backends.backend_pdf	import PdfPages

from eirene.get_triangles_on_line		import get_triangles_on_line
from files.load_neighbors				import load_neighbors
from routines.h5_routines				import h5_read
from routines.intersect_contour			import intersect_2contours
from routines.utils_walls				import get_dmax_points_walls
from mesh.get_rz_core_sep				import get_rz_core_sep
from mesh.get_rho_in_out_core_sep		import get_rho_in_out_core_sep
from files.load_soledge_mesh_file		import load_soledge_mesh_file
from files.load_exp_data				import load_exp_data
from files.load_ions_list				import load_ions_list
from files.load_text_data				import load_text_data
from files.load_refpar_file				import load_refpar_file



#==============================================================================
# This routine plots ne/ni and Te/Ti ionization and gas pressure on eirene mesh
#==============================================================================

def plot1d_tri(path="./", evolution=[], rz0_line = [2.,0.], theta_line=5., mod_file="", exp="", shot=0, time=0, log_scale=0, rho_scale=0, plot_all=0, one_plot=0, save="none"):


	matplotlib_ver = matplotlib.__version__
	
	print("plot1d_tri")

	i_plot_file = 0

#	Read reference parameters

	RefPar = load_refpar_file(path+"Results/")

	ions = load_ions_list(path)


#	Read mesh

	Config = load_soledge_mesh_file(path+"mesh.h5")

#	Find all triangles along line

	if(len(rz0_line) == 0):
		Rcore, Zcore, CoreMegazone = get_rz_core_sep(Config, core_and_sep = False)
		rz0_line = [0.5*(Rcore.min() + Rcore.max()), 0.]

	rMax		= 4*get_dmax_points_walls(Config, rz0_line[0], z0_line[1], plasma_wall=True, eirene_wall=False, extra_wall=False)
	theta_line	= theta_line*np.pi/180.
	RLine		= np.array([rz0_line[0], rz0_line[0]+rMax*np.cos(theta_line)])
	ZLine		= np.array([rz0_line[1], rz0_line[1]+rMax*np.sin(theta_line)])
	IntTriangles, IntSides, dist, IntR, IntZ, IntKnots, IntWKnots = get_triangles_on_line(path, RLine, ZLine)

	Rho, In_Sep, Out_Sep, RZcore, RZsep = get_rho_in_out_core_sep(Config, IntR, IntZ)

	if(len(In_Sep) > 0):
		if(Out_Sep[-1] < In_Sep[0]):
			Out_Sep = np.append(Out_Sep, In_Sep[0])
		else:
			Out_Sep = np.append(In_Sep[-1], Out_Sep)

	if(rho_scale != 0):
		dist = Rho
	else:
		Ri, Zi, is1, is2  = intersect_2contours(RZsep[:,0], RZsep[:,1], IntR, IntZ)
		dsep = sqrt((IntR[0] - Ri[0])**2 + (IntZ[0] - Zi[0])**2)
		dist = -(dist-dsep)

	exp_data_ok = False
	if((len(exp) > 0) and (shot > 0) and (time > 0)):
		exp_data = load_exp_data(exp, shot, time)
	
		if(exp_data.ts.ok and (rho_scale != 0)):
			ii = np.where((exp_data.ts.rho >= dist.min()) & (exp_data.ts.rho <= dist.max()) & \
						  ((exp_data.ts.te > 0.) | (exp_data.ts.ne > 0.))); ii = ii[0]

			if(len(ii) >= 0): exp_data_ok = True

		if(exp_data_ok):
			exp_dist = exp_data.ts.rho[ii]
			exp_te	 = exp_data.ts.te[ii]
			exp_dte	 = exp_data.ts.dte[ii]
			exp_ne	 = exp_data.ts.ne[ii]
			exp_dne	 = exp_data.ts.dne[ii]

	mod_data_ok = False
	if(len(mod_file) > 0):
		Headers, TextData = load_text_data(mod_file)

		iRho = indexes(Headers, "rho_psi")
		if(len(iRho) != 1): iRho = indexes(Headers, "rho")

		il = indexes(Headers, "x")
		if(len(il) != 1): il = indexes(Headers, "l")
		if(len(il) != 1): il = indexes(Headers, "d")

		iNe = indexes(Headers, "Ne")
		iNi = indexes(Headers, "Ni")
		if(len(iNi) == 0): iNi = iNe
		if(len(iNe) == 0): iNe = iNi
		iTe = indexes(Headers, "Te")
		iTi = indexes(Headers, "Ti")

		
		if(((len(il) == 1) or (len(iRho) == 1)) and (len(iNe) == 1) and (len(iTe) == 1) and (len(iTi) == 1)):
			if(len(iRho) == 1):
				InRange = np.where((TextData[:,iRho[0]]  >= Rho.min()) &  (TextData[:,iRho[0]]  <= Rho.max()))[0]
			else:
				InRange = np.where((TextData[:,il[0]]  >= dist.min()) &  (TextData[:,il[0]]  <= dist.max()))[0]

		if(len(InRange) > 0):
				mod_data_ok = True

				if((len(iRho) == 1) and (rho_scale != 0)):
					mod_dist	= TextData[InRange,iRho[0]]													#Rho to rho
				elif((len(il) == 1) and (rho_scale == 0)):
					mod_dist	= TextData[InRange,il[0]]													#l to l
				elif((len(iRho) == 1) and (rho_scale == 0)):
					if(Rho[-1] > Rho[0]):
						mod_dist	= np.interp(TextData[InRange,iRho[0]], Rho, dist)						#Rho to dist
					else:
						mod_dist	= np.interp(TextData[InRange,iRho[0]], Rho[::-1], dist[::-1])			#Rho to dist
				elif((len(il) == 1) and (rho_scale != 0)):
					mod_dist	= np.interp(TextData[InRange,il[0]], dist, Rho)								#dist to Rho

				mod_density     = [TextData[InRange,iNe[0]], TextData[InRange,iNi[0]]]
				mod_temperature = [TextData[InRange,iTe[0]], TextData[InRange,iTi[0]]]
		


#	Read electron end main specie Ion temperature

	if(len(evolution) == 0):
		base_plasma_name = [path+"Results/"]
	else:
		base_plasma_name = []
		for ev  in evolution: base_plasma_name.append(path+"Evolution/{:d}_".format(ev))

	if(save == "csv"):
		file_data = np.empty((len(dist),len(ions)*5+4), dtype='f4')

	if(save == "pdf"):	pdf = PdfPages("plot1d_tri_t={:.3f}.".format(RefPar.time)+save)   #pdf in one file only

	colors = ['b','g','r','c','m','y','b','g','r','c','m','y']
#	Read parameters to plot
	Te	= []
	Ti	= []
	for iPlasma in range(len(ions)):
		temperature = []
		density		= []
		M			= []
		for iEv in range(len(base_plasma_name)):
			if(iPlasma == 0):
				if_plasma	= h5py.File(base_plasma_name[iEv]+"plasma_0", "r")
				Te.append(h5_read(if_plasma,"triangles/temperature")*RefPar.T0eV*1e-3)
				if_plasma.close()

				if_plasma	= h5py.File(base_plasma_name[iEv]+"plasma_1", "r")
				Ti.append(h5_read(if_plasma,"triangles/temperature")*RefPar.T0eV*1e-3)
				try:
					mass	= h5_read(if_plasma, 'mass', keep_array=False)
				except:
					mass	= 2.
				if_plasma.close()

				Te[iEv] 	= Te[iEv][IntKnots[:,0]]*IntWKnots[:,0]	+ Te[iEv][IntKnots[:,1]]*IntWKnots[:,1]
				Ti[iEv] 	= Ti[iEv][IntKnots[:,0]]*IntWKnots[:,0]	+ Ti[iEv][IntKnots[:,1]]*IntWKnots[:,1]

			try:
				if_plasma	= h5py.File(base_plasma_name[iEv]+"plasma_"+str(iPlasma), "r")
			except:
				break

			temperature.append(h5_read(if_plasma,"triangles/temperature")*RefPar.T0eV*1e-3)
			density.append(h5_read(if_plasma,"triangles/density")*RefPar.n0*1e-19)
			velocity	= h5_read(if_plasma,"triangles/velocity")*RefPar.c0

			temperature[iEv] = temperature[iEv][IntKnots[:,0]]*IntWKnots[:,0]	+ temperature[iEv][IntKnots[:,1]]*IntWKnots[:,1]
			density[iEv] 	 = density[iEv][IntKnots[:,0]]*IntWKnots[:,0]		+ density[iEv][IntKnots[:,1]]*IntWKnots[:,1]
			velocity 		 = velocity[IntKnots[:,0]]*IntWKnots[:,0]	+ velocity[IntKnots[:,1]]*IntWKnots[:,1]

			M.append(velocity/np.sqrt(1.e3*(Te[iEv]+Ti[iEv])/mass)*(sqrt(RefPar.T0eV)/RefPar.c0))

			if(iPlasma > 0):
				try:
					Nrad_tri	= -h5_read(if_plasma,"triangles/NRad")/RefPar.tau0*1e-3
					if(Nrad_tri.max() - Nrad_tri.min() > 0.):
						Nrad_tri = set_min_positive(Nrad_tri)
						Nrad.append(Nrad_tri[IntTriangles])
#						print("Nrad.min(), Nrad.max()",Nrad.min(), Nrad.max())
					else: Nrad = []
					Nrad_tri = 0
				except:
#					print("plot1d_tri: Nrad is not available in "+base_plasma_name[iEv]+"plasma_"+str(iPlasma))
					Nrad = []

				try:
					Sn_tri	= h5_read(if_plasma,"triangles/Sn")*RefPar.n0/RefPar.tau0
					if(Sn_tri.max() - Sn_tri.min() > 0.):
						Sn_tri = set_min_positive(Sn_tri)
						Sn.append(Sn_tri[IntTriangles])
	#					print("Sn.min(), Sn.max()",Sn.min(), Sn.max())
					else: Sn = []
					Sn_tri = 0
				except:
#					print("plot1d_tri: Sn is not available in "+base_plasma_name[iEv]+"plasma_"+str(iPlasma))
					Sn = []

				try:
					Nn	= h5_read(if_plasma,"triangles/Nn")
					Nm	= h5_read(if_plasma,"triangles/Nm")
					Tn	= h5_read(if_plasma,"triangles/Tn")
					Tm	= h5_read(if_plasma,"triangles/Tm")
					
					Nn = Nn[IntTriangles]
					Nm = Nm[IntTriangles]
					Tn = Tn[IntTriangles]
					Tm = Tm[IntTriangles]

					Tn = set_min_positive(Tn)
					Tm = set_min_positive(Tm)

					Pn.append((Nn*Tn+Nm*Tm)*1.6e-19*RefPar.n0*1e-19)
					Nn	= 0.; Nm = 0.; Tn = 0.; Tm	= 0.
					if(Pn[-1].max() - Pn[-1].min() < 0.): Pn = []
				except:
#					print("plot1d_tri: Nn is not available in "+base_plasma_name[iEv]+"plasma_"+str(iPlasma))
					Pn = []
			else:
				Nrad = []; Sn= []; Pn = []
			if_plasma.close()

	#		skip phi

			phi = []
			"""
			try:
				if_ef	= h5py.File(str(path)+"/Results/electric_fields", "r")
				phi.append(h5_read(if_ef,"/triangles/phi"))
				if_ef.close()
				phi[iEv] = phi[iEv][IntKnots[:,0]]*IntWKnots[:,0]	+ phi[iEv][IntKnots[:,1]]*IntWKnots[:,1]
			except:
				phi		= []
			"""

#		Prepare for plotting
		Fig = []
		Ax  = []
		if(rho_scale == 0):
			xLabels    = ["$R\ (m)$", "$l\ (m)$",	 "$l\ (m)$",  "$l\ (m)$"]
			file_header = '"l (m)"'
		else:
			if(matplotlib_ver[:3] == "1.3"):
				xLabels    = ["$R\ (m)$", "$\\rho$",	 "$\\rho$",  "$\\rho$"]
			else:
				xLabels    = ["$R\ (m)$", "$rho$",	 "$rho$",  "$rho$"]

		yLabels     = ["$Z\ (m)$", "$n\/(*10^{19}\/m^{-3})$",	  "$T\/(keV)$",	 "$Mach\/number$"]
		yLogLabels  = ["$Z\ (m)$", "$Ln(n)\/(*10^{19}\/m^{-3})$",  "$Ln(T)\/(keV)$", "$Mach\/number$"]
		file_header = '"rho (m)", "n (*10^19 m-3)", "T (keV)", "Mach"'

		if(len(Nrad) > 0):
			xLabels.append(xLabels[1])
			yLabels.append("$Rad\/("+ions[iPlasma][0]+")\ (kW)$")
			yLogLabels.append("$Ln(Rad)$")
			file_header = file_header + ', "Rad_'+ions[iPlasma]+'"'

		if(len(Sn) > 0):
			xLabels.append(xLabels[1])
			yLabels.append("$Ionization\/("+ions[iPlasma][0]+")$")
			yLogLabels.append("$Ln(Ion\ )$")
			file_header = file_header + ', "Ion_'+ions[iPlasma]+'"'

		if(len(Pn) > 0):
			xLabels.append(xLabels[1])
			yLabels.append("$P\/("+ions[iPlasma][0]+")$")
			yLogLabels.append("$Ln(P)$")
			file_header = file_header + ', "P_'+ions[iPlasma]+'"'

		if(len(phi) > 0):
			xLabels.append(xLabels[1])
			yLabels.append("Phi")
			yLogLabels.append("Ln(Phi)")
			file_header = file_header + ', "Phi_'+ions[iPlasma]+'"'

		if(one_plot != 1):
			Fig.append(pyp.figure())
			for i in range(4): Ax.append(Fig[-1].add_subplot(2,2,i+1))
			Fig[-1].tight_layout(pad=2., w_pad=3., h_pad=3.)

			if(len(yLabels) > 4):
				Fig.append(pyp.figure())
				for i in range(len(yLabels)-4): Ax.append(Fig[-1].add_subplot(2,2,i+1))
				Fig[-1].tight_layout(pad=2., w_pad=3., h_pad=3.)
		else:
			for i in range(len(yLabels)):
				Fig.append(pyp.figure())
				Ax.append(Fig[i].add_subplot(111))
			
		Ax[0].set_title(os.path.basename(os.path.abspath(path))+" @ t={:.3f} s".format(RefPar.time))

		for i in range(len(xLabels)):
			Ax[i].autoscale(enable=True, axis='both', tight=True)
			Ax[i].set_xlabel(xLabels[i])
			if((log_scale == 0) or (i == 0) or (i == 3)):
				Ax[i].set_ylabel(yLabels[i])
				Ax[i].set_yscale('linear')
			else:
				Ax[i].set_ylabel(yLogLabels[i])
				Ax[i].set_yscale('log')
			if(i == 0): Ax[i].set_aspect(1.)

#			Ax[i].set_axis_bgcolor("#bdb76b")

#		Plot parameters

		pl = 0
		plot2d_walls(Ax[pl], Config.Walls, plasma_wall=True, eirene_wall=False, extra_wall=True):
		Ax[pl].plot(RZcore[:,0],  RZcore[:,1],  'b-')
		Ax[pl].plot(RZsep[:,0],   RZsep[:,1],  'g-')

		Ax[pl].plot(IntR[Out_Sep], IntZ[Out_Sep], 'g-')
		Ax[pl].plot(IntR[In_Sep],  IntZ[In_Sep],  'b-')

#		Density

		pl += 1
		if(len(evolution) < 2):
			Ax[pl].plot(dist[In_Sep],  density[0][In_Sep],  'b-')
			Ax[pl].plot(dist[Out_Sep], density[0][Out_Sep], 'g-')
		else:
			for iEv in range(len(evolution)):
				Ax[pl].plot(dist[In_Sep],  density[iEv][In_Sep],  '-', color = colors[iEv], label="{:d}".format(evolution[iEv]))
				Ax[pl].plot(dist[Out_Sep], density[iEv][Out_Sep], '-', color = colors[iEv])

		if(exp_data_ok):
			Ax[pl].errorbar(exp_dist, exp_ne*1e-19, yerr=exp_dne*1e-19, fmt = 'o', color='r')

		if(mod_data_ok):
			Ax[pl].plot(mod_dist, mod_density[iPlasma]*1e-19, 'ro')

		if(len(evolution) != 1):
			Ax[pl].set_title(ions[iPlasma])
		else:
			Ax[pl].set_title(ions[iPlasma]+" evol. {:d}".format(evolution[0]))

#		Temperature

		pl += 1
		if(len(evolution) < 2):
			Ax[pl].plot(dist[In_Sep],  temperature[0][In_Sep],  'b-')
			Ax[pl].plot(dist[Out_Sep], temperature[0][Out_Sep], 'g-')
		else:
			for iEv in range(len(evolution)):
				Ax[pl].plot(dist[In_Sep],  temperature[iEv][In_Sep],  '-', color = colors[iEv], linewidth = 2, label="{:d}".format(evolution[iEv]))
				Ax[pl].plot(dist[Out_Sep], temperature[iEv][Out_Sep], '-', color = colors[iEv])

		if(exp_data_ok):
			Ax[pl].errorbar(exp_dist, exp_te*1e-3, yerr=exp_dte*1e-3, fmt = 'o', color='r')
			Ax[pl].set_title(exp+"# {:d}@{:.3f}".format(exp_data.shot, exp_data.ts.time))
		if(mod_data_ok):
			Ax[pl].plot(mod_dist, mod_temperature[iPlasma]*1e-3, 'ro')

#		M
				
		pl += 1
		if(len(evolution) < 2):
			Ax[pl].plot(dist[In_Sep], M[0][In_Sep],  'b-')
			Ax[pl].plot(dist[Out_Sep], M[0][Out_Sep], 'g-')
		else:
			for iEv in range(len(evolution)):
				Ax[pl].plot(dist[In_Sep], M[iEv][In_Sep],   '-', color = colors[iEv], linewidth = 2, label="{:d}".format(evolution[iEv]))
				Ax[pl].plot(dist[Out_Sep], M[iEv][Out_Sep], '-', color = colors[iEv])


		if(save == "csv"):
			if(iPlasma == 0):
				file_data[:, 0] = dist
				n_file_data		= 1

			file_data[:, n_file_data] = density[0];			n_file_data	+= 1	
			file_data[:, n_file_data] = temperature[0];		n_file_data	+= 1
			file_data[:, n_file_data] = M[0];				n_file_data	+= 1

#		Rad

		if(len(Nrad) > 0):
#			Nrad = np.where(Nrad > 0.1,Nrad, 0.1)
			pl += 1
			if(len(evolution) < 2):
				Ax[pl].plot(dist[In_Sep],  Nrad[0][In_Sep]*1e-3,  'b-')
				Ax[pl].plot(dist[Out_Sep], Nrad[0][Out_Sep]*1e-3, 'g-')
			else:
				for iEv in range(len(evolution)):
					Ax[pl].plot(dist[In_Sep],  Nrad[iEv][In_Sep]*1e-3,  '-', color = colors[iEv], linewidth = 2, label="{:d}".format(evolution[iEv]))
					Ax[pl].plot(dist[Out_Sep], Nrad[iEv][Out_Sep]*1e-3, '-', color = colors[iEv])

			if(save == "csv"):
				file_data[:, n_file_data]  = Nrad[0];	n_file_data	+= 1

#		Ionization

		if(len(Sn) > 0):
#			Sn = np.where(Sn > 0.1,Sn, 0.1)
			pl += 1
			if(len(evolution) < 2):
				Ax[pl].plot(dist[In_Sep],  Sn[0][In_Sep]*1e-19,  'b-')
				Ax[pl].plot(dist[Out_Sep], Sn[0][Out_Sep]*1e-19, 'g-')
			else:
				for iEv in range(len(evolution)):
					Ax[pl].plot(dist[In_Sep],  Sn[iEv][In_Sep]*1e-19,  '-', color = colors[iEv], linewidth = 2, label="{:d}".format(evolution[iEv]))
					Ax[pl].plot(dist[Out_Sep], Sn[iEv][Out_Sep]*1e-19, '-', color = colors[iEv])

			if(save == "csv"):
				file_data[:, n_file_data]  = Sn[0];	n_file_data	+= 1

#		Pressure

		if(len(Pn) > 0):
#			Pn = np.where(Pn > 0.1,Pn, 0.1)
			pl += 1
			if(len(evolution) < 2):
				Ax[pl].plot(dist[In_Sep],  Pn[0][In_Sep],  'b-')
				Ax[pl].plot(dist[Out_Sep], Pn[0][Out_Sep], 'g-')
			else:
				for iEv in range(len(evolution)):
					Ax[pl].plot(dist[In_Sep],  Pn[iEv][In_Sep],  '-', color = colors[iEv], linewidth = 2, label="{:d}".format(evolution[iEv]))
					Ax[pl].plot(dist[Out_Sep], Pn[iEv][Out_Sep], '-', color = colors[iEv])

			if(save == "csv"):
				file_data[:, n_file_data] = Pn[0];		n_file_data	+= 1

#		Phi

		if(len(phi) > 0):
			pl += 1
			if(len(evolution) < 2):
				pl.append(Ax[pl].plot(dist[In_Sep],  phi[0][In_Sep],  'b-'))
				Ax[pl].plot(dist[Out_Sep], phi[0][Out_Sep], 'g-')
			else:
				for iEv in range(len(evolution)):
					Ax[pl].plot(dist[In_Sep],  phi[iEv][In_Sep],  '-', color = colors[iEv], linewidth = 2, label="{:d}".format(evolution[iEv]))
					Ax[pl].plot(dist[Out_Sep], phi[iEv][Out_Sep], '-', color = colors[iEv])

			if(save == "csv"):
				file_data[:, n_file_data]  = phi;	n_file_data	+= 1

		if(rho_scale != 0): xv = 1.
		else:				xv = 0.
		for i in range(1,len(Ax)):	Ax[i].axvline(x=xv, color='k', linestyle='dashed')

		if(len(evolution) > 1):
			for k in range(1, len(Ax)):
				Ax[k].legend(fontsize='small', loc='lower left')
#				pyp.setp(Ax[k].yaxis.get_ticklabels(), rotation=90.)

		if(save != "none"):
			if(save == "csv"):
				np.savetxt("line_data_"+os.path.basename(os.path.abspath(path))+".csv", file_data[:,:n_file_data], delimiter=", ", fmt="%.4e", header= file_header, comments='')
				
			else:
				for i in range(len(Fig)):
					i_plot_file += 1
					if(one_plot != 1): Fig[i].set_size_inches(10.05,7.44)
					if(save == "pdf"):
						pdf.savefig(Fig[i])
					else:
						Fig[i].savefig("plot1d_tri_t={:.3f}_{:d}.".format(RefPar.time,i_plot_file)+save)

			pyp.show(block=False)
			pyp.close()
		else:
			pyp.show()

	if(save == "pdf"):	pdf.close()

	print("plot1d_tri: Completed")

	return


#
# This routine set to the minimum positive value all negative values
#

def set_min_positive(Values):
	index_n = np.where(Values < 0.)
	if(len(index_n[0]) > 0):
		index_p = np.where(Values >= 0.)
		min_p   = np.min(Values[index_p])
		Values[index_n] = min_p

	return Values


def indexes(list_strs, str):
	return  [i for i, j in enumerate(list_strs) if j == str]