#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot1d_wall_ntmpi 		import plot1d_wall_ntmpi
	from routines.cli_routines	import *

	if(cli_present("-h",sys.argv)):
		print("\nThis script plot plasma parameters on wall\n")
		print("plot1d_wall_ntmpi options")
		print("\t-path        Directory with simulation [d='./']")
		print("\t-evolution   If > 0 plot evolution file [d=0]")
		print("\t-rz0_line    Core coordinate direction line [D=[[0.,0.]]")
		print("\t-theta_line  Angular (deg) direction line [D=[0.]")
		print("\t-exp         Experiment to compare data [d='']")
		print("\t-shot        Shot to compare [d=0]")
		print("\t-time        Time to compare [d=0.]")
		print("\t-log_scale   Use log scale for y axis [d=false]")
		print("\t-one_plot    One plot on each figure[d=false]")
		print("\t-save        Save figures/data on files, values=none/png/ps/eps/pdf/csv [D='none']")
		print()
		exit()

	path	 	= cli_get_value("-path",			sys.argv,"./")
	evolution	= cli_get_value("-evolution",		sys.argv,	0)
	rz0_line	= cli_get_value("-rz0_line",		sys.argv,	[0.,0.])
	theta_line	= cli_get_value("-theta_line",	sys.argv,	0.)
	exp			= cli_get_value("-exp",			sys.argv,"")
	shot		= cli_get_value("-shot",			sys.argv,	0)
	time		= cli_get_value("-time",			sys.argv,	0.)
	save		= cli_get_value("-save",			sys.argv, "none")
	log_scale	= cli_present("-log_scale",		sys.argv)
	one_plot	= cli_present("-one_plot",		sys.argv)
 
 
	##########################
	# SETTINGS OVERRIDES... CLI CAPTURE BROKEN
	##########################
	# save = "csv"
	save = "none"

 
	plot1d_wall_ntmpi(path=path, evolution=evolution, rz0_line = rz0_line, theta_line=theta_line, exp=exp, shot=shot, time=time, log_scale=log_scale, one_plot=one_plot, save=save)
	exit()

#=======================================

# Function definition is here

import os
import h5py
import numpy									as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot						as pyp
from matplotlib.backends.backend_pdf			import PdfPages
from math										import sqrt
from eirene.get_wall_triangle					import get_wall_triangle
from routines.h5_routines						import h5_read
from routines.get_strike_points					import get_strike_points
from files.load_exp_data						import load_exp_data
from files.load_ions_list						import load_ions_list
from files.load_refpar_file						import load_refpar_file
from files.load_eirene_triangles				import load_eirene_triangles


#==============================================================================
# This routine plots ne/ni and Te/Ti ionization and gas pressure on eirene mesh
#==============================================================================

def plot1d_wall_ntmpi(path="./", evolution=0, rz0_line = [0.,0.], theta_line=0., exp="", shot=0, time=0, log_scale=0, one_plot=0, save="none"):

	print("plot1d_wall_ntmpi")

	Eirene = load_eirene_triangles(path+"triangles.h5")
	eirene_neutrals = h5py.File(os.path.join(path, "Results", "eirene_neutrals"), 'r')

	i_plot_file = 0

	exp_data_ok = False
	if((len(exp) > 0) and (shot > 0) and (time > 0)):
		exp_data = load_exp_data(exp, shot, time)
		if(exp_data.lp.ok): exp_data_ok = True

#	Find wall intersection along line

	ZeroTriangle, ZeroSide, Ri, Zi, RWallTriangles, ZWallTriangles, iWallTriangles, iWallSide, iWallKnots  = \
				get_wall_triangle(Eirene, rz0_line=rz0_line, theta_line=theta_line, no_plot=1, no_print=1, no_triangles=0)

	RWallTriangles = np.append(RWallTriangles, RWallTriangles[0])
	ZWallTriangles = np.append(ZWallTriangles, ZWallTriangles[0])
	iWallKnots	   = np.append(iWallKnots, iWallKnots[0])

	WalldL  	   = np.sqrt((RWallTriangles[1:]-RWallTriangles[:-1])**2 + (ZWallTriangles[1:]-ZWallTriangles[:-1])**2)
	DistKnots	   = np.cumsum(np.append(0., WalldL))
	DistTriangles  = 0.5*(DistKnots[:-1]+DistKnots[1:])
		

#	Read references data

	RefPar = load_refpar_file(path+"Results/")

	if(evolution == 0):
		base_plasma_name = path+"/Results/"
	else:
		base_plasma_name = path+"/Evolution/{:d}_".format(evolution)

	if_plasma	= h5py.File(base_plasma_name+"plasma_0", "r")
	Te			= h5_read(if_plasma,"triangles/temperature")*RefPar.T0eV
	if_plasma.close()

	if_plasma	= h5py.File(base_plasma_name+"plasma_1", "r")
	Ti			= h5_read(if_plasma,"triangles/temperature")*RefPar.T0eV
	if_plasma.close()
	print(["*"]*12)	
	print(iWallKnots)

	Te 	= Te[iWallKnots]
	Ti 	= Ti[iWallKnots]

	ions = load_ions_list(path)

#	Prepare for plotting

	if(save == "pdf"):	pdf = PdfPages("plot1d_wall_ntmpi_t={:.3f}.".format(RefPar.time)+save)   #pdf in one file only

#	Read parameters to plot

	for iPlasma in range(len(ions)):
		try:
			if_plasma	= h5py.File(base_plasma_name+"plasma_"+str(iPlasma), "r")
		except:
			break

		temperature	= h5_read(if_plasma,"triangles/temperature")*RefPar.T0eV
		density		= h5_read(if_plasma,"triangles/density")*RefPar.n0
		velocity	= h5_read(if_plasma,"triangles/velocity")*RefPar.c0

		temperature = temperature[iWallKnots]
		density 	= density[iWallKnots]
		velocity 	= velocity[iWallKnots]

		Jsat		= np.abs(1.6022e-19*velocity*density)*1e-3											#eletronic charge

		M			= velocity/np.sqrt((Te+Ti)/2)*(np.sqrt(RefPar.T0eV)/RefPar.c0)										#warning m_i = 2 for deuterium
		
		if(iPlasma > 0):
			try:
				Sn_tri	= h5_read(if_plasma,"triangles/Sn")*RefPar.n0/RefPar.tau0
				if(Sn_tri.max() - Sn_tri.min() > 0.):
					Sn_tri = set_min_positive(Sn_tri)
					Sn	   = Sn_tri[iWallTriangles]
					Sn	   = 0.5*np.append(np.append(Sn[0]+Sn[-1], Sn[:-1]+Sn[1:]), Sn[0]+Sn[-1])
#					print("Sn.min(), Sn.max()",Sn.min(), Sn.max())
				else: Sn = 0.
				Sn_tri = 0
#				Sn	= -h5_read(if_plasma,"triangles/Sn")*RefPar.n0/RefPar.tau0									#recombination
			except:
				print("plot1d_wall_ntmpi: Sn is not available")
				Sn  = 0.

			try:
				Nn = eirene_neutrals["atomic_species"]["dens_1"][:]
				Tn = eirene_neutrals["atomic_species"]["T_1"][:]
				Nm = eirene_neutrals["molecular_species"]["dens_1"][:]
				Tm = eirene_neutrals["molecular_species"]["T_1"][:]
				
				Nn = Nn[iWallTriangles]
				Nm = Nm[iWallTriangles]
				Tn = Tn[iWallTriangles]
				Tm = Tm[iWallTriangles]

				Tn = set_min_positive(Tn)
				Tm = set_min_positive(Tm)
				Nn	= 0.5*np.append(np.append(Nn[0]+Nn[-1], Nn[:-1]+Nn[1:]), Nn[0]+Nn[-1])*RefPar.n0
				Nm	= 0.5*np.append(np.append(Nm[0]+Nm[-1], Nm[:-1]+Nm[1:]), Nm[0]+Nm[-1])*RefPar.n0
				Tn	= 0.5*np.append(np.append(Tn[0]+Tn[-1], Tn[:-1]+Tn[1:]), Tn[0]+Tn[-1])
				Tm	= 0.5*np.append(np.append(Tm[0]+Tm[-1], Tm[:-1]+Tm[1:]), Tm[0]+Tm[-1])

				Pn	= (Nn*Tn+Nm*Tm)*1.6e-19
#				Pn	= 0.5*np.append(np.append(Pn[0]+Pn[-1], Pn[:-1]+Pn[1:]), Pn[0]+Pn[-1])

#				Nn	= 0.; Nm = 0.; Tn = 0.; Tm	= 0.
				if(Pn.max() - Pn.min() < 0.): Pn = 0.
			except:
				print("plot1d_wall_ntmpi: Nn is not available")
				Pn  = 0.
		else:
			Sn	= 0.; Pn	= 0.

		if_plasma.close()

#		Skip Phi
		phi = 0.
		"""
		try:
			if_ef	= h5py.File(str(path)+"/Results/electric_fields", "r")
			phi		= h5_read(if_ef,"/triangles/phi")
			if_ef.close()
			phi = phi[iWallKnots]
		except:
			phi		= 0.
		"""
		Fig = []
		Ax  = []

		yLabels     = ["$n\/(*10^{19}\/m^{-3})$",	  "$T\/(eV)$",	 "$J_{sat}\/(kAm^{-2})$",	"$Mach\/number$"]
		yLogLabels  = ["$Ln(n)\/(*10^{19}\/m^{-3})$",  "$Ln(T)\/(eV)$", "$Ln(J_{sat})\/(kAm^{-2})$", "$Mach\/number$"]

		if(isinstance(Sn, np.ndarray)):
			yLabels.append("$Ioniz.\/("+ions[iPlasma][:-2]+")$")
			yLogLabels.append("$Ln(Ioniz.)$")

		if(isinstance(Pn, np.ndarray)):
			yLabels.append("$n_A\ &\ n_M\ (*10^{19}\ m^{-3})$")
			yLogLabels.append("$Ln(n)$")
			yLabels.append("$T_A\ &\ T_M\ (K)$")
			yLogLabels.append("$Ln(T)$")
			yLabels.append("$P_A+P_m\ (Pa)$")
			yLogLabels.append("$Ln(P)$")
		if(isinstance(phi, np.ndarray)):
			yLabels.append("Phi")
			yLogLabels.append("Ln(Phi)")

		if(one_plot != 1):
			Fig.append(pyp.figure())
			for i in range(4): Ax.append(Fig[0].add_subplot(2,2,i+1))
			Fig[-1].tight_layout()

			if(len(yLabels) > 4):
				Fig.append(pyp.figure())
				for i in range(len(yLabels)-4): Ax.append(Fig[1].add_subplot(2,2,i+1))
				Fig[-1].tight_layout()
		else:
			for i in range(len(yLabels)):
				Fig.append(pyp.figure())
				Ax.append(Fig[i].add_subplot(111))

		for figure in Fig:  figure.patch.set_facecolor('white')
		
		Ax[0].set_title(os.path.basename(os.path.abspath(path))+" @ t={:.3f} s".format(RefPar.time))
		Ax[1].set_title(ions[iPlasma])

		for i in range(len(yLabels)):
			Ax[i].autoscale(enable=True, axis='both', tight=True)
			Ax[i].set_xlabel("$l\ (m)$")
			if(log_scale == 0):
				Ax[i].set_ylabel(yLabels[i])
				Ax[i].set_yscale('linear')
			else:
				Ax[i].set_ylabel(yLogLabels[i])
				Ax[i].set_yscale('log')

#			Ax[i].set_axis_bgcolor("#bdb76b")

	#	Plot parameters

		pl = []

		pl.append(Ax[len(pl)].plot(DistKnots,  density*1e-19,  'b-'))

		pl.append(Ax[len(pl)].plot(DistKnots,  temperature,  'b-'))	

		pl.append(Ax[len(pl)].plot(DistKnots,  Jsat,  'b-'))

		pl.append(Ax[len(pl)].plot(DistKnots, M,  'b-'))
		if(exp_data_ok and (iPlasma == 0)):
			Ax[0].plot(exp_data.lp.pos_l,  exp_data.lp.dens*1e-19, 'ro')
			Ax[1].plot(exp_data.lp.pos_l,  exp_data.lp.te,   'ro')
			Ax[2].plot(exp_data.lp.pos_l,  exp_data.lp.jsat*1e-3, 'ro')

		if(save == "csv"):
			if(iPlasma == 0): 
				file_header = "l (m)" 
				file_data	= [DistKnots]

			file_data.extend([density*1e-19, temperature, Jsat*1e-3, M])
			if(iPlasma == 0):
				file_header = file_header + ", ne (*10^19 m-3), Te (eV), Jsat_e (kA/m^2), Mach_e"
			else:
				file_header = file_header + ", ni (*10^19 m-3), Ti (eV), Jsat_i (kA/m^2), Mach_i"

		if(isinstance(Sn, np.ndarray)):
			Sn = np.where(Sn > 0.1,Sn, 0.1)
			pl.append(Ax[len(pl)].plot(DistKnots,  Sn,  'b-'))
			if(save == "csv"):
				file_data.append(Sn)
				file_header = file_header + ', Ioniz_'+ions[iPlasma][:-2]

		if(isinstance(Pn, np.ndarray)):
			Pn = np.where(Pn > 0.1,Pn, 0.1)
			pl.append(Ax[len(pl)].plot(DistKnots,  Nn*1e-19,  'b-'))
			Ax[len(pl)-1].plot(DistKnots,  Nm*1e-19,  'g-')
			pl.append(Ax[len(pl)].plot(DistKnots,  Tn*1.1604e4,  'b-'))
			Ax[len(pl)-1].plot(DistKnots,  Tm*1.1604e4,  'g-')
			pl.append(Ax[len(pl)].plot(DistKnots,  Pn,  'b-'))
			if(save == "csv"):
				file_data.extend([Nn, Tn, Pn])
				file_header = file_header + ', n_'+ions[iPlasma][:-2]+ ' (*10^19 m-3), "T_'+ions[iPlasma][:-2] + ' (eV), P_'+ions[iPlasma][:-2]+" (Pa)"

		if(isinstance(phi, np.ndarray)):
			pl.append(Ax[len(pl)].plot(DistKnots,  phi,  'b-'))
			if(save == "csv"):
				file_data.append(Phi)
				file_header = file_header + ', Phi_'+ions[iPlasma][:-2]

		if(save != "none"):
			if(save == "csv"):
				np.savetxt("wall_ntmpi_{:s}_t={:.3f}.csv".format(os.path.basename(os.path.abspath(path)),RefPar.time), np.array(file_data).T, delimiter=", ", fmt="%.4e", header= file_header, comments='')
				
			else:
				for i in range(len(Fig)):
					i_plot_file += 1
					if(one_plot != 1): Fig[i].set_size_inches(20.,15.)
					if(save == "pdf"):
						pdf.savefig(Fig[i])
					else:
						Fig[i].savefig("plot1d_wall_ntmpi_t={:.3f}_{:d}.".format(RefPar.time,i_plot_file)+save)

			pyp.show(block=False)
			pyp.close()
		else:
			pyp.show()

	if(save == "pdf"):	pdf.close()

	print("plot1d_wall_ntmpi: Completed")

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

