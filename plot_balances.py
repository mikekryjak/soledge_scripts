#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot_balances import plot_balances
	from routines.cli_routines import *

	if(cli_present("-h",sys.argv)):
		print("\nThis script plot balances and residuals files\n")
		print("plot_balances options")
		print("\t-path        Directory with simulation [d='./']")
		print("\t-last        If <= 0 ignored, <= 10 plot last time interval, > 10 plot last points [D=0]")
		print("\t-y_ranges    If different from zero ranges of yscale +-yRange (Ex. [1e22,4.e7,0.1]) [D=[]]")
		print("\t-use_time    Plot versus time [d=false]")
		print("\t-plot_all    Plot all ions non only e & D and Atoms [d=false]")
		print("\t-log_scale   Use log scale for y axis [d=false]")
		print("\t-one_plot    One plot on each figure [d=false]")
		print("\t-no_samex    No same x scale [d=false]")
		print("\t-use_pts     Use points instead than lines [d=false]")
		print("\t-save        Save figures on files, values=none/png/ps/eps/pdf [D='none']")
		print()
		exit()

	path	 	= cli_get_value("-path",			sys.argv,"")
	last		= cli_get_value("-last",			sys.argv,	0.)
	y_ranges	= cli_get_value("-y_ranges",		sys.argv,	[])
	save		= cli_get_value("-save",			sys.argv, "none")
	use_time	= cli_present("-use_time",		sys.argv)
	plot_all	= cli_present("-plot_all",		sys.argv)
	log_scale	= cli_present("-log_scale",		sys.argv)
	one_plot	= cli_present("-one_plot",		sys.argv)
	no_samex	= cli_present("-no_samex",		sys.argv)
	use_pts    	= cli_present("-use_pts",		sys.argv)
	plot_balances(path=path, last=last, y_ranges=y_ranges, use_time=use_time, use_pts=use_pts ,log_scale=log_scale, plot_all=plot_all, one_plot=one_plot, no_samex=no_samex, save=save)
	exit()

#=======================================

# Function definition is here

import h5py
import os
import numpy								as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot					as pyp
from matplotlib.backends.backend_pdf		import PdfPages
from routines.h5_routines					import h5_read
from files.load_ions_list				import load_ions_list
import pandas as pd

#==============================================================================
# This routine plots ne/ni and Te/Ti ionization and gas pressure on eirene mesh
#==============================================================================

def plot_balances(path="", last=0., y_ranges=[], use_time=0, use_pts=0, log_scale=0, plot_all=0, one_plot=0, no_samex=0, no_plot = 0, save="none"):

	print("plot_balances")

	if((len(path) > 0) and (path[-1] != "/")): path = path + "/"

#	Read simulation time

	try:
		if_glob = h5py.File(path+"Results/globals", "r")
		tempus	= h5_read(if_glob, "tempus",keep_array=False)
		if_glob.close()
	except:
		tempus	= 0.

	try:
		residuals	 = np.loadtxt(path+"residuals", dtype='f8')
	except:
		print("\tError: Not found or unable to read "+ path+"residuals")
		exit()
	
	LenBalRes = residuals.shape[0]

	ions	= load_ions_list(path)										#read ions list

	nPlasmas = len(ions)
	Atoms  = [ions[1][:ions[1].find("+")-1]]
	lAtoms = [len(Atoms[-1])]
	for i in range(1,nPlasmas):
		if(Atoms[-1] != ions[i][:lAtoms[-1]]):
			Atoms.append(ions[i][:ions[i].find("+")-1])
			lAtoms.append(len(Atoms[-1]))
			
	nAtoms = len(Atoms)

	balances = []
	for i in range(nPlasmas):
		balances.append(np.loadtxt(path+"balances_{:d}".format(i), dtype='f8'))
		LenBalRes = min(balances[0].shape[0], LenBalRes)

	t_balances     = balances[0][:,-1]						#Time array

	if((last > 0.) and (last <= 10.)):
		i_after = np.where(t_balances >= last)[0]
		if(len(i_after) > 0):
			last = LenBalRes - np.min(i_after)
		else:
			last = LenBalRes
	elif(last > 10.):
		last = min(LenBalRes, int(last))
	else:
		last = LenBalRes

#	0 = Residual Ne
#	1 = Residual ND
#	2 = Residual Ge
#	3 = Residual GD
#	4 = Residual Te
#	5 = Residual TD

	R_iValues = [0,1,nPlasmas,nPlasmas+1,2*nPlasmas,2*nPlasmas+1]
	if(use_pts == 0):	R_lines = ["-b","--b","-r","-c","--r","-y","--y"]
	else:				R_lines = [".b","ob", ".r",".c","or", ".y","oy"]
	R_labels  = ["$Res_Ne*10^{-3}$", "$Res_ND*10^{-3}$", "$Res_Ge*10^{-3}$", "$Res_GD*10^{-3}$","$Res_{Te}*10^{-3}$", "$Res_{TD}*10^{-3}$"]
	R_facts   = [1e3,1e3,1e3,1e3,1e3,1e3]

#	3 = Energy input
#	4 = Energy output
#	5 = Energy source
#	6 = Energy radiation
#	7 = Energy variation 
#	8 = Energy stored

	E_iValues  = [3,4,5,6,7,8]

	if(use_pts == 0): E_lines	= ['-b',  '-g',   '-r',  '-c', '-m',  '-y', '-k',"--m", '--k']
	else:			  E_lines	= ['.b',  '.g',   '.r',  '.c', '.m',  '.y', 'k.',"om",  'ok']
	E_labels = ["$P_{in}$", "$P_{out}$", "$P_{sou}$", "$P_{rad}$", "$P_{var}$", "$E_{stored}$", "$P_{DIFF}$"]
	E_facts = [1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6,1e-6,1e-6]

#	0  = Particle input
#	1  = Particle output
#	2  = Particle source
#	9  = Particle variation
#	10 = Particle stored

	N_iValues  = [0,1,2,9,10]
	if(use_pts == 0):	N_lines  = ["-b","-g","-r","-m","-y"]
	else:				N_lines  = [".b",".g",".r",".m",".y"]
	N_labels = ["$\Gamma_{in}*10^{21}$", "$\Gamma_{out}*10^{24}$", "$\Gamma_{sou}*10^{24}$", "$N_{Var}*10^{22}$", "$N_{Stored}*10^{21}$"]
	N_facts = [1e-21,1e-24,1e-24,1e-22,1e-21]

	i_plot_file = 0

	if(save == "pdf"):	pdf = PdfPages("plot_balances_t={:.3f}.".format(tempus)+save)   #pdf in one file only

	nRows		= 4
	nCols		= 1
	PlotPerFig	= nRows*nCols

	nPLots		= 2*(nAtoms)
	if(plot_all != 0): nPLots += 2*(nPlasmas-2)

	nFigs		= int(nPLots/PlotPerFig)
	if(nFigs*PlotPerFig < nPLots): nFigs += 1

	if(use_time == 1):
		x_balances = t_balances[-last:]
		x_Labels   = "$t\ (s)$"
	else:
		x_balances = np.arange(last,dtype='f8')
		x_Labels    = "$n_{Iter}\ (*10^3)$"

	y_Labels	= ["$Residuals$",  "$(MW)&(MJ)$", "$(s^{-1})$", "$(MW)&(MJ)$", "$(MW)&(MJ)$"]

	for i in range(int(nPLots/2)-1):
		y_Labels.append("$(MW)&(MJ)$")
		y_Labels.append("$(s^{-1})$")

	if no_plot == 0:
		Fig = []
		Ax  = []
		if(one_plot != 1):
			Fig.append(pyp.figure())
			for i in range(3):
				if((i == 0) or (no_samex != 0)):	Ax.append(Fig[-1].add_subplot(3,1,i+1))
				else:								Ax.append(Fig[-1].add_subplot(3,1,i+1, sharex = Ax[0]))
			Fig[-1].tight_layout()
			for i in range(nFigs):	
				Fig.append(pyp.figure())
				for k in range(min(PlotPerFig,nPLots-i*PlotPerFig)):
					if(no_samex != 0):	Ax.append(Fig[-1].add_subplot(nRows,nCols,k+1))
					else:				Ax.append(Fig[-1].add_subplot(nRows,nCols,k+1, sharex = Ax[0]))
				Fig[-1].tight_layout()
		else:
			for i in range(nPLots+3):
				Fig.append(pyp.figure())
				if((i == 0) or (no_samex != 0)):	Ax.append(Fig[i].add_subplot(111))
				else:								Ax.append(Fig[i].add_subplot(111, sharex = Ax[0]))

		nFigs  += 1													#Add residual figure
		nPLots += 3													#Add first page plots

		for i in range(len(y_Labels)):
			Ax[i].set_xlabel(x_Labels)
			Ax[i].set_ylabel(y_Labels[i])
			if(log_scale == 0):
				Ax[i].set_yscale('linear')
			else:
				Ax[i].set_yscale('log')

			k = i
			if(i > 0):  k = int((i-1)/2) + 1
			if(i == 3): k = 1 
			if((len(y_ranges) > i) and (y_ranges[k] > 0.)):
				Ax[i].set_ylim((-y_ranges[k], y_ranges[k]))
				Ax[i].autoscale(enable=True,  axis='x', tight=True)
				Ax[i].autoscale(enable=False, axis='y', tight=True)
			else:
				Ax[i].autoscale(enable=True, axis='both', tight=True)
			Ax[i].locator_params(axis='y',nbins=4)

#	Prepare sum of energies

	E_summ = []
	for k in range(len(E_iValues)):
		E_summ.append(np.zeros(last, dtype='f8'))
	for i in range(nPlasmas):
		for k in range(len(E_iValues)):
			E_summ[k] += balances[i][-last:,E_iValues[k]]

	E_summ[3] -= balances[0][-last:,E_iValues[3]]							#remove electron rad which alread contains all rads

	E_summ.append(E_summ[0] - E_summ[1] + E_summ[2]) 	#Ediff
	if(use_time == 1):
		E_summ.append(np.append(0, 0.5*np.cumsum((E_summ[4][1:] + E_summ[4][:-1])*(x_balances[1:]-x_balances[:-1]))))	#summ of 	energy variation
		E_summ[-1] += E_summ[5][0]
		E_labels.append("$Es_{var}$")
		E_summ.append(np.append(0, 0.5*np.cumsum((E_summ[6][1:] + E_summ[6][:-1])*(x_balances[1:]-x_balances[:-1]))))	#summ of 	energy variation
		E_summ[-1] += E_summ[5][0]
		E_labels.append("$Es_{diff}$")

	if no_plot == 0:
		Ax[0].set_title(os.path.basename(os.path.abspath(path))+" @ t={:.3f} s".format(tempus))

		for i in range(len(R_iValues)):
			Ax[0].plot(x_balances, residuals[-last:,R_iValues[i]]*R_facts[i], R_lines[i], label=R_labels[i])

		all_ions = "e"
		for Atom in Atoms: all_ions = all_ions + " - " + Atom
		Ax[1].set_title(all_ions)
		for i in range(len(E_summ)):
			Ax[1].plot(x_balances, E_summ[i]*E_facts[i], E_lines[i], label=E_labels[i])

		Ax[2].set_title(ions[1])
		for i in range(len(N_iValues)):
			Ax[2].plot(x_balances, balances[1][-last:,N_iValues[i]]*N_facts[i], N_lines[i], label=N_labels[i])

		Ax[3].set_title(ions[0])
		Ax[4].set_title(ions[1])
		for i in range(len(E_iValues)):
			Ax[3].plot(x_balances,  balances[0][-last:,E_iValues[i]]*E_facts[i], E_lines[i], label=E_labels[i])
			Ax[4].plot(x_balances,  balances[1][-last:,E_iValues[i]]*E_facts[i], E_lines[i], label=E_labels[i])

		dfs = dict()
	
		ip = 4
		# Only active if multiple ion species
		if(nAtoms > 1):
			E_lines  = E_lines[:-1]
			E_labels = E_labels[:-1]
			iPlasma  = 2
			for iAtom in range(1,nAtoms):
				E_summ   = []
				for k in range(len(E_iValues)):
					E_summ.append(np.zeros(last, dtype='f8'))

				N_summ   = []
				for k in range(len(N_iValues)):
					N_summ.append(np.zeros(last, dtype='f8'))

				while((iPlasma < nPlasmas) and (Atoms[iAtom] == ions[iPlasma][:lAtoms[iAtom]])):
					for i in range(len(E_iValues)): E_summ[i] += balances[iPlasma][-last:,E_iValues[i]]
					for i in range(len(N_iValues)): N_summ[i] += balances[iPlasma][-last:,N_iValues[i]]
					iPlasma += 1

				ip += 1
				Ax[ip].set_title(Atoms[iAtom])
				for i in range(len(E_iValues)):
					Ax[ip].plot(x_balances, E_summ[i]*E_facts[i], E_lines[i], label=E_labels[i])

				ip += 1
				Ax[ip].set_title(Atoms[iAtom])
				for i in range(len(N_iValues)):
					Ax[ip].plot(x_balances, N_summ[i]*N_facts[i], N_lines[i], label=N_labels[i])

		if(plot_all != 0):
			# Iterate through all species
			print(ions)
			for iPlasma in range(2,nPlasmas):
				ip += 1
				Ax[ip].set_title(ions[iPlasma])
				df = pd.DataFrame()
				for i in range(len(E_iValues)):
					Ax[ip].plot(x_balances,  balances[iPlasma][-last:,E_iValues[i]]*E_facts[i], E_lines[i], label=E_labels[i])
					df[E_labels[i]] = balances[iPlasma][-last:,E_iValues[i]]

				ip += 1
				Ax[ip].set_title(ions[iPlasma])
				for i in range(len(N_iValues)):
					Ax[ip].plot(x_balances,  balances[iPlasma][-last:,N_iValues[i]]*N_facts[i], N_lines[i], label=N_labels[i])
				

		for k in range(len(Ax)):
			Ax[k].legend(fontsize='small', loc='upper left')
			pyp.setp(Ax[k].yaxis.get_ticklabels(), rotation=90.)

		for figure in Fig:  figure.patch.set_facecolor('white')

		if(save != "none"):
			for i in range(len(Fig)):
				i_plot_file += 1
				if(save == "pdf"):
					Fig[i].set_size_inches(11.7,8.26)
					pdf.savefig(Fig[i])
				else:
					if(one_plot != 1): Fig[i].set_size_inches(20.,15.)
					Fig[i].savefig("plot_balances_t={:.3f}_{:d}.".format(tempus,i_plot_file)+save)

			pyp.show(block=False)
			pyp.close()
		else:
			pyp.show()

	if(save == "pdf"):	pdf.close()
 
	# MK tabular output
 
	# Collect all plotted variables to Pandas dataframe (you can convert it to a dict if you don't like pandas)
	df = pd.DataFrame()
 
	# Residuals
	for i in range(len(R_iValues)):
		df[R_labels[i]] = residuals[-last:,R_iValues[i]]*R_facts[i]
	
	# Density
	for i in range(len(N_iValues)):
		df[N_labels[i]] = balances[1][-last:,N_iValues[i]] # Do not multiply by N_facts to maintain SI
		
	# Energy (total)
	for i in range(len(E_summ)):
		df[E_labels[i]+"_tot"] = E_summ[i] * E_facts[i]   # Keep factors to stay in [MW]
		
	# Energy (per species)
	for i in range(len(E_iValues)):
		df[E_labels[i]+"_e"] = balances[0][-last:,E_iValues[i]]*E_facts[i] # Keep factors to stay in [MW]
		df[E_labels[i]+"_i"] = balances[1][-last:,E_iValues[i]]*E_facts[i] # Keep factors to stay in [MW]

	# Now rename the labels to drop latex syntax and unit scaling
	# Manually specify label mapping for most variables
	label_map = {
		"$Res_Ne*10^{-3}$" : "Res_Ne",
		"$Res_ND*10^{-3}$" : "Res_ND",
		"$Res_Ge*10^{-3}$" : "Res_Ge",
		"$Res_GD*10^{-3}$" : "Res_GD",
		"$Res_{Te}*10^{-3}$" : "Res_Te",
		"$Res_{TD}*10^{-3}$" : "Res_Td",
		"$\Gamma_{in}*10^{21}$" : "Gamma_in",
		"$\Gamma_{out}*10^{24}$" : "Gamma_out",
		"$\Gamma_{sou}*10^{24}$" : "Gamma_sou",
		"$N_{Var}*10^{22}$" : "N_Var",
		"$N_{Stored}*10^{21}$" : "N_stored",
	}

	# Construct label mapping for the energy variables
	e_map = dict()
	for label in E_labels:
		new_label = label.split("$")[1].replace("{","").replace("}","")
		
		for suffix in ["_tot", "_i", "_e"]:
			e_map[label+suffix] = new_label+suffix
			
	# Combine label map
	label_map = {**label_map, **e_map}
	df = df.rename(label_map, axis = 1)

	print("plot_balances: Completed")

	return df

