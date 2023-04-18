#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot_target 					import plot_target
	from routines.cli_routines		import *

	if(cli_present("-h",sys.argv) and cli_present("-exp_file",sys.argv) ):
		from routines.exp_data_routines			import get_exp_data_point_help
		from files.load_exp_descr						import load_exp_descr_help
		load_exp_descr_help()
		get_exp_data_point_help()
		exit()

	if(cli_present("-h",sys.argv)):
		print("\nThis script plot plasma parameters on soledge mesh along a line starting from plasma core\n")
		print("plot_target options")
		print("\t-path           [BROKEN] Directory with simulation [d='./']")
		print("\t-r_bias         Bias of R coordinate[m] [d='0.00']")
		print()
		exit()

	path = cli_get_value("-path",				sys.argv,   [""])
	r_bias = cli_get_value("-r_bias",		sys.argv,	0.00)
	plot_target(path=path, r_bias=r_bias)
	exit()

#=======================================

# Function definition is here

import types
import os
from os import listdir
from os.path import isfile, join
import h5py
from math										import sqrt, exp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import read_exp_data as red
import csv
import sys
from scipy 										import interpolate
from scipy.interpolate 						import UnivariateSpline

#==============================================================================
# This routine plots ne, Te e transport parameters
#==============================================================================

def plot_target(path=[], r_bias=0.0):

	file = h5py.File('data55064.h5', 'r')
#	print(np.array(file.keys()))
#	print(np.array(file['OMP']))
	
#	dict_group_load = file['OMP']
#	dict_group_keys = dict_group_load.keys()
	
#	hf= {}
#	for k in dict_group_keys:
#		hf[k]= dict_group_load[k][:]
	
#	R_sep = np.array(hf['Rsep'])
	
#	print(np.array(file.keys()))

# ################## #
#	Lettura dati sperimentali langmuir outer target
# ################## #
	dict_group_load = file['LPo']
	dict_group_keys = dict_group_load.keys()
	
	hf= {}
	for k in dict_group_keys:
		hf[k]= dict_group_load[k][:]
	#print(hf.keys())
	R = np.array(hf['R'])
	Te = np.array(hf['Te'])
	ne = np.array(hf['ne'])
	qperpLP = np.array(hf['qperp'])
	JparLP = np.array(hf['Jpar'])
	Terms = np.array(hf['Terms'])
	nerms = np.array(hf['nerms'])
#	idx=np.argmin(abs(qperpLP-max(qperpLP)))
#	R_sep=R[idx]																	#dal massimo di q delle langmuir 
	R_sep=2.24122																#dalla simulazione
	print("Rsep_ot=",R_sep)

#	Lettura dati sperimentali IR
	dict_group_load = file['IRo']
	dict_group_keys = dict_group_load.keys()

	hf= {}
	for k in dict_group_keys:
		hf[k]= dict_group_load[k][:]
	Rir= np.array(hf['R'])
	qparIR= np.array(hf['qpar'])
	qperpIR= np.array(hf['qperp'])

# ################## #
#	Lettura dati sperimentali langmuir outer target
# ################## #
	dict_group_load = file['LPi']
	dict_group_keys = dict_group_load.keys()
	
	hf= {}
	for k in dict_group_keys:
		hf[k]= dict_group_load[k][:]
	#print(hf.keys())
	R_it = np.array(hf['R'])
	Te_it = np.array(hf['Te'])
	ne_it = np.array(hf['ne'])
	qperpLP_it = np.array(hf['qperp'])
	JparLP_it = np.array(hf['Jpar'])
	Terms_it = np.array(hf['Terms'])
	nerms_it = np.array(hf['nerms'])
	R_sep_it=2.13515															#dalla simulazione
	print("Rsep_it=",R_sep_it)

	file.close()

	font = {'family' : 'normal',
		'weight' : 'bold',
		'size'   : 15}

	mpl.rc('font', **font)

	#some physical constants
	c_vac=2.99792458e8
	h_planck=6.626070040e-34
	kB=1.3806488E-23
	qe=1.60217657E-19
	mproton=1.67262178e-27
	melectron=9.10938291e-31 
	gammai = 2.5
	gammae  = 4.5
	Eiondiss=15.6
	amass = 2
	mi = amass*mproton

	folderList =["./"]
	labelc = [""]
	markerc = [".",".","v","p","x","x"]
	Zorderc = [10,10]
	colorc = ["green","orange","red","grey","blue"]


	fig,ax  = plt.subplots(1,4)#,sharex = True);
	fig.patch.set_facecolor("white")
	fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)


	for i,el in enumerate(folderList):
		onlyfiles = [f for f in listdir(el) if "meter_data.csv" in f]
		Dict_case= red.read_csv(el+onlyfiles[0],flag = False)
		MM = Dict_case["DIST"]
		onlyfiles = [f for f in listdir(el) if "rho_data.csv" in f]
		Dict_case= red.read_csv(el+onlyfiles[0],flag = False)
		print(Dict_case)
		rho_rho = np.array(Dict_case['Psi_pol'])-1
		spl = UnivariateSpline(rho_rho,MM)
		onlyfiles = [f for f in listdir(el) if "target_meter" in f]
		Dict_case= red.read_csv(el+onlyfiles[0],flag = True)
		
		#RRRR = (np.array(spl(np.array(Dict_case['dist_sp2'])))) #se carichi i dati 
		RRRR = np.array(Dict_case['dist_sp2'])-R_sep+2.24122 #se carichi i dati 
		
		ax[0].plot(RRRR,np.array(Dict_case['n_e_sp2 (*10^{19}m^{-3})']),color = colorc[i], linewidth = 5,zorder =Zorderc[i])
		ax[1].plot(RRRR,np.array(Dict_case['T_e_sp2 (eV)']),color = colorc[i],linewidth = 10,zorder =Zorderc[i],label = labelc[i])
		#ax[0].plot(R-R_sep+0.71,ne*1e-19,"ko")
		ax[0].errorbar((R-R_sep+r_bias)*1.078,ne*1e-19, yerr=nerms*1e-19,marker = "o",fmt='.',color = "grey",ms = 5)															#Li metto nella stessa posizione e moltiplico per il coeff. angolare
		ax[1].errorbar((R-R_sep+r_bias)*1.078,Te, yerr=Terms,marker = "o",fmt='.',color = "grey",ms = 5)
		ax[2].plot(RRRR,np.array(Dict_case['\\Gamma_{Etot}_i_sp2 (MW/m^{2})']),color = colorc[i],linewidth = 10,zorder =Zorderc[i],label = labelc[i])
		ax[3].plot(RRRR,np.array(Dict_case['J_{sat-par}_i_sp2 (kA/m^{2})']),color = colorc[i],linewidth = 10,zorder =Zorderc[i],label = labelc[i])
		ax[2].plot((Rir-R_sep+r_bias)*1.078,qperpIR*1e-6)
		ax[2].plot((R-R_sep+r_bias)*1.078,qperpLP*1e-6)
		ax[3].plot((R-R_sep+r_bias)*1.078,JparLP*1e-3)

		# ################## #
		#Calcolo integrale della j_sat simulata
		jsat_sim=np.array(Dict_case['J_{sat-par}_i_sp2 (kA/m^{2})'])
		idx=np.argmin(jsat_sim)
		R_sim_min=RRRR[idx]
		J_sat_sim_int=0.0
		print("idx=",idx,"  j_sat_mix=",jsat_sim[idx], "R_sim_min=",R_sim_min)
		for j in range(idx,len(RRRR)-1,1):
			J_sat_sim_int=J_sat_sim_int+2*3.1416*(RRRR[j]+R_sep)*(RRRR[j+1]-RRRR[j])*jsat_sim[j]
		# ################## #
		#Calcolo integrale della j_sat sperimentale
#		idx_exp=np.argmin(abs(R-R_sep+r_bias-R_sim_min))
#		print("idx_exp=",idx_exp,"R_exp_min=",R[idx_exp]-R_sep+r_bias)
		J_sat_exp_int=0.0
		qperp_exp_LP_ot_int=0.0
		qperp_exp_LP_it_int=0.0
		j_sat_exp=JparLP*1e-3																								#jpar in kA/m2
		qperp_exp_LP_ot=qperpLP*1e-6																				#qper in MW/m2
		qperp_exp_LP_it=qperpLP_it*1e-6																				#qper in MW/m2
		for j in range(len(R)-1):
			J_sat_exp_int=J_sat_exp_int+2*3.1416*(R[j]+r_bias)*(R[j+1]-R[j])*j_sat_exp[j]*1.078
			qperp_exp_LP_ot_int=qperp_exp_LP_ot_int+2*3.1416*(R[j]+r_bias)*(R[j+1]-R[j])*qperp_exp_LP_ot[j]*1.078
		for j in range(len(R_it)-1):
			qperp_exp_LP_it_int=qperp_exp_LP_it_int+2*3.1416*(R_it[j]+r_bias)*(R_it[j+1]-R_it[j])*qperp_exp_LP_it[j]*1.078
			
		print("J_sat_sim_int=",J_sat_sim_int,"kA","   J_sat_exp_int=",J_sat_exp_int,"kA")
		print("qperp_exp_LP_int outer target=",qperp_exp_LP_ot_int )
		print("qperp_exp_LP_int inner target=",qperp_exp_LP_it_int )
	
	ax[1].legend(loc = "best")
	
	
	ax[0].axvline(x = 0.0,color = "black", linestyle = "--")
	ax[1].axvline(x = 0.0,color = "black", linestyle = "--")
	
	
	ax[0].set_title(r"$n_e(10^{20}m^{-3})$", y=1.022,size =20)
	ax[1].set_title("$Te(eV)$",y=1.022,size =20)
	#ax[2].set_title(r"$\Gamma_{Etot}_i_sp2 (MW/m^{2})$", y=1.022,size =20)
	ax[3].set_title("$J_{satpar} (kA/m^{2})$",y=1.022,size =20)
	ax[0].set_xlabel(r"$r-r_{sep} (mm)$",size =20)
	ax[1].set_xlabel(r"$r-r_{sep} (mm)$",size =20)
	
	
	plt.show()

	return