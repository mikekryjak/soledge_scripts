#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot_OMP 						import plot_OMP
	from routines.cli_routines		import *

	if(cli_present("-h",sys.argv) and cli_present("-exp_file",sys.argv) ):
		from routines.exp_data_routines			import get_exp_data_point_help
		from files.load_exp_descr						import load_exp_descr_help
		load_exp_descr_help()
		get_exp_data_point_help()
		exit()

	if(cli_present("-h",sys.argv)):
		print("\nThis script plot plasma parameters on soledge mesh along a line starting from plasma core\n")
		print("plot_OMP options")
		print("\t-path           [BROKEN] Directory with simulation [d='./']")
		print("\t-r_bias         Bias of R coordinate[m] [d='0.00']")
		print()
		exit()

	path = cli_get_value("-path",			sys.argv,   [""])
	r_bias = cli_get_value("-r_bias",		sys.argv,	0.00)
	plot_OMP(path=path, r_bias=r_bias)
	exit()

#=======================================

# Function definition is here

import types
import os
from os import listdir
from os.path import isfile, join
import h5py
from math										import sqrt, exp
import numpy 								as np
import matplotlib.pyplot 				as plt
import csv

import matplotlib as mpl
import read_exp_data as red
import sys
from scipy import interpolate


#==============================================================================
# This routine plots ne, Te e transport parameters
#==============================================================================

def plot_OMP(path=[], r_bias=0.0):
	
	print("plot_OMP")

	file = h5py.File('data55064.h5', 'r')
	print(np.array(file.keys()))
	print(np.array(file['OMP']))
	dict_group_load = file['OMP']
	dict_group_keys = dict_group_load.keys()
	
	hf= {}
	for k in dict_group_keys:
		hf[k]= dict_group_load[k][:]
	
	R = np.array(hf['R'])+r_bias
	R_sep = np.array(hf['Rsep'])
	ne_iterf = np.array(hf['ne_iterf'])
	ne_refl =  np.array(hf['ne_refl'])
	ne_iterf_rms = np.array(hf['nerms_interf'])
	ne_refl_rms=  np.array(hf['nerms_refl'])
	file.close()



	font = {'family' : 'normal',
			'weight' : 'bold',
			'size'   : 20}

	mpl.rc('font', **font)

	controlwidth = 4



	folderList =["./"]
	labelc = ["test data west"]
	markerc = [".",".","v","p","x","x"]
	Zorderc = [10]
	colorc = ["green","orange","red","grey","blue"]


	fig,ax  = plt.subplots(1,3)#,sharex = True);
	fig.patch.set_facecolor("white")
	fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)



	for i,el in enumerate(folderList):
		onlyfiles = [f for f in listdir(el) if "meter_data" in f]
		print(el,onlyfiles)
		Dict_not_cost= red.read_csv(el+onlyfiles[0])
		ax[0].plot((np.array(Dict_not_cost['DIST'])),np.array(Dict_not_cost['Dense']),linewidth = controlwidth,color = colorc[i],label = labelc[i],zorder =Zorderc[i])
		ax[1].plot((np.array(Dict_not_cost['DIST'])),np.array(Dict_not_cost['Tempe']),linewidth = controlwidth,color = colorc[i],zorder =Zorderc[i])
		
		ax[2].plot((np.array(Dict_not_cost['DIST'])),np.array(Dict_not_cost['D']),color = "blue",label = "D",linewidth = controlwidth,zorder =Zorderc[i])
		ax[2].plot((np.array(Dict_not_cost['DIST'])),np.array(Dict_not_cost['Chie']),color = "red",label = "chi",linewidth = controlwidth,zorder =Zorderc[i])
		#exp_data
		ax[0].errorbar(R-R_sep,ne_refl, yerr=ne_refl_rms,marker = "o",fmt='.',color = "grey",ms = 2)

	ax[0].legend(loc= "best")
	ax[2].legend(loc= "best")
	#ax[0].set_xlabel(r"$\psi_{norm}$",size =39)
	#ax[1].set_xlabel(r"$\psi_{norm}$",size =39)
	#ax[2].set_xlabel(r"$\psi_{norm}$",size =39)
	ax[0].set_title(r"$n(m^{-3})$", y=1.022,size =50)
	ax[1].set_title(r"$Te(eV)$", y=1.022,size =50)
	ax[2].set_title(r"$D(m^{2}/s)$,$\chi(m^{2}/s)$", y=1.022,size =50)
	ax[2].set_ylim(0,1.5)
	for i in range(0,3):
		ax[i].axvline(x = 0.0,color = "black", linestyle = "--",lw= 2)
		ax[i].grid(which='minor', alpha=0.2)
		ax[i].grid(which='major', alpha=0.5)
		#ax[i].set_xlim(-0.2,0.3)
		ax[i].set_xlabel(r"$R-R_{sep} (m)$",size =39)



	plt.show()

	return