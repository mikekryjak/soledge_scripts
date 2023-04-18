#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from routines.cli_routines import cli_present, cli_get_value, cli_get_value
	from evolution_to_results import evolution_to_results
	 
	if(cli_present("-h",sys.argv)):
		print("\nThis script copy Evolution plasma files to Results and fix balances and residuals files\n")
		print("evolution_to_results options")
		print("\t-evolution <value>  Version to copy in Results, 100=Last, -1=Progress, 0=Results [d=100]")
		print()
		exit()

	evolution	 = cli_get_value("-evolution",	sys.argv, 100)
	evolution_to_results(evolution=evolution)
	exit()

#=======================================

# Function definition is here

import os
import glob
import shutil
import numpy as np
import h5py
from files.load_ions_list				import load_ions_list
from routines.h5_routines				import h5_read, h5_write
from routines.globals					import DEBUG

def evolution_to_results(evolution=100):

#	Move plasma files from Evolution/Progress to Results

	if(DEBUG > 0): print("evolution_to_results")

	Results_exists		= os.path.isdir("Results")
	Evolution_exists	= os.path.isdir("Evolution")
	Progress_exists		= os.path.isdir("Progress")

	if(not Results_exists):
		print("\tERROR: Results dose not existS")
		return

	ions = load_ions_list("./")

#	Find maximum evolution file

	if(evolution == 100):
		evo_times = []
		evo_numbers = []
		if(Progress_exists):
			if_global = h5py.File("Progress/globals", "r")
			evo_times.append(h5_read(if_global, "tempus", keep_array=False))								#read time
			if_global.close()
			evo_numbers.append(-1)

		if(Evolution_exists):			
			plasma_0 = glob.glob("Evolution/*plasma_0")
			for name in plasma_0: 
				if_plasma = h5py.File(name, "r")
				evo_times.append(h5_read(if_plasma, "tempus", keep_array=False))							#read time
				if_plasma.close()
				evo_numbers.append(eval(name[10:name.index("_")]))

		it_max	  = np.argmax(np.array(evo_times))
		evo_time  = evo_times[it_max]
		evolution = evo_numbers[it_max]

	if(evolution == 0):
		if_global = h5py.File("Results/globals", "r")
		evo_time = h5_read(if_global, "tempus", keep_array=False)									#read time
		if_global.close()

		if(DEBUG > 1): print("\tKeep Results and cut balances and residuals at time {:.5f} (s)".format(evo_time))

	elif(evolution == -1):
		try:	
			if_global = h5py.File("Progress/globals", "r")
			evo_time = h5_read(if_global, "tempus", keep_array=False)								#read time
			if_global.close()
		except:
			if(DEBUG > 0): print("\tERRROR: Unable to open Progress/globals")
			return

		if(DEBUG > 1): print("\tUse Progress with time {:.5f} (s)".format(evo_time))

		for ion in range(len(ions)):
			plasma_name = "plasma_{:d}".format(ion)
			fluxes_name = "fluxes_{:d}".format(ion)
			shutil.copy2("Progress/"+plasma_name, "Results/"+plasma_name)
			shutil.copy2("Progress/"+fluxes_name, "Results/"+fluxes_name)

		shutil.copy2("Progress/globals", "Results/globals")
		shutil.copy2("Progress/eirene_neutrals", "Results/eirene_neutrals")

	elif(evolution > 0):
		try:	
			if_plasma = h5py.File("Evolution/{:d}_plasma_0".format(evolution), "r")
			evo_time = h5_read(if_plasma, "tempus", keep_array=False)								#read time
			if_plasma.close()
		except:
			if(DEBUG > 0): print("\tERRROR: Unable to open Evolution/{:d}_plasma_0".format(evolution))
			return

		if(DEBUG > 1): print("\tUse evolution n. {:d} with time {:.5f} (s)".format(evolution,evo_time))

		for ion in range(len(ions)):
			plasma_name = "Results/plasma_{:d}".format(ion)
			evo_name	= "Evolution/{:d}_plasma_{:d}".format(evolution,ion)
			shutil.copy2(evo_name, plasma_name)

	#	set time in globals

		if_global = h5py.File("Results/globals", "w")
		h5_write(if_global, "tempus", evo_time)
		if_global.close()

	else:
		print("\tERROR: Invalid evolution number {:d}".format(evolution))
		return

#	cut balances_## and residual files to evo file time

	FirstFile = True
	for ion in range(len(ions)):
		FileData = np.loadtxt("balances_{:d}".format(ion), dtype='f8')
		if(FirstFile):
			FirstFile = False
			times     = FileData[:,-1]
			nOkTimes  = np.max(np.where(times <= evo_time)[0])+1


		FileData = FileData[:nOkTimes,:]
		np.savetxt("balances_{:d}".format(ion), FileData, fmt="%15.7e", delimiter="")

	FileData = np.loadtxt("residuals", dtype='f8')
	FileData = FileData[:nOkTimes,:]
	np.savetxt("residuals", FileData, fmt="%15.7e", delimiter="")

	if(DEBUG > 0): print("evolution_to_results: Completed")

	
	return
	
