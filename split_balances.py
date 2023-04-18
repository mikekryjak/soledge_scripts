#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from routines.cli_routines import cli_present, cli_get_value, cli_get_value
	from split_balances import split_balances
	 
	if(cli_present("-h",sys.argv)):
		print("\nThis script split balances and residuals files to speed-up reading\n")
		print("split_balances options")
		print("\t-path        Directory with simulation [d='./']")
		print("\t-keep <value>  numbero of points to keep [d=2000]")
		print()
		exit()

	path	 	= cli_get_value("-path", sys.argv, "")
	keep		= cli_get_value("-keep", sys.argv, 2000)
	split_balances(path=path, keep=keep)
	exit()

#=======================================

# Function definition is here

import os
import numpy as np
from routines.h5_routines				import h5_read, h5_write
from routines.globals					import DEBUG
from files.load_ions_list				import load_ions_list

def split_balances(path="", keep=2000):

#	Split balances and residual files

	if(DEBUG > 0): print("split_balances")

	if((len(path) > 0) and (path[-1] != "/")): path = path + "/"
	ions	= load_ions_list(path)										#read ions list

#	cut balances_## and residual files 

	for ion in range(len(ions)+1):
		if(ion == len(ions)):
			file_name	  = path+"residuals"
			pre_file_name = path+"pre_residuals"
		else:
			file_name	  = path+"balances_{:d}".format(ion)
			pre_file_name = path+"pre_balances_{:d}".format(ion)

		FileData = np.loadtxt(file_name, dtype='f8')

		if(FileData.shape[0] > keep):
			if(os.path.exists(pre_file_name+".npy")):
				Pre_FileData = np.load(pre_file_name+".npy")
				Pre_FileData = np.append(Pre_FileData, np.float32(FileData[:-keep,:]), axis =0)
				np.save(pre_file_name+".npy",Pre_FileData)
				Pre_FileData = 0
			else:
				np.save(pre_file_name+".npy",np.float32( FileData[:-keep,:]))

			np.savetxt(file_name, FileData[-keep:,:], fmt="%15.7e", delimiter="")

	if(DEBUG > 0): print("split_balances: Completed")
	
	return
	
