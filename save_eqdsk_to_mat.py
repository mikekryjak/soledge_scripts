#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from save_eqdsk_to_mat			import save_eqdsk_to_mat
	from routines.cli_routines		import *

	if(cli_present("-h",sys.argv)):
		print("\nThis script save in equilibrium fields from eqdsk file to matlab format for matalb grid-gen\n")
		print("save_eqdsk_to_mat options")
		print("\t-file <name>   input eqdsk filename [mag.eqdsk]")
		print("\t-mag_file   <name>     output fields matlab filename [mag.mat]")
		print("\t-wall_file  <name>    output wall matlab filename [wall.mat]")
		print("\t-nr         <value>   if > than in eqdsk grid sinze in r [0]")
		print("\t-nz         <value>   if > than in eqdsk grid sinze in z [0]")
		print("\t-smooth     <value>   smooth factor (>= 0.5) used if nr or nz > grid size [1]")
		print("\t-out_length <value>   if > -1. length to cut grid outside wall [D=-1]")
		print()
		exit()

	eqdsk_file		= cli_get_value("-file",			sys.argv, "mag.eqdsk")
	mag_file		= cli_get_value("-mag_file",		sys.argv, "mag.mat")
	wall_file		= cli_get_value("-wall_file",	sys.argv, "wall.mat")
	nr				= cli_get_value("-nr",				sys.argv, 0)
	nz				= cli_get_value("-nz",				sys.argv, 0)
	smooth			= cli_get_value("-smooth",		sys.argv, 1.)
	out_length		= cli_get_value("-out_length",	sys.argv, -1.)
	save_eqdsk_to_mat(eqdsk_file, mag_file, wall_file, nr=nr, nz=nz, smooth=smooth, out_length=out_length)
	exit()

#=======================================

# Function definition is here

import scipy.io
import numpy 					as np
from files.eqdsk_routines		import load_eqdsk_file, eqdsk_compute_fields
from routines.globals			import DEBUG

#=========================================================
# This routine write edsk file to mat file
#=========================================================
#

def save_eqdsk_to_mat(eqdsk_file="mag.eqdsk", mag_file="mag.mat", wall_file="wall.mat", nr=0, nz=0, smooth=1., out_length=-1.):

	if(DEBUG > 0):	print("save_eqdsk_to_mat")
	if(DEBUG > 1):	print("\treading magnetic file")

	eqdsk = load_eqdsk_file(eqdsk_file)

	Wall  = {}
	Wall['Rwall'] = np.reshape(eqdsk.rlim,(eqdsk.rlim.shape[0],1))
	Wall['Zwall'] = np.reshape(eqdsk.zlim,(eqdsk.rlim.shape[0],1))

	if(DEBUG > 1):	print("\tcomputing magnetic file")

	rz_b_psi  = eqdsk_compute_fields(eqdsk, nr=nr, nz=nz, smooth=smooth, out_length=out_length)

	Mag = {}
	Mag['r2D']			= np.copy(rz_b_psi[:,:,0])
	Mag['z2D']			= np.copy(rz_b_psi[:,:,1])
	Mag['Br2D']		= np.copy(rz_b_psi[:,:,2])
	Mag['Bz2D']		= np.copy(rz_b_psi[:,:,3])
	Mag['Bphi2D']		= np.copy(rz_b_psi[:,:,4])
	if(eqdsk.simag > eqdsk.sibry):
		Mag['flux2D']	= np.copy(-rz_b_psi[:,:,5])
	else:
		Mag['flux2D']	= np.copy(rz_b_psi[:,:,5])


	if(DEBUG > 1):	print("\twriting magnetic file")

	scipy.io.savemat(mag_file, Mag)

	if(DEBUG > 1):	print("\twriting wall file")

	scipy.io.savemat(wall_file, Wall)

	if(DEBUG > 0):	print("save_eqdsk_to_mat: Completed")

