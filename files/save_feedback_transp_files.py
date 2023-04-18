#!/usr/bin/python

# Function definition is here

import numpy				as np
import os
from routines.globals		import DEBUG

#=========================================================
# This routine to write data to the H5DF mesh file
#=========================================================

def save_feedback_transp_files(FeedbackFilePath, FeedbackTransp): 

	if(DEBUG > 0): print("save_feedback_transp_files: Saving to ",FeedbackFilePath)

	try:
		os.mkdir(FeedbackFilePath)
	except OSError:
		pass

#	FeedbackTransp profiles input file

	"""
#	Old versione

	FeedProfiles = np.zeros((FeedbackTransp.Data.Profiles.shape[0]+2, FeedbackTransp.Data.Profiles.shape[1]),dtype=FeedbackTransp.Data.Profiles.dtype)
	FeedProfiles[3,:]   = FeedbackTransp.Data.Profiles[0,:]							#Ne
	FeedProfiles[4:6,:] = FeedbackTransp.Data.Profiles[2:4,:]							#Te & Ti
	FeedProfiles[6,:]   = FeedbackTransp.Data.Profiles[5,:]							#Chi
	FeedProfiles[7,:]   = FeedbackTransp.Data.Profiles[4,:]							#D
	FeedProfiles[8,:]   = FeedbackTransp.Data.Profiles[-1,:]							#Length, flux, flux norm or Rho in profile definition
	"""

	FeedProfiles = np.zeros((FeedbackTransp.Data.Profiles.shape[0]+2, FeedbackTransp.Data.Profiles.shape[1]+1),dtype=FeedbackTransp.Data.Profiles.dtype)
	FeedProfiles[3,1:]   = FeedbackTransp.Data.Profiles[0,:]							#Ne
	FeedProfiles[4:6,1:] = FeedbackTransp.Data.Profiles[2:4,:]							#Te & Ti
	FeedProfiles[6,1:]   = FeedbackTransp.Data.Profiles[5,:]							#Chi
	FeedProfiles[7,1:]   = FeedbackTransp.Data.Profiles[4,:]							#D
	FeedProfiles[8,1:]   = FeedbackTransp.Data.Profiles[-1,:]							#Length, flux, flux norm or Rho in profile definition

	FeedProfiles[:,0]	 = FeedProfiles[:,1]									#Auto set boundary values as the internal one						  
	np.savetxt(FeedbackFilePath+"/radialFeedbackProfiles.txt",	FeedProfiles.T, fmt="%15.7e")
	FeedProfiles = 0

#	FeedbackTransp option input file
	
	fid = open(FeedbackFilePath+"/radialFeedback.txt",'w')

	fid.write("#########################################################################################\n")
	fid.write("####   File containing informations for automatic adjustment of diffusion        ########\n")
	fid.write("####   coefficients to match prescribed density and temperature radial profiles  ########\n")
	fid.write("#########################################################################################\n")
	fid.write("\n")
	fid.write("Position of the profile to fit:\n")
	fid.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
	fid.write("\n")
	fid.write("Number of zones involved:\n")
	fid.write("Nzones = {:d}\n".format(len(FeedbackTransp.Data.nZ)))
	fid.write("\n")
	
	fid.write("List of zones:\n")
	fid.write("zones = ")
	for k in range(len(FeedbackTransp.Data.nZ)-1): fid.write("{:d}, ".format(FeedbackTransp.Data.nZ[k]+1))
	fid.write("{:d}\n".format(FeedbackTransp.Data.nZ[-1]+1))

	fid.write("\n")
	fid.write("Poloidal position:\n")
	fid.write("Nz = {:d}\n".format(FeedbackTransp.Data.jZ+1))
	fid.write("\n")
	fid.write("\n")
	fid.write("Thresholds:\n")
	fid.write("~~~~~~~~~~~\n")
	fid.write("\n")
	fid.write("Minimal diffusivity:\n")
	fid.write("Dmin = {:0.2e}\n".format(FeedbackTransp.Data.Dmin))
	fid.write("\n")
	fid.write("Maximal diffusivity:\n")
	fid.write("Dmax = {:0.2e}\n".format(FeedbackTransp.Data.Dmax))
	fid.write("\n")
	fid.write("Iteration window to evolve D:\n")
	fid.write("keep = {:0.2e}\n".format(FeedbackTransp.Data.Keep))
	fid.write("\n")
	fid.write("Gain for feedback loop (0 no FB):\n")
	fid.write("gain = {:0.2e}\n".format(FeedbackTransp.Data.Gain))
	fid.write("\n")
	fid.write("Gradient Gain for feedback loop (0 no FB):\n")
	fid.write("gainG = {:0.2e}\n".format(FeedbackTransp.Data.GainG))
	fid.write("\n")
	fid.write("Decay length in private:\n")
	fid.write("lambda = {:0.2e}\n".format(FeedbackTransp.Data.Lambda))
	fid.write("\n")
	fid.write("Use ballooning pattern stored in mesh.h5:\n")
	if(FeedbackTransp.Data.UsePattern == "No"):
		fid.write("use_pattern = F\n")
	else:
		fid.write("use_pattern = T\n")
		
	fid.close()

	return