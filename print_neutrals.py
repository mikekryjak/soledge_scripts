#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from print_neutrals 		import print_neutrals
	from routines.cli_routines	import *

	if(cli_present("-h",sys.argv)):
		print("\nThis neutral parameters on point or averaged\n")
		print("print_neutrals options")
		print("\t-path        Directory with simulation or list of directories[d='./']")
		print("\t-rp          R coordinate point [d=0.]")
		print("\t-zp          Z coordinate point [d=0.]")
		print("\t-deltar      Amplitude in R rectangle to average around rp [d=0.]")
		print("\t-deltaz      Amplitude in Z rectangle to average around zp [d=0.]")
		print("\t-radius      Radius circle to average around ((rp,zp) [d=0.]")
		print()
		exit()

	path	 	= cli_get_value("-path",				sys.argv,   "")
	Rp			= cli_get_value("-rp",				sys.argv,	0.)
	Zp			= cli_get_value("-zp",				sys.argv,	0.)
	DeltaR		= cli_get_value("-deltar",		sys.argv,	0.)
	DeltaZ		= cli_get_value("-deltaz",		sys.argv,	0.)
	Radius		= cli_get_value("-radius",		sys.argv,	0.)

	print_neutrals(path=path,Rp=Rp,Zp=Zp,DeltaR=DeltaR,DeltaZ=DeltaZ,Radius=Radius)
	exit()

#=======================================

# Function definition is here

import types
import os
import h5py
from math								import sqrt
import numpy							as np
import matplotlib.pyplot				as pyp
import matplotlib.tri					as tri
from matplotlib.backends.backend_pdf	import PdfPages
from matplotlib.colors 					import LogNorm

from routines.h5_routines				import h5_read
from routines.utils_walls 				import plot2d_walls
from routines.globals					import DEBUG, KB, BALLOONING_NAMES

from mesh.get_rz_core_sep				import get_rz_core_sep
from files.load_soledge_mesh_file		import load_soledge_mesh_file
from files.load_plasma_files			import load_plasma_files


#==============================================================================
# This routine prints useful informations
#==============================================================================

def print_neutrals(path="",Rp=0.,Zp=0.,DeltaR=0.,DeltaZ=0.,Radius=0.):


	print("print_neutrals")

	if((len(path) > 0) and (path[-1] != "/")): path = path + "/"

#	Read mesh

	Config = load_soledge_mesh_file(path+"mesh.h5")
	nZones  = len(Config.Zones)

#	load Eirene mesh

	if_tri	 = h5py.File(path+"triangles.h5", "r")

	TriKnots = h5_read(if_tri,"triangles/tri_knots")
	TriKnots = TriKnots - 1
	R = h5_read(if_tri,"knots/R")*0.01
	Z = h5_read(if_tri,"knots/Z")*0.01
	if_tri.close()

# Load plasma results

	Plasmas = load_plasma_files(path, nZones=nZones, ToKnodes = 0)

# Define paramters to compute and print

	Parameters = ["Nni","Tni","Nmi","Tmi","Pni"]
	Names=["Atoms density (m-3)","Atoms temperature (eV)","Molecules density (m-3)","Molecules temperature (eV)", \
					"Neutral pressure (Pa)"]

	iPars=[]
	for par in Parameters:
		try:
			iPars.append(Plasmas[1][0].Triangles.VNames.index(par))
		except:
			print("\tERROR: not found parameter ",par)
			exit()

#	ATTENTION: averaging is not weighted on triangle size!
#	---------------------------------------------------------------------

	Values=[]
	if((DeltaR > 0.) and (DeltaZ > 0.) or (Radius > 0.)):
		if((DeltaR > 0.) and (DeltaZ > 0.)):																					#Average on rectangle
			ipts=np.where((np.abs(R-Rp) < 0.5*DeltaR) & (np.abs(Z-Zp) < 0.5*DeltaZ))[0]
		else:																																	#Average on circle
			d=np.sqrt((R-Rp)**2+(Z-Zp)**2)
			ipts=np.where(d<Radius)[0]

		print("R[ipts]=",R[ipts])
		print("Z[ipts]=",Z[ipts])
		Re=np.mean(R[ipts])
		Ze=np.mean(Z[ipts])
		for iPar in iPars:
			Values.append(np.mean(Plasmas[1][0].Triangles.Values[iPar][ipts]))
	else:
		d=np.sqrt((R-Rp)**2+(Z-Zp)**2)
		ip=np.argmin(d)
		Re=R[ip]
		Ze=Z[ip]
		for iPar in iPars:
			Values.append(np.mean(Plasmas[1][0].Triangles.Values[iPar][ip]))

	print("{:30} = {:8.4f}".format("R_eirene (m)",Re))
	print("{:30} = {:8.4f}".format("Z_eirene (m)",Ze))
	for i in range(len(iPars)):
		print("{:30} = {:8.4e}".format(Names[i],Values[i]))

# Plot triangles to check position

	Fig = pyp.figure()
	Fig.patch.set_facecolor('white')
	Ax = Fig.add_subplot(111)

	Ax.set_xlabel("$R\ (m)$")
	Ax.set_ylabel("$Z\ (m)$")
	Ax.set_title(os.path.basename(os.path.abspath(path)))

	Ax.tick_params(axis='both', which='major', labelsize=8)
	Ax.autoscale(axis='both', tight=True)

	Ax.triplot(R, Z, TriKnots, 'b-')
	plot2d_walls(Ax, Config.Walls)

	if((DeltaR > 0.) and (DeltaZ > 0.)):
		Ax.plot(R[ipts], Z[ipts], 'o', markersize=10,color="g")
		Ax.plot([Rp-DeltaR/2,Rp+DeltaR/2,Rp+DeltaR/2,Rp-DeltaR/2,Rp-DeltaR/2],
					 [Zp-DeltaZ/2,Zp-DeltaZ/2,Zp+DeltaZ/2,Zp+DeltaZ/2,Zp-DeltaZ/2], linestyle='-', color='g')
	elif(Radius > 0.):
		Ax.plot(R[ipts], Z[ipts], 'o', markersize=10,color="g")
		theta = np.linspace(0, 2*np.pi, 100)
		Rc = Radius*np.cos(theta)+Rp
		Zc = Radius*np.sin(theta)+Zp
		Ax.plot(Rc, Zc, linestyle='-', color='g')
	else:
		Ax.plot(Re, Ze, 'o', markersize=10,color="g")


	Ax.plot(Rp, Zp, 'x', markersize=10,color="r")

	Ax.set_aspect(1.)
	pyp.show()

	print("print_neutrals: Completed")

	return

