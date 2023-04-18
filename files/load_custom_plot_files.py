# Function definition is here

import numpy as np
import types
import os
from math					 import sqrt
from files.load_refpar_file	 import load_refpar_file
from files.load_ions_list	 import load_ions_list
from routines.globals		 import DEBUG, EV, KB

#=========================================================
# This routine read custum plot files
#=========================================================

def load_custom_plot_temporal_files(Custom_name, Path="", iCustoms = []): 

	if((len(Path) > 0) and (Path[-1] != "/")): Path = Path + "/"

	try:
		ions = load_ions_list(Path+"../")
	except:
		ions = ["e","D+"]

	if(len(iCustoms) == 0):	iCustoms = [k for k in range(len(ions))]

#	Custom files

	Customs = []
	for i in range(len(iCustoms)):
		Customs.append(types.SimpleNamespace())
		fid = open(Path+"custom_plots/"+Custom_name+"_{:d}.txt".format(iCustoms[i]),'r')

		Customs[i].Type=int(eval(fid.readline()[:-1].split("=")[1]))
		if(Customs[i].Type != 3):
			print("\tInvalid file type for: ",Path+"custom_plots/"+Custom_name+"_{:d}.txt".format(iCustoms[i]))
			fid.close
			continue

		RZgeom = fid.readline()
		header = fid.readline()
		fid.close

		RZgeom = eval(RZgeom[:-1].split("=")[1])
		Customs[i].Rgeom=RZgeom[0]
		Customs[i].Zgeom=RZgeom[1]

		Customs[i].Names=header[:-1].split(",")

		Customs[i].Values=np.loadtxt(Path+"custom_plots/"+Custom_name+"_{:d}.txt".format(iCustoms[i]), delimiter=",",skiprows=3)

	return Customs


def load_custom_plot_line_files(Custom_name, Path="", iCustoms = []): 

	if((len(Path) > 0) and (Path[-1] != "/")): Path = Path + "/"

	try:
		ions = load_ions_list(Path+"../")
	except:
		ions = ["e","D+"]

	if(len(iCustoms) == 0):	iCustoms = [k for k in range(len(ions))]

#	Custom files

	Customs = []
	for i in range(len(iCustoms)):
		Customs.append(types.SimpleNamespace())
		fid = open(Path+"custom_plots/"+Custom_name+"_{:d}.txt".format(iCustoms[i]),'r')

		Customs[i].Type=int(eval(fid.readline()[:-1].split("=")[1]))
		if((Customs[i].Type < 1) or (Customs[i].Type > 2)):
			print("\tInvalid file type for: ",Path+"custom_plots/"+Custom_name+"_{:d}.txt".format(iCustoms[i]))
			fid.close
			continue

		PointsVar = eval(fid.readline()[:-1].split("=")[1])
		Customs[i].nPoints=int(PointsVar[0])
		Customs[i].nVar=PointsVar[1]
		Customs[i].Rgeom=np.empty(Customs[i].nPoints, dtype='f8')
		Customs[i].Zgeom=np.empty(Customs[i].nPoints, dtype='f8')
		nLoops=int(Customs[i].nPoints/Customs[i].nVar)

		if(nLoops*Customs[i].nVar < Customs[i].nPoints): nLoops += 1
		k1=0
		for k in range(nLoops):
			k2=min(k1+Customs[i].nVar,Customs[i].nPoints)
			Customs[i].Rgeom[k1:k2] = eval(fid.readline()[:-1])
			Customs[i].Zgeom[k1:k2] = eval(fid.readline()[:-1])
			k1=k2

		header = fid.readline()
		Names=header[:-1].split(",")
		Customs[i].Names=[]
		for Name in Names: Customs[i].Names.append(Name.strip())
		
		fid.close

		Data=np.loadtxt(Path+"custom_plots/"+Custom_name+"_{:d}.txt".format(iCustoms[i]), delimiter=",",skiprows=2*nLoops+3)
		Customs[i].nTimes	= int(Data.shape[0]/(Customs[i].nPoints+1))
		Data				= Data.reshape(Customs[i].nTimes,Customs[i].nPoints+1,Customs[i].nVar)
		Customs[i].Times	= np.copy(Data[:,0,0])

		NoZero=np.where(Customs[i].Rgeom > 0.)[0]
		if(len(NoZero) < Customs[i].nPoints):
			Customs[i].Rgeom=Customs[i].Rgeom[NoZero]
			Customs[i].Zgeom=Customs[i].Zgeom[NoZero]
			Customs[i].Values	= np.copy(Data[:,NoZero+1,:])
			Customs[i].nPoints	= len(NoZero)
		else:
			Customs[i].Values	= np.copy(Data[:,1:,:])

	return Customs
