
import numpy as np
import os
import types
from routines.globals			import DEBUG

#=================================================
# This routine write soledge2D.elemente file
#=================================================

def load_custom_save_file(CustomPath):

	if(DEBUG > 0):	print("load_custom_save_file")

	try:
		fid = open(CustomPath+"custom_save.txt",'r')
	except:
		return []

	CustomPlots=[[],[],[]]
	nPlots = eval(fid.readline()[:-1].split("=")[1])
	fid.readline()
	for iPlot in range(nPlots):
		fid.readline()
		Name	=fid.readline()[:-1].split("=")[1].strip()
		iType		= eval(fid.readline()[:-1].split("=")[1])-1

		CustomPlots[iType].append(types.SimpleNamespace())
		CustomPlots[iType][-1].Name = Name
		nZones=eval(fid.readline()[:-1].split("=")[1])
		CustomPlots[iType][-1].Coords = np.zeros(3,dtype='i4')
		if(iType != 2):
			CustomPlots[iType][-1].kZones = np.array(eval(fid.readline()[:-1].split("=")[1]))-1
			CustomPlots[iType][-1].Coords[0] = CustomPlots[iType][-1].kZones[0]
			CustomPlots[iType][-1].Coords[iType+1] = eval(fid.readline()[:-1].split("=")[1])-1
		else:
			CustomPlots[iType][-1].Coords[0] = eval(fid.readline()[:-1].split("=")[1])-1
			CustomPlots[iType][-1].Coords[1] = eval(fid.readline()[:-1].split("=")[1])-1
			CustomPlots[iType][-1].Coords[2] = eval(fid.readline()[:-1].split("=")[1])-1

		CustomPlots[iType][-1].SubPlots	= []
		CustomPlots[iType][-1].nRows		= 1
		CustomPlots[iType][-1].nCols			= 1
		CustomPlots[iType][-1].SameX		= True


		fid.readline()

	fid.close()

	if(DEBUG > 0):	print("load_custom_save_file: Completed")

	return CustomPlots

