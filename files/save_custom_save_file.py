
import numpy as np
import os
from routines.globals			import DEBUG

#=================================================
# This routine write soledge2D.elemente file
#=================================================

def save_custom_save_file(CustomPath, CustomPlots):

	if(DEBUG > 0):	print("save_custom_save_file")

	nPlots = 0
	for iType in range(len(CustomPlots)): nPlots += len(CustomPlots[iType])

	fid = open(CustomPath+"custom_save.txt",'w')
	fid.write("nPlots = {:d}\n".format(nPlots))
	fid.write("Type 1 = parallel ; Type 2 = perp ; Type 3 = temporal\n")

	nPlots = 0
	for iType in range(len(CustomPlots)):
		if(len(CustomPlots[iType]) > 0):
			for iPlot in range(len(CustomPlots[iType])):
				nPlots += 1
				fid.write("########## plot {:d} ############\n".format(nPlots))
				fid.write("Name = {:s}\n".format(CustomPlots[iType][iPlot].Name))
				fid.write("Type = {:d}\n".format(iType+1))
				fid.write("nZones = {:d}\n".format(len(CustomPlots[iType][iPlot].kZones)))
				if(iType != 2):
					fid.write("Zones = ") 
					for iz in range(len(CustomPlots[iType][iPlot].kZones)-1): fid.write("{:d},".format(CustomPlots[iType][iPlot].kZones[iz]+1))
					fid.write("{:d}\n".format(CustomPlots[iType][iPlot].kZones[-1]+1))
					fid.write("Coord = {:d}\n".format(CustomPlots[iType][iPlot].Coords[iType+1]+1))
				else:
					fid.write("Zones = {:d}\n".format(CustomPlots[iType][iPlot].Coords[0]+1)) 
					fid.write("Coord = {:d}\n".format(CustomPlots[iType][iPlot].Coords[1]+1))
					fid.write("Coord2 = {:d}\n".format(CustomPlots[iType][iPlot].Coords[2]+1))

				fid.write("\n")

	fid.close()

	if(DEBUG > 0):	print("save_custom_save_file: Completed")
	
