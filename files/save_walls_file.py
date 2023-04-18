
import numpy as np
import os
import types
from routines.globals			import DEBUG

#=================================================
# 
#=================================================

def save_walls_file(WallFile, Walls):

	if(DEBUG > 0):	print("save_walls_file")

	fid = open(WallFile,'w')

	for iWall in range(len(Walls)):
		Wall = Walls[iWall]
		fid.write("Wall = {:d}\n".format(iWall+1))
		line = "Type = {:d} LineType = {:d}, {:d}, {:d}\n".format(Wall.Type,Wall.LineType[0],Wall.LineType[1],Wall.LineType[2])
		fid.write(line)
		for k in range(len(Wall.Rwall)):
			fid.write("{:15.7e}, {:15.7e}\n".format(Wall.Rwall[k],Wall.Zwall[k]))

	fid.close()

	if(DEBUG > 0):	print("save_walls_file: Completed")
	
	return Walls