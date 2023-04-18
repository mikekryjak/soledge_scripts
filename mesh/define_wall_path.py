from matplotlib.path	import Path
import numpy			as np
from routines.globals	import *

def non_necessaria_define_wall_path(Config):

#	define boundary path

	Config.iPwalls = []
	Config.iEwalls = []
	for iWall in range(len(Config.Walls)):
		DefShape = False
		Wall = Config.Walls[iWall]
		if((Wall.Type == EXTERNAL_PLASMA_WALL) or (Wall.Type == INTERNAL_PLASMA_WALL)):
			Config.iPwalls.append(iWall)
			DefShape = True
		if((Wall.Type == EXTERNAL_EIRENE_WALL) or (Wall.Type == INTERNAL_EIRENE_WALL)):
			Config.iEwalls.append(iWall)
			DefShape = True

		if(DefShape):
			nWall = Wall.Rwall.shape[0]
			if((Wall.Rwall[0] == Wall.Rwall[-1]) and (Wall.Zwall[0] == Wall.Zwall[-1])):
				nWallPath = nWall
			else:
				nWallPath = nWall + 1

			WallVerts = np.empty((nWallPath,2),dtype='f8')

			if((Wall.Rwall[0] == Wall.Rwall[-1]) and (Wall.Zwall[0] == Wall.Zwall[-1])):
				WallVerts[:,0] = Wall.Rwall
				WallVerts[:,1] = Wall.Zwall	
			else:
				WallVerts[:-1,0] = Wall.Rwall
				WallVerts[:-1,1] = Wall.Zwall
				WallVerts[-1,:]  = WallVerts[0,:]

			Wall.WallPath = Path(WallVerts, closed=True)

		if((Wall.Rwall[0] == Wall.Rwall[-1]) and (Wall.Zwall[0] == Wall.Zwall[-1])): Config.Walls[iWall].Closed = 1
		else:																		 Config.Walls[iWall].Closed = -1

		Area = np.sum((Wall.Zwall[1:]+Wall.Zwall[:-1])*(Wall.Rwall[1:]-Wall.Rwall[0:-1]))
		if(Area > 0.):	Config.Walls[iWall].Clockwise = 1
		else:			Config.Walls[iWall].Clockwise = -1
