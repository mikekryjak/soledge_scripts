import numpy as np
from math import floor
from routines.globals		import DEBUG, EXTERNAL_PLASMA_WALL, INTERNAL_PLASMA_WALL

def clean_chis(Config):

	if(DEBUG > 0):	print("clean_chis")

	Zones 		= Config.Zones
	MagZones	= Config.MagZones
	X_points	= Config.X_points
	nZones		= len(Zones)

	if(nZones != len(MagZones)):
		print("\tThere is an error in grid-gen!")
		print("\t\tclean_chis must be called before zones plitting!")
		exit()

	nPwall = len(Config.iPwalls)
	for iPwall in range(nPwall):
		Wall = Config.Walls[Config.iPwalls[iPwall]]
		WallPath = Wall.WallPath
		if(Wall.Type == EXTERNAL_PLASMA_WALL):							#Outer wall
			for k in range(nZones):
				if((MagZones[k].pA.coord[2]  == -1) and (MagZones[k].pA.coord[3] == 1)):
					if(not WallPath.contains_point([X_points[MagZones[k].pA.coord[1]].R,X_points[MagZones[k].pA.coord[1]].Z])):		Zones[k].Chi[0,0] = 1

				if((MagZones[k].pB.coord[2] == -1) and (MagZones[k].pB.coord[3] == 1)):
					if(not WallPath.contains_point([X_points[MagZones[k].pB.coord[1]].R,X_points[MagZones[k].pB.coord[1]].Z])):		Zones[k].Chi[-1,0] = 1

				if((MagZones[k].pC.coord[2] == -1) and (MagZones[k].pC.coord[3] == 1)):
					if(not WallPath.contains_point([X_points[MagZones[k].pC.coord[1]].R,X_points[MagZones[k].pC.coord[1]].Z])):	Zones[k].Chi[-1,-1] = 1

				if((MagZones[k].pD.coord[2] == -1) and (MagZones[k].pD.coord[3] == 1)):
					if(not WallPath.contains_point([X_points[MagZones[k].pD.coord[1]].R,X_points[MagZones[k].pD.coord[1]].Z])):		Zones[k].Chi[0,-1] = 1

		elif(Wall.Type == INTERNAL_PLASMA_WALL):							#Inner component
			for k in range(nZones):
				if((MagZones[k].pA.coord[2]  == -1) and (MagZones[k].pA.coord[3] == 1)):
					if(WallPath.contains_point([X_points[MagZones[k].pA.coord[1]].R,X_points[MagZones[k].pA.coord[1]].Z])):		Zones[k].Chi[0,0] = 1

				if((MagZones[k].pB.coord[2] == -1) and (MagZones[k].pB.coord[3] == 1)):
					if(WallPath.contains_point([X_points[MagZones[k].pB.coord[1]].R,X_points[MagZones[k].pB.coord[1]].Z])):		Zones[k].Chi[-1,0] = 1

				if((MagZones[k].pC.coord[2] == -1) and (MagZones[k].pC.coord[3] == 1)):
					if(WallPath.contains_point([X_points[MagZones[k].pC.coord[1]].R,X_points[MagZones[k].pC.coord[1]].Z])):	Zones[k].Chi[-1,-1] = 1

				if((MagZones[k].pD.coord[2] == -1) and (MagZones[k].pD.coord[3] == 1)):
					if(WallPath.contains_point([X_points[MagZones[k].pD.coord[1]].R,X_points[MagZones[k].pD.coord[1]].Z])):		Zones[k].Chi[0,-1] = 1

	for k in range(nZones):
		Chi  = Zones[k].Chi
		m	 = Chi.shape[0]
		p	 = Chi.shape[1]

#	Scan in poloidal direction to avoid less than 3 points with chi=0

		West = Zones[k].Neighbour.west
		East = Zones[k].Neighbour.east
		for i in range(m):
			l = 0

#			Consider last two point at west

			if(West > -1):
				if(Zones[West].Chi[i,-1] == 0):
					l +=1
					if(Zones[West].Chi[i,-2] == 0): l +=1

			for j in range(p):
				if(Chi[i,j] == 1):
					if(l  <  3):
						for ll in range(1,l+1): Chi[i,j-ll] = 1
					l = 0
				elif(Chi[i,j] == 0):
					l += 1

#			Check first two point at east

			if(East > -1):
				if(Zones[East].Chi[i,0] == 1) :
					if(l  <  3):
						for ll in range(1,l+1): Chi[i,-ll] = 1
					l = 0
				elif(Zones[East].Chi[i,0] == 0):
					l += 1
	
				if((Zones[East].Chi[i,1] == 1)  and (l < 3)): 
					Chi[i,-1] = 1
			else:
				if((Chi[i,-1] == 0) and (Chi[i,-2] == 1)): 
					Chi[i,-1] = 1
					l = 0
		
		for i in range(m):
			for j in range(1,p-1):
				if((Chi[i,j-1] == 0) and (Chi[i,j] == 1) and (Chi[i,j+1] == 0)):	   Chi[i,j] = 0

			for j in range(1,p-2):
				if(( Chi[i,j-1] == 0) and (Chi[i,j] == 1) and (Chi[i,j+1] == 1) and (Chi[i,j+2] == 0)):
					Chi[i,j]   = 0
					Chi[i,j+1] = 0

#		Remove single points with chi=1 in poloidal direction

		for i in range(m):
			for j in range(1,p-1):
				if((Chi[i,j-1] == 0) and (Chi[i,j] == 1) and (Chi[i,j+1] == 0)):	   Chi[i,j] = 0

			for j in range(1,p-2):
				if(( Chi[i,j-1] == 0) and (Chi[i,j] == 1) and (Chi[i,j+1] == 1) and (Chi[i,j+2] == 0)):
					Chi[i,j]   = 0
					Chi[i,j+1] = 0

#	Scan in radial direction to avoid less than 3 points with chi=0

		South = Zones[k].Neighbour.south
		North = Zones[k].Neighbour.north
		for j in range(p):
			l = 0

#			Consider last two point at south

			if(South > -1):
				if(Zones[South].Chi[-1,j] == 0):
					l +=1
					if(Zones[South].Chi[-2,j] == 0): l +=1

			for i in range(m):
				if(Chi[i,j] == 1):
					if(l  <  3):
						for ll in range(1,l+1): Chi[i-ll,j] = 1
					l = 0
				elif(Chi[i,j] == 0):
					l += 1

#			Check first two point at north

			if(North > -1):
				if(Zones[North].Chi[0,j] == 1) :
					if(l  <  3):
						for ll in range(1,l+1): Chi[-ll,j] = 1
					l = 0
				elif(Zones[North].Chi[0,j] == 0):
					l += 1
	
				if((Zones[North].Chi[1,j] == 1)  and (l < 3)): 
					Chi[-1,j] = 1
			else:
				if((Chi[-1,j] == 0) and (Chi[-2,j] == 1)): 
					Chi[-1,j] = 1
					l = 0

#		Remove single points with chi=1 in radial direction

		for j in range(p):
			for i in range(1,m-1):
				if((Chi[i-1,j] == 0) and (Chi[i,j] == 1) and (Chi[i+1,j] == 0)):	   Chi[i,j] = 0

			for i in range(1,m-2):
				if(( Chi[i-1,j] == 0) and (Chi[i,j] == 1) and (Chi[i+1,j] == 1) and (Chi[i+2,j] == 0)):
					Chi[i,j]   = 0
					Chi[i+1,j] = 0

		Zones[k].Chi = np.copy(Chi)


	for k in range(len(Zones)):
		Chi 			= Zones[k].Chi
		Chi_			= np.zeros((Chi.shape[0]+2,Chi.shape[1]+2), dtype='f8')
		Chi_[1:-1,1:-1]	= Chi
		North			= Zones[k].Neighbour.north
		South			= Zones[k].Neighbour.south
		East			= Zones[k].Neighbour.east
		West			= Zones[k].Neighbour.west
		if(West > -1):
			Chi_[1:-1,0] = Zones[West].Chi[:,-1]
		else:
			Chi_[1:-1,0] = Chi_[1:-1,1]

		if(East > -1):
			Chi_[1:-1,-1] = Zones[East].Chi[:,0]
		else:
			Chi_[1:-1,-1] = Chi_[1:-1,-2]

		if(North > -1):
			Chi_[-1,1:-1] = Zones[North].Chi[0,:]
		else:
			Chi_[-1,1:-1] = Chi_[-2,1:-1]

		if(South > -1):
			Chi_[0,1:-1] = Zones[South].Chi[-1,:]
		else:
			Chi_[0,1:-1] = Chi_[1,1:-1]

		ii, jj = np.where(Chi == 0)
		for l in range(len(ii)):
			if(np.sum(Chi_[ii[l]:ii[l]+3,jj[l]:jj[l]+3]) >= 7):
				Chi[ii[l],jj[l]] = 1
			elif(np.sum(Chi_[ii[l]:ii[l]+3,jj[l]:jj[l]+3]) >= 6):
				if(((Chi_[ii[l], jj[l]+1] == 0) and (Chi_[ii[l]+2, jj[l]+1]==0)) or
				   ((Chi_[ii[l]+1, jj[l]] == 0) and (Chi_[ii[l]+1, jj[l]+2]==0))): Chi[ii[l], jj[l]] = 1

#		update Chi_
		Chi_[1:-1,1:-1] = Chi

#		remove spikes of 1 surrounded by zeros
		ii, jj = np.where((Chi_[1:-1,1:-1] == 1) & (Chi_[1:-1,:-2] == 0) & (Chi_[1:-1,2:] == 0))		#poloidal direction
		Chi[ii,jj] = 0

		ii, jj = np.where((Chi_[1:-1,1:-1] == 1) & (Chi_[:-2,1:-1] == 0) & (Chi_[2:,1:-1] == 0))		#radial direction
		Chi[ii,jj] = 0

#		update Chi_
		Chi_[1:-1,1:-1] = Chi

#		remove hole of zero
		ii, jj = np.where((Chi_[1:-1,1:-1] == 0) & (Chi_[1:-1,:-2] == 1) & (Chi_[1:-1,2:] == 1))		#polidal direction
		Chi[ii,jj] = 1

		ii, jj = np.where((Chi_[1:-1,1:-1] == 0) & (Chi_[:-2,1:-1] == 1) & (Chi_[2:,1:-1] == 1))		#radial direction
		Chi[ii,jj] = 1

#		update Chi_
		Chi_[1:-1,1:-1] = Chi
		
		Zones[k].Chi = np.copy(Chi)

	if(DEBUG > 0):	print("clean_chis: Completed")
