from tkinter 						import messagebox
from tkinter.messagebox 			import Message
import numpy as np
from matplotlib.path 				import Path
from routines.find_closest_segment	import find_closest_segment
from routines.globals				import DEBUG, DRAWING_WALL, EXTERNAL_PLASMA_WALL, INTERNAL_PLASMA_WALL, EXTERNAL_EIRENE_WALL, INTERNAL_EIRENE_WALL

def find_inplasma(Config, Zones, RKnots, ZKnots):

	if(DEBUG > 0):	print("find_inplasma")


	InPWallsKnots = []

#	Compute quantities In/Out of plasma walls

	nPwall = len(Config.iPwalls)
	for iPwall in range(nPwall):
		Wall = Config.Walls[Config.iPwalls[iPwall]]
		WallPath = Wall.WallPath

		PKnots = np.empty((len(RKnots),2), dtype='f8')
		PKnots[:,0] = RKnots
		PKnots[:,1] = ZKnots
		InWallKnots = np.where(WallPath.contains_points(PKnots), 1, 0)
		
		Radius = 1.e-4
		if(WallPath.contains_point([Wall.Rwall[0],Wall.Zwall[0]], radius= -Radius)): Radius = -Radius				#reverse radius if needed

		NearWallKnots = np.where((np.where(WallPath.contains_points(PKnots, radius= Radius), 1, 0) - \
								  np.where(WallPath.contains_points(PKnots, radius=-Radius), 1, 0)) > 0)[0]


		nOnWallKnots  = 0
		if(len(NearWallKnots) > 0):
			OnWallKnots	  = np.zeros(len(NearWallKnots), dtype = 'i4')
			for i in range(len(NearWallKnots)):
				iKnot = NearWallKnots[i]
				d, iSeg = find_closest_segment(RKnots[iKnot], ZKnots[iKnot], Wall.Rwall, Wall.Zwall)
				if(d < 1e-6):
					OnWallKnots[nOnWallKnots] = iKnot
					nOnWallKnots			 += 1
		else:
			print("\tATTENTION!! COULD BE SOMETHING WRONG IN SEARCHING NODES CLOSE TO THE WALL")
			print("\tATTENTION!! COULD BE SOMETHING WRONG IN SEARCHING NODES CLOSE TO THE WALL")
			print("\tATTENTION!! COULD BE SOMETHING WRONG IN SEARCHING NODES CLOSE TO THE WALL")

		if(nOnWallKnots > 0):
			OnWallKnots = OnWallKnots[:nOnWallKnots]
			InWallKnots[OnWallKnots] = 0							#at the begining set "on wall nodes" as were "out of wall nodes" to compute IsCrossed checking a sum
			print("\tNumber of point on wall  : ",len(OnWallKnots))

		if(Wall.Type == INTERNAL_PLASMA_WALL): InWallKnots = 1 - InWallKnots

		nZones = len(Zones)	

		for k in range(nZones):
			
			Nx = Zones[k].Nx
			Nz = Zones[k].Nz

			if(iPwall == 0):
				Zones[k].IsCrossed  = np.zeros((Nx, Nz, 4, nPwall), dtype='i4')
				Zones[k].InPlasmaRZ = np.zeros((Nx, Nz), dtype = 'i4')
				Zones[k].InPlasmaA  = np.zeros((Nx, Nz), dtype = 'i4')
				Zones[k].InPlasmaB  = np.zeros((Nx, Nz), dtype = 'i4')
				Zones[k].InPlasmaC  = np.zeros((Nx, Nz), dtype = 'i4')
				Zones[k].InPlasmaD  = np.zeros((Nx, Nz), dtype = 'i4')

			InPlasmaA = InWallKnots[Zones[k].KnotA]
			InPlasmaB = InWallKnots[Zones[k].KnotB]
			InPlasmaC = InWallKnots[Zones[k].KnotC]
			InPlasmaD = InWallKnots[Zones[k].KnotD]	

			Zones[k].InPlasmaA += InPlasmaA
			Zones[k].InPlasmaB += InPlasmaB
			Zones[k].InPlasmaC += InPlasmaC
			Zones[k].InPlasmaD += InPlasmaD

			Zones[k].IsCrossed[:,:,0,iPwall] = np.where((InPlasmaA + InPlasmaB) == 1, 1, 0) 
			Zones[k].IsCrossed[:,:,1,iPwall] = np.where((InPlasmaB + InPlasmaC) == 1, 1, 0) 
			Zones[k].IsCrossed[:,:,2,iPwall] = np.where((InPlasmaC + InPlasmaD) == 1, 1, 0) 
			Zones[k].IsCrossed[:,:,3,iPwall] = np.where((InPlasmaD + InPlasmaA) == 1, 1, 0) 

			RZ = np.array([Zones[k].gridRc.reshape(Zones[k].gridRc.size), Zones[k].gridZc.reshape(Zones[k].gridZc.size)]).T
			if(Wall.Type == EXTERNAL_PLASMA_WALL):   Zones[k].InPlasmaRZ += np.where(WallPath.contains_points(RZ), 1, 0).reshape(Zones[k].gridRc.shape)
			elif(Wall.Type == INTERNAL_PLASMA_WALL): Zones[k].InPlasmaRZ += np.where(WallPath.contains_points(RZ), 0, 1).reshape(Zones[k].gridRc.shape)

		if((Wall.Type == EXTERNAL_PLASMA_WALL) and (nOnWallKnots > 0)): InWallKnots[OnWallKnots] = 2				#now I can correctly set on wall nodes

		InPWallsKnots.append(InWallKnots)

	for k in range(nZones):
		Zones[k].InPlasma = np.where((Zones[k].InPlasmaA == 1) | (Zones[k].InPlasmaB == 1) | (Zones[k].InPlasmaC == 1) | (Zones[k].InPlasmaD == 1), 1, 0) 

	InWallKnots = InPWallsKnots[0]
	for iPwall in range(1,nPwall):
		InWallKnots = np.where(InWallKnots != 0, InWallKnots, InPWallsKnots[iPwall])

#	replace previous InPlasma values with walues (0,1,2)

	for k in range(nZones):
		Zones[k].InPlasmaA = InWallKnots[Zones[k].KnotA]
		Zones[k].InPlasmaB = InWallKnots[Zones[k].KnotB]
		Zones[k].InPlasmaC = InWallKnots[Zones[k].KnotC]
		Zones[k].InPlasmaD = InWallKnots[Zones[k].KnotD]

#	Check for overllaping wall

	for k in range(nZones):
		if(np.max(Zones[k].InPlasmaRZ) > 1):
			messagebox.showwarning("Find in Plasma", "ATTENTION: some wall is overlapping")

	return InWallKnots