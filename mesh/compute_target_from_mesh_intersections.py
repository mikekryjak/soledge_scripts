from math								import sqrt, asin
import numpy							as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot				as pyp

from routines.get_in_out_walls			import get_in_out_walls
from routines.intersect_contour			import intersect_2contours

from mesh.get_rho_in_out_core_sep		import get_rho_in_out_core_sep
from mesh.find_zones_intersections		import find_zones_intersections
from mesh.compute_mesh_intersections	import compute_mesh_intersections


def compute_target_from_mesh_intersections(Config, Eirene, RZLine, d_from_sep=0.01, east_west=0, l_pol=1):

	Zones	= Config.Zones

	WallTriangles = Eirene.WallTriangles

	Cut = find_zones_intersections(Config, RZLine)
	Lengths, IntRZ, IntCEll = compute_mesh_intersections(Config, Cut, also_pos=True, use_mag_zones=False)
	in_wall, out_wall = get_in_out_walls(Config, IntRZ[:,0], IntRZ[:,1])
	Lengths = Lengths[in_wall,:]
	IntRZ	= IntRZ[in_wall,:]
	IntCEll	= IntCEll[in_wall,:]

	dist = Lengths[:,0] - Lengths[0,0]
	Rho, In_Sep, Out_Sep, RZcore, RZsep = get_rho_in_out_core_sep(Config, IntRZ[:,0], IntRZ[:,1])
	
	if(len(In_Sep) > 0):
		if(Out_Sep[-1] < In_Sep[0]):
			Out_Sep = np.append(Out_Sep, In_Sep[0])
		else:
			Out_Sep = np.append(In_Sep[-1], Out_Sep)

	Ri, Zi, is1, is2  = intersect_2contours(RZsep[:,0], RZsep[:,1], IntRZ[:,0], IntRZ[:,1])
	if(len(Ri)==0):
		pyp.plot(RZsep[:,0], RZsep[:,1],'k-')
		pyp.plot(IntRZ[:,0], IntRZ[:,1],'r-')
		pyp.show()

	dsep = sqrt((IntRZ[0,0] - Ri[0])**2 + (IntRZ[0,1] - Zi[0])**2)
	dist -= dsep

	if(d_from_sep == 0): d_from_sep = 2*dist.max()
	OutSep	= np.where((dist > 0) & (dist < d_from_sep))[0]
	dist	= dist[OutSep]
	IntRZ	= IntRZ[OutSep,:]
	IntCEll = IntCEll[OutSep,:]

	SpCell = np.empty((IntCEll.shape[0],4), dtype = 'i4')
	SpPos  = np.empty((IntCEll.shape[0],3), dtype = 'f8')

	for iCell in range(IntCEll.shape[0]):											#find nearest S
		iZone 	= IntCEll[iCell,0]
		ix		= IntCEll[iCell,1]
		iTheta	= IntCEll[iCell,2]

#		Find zones along poloidal coordinate

		if(Zones[iZone].Chi[ix,-1] == 1):
			iThEast = np.array([np.min(np.where(Zones[iZone].Chi[ix,iTheta:] == 1)[0])+iTheta])
			East = -1
		else:
			iThEast = np.array([Zones[iZone].Chi.shape[1]])
			East = Zones[iZone].Neighbour.east
		
		nThetaPts = iThEast[0] - iTheta

		if(Zones[iZone].Chi[ix,0] == 1):
			iThWest = np.array([np.max(np.where(Zones[iZone].Chi[ix,:iTheta] == 1)[0])])
			West = -1
		else:
			iThWest = np.array([0])
			West = Zones[iZone].Neighbour.west

		iThetaOff  = iTheta - iThWest[0]
		nThetaPts  = iThEast[0] - iThWest[0]
		iZones	   = np.array([iZone])

#		Look East

		while (East > -1):
			iZones = np.append(iZones, East)
			iThWest = np.append(iThWest,0)
			if(Zones[East].Chi[ix,-1] == 1):
				iThEast = np.append(iThEast, np.min(np.where(Zones[East].Chi[ix,:] == 1)[0]))
				if(iThEast[-1] == 0):
					iZones  = iZones[:-1]
					iThWest = iThWest[:-1]
					iThEast = iThEast[:-1]
					East	= -2
				else: East = -1
			else:
				iThEast = np.append(iThEast, Zones[East].Chi.shape[1])
				East 	 = Zones[East].Neighbour.east

			if(East > -2):
				nThetaPts += iThEast[-1]

#		Look West

		while (West > -1):
			iZones = np.append(West, iZones)
			iThEast = np.append(Zones[West].Chi.shape[1], iThEast)
			if(Zones[West].Chi[ix,0] == 1):
				iThWest = np.append(np.max(np.where(Zones[West].Chi[ix,:] == 1)[0])+1, iThWest)
				if(iThWest[0] == Zones[West].Chi.shape[1]):
					iThWest = iThWest[1:]
					iZones  = iZones[1:]
					iThEast = iThEast[1:]
					West = -2
				else: West = -1
			else:
				iThWest = np.append(0, iThWest)
				West = Zones[West].Neighbour.west

			if(West > -2):
				nThetaPts += iThEast[0] - iThWest[0]
				iThetaOff += iThEast[0] - iThWest[0]

		Rpol  = np.empty((nThetaPts), dtype = 'f8')
		Zpol  = np.empty((nThetaPts), dtype = 'f8')
#		Bfact = np.empty((nThetaPts), dtype = 'f8')
		jOff = 0
		for k in range(len(iZones)):
			Rpol[jOff: jOff + iThEast[k] - iThWest[k]] = Zones[iZones[k]].gridRc[ix, iThWest[k]:iThEast[k]]
			Zpol[jOff: jOff + iThEast[k] - iThWest[k]] = Zones[iZones[k]].gridZc[ix, iThWest[k]:iThEast[k]]
#			Bfact[jOff: jOff + iThEast[k] - iThWest[k]] = np.sqrt(Zones[iZones[k]].Bphi[ix, iThWest[k]:iThEast[k]]**2 + Zones[iZones[k]].Br[ix, iThWest[k]:iThEast[k]]**2 + Zones[iZones[k]].Bz[ix, iThWest[k]:iThEast[k]]**2)/ \
#														  np.sqrt(Zones[iZones[k]].Br[ix, iThWest[k]:iThEast[k]]**2 + Zones[iZones[k]].Bz[ix, iThWest[k]:iThEast[k]]**2)

			jOff += iThEast[k] - iThWest[k]

		if(l_pol == 0):
#			Lpara1 = np.append(0., np.cumsum(np.sqrt((Rpol[1:]-Rpol[:-1])**2 + (Zpol[1:]-Zpol[:-1])**2)*(Bfact[1:]+Bfact[:-1])*0.5))

			dl = np.empty((nThetaPts), dtype = 'f8')
			jOff = 0
			for k in range(len(iZones)):
				dtheta = (Zones[iZones[k]].zb[ix, 1:] - Zones[iZones[k]].zb[ix, :-1])*2.*np.pi
				dlZone	   = 2.*dtheta/Gmet[iZones[k]][ix+1, 1:-1]
				
				dl[jOff: jOff + iThEast[k] - iThWest[k]] = dlZone[iThWest[k]:iThEast[k]]
				jOff += iThEast[k] - iThWest[k]
			dlZone = 0
			Lpara = np.cumsum(dl)
			dl	  = 0
		else:
			Lpara = np.append(0., np.cumsum(np.sqrt((Rpol[1:]-Rpol[:-1])**2 + (Zpol[1:]-Zpol[:-1])**2)))

		Lpara = Lpara - Lpara[iThetaOff]

		if(east_west == 0):
			if(Lpara[-1] >= -Lpara[0]): east_west = -1
			else:						east_west = 1

		if(east_west == -1):
			SpCell[iCell,0] = iZones[-1]						#SpCell[k_zone, i_grid, j_grid, ind_wall_tri]
			SpCell[iCell,1] = ix
			SpCell[iCell,2] = iThEast[-1]-1
			SpPos[iCell,2]	= -Lpara[-2]						#SpPos[R_target, Z_target,Line_lenght]
		else:
			SpCell[iCell,0] = iZones[0]
			SpCell[iCell,1] = ix
			SpCell[iCell,2] = iThWest[0]
			SpPos[iCell,2]	= Lpara[1]


		"""
		gridR = SpCell[iCell,0]].gridR
		gridZ = SpCell[iCell,0]].gridZ
		jsp = SpCell[iCell,2]
		R1	= gridR[ix, jsp]; R2 = gridR[ix, jsp+1]; R3 = gridR[ix+1, jsp+1]; R4 = gridR[ix+1, jsp]
		Z1	= gridZ[ix, jsp]; Z2 = gridZ[ix, jsp+1]; Z3 = gridZ[ix+1, jsp+1]; Z4 = gridZ[ix+1, jsp]
		DeltaRTheta = 0.5*((R2-R1) + (R3-R4))
		DeltaZTheta = 0.5*((Z2-Z1) + (Z3-Z4))
		DeltaLTheta = sqrt(DeltaRTheta**2 + DeltaRTheta**2)
		DeltaRx		= 0.5*((R4-R1) + (R3-R2))
		DeltaZx		= 0.5*((Z4-Z1) + (Z3-Z2))
		DeltaLPhi	= sqrt(DeltaRx**2 + DeltaRx**2)
		DeltaLPhi	= (DeltaRx*DeltaRTheta + DeltaZx*DeltaZTheta)/DeltaLTheta
		"""

		SpPos[iCell,0]	= Zones[SpCell[iCell,0]].gridRc[ix, SpCell[iCell,2]]
		SpPos[iCell,1]	= Zones[SpCell[iCell,0]].gridZc[ix, SpCell[iCell,2]]

		iZoneEnd	= SpCell[iCell,0]
		iThetaEnd	= SpCell[iCell,2]
		for iOff in range(5):
			if(iOff == 0):
				aa = np.where((WallTriangles.k == iZoneEnd) & (WallTriangles.i == ix) & (WallTriangles.j == iThetaEnd))[0]
			else:
				if(iThetaEnd+iOff < Zones[iZoneEnd].gridZc.shape[1]):
					aa = np.where((WallTriangles.k == iZoneEnd) & (WallTriangles.i == ix) & (WallTriangles.j == iThetaEnd+iOff))[0]
				else:
					aa = np.where((WallTriangles.k == Zones[iZoneEnd].Neighbour.east) & (WallTriangles.i == ix) & (WallTriangles.j == iThetaEnd+iOff-Zones[iZoneEnd].gridZc.shape[1]))[0]
				if(len(aa) > 0): break

				if(iThetaEnd-iOff >= 0):
					aa = np.where((WallTriangles.k == iZoneEnd) & (WallTriangles.i == ix) & (WallTriangles.j == iThetaEnd-iOff))[0]
				else:
					aa = np.where((WallTriangles.k == Zones[iZoneEnd].Neighbour.west) & (WallTriangles.i == ix) & (WallTriangles.j == Zones[Zones[iZoneEnd].Neighbour.west].gridZc.shape[1]+iThetaEnd-iOff))[0]
				
			if(len(aa) > 0): break

		if(len(aa) > 1):
			print("Attention more than one wall triangle in soledge quadrangle: triangles=",aa+1)

		SpCell[iCell,3] = aa[0]

	return dist, SpCell, SpPos
