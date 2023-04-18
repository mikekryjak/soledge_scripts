import sys
import types
from math					import sqrt, atan2
from math					import pi as PI
TWOPI = 2*PI

import numpy as np
from matplotlib.path		import Path
from scipy.spatial			import Delaunay

from eirene.direct_or_indirect		import triangles_direct_or_indirect
from eirene.triangles_routines		import triangles_reshape
from eirene.get_ext_plasma_triseq	import get_ext_plasma_triseq

from routines.globals				import *
from routines.find_closest_point	import get_dist_from_segment
from routines.utils_routines		import min_and_arg, get_clockwise
from routines.intersect_segments	import intersect_segments


def generate_vacuum_triangles(Root, Config, EireneWall, Triangles, WallTriangles, RKnotsP, ZKnotsP):

	print("generate_vacuum_triangles")

	CWalls		 = Config.Walls
	iIntEWalls	 = np.array(Config.iIntEWalls)
	EWalls		 = EireneWall.EWalls
	TriSequences = EireneWall.TriSequences
	TriFaces	 = EireneWall.TriFaces
	R12			 = EireneWall.R12
	Z12			 = EireneWall.Z12

	iTriPlasma = get_ext_plasma_triseq(Config, EireneWall)

	RwallEP = R12[iTriPlasma][:,0]												#Coordinate Eirene external plasma wall
	ZwallEP = Z12[iTriPlasma][:,0]
	ClockwiseEP = get_clockwise(RwallEP, ZwallEP) 

	nTri   = np.copy(WallTriangles.ntri[TriSequences[iTriPlasma]])
	TriKnots	= np.array([Triangles.p1,Triangles.p2,Triangles.p3]).T
	nwallEP 	= TriKnots[nTri,TriFaces[iTriPlasma]]							#Get number of external plasma wall knots
	DirectEP	= triangles_direct_or_indirect(TriKnots, RKnotsP, ZKnotsP)		#get point orientation:  clockwise/counterclockwise
	del nTri, TriKnots

	
#	define one point for each internal EIRENE walls
#	==========================

	if(len(iIntEWalls) > 0):
		RIntEWalls = np.empty(len(iIntEWalls), dtype='f8')
		ZIntEWalls = np.empty(len(iIntEWalls), dtype='f8')
		for iIntEWall in range(len(iIntEWalls)):
			RIntEWalls[iIntEWall] = CWalls[iIntEWalls[iIntEWall]].Rwall[0]
			ZIntEWalls[iIntEWall] = CWalls[iIntEWalls[iIntEWall]].Zwall[0]

		RZIntEWalls	= np.array([RIntEWalls,ZIntEWalls]).T
		del RIntEWalls, ZIntEWalls

#	Sort SubEzones from the smaller dl_knots to the largests
#	=================================

	for iEWall in range(len(EWalls)):							
		EWall		= EWalls[iEWall]
		ESubZones	= EWall.ESubZones

		dl_knots = np.empty(len(ESubZones), dtype='f8')
		for iESubZone in range(len(ESubZones)): dl_knots[iESubZone] = np.sqrt(ESubZones[iESubZone].dr_knots**2 + ESubZones[iESubZone].dz_knots**2)
		EWall.dlKnotsSort = np.argsort(dl_knots)

#	Define  SubZones border
#	===============

	for iEWall in range(len(EWalls)):
		EWall		  = EWalls[iEWall]
		CWall		  = EWall.Wall
		EWall.RWbord  = np.copy(CWall.Rwall)
		EWall.ZWbord  = np.copy(CWall.Zwall)
		EWall.RKnotsE = np.empty(0, dtype='f8')									#Array to store eirene wall knots where external plasma wall is is attached to external eirene wall
		EWall.ZKnotsE = np.empty(0, dtype='f8')
		EWall.nKnotsE = np.empty(0, dtype='i4')
		for iESubZone in range(len(EWall.ESubZones)):
			ESubZone = EWall.ESubZones[iESubZone]
			ESubZone.RKbord  = np.copy(ESubZone.Rbord)								#knots SubZones boundaries
			ESubZone.ZKbord  = np.copy(ESubZone.Zbord)
			ESubZone.IsKnot  = np.ones(len(ESubZone.Rbord), dtype='i4')

#	Find ESubZones common borders
#	====================

	d2Max=1.e-8
	for iEWall1 in range(-1,len(EWalls)):
		if(iEWall1 == -1):														#Special treatment for wirst wall wich is External plasma wall
			EWall1 = CWalls[Config.iExtWalls[0]]	
			iESubZones1 = [0]
		else:
			EWall1 = EWalls[iEWall1]	
			iESubZones1 = np.arange(len(EWall1.dlKnotsSort)-1,dtype='i4')

		for iESubZone1 in iESubZones1:
			if(iEWall1 == -1):													#Special treatment for wirst wall wich is External plasma wall
				RKbord1 = EWall1.Rwall
				ZKbord1 = EWall1.Zwall
			else:			
				ESubZone1 = EWall1.ESubZones[EWall1.dlKnotsSort[iESubZone1]]
				RKbord1  = ESubZone1.RKbord
				ZKbord1  = ESubZone1.ZKbord

			Clockwise1 = get_clockwise(RKbord1, ZKbord1) 
			Len1  	   = len(RKbord1)

			if(iEWall1 == -1): iEWalls2 = np.arange(len(EWalls))
			else:			   iEWalls2 = [iEWall1]
			for iEWall2 in iEWalls2:
				EWall2 = EWalls[iEWall2]

				if(iEWall1 != iEWall2):
					iESubZones2 = np.arange(len(EWall2.ESubZones))
					ClockwiseW2 = get_clockwise(EWall2.RWbord, EWall2.ZWbord) 
				else:
					iESubZones2 = EWall1.dlKnotsSort[iESubZone1+1:]

				for iESubZone2 in iESubZones2:
					ESubZone2  = EWall2.ESubZones[iESubZone2]
					RKbord2	   = ESubZone2.RKbord
					ZKbord2	   = ESubZone2.ZKbord
					IsKnot2	   = ESubZone2.IsKnot
					Clockwise2 = get_clockwise(RKbord2, ZKbord2) 

					Len2   = len(RKbord2)
					iPt1   = 0
					while (iPt1 < Len1-1):																#for closed wall last point is  equal to the first one
						iPt2 = np.argmin((RKbord2-RKbord1[iPt1])**2 + (ZKbord2-ZKbord1[iPt1])**2)
						if((RKbord2[iPt2]-RKbord1[iPt1])**2 + (ZKbord2[iPt2]-ZKbord1[iPt1])**2 < d2Max):
							forward2 = 0
							if(Clockwise1 != Clockwise2):
								if(((RKbord2[(iPt2+1)%Len2]-RKbord1[(iPt1+1)%Len1])**2 + (ZKbord2[(iPt2+1)%Len2]-ZKbord1[(iPt1+1)%Len1])**2) < d2Max): forward2 = +1
							else:
								if(iPt2 > 0):
									if(((RKbord2[iPt2-1]-RKbord1[(iPt1+1)%Len1])**2 + (ZKbord2[iPt2-1]-ZKbord1[(iPt1+1)%Len1])**2) < d2Max): forward2 = -1
								else:
									if(((RKbord2[-2]-RKbord1[(iPt1+1)%Len1])**2 + (ZKbord2[-2]-ZKbord1[(iPt1+1)%Len1])**2) < d2Max): forward2 = -1			#for closed wall last point is  equal to the first one

							if(forward2 != 0):
								iPt1s = iPt1
								iPt2s = iPt2
								iPt1  = (iPt1 + 2)
								if((forward2 < 0) and (iPt2 == 0)):	iPt2  = -3
								else:								iPt2  = (iPt2 + 2*forward2)
								
								iPt1c = iPt1%Len1	
								iPt2c = iPt2%Len2				
								if(iPt1 < 2*Len1):
									while (((RKbord2[iPt2c]-RKbord1[iPt1c])**2 + (ZKbord2[iPt2c]-ZKbord1[iPt1c])**2) < d2Max):
										iPt1 += 1
										iPt2 += forward2
										iPt1c = iPt1%Len1	
										iPt2c = iPt2%Len2				
										
								iPt1e = (iPt1 - 1)%Len1
								iPt2e = (iPt2 - forward2)%Len2
								if((iPt1s > 0) or (((RKbord2[(iPt2s -forward2)%Len2]-RKbord1[iPt1s-1])**2 + (ZKbord2[(iPt2s -forward2)%Len2]-ZKbord1[iPt1s-1])**2) >= d2Max)): #skip  if between end and begin, it will be found later				

									if(iPt2s < iPt2e):	iPt2l, iPt2h = iPt2s, iPt2e
									else:				iPt2l, iPt2h = iPt2e, iPt2s
									if(iEWall1 == -1):																	#Now for external plasma wall fix eirene knots

#										Fix SubZone border
#										===========

										diPtEs, iPtEs = min_and_arg((RwallEP-RKbord1[iPt1s])**2 + (ZwallEP-ZKbord1[iPt1s])**2)
										diPtEe, iPtEe = min_and_arg((RwallEP-RKbord1[iPt1e])**2 + (ZwallEP-ZKbord1[iPt1e])**2)
										
										if(iPtEs < iPtEe):	iPtEl, iPtEh = iPtEs,iPtEe
										else:				iPtEl, iPtEh = iPtEe,iPtEs
										if(((ClockwiseEP == Clockwise1) and (iPtEs < iPtEe)) or ((ClockwiseEP != Clockwise1) and (iPtEs > iPtEe))):
											RKnotsE = RwallEP[iPtEl:iPtEh+1]
											ZKnotsE = ZwallEP[iPtEl:iPtEh+1]
											nKnotsE = nwallEP[iPtEl:iPtEh+1]
										else:
											RKnotsE	 = np.append(RwallEP[iPtEh:],RwallEP[:iPtEl+1])
											ZKnotsE	 = np.append(ZwallEP[iPtEh:],ZwallEP[:iPtEl+1])
											nKnotsE	 = np.append(nwallEP[iPtEh:],nwallEP[:iPtEl+1])

										if(ClockwiseEP == Clockwise2):
											RKnotsE = RKnotsE[::-1]
											ZKnotsE = ZKnotsE[::-1]
											nKnotsE = nKnotsE[::-1]

										EWall2.RKnotsE = np.append(EWall2.RKnotsE, RKnotsE)
										EWall2.ZKnotsE = np.append(EWall2.ZKnotsE, ZKnotsE)
										EWall2.nKnotsE = np.append(EWall2.nKnotsE, nKnotsE)

										IsKnotE		   = np.zeros(len(RKnotsE), dtype='i4')

										if(iPt2s < iPt2e):	diPt2l, diPt2h = diPtEs, diPtEe
										else:				diPt2l, diPt2h = diPtEe, diPtEs
										LenSeg2l = sqrt((RKnotsE[1] -RKnotsE[0])**2 +(ZKnotsE[1] -ZKnotsE[0])**2)
										LenSeg2h = sqrt((RKnotsE[-1]-RKnotsE[-1])**2+(ZKnotsE[-1]-ZKnotsE[-2])**2)
										diPt2l = get_dist_from_segment(RKnotsE[0],  ZKnotsE[0],  RKnotsE[1],  ZKnotsE[1],  RKbord2[iPt2l], ZKbord2[iPt2l])
										diPt2h = get_dist_from_segment(RKnotsE[-2], ZKnotsE[-2], RKnotsE[-1], ZKnotsE[-1], RKbord2[iPt2h], ZKbord2[iPt2h])
										if(((Clockwise2 != Clockwise1) and (iPt2s < iPt2e)) or ((Clockwise2 == Clockwise1) and (iPt2s > iPt2e))):
											if(diPt2l < LenSeg2l/2): iPt2l -= 1
											if(diPt2h < LenSeg2h/2): iPt2h += 1
#											if(sqrt(diPt2l) < 1e-3): iPt2l -= 1
#											if(sqrt(diPt2h) < 1e-3): iPt2h += 1
											RPt2l, ZPt2l   = RKbord2[iPt2l], ZKbord2[iPt2l]
											RPt2h, ZPt2h   = RKbord2[iPt2h], ZKbord2[iPt2h]
											RKbord2 = np.concatenate((RKbord2[:iPt2l+1], RKnotsE, RKbord2[iPt2h:]))		#border inside array
											ZKbord2 = np.concatenate((ZKbord2[:iPt2l+1], ZKnotsE, ZKbord2[iPt2h:])) 
											IsKnot2 = np.concatenate((IsKnot2[:iPt2l+1], IsKnotE, IsKnot2[iPt2h:]))
										else:								
											if(diPt2l < LenSeg2l/2): iPt2l += 1
											if(diPt2h < LenSeg2h/2): iPt2h -= 1
#											if(sqrt(diPt2l) < 1e-3): iPt2l += 1
#											if(sqrt(diPt2h) < 1e-3): iPt2h -= 1
											RPt2l, ZPt2l   = RKbord2[iPt2l], ZKbord2[iPt2l]
											RPt2h, ZPt2h   = RKbord2[iPt2h], ZKbord2[iPt2h]
											RKbord2 = np.concatenate((RKbord2[iPt2l:iPt2h+1], RKnotsE))							#border crossing array end
											ZKbord2 = np.concatenate((ZKbord2[iPt2l:iPt2h+1], ZKnotsE)) 
											IsKnot2 = np.concatenate((IsKnot2[iPt2l:iPt2h+1], IsKnotE))

										ESubZone2.RKbord = RKbord2 
										ESubZone2.ZKbord = ZKbord2
										ESubZone2.IsKnot = IsKnot2

#										Now fix wall border
#										===========

										RWbord = EWall2.RWbord
										ZWbord = EWall2.ZWbord

										iPtWl = np.argmin((RWbord-RPt2l)**2 + (ZWbord-ZPt2l)**2)
										iPtWh = np.argmin((RWbord-RPt2h)**2 + (ZWbord-ZPt2h)**2)
										if(iPt2s < iPt2e):	iPtWs, iPtWe = iPtWl, iPtWh
										else:				iPtWs, iPtWe = iPtWh, iPtWl

										if(iPtWs < iPtWe): iPtWl, iPtWh = iPtWs, iPtWe
										else:			   iPtWl, iPtWh = iPtWe, iPtWs
										if(((ClockwiseEP != ClockwiseW2) and (iPtWs < iPtWe)) or ((ClockwiseEP == ClockwiseW2) and (iPtWs > iPtWe))):
											RWbord = np.concatenate((RWbord[:iPtWl+1], RKnotsE, RWbord[iPtWh:]))		#border inside array
											ZWbord = np.concatenate((ZWbord[:iPtWl+1], ZKnotsE, ZWbord[iPtWh:]))
										else:								
											RWbord = np.concatenate((RWbord[iPtWl:iPtWh+1], RKnotsE, np.array([RWbord[iPtWl]])))	#border crossing array end
											ZWbord = np.concatenate((ZWbord[iPtWl:iPtWh+1], ZKnotsE, np.array([ZWbord[iPtWl]])))	#Last point to close border

										EWall2.RWbord = np.copy(RWbord)
										EWall2.ZWbord = np.copy(ZWbord)
									else:
										if(((forward2 == 1) and (iPt2s < iPt2e)) or ((forward2 == -1) and (iPt2s > iPt2e))):
											IsKnot2[iPt2l:iPt2h+1] = 0								#Set nearest zone  node to false
										else:
											IsKnot2[:iPt2l+1] = 0									#border crossing array end
											IsKnot2[iPt2h:] = 0
						iPt1 +=1

#	Build triangles
#	====================

	for iEWall in range(len(EWalls)):
		RKnotExtB		= np.array([], dtype='f8')				#Knots external wall border
		ZKnotExtB	 	= np.array([], dtype='f8')
		RKnotIntB	 	= np.array([], dtype='f8')				#Knots internal walls border
		ZKnotIntB	 	= np.array([], dtype='f8')
		nKnotIntBStart	= np.array([], dtype='i4')				#number first knot internal wall border
		RKnotIn	 	= np.array([], dtype='f8')					#knots internal zones
		ZKnotIn	 	= np.array([], dtype='f8')
		EWall	  	= EWalls[iEWall]
		ESubZones 	= EWall.ESubZones
		AllIntEWalls = np.array([], dtype='i4')
		for iESubZone in EWall.dlKnotsSort:
			ESubZone = ESubZones[iESubZone]
			dr_knots = ESubZone.dr_knots
			dz_knots = ESubZone.dz_knots
			RKbord	 = ESubZone.RKbord
			ZKbord	 = ESubZone.ZKbord
			IsKnot	 = ESubZone.IsKnot

			Area	= np.sum((ZKbord[1:]+ZKbord[:-1])*(RKbord[1:]-RKbord[:-1]))+(ZKbord[0]+ZKbord[-1])*(RKbord[0]-RKbord[-1])
			if(Area > 0.):
				RKbord = RKbord[::-1]			
				ZKbord = ZKbord[::-1]	
			SubPath = Path(np.array([RKbord,ZKbord]).T, closed=True)

#			define knots on zone borders
#			==================

			LKbord	 = len(ZKbord)-1										#border are closed and so last point=start point
			for i1 in range(LKbord):										#define border knots 
				if(IsKnot[i1] > 0):
					RKnotExtB = np.append(RKnotExtB,RKbord[i1])
					ZKnotExtB = np.append(ZKnotExtB,ZKbord[i1])
				i2 = i1+1
				if((IsKnot[i1] > 0) or (IsKnot[i2] > 0)):
					npt = max(abs((RKbord[i2]-RKbord[i1])/dr_knots), abs((ZKbord[i2]-ZKbord[i1])/dz_knots))
					if(npt > 1.):
						npt = int(npt+1)
						dR  = (RKbord[i2] - RKbord[i1])/npt
						dZ  = (ZKbord[i2] - ZKbord[i1])/npt
						for k in range(1,npt):
							RKnotExtB = np.append(RKnotExtB, RKbord[i1] +  dR*k)
							ZKnotExtB = np.append(ZKnotExtB, ZKbord[i1] +  dZ*k)

#			define knots for internal EIRENE wall
#			=====================

			if(len(iIntEWalls) > 0):
				IsIn = np.where(SubPath.contains_points(RZIntEWalls))[0]					#find internal EIRENE wall inside sub-zone
				if(len(IsIn) > 0):
					InIntEWalls = iIntEWalls[IsIn]
					for iInEWall in InIntEWalls:
						nKnotIntBStart = np.append(nKnotIntBStart, len(RKnotIntB))
						RIntbord = CWalls[iInEWall].Rwall
						ZIntbord = CWalls[iInEWall].Zwall
						for i1 in range(len(RIntbord)-1):									#define border knots for internal EIRENE walls
							RKnotIntB = np.append(RKnotIntB,RIntbord[i1])
							ZKnotIntB = np.append(ZKnotIntB,ZIntbord[i1])
							i2 = i1+1
							npt = max(abs((RIntbord[i2]-RIntbord[i1])/dr_knots), abs((ZIntbord[i2]-ZIntbord[i1])/dz_knots))
							if(npt > 1.):
								npt = int(npt+1)
								dR  = (RIntbord[i2] - RIntbord[i1])/npt
								dZ  = (ZIntbord[i2] - ZIntbord[i1])/npt
								for k in range(1,npt):
									RKnotIntB = np.append(RKnotIntB, RIntbord[i1] +  dR*k)
									ZKnotIntB = np.append(ZKnotIntB, ZIntbord[i1] +  dZ*k)
			
					nKnotIntBStart = np.append(nKnotIntBStart, len(RKnotIntB))
			else:	InIntEWalls = np.array([], dtype='i4')
			AllIntEWalls = np.append(AllIntEWalls, InIntEWalls)

			Rmin	 = np.min(RKbord)
			Rmax	 = np.max(RKbord)
			Zmin	 = np.min(ZKbord)
			Zmax	 = np.max(ZKbord)
			nPtR	 = int((Rmax - Rmin)/dr_knots) + 1
			nPtZ	 = int((Zmax - Zmin)/dz_knots) + 1
			dR 		 = (Rmax - Rmin)/nPtR
			dZ		 = (Zmax - Zmin)/nPtZ

			Rknots1D = np.linspace(Rmin+dR, Rmax-dR, nPtR-1)
			Zknots1D = np.linspace(Zmin+dZ, Zmax-dZ, nPtZ-1)
			Rknots2D,Zknots2D = np.meshgrid(Rknots1D, Zknots1D)			#coordinate matrix knots
			Rknots2D = Rknots2D.reshape(-1)
			Zknots2D = Zknots2D.reshape(-1)

#			Remove matrix knots inside  internal EIRENE walls
#			==============================

			dMin = min(dr_knots, dz_knots)*0.9
			IsIn = np.where(SubPath.contains_points(np.array([Rknots2D,Zknots2D]).T, radius=-dMin))[0]				#knots inside sub-zone border
			Rknots2D = Rknots2D[IsIn]
			Zknots2D = Zknots2D[IsIn]
		
			if(len(InIntEWalls) > 0):
				for iInEWall in InIntEWalls:
					IntPath  = CWalls[iInEWall].WallPath
					IsOut = np.where(~IntPath.contains_points(np.array([Rknots2D,Zknots2D]).T, radius=dMin))[0]		#knots inside sub-zone border
					Rknots2D = Rknots2D[IsOut]
					Zknots2D = Zknots2D[IsOut]

			RKnotIn   = np.append(RKnotIn, Rknots2D)				#coordinate all inner knots
			ZKnotIn   = np.append(ZKnotIn, Zknots2D)

		RKnotAll	= np.concatenate((EWall.RKnotsE, RKnotExtB, RKnotIntB, RKnotIn))
		ZKnotAll	= np.concatenate((EWall.ZKnotsE, ZKnotExtB, ZKnotIntB, ZKnotIn))
		DelTri		= Delaunay(np.array([RKnotAll,ZKnotAll]).T)
		VTri		= np.copy(DelTri.simplices)

#		Remove triangles outside external wall
#		======================

		RWbord = np.copy(EWall.RWbord)
		ZWbord = np.copy(EWall.ZWbord)
		Area	= np.sum((ZWbord[1:]+ZWbord[:-1])*(RWbord[1:]-RWbord[:-1]))
		if(Area > 0.):
			RWbord = RWbord[::-1]
			ZWbord = ZWbord[::-1]
		BorderPath = Path(np.array([RWbord, ZWbord]).T, closed=True)

#		Find seguence of knots in the external wall

		RKnotExtBandE = np.append(EWall.RKnotsE, RKnotExtB)
		ZKnotExtBandE = np.append(EWall.ZKnotsE, ZKnotExtB)
		IsOnW = np.where(~BorderPath.contains_points(np.array([RKnotExtBandE, ZKnotExtBandE]).T,radius=-1.e-5))[0]		#Find external knots on external wall
		RKnotExtBandEW = RKnotExtBandE[IsOnW]
		ZKnotExtBandEW = ZKnotExtBandE[IsOnW]
		iSeqKnotsW = get_wall_knots_sequence(RWbord, ZWbord, RKnotExtBandEW, ZKnotExtBandEW)
		iSeqKnotsW = IsOnW[iSeqKnotsW]																					#indexes of external wall knots ordered
#		del RKnotExtBandE, ZKnotExtBandE, RKnotExtBandEW, ZKnotExtBandEW, IsOnW					#TO BE REACTIVATED

		if(False):
			print("\n\nAfter triangles generation for iEWall=",iEWall)
			Root.eirene_gen.Ax.triplot(RKnotAll, ZKnotAll, VTri, 'c-')
			Root.eirene_gen.Ax.plot(EWall.RWbord, EWall.ZWbord, color="blue", linestyle="-", linewidth=1, marker=".", markersize = 6)
			Rc = (RKnotAll[VTri[:,0]] + RKnotAll[VTri[:,1]] + RKnotAll[VTri[:,2]])/3.
			Zc = (ZKnotAll[VTri[:,0]] + ZKnotAll[VTri[:,1]] + ZKnotAll[VTri[:,2]])/3.
			Root.eirene_gen.Ax.plot(Rc, Zc, 'r.')
			Root.eirene_gen.Ax.plot(RKnotExtBandE[iSeqKnotsW], ZKnotExtBandE[iSeqKnotsW], color="green", linestyle="-", linewidth=1, marker="s", markersize = 4)
			Root.eirene_gen.Ax.plot(Rknots2D, Zknots2D, color="blue", linestyle="none", linewidth=1, marker="s", markersize = 4)
			Root.eirene_gen.Fig.canvas.draw()
			return

#		Find triangles out of external wall
#		===================

		VTriExt = np.copy(VTri)
		niSeqKnotsW = len(iSeqKnotsW)
		IsKnotWCon	= np.zeros(len(iSeqKnotsW), dtype='i4')								#Flag wall knots connected
		IsTriOn		= np.ones(VTri.shape[0], dtype='i4')									#Flag tringle is inside wall
		for kKnot in range(len(iSeqKnotsW)):
			nKnotW   	= iSeqKnotsW[kKnot]
			nKnotT   	= iSeqKnotsW[kKnot-1]
			ThetaWPrev	= atan2(ZKnotAll[nKnotT]-ZKnotAll[nKnotW],RKnotAll[nKnotT]-RKnotAll[nKnotW])
			if(ThetaWPrev < 0.): ThetaWPrev += TWOPI

			nKnotT   	= iSeqKnotsW[(kKnot+1)%niSeqKnotsW]
			ThetaWNext	= atan2(ZKnotAll[nKnotT]-ZKnotAll[nKnotW],RKnotAll[nKnotT]-RKnotAll[nKnotW])
			if(ThetaWNext < 0.): 		ThetaWNext += TWOPI
			if(ThetaWNext < ThetaWPrev):	ThetaWNext += TWOPI
			DThetaW	= ThetaWNext - ThetaWPrev

#			print("\R={:f}, Z={:f},ThetaWPrev={:f},ThetaWNext={:f},DTheta={:f}".format(RKnotAll[nKnotW],ZKnotAll[nKnotW],ThetaWPrev*PI_TO_DEG,ThetaWNext*PI_TO_DEG,DThetaW*PI_TO_DEG))
			iTriKnot, iVertKnot = np.where(VTri == nKnotW)

			VTriExt[iTriKnot,iVertKnot] = -(nKnotW + 1)
			for kTri in range(len(iVertKnot)):												#Search segments of triangles starting from nknotW out of wall
				for kVert in range(3):
					if((kVert != iVertKnot[kTri]) and (VTriExt[iTriKnot[kTri],kVert] >= 0)):
						nKnotT = VTriExt[iTriKnot[kTri],kVert]								#Number of knot on Triangle
						if(nKnotT == iSeqKnotsW[kKnot-1]):									#connected to previous wall knot
							IsKnotWCon[kKnot-1] = 1											#this should happen only for kknot=0
						elif(nKnotT == iSeqKnotsW[(kKnot+1)%niSeqKnotsW]):					#connected to next wall knot
							IsKnotWCon[kKnot] = 1
						else:
							ThetaKnotT	= atan2(ZKnotAll[nKnotT]-ZKnotAll[nKnotW],RKnotAll[nKnotT]-RKnotAll[nKnotW])
							if(ThetaKnotT < 0.): 			ThetaKnotT += TWOPI
							if(ThetaKnotT < ThetaWPrev):	ThetaKnotT += TWOPI
							if(ThetaKnotT-ThetaWPrev <= DThetaW+1e-10):						#Triangle out of wall (+1e-10 to account alligned points)
								IsTriOn[iTriKnot[kTri]]   = 0								#set ready to remove triangle
								VTriExt[iTriKnot[kTri],:] = -1								#set ready to remove triangle
								print("\tOff: R={:f}, Z={:f},Theta={:f},DTheta={:f}".format(RKnotAll[nKnotT],ZKnotAll[nKnotT],ThetaKnotT*PI_TO_DEG,(ThetaKnotT-ThetaWPrev)*PI_TO_DEG))
#
		IsIn    = np.where(IsTriOn == 1)[0]
		VTri   	= VTri[IsIn,:]																#Remove external triangles
		IsTriOn	= np.ones(VTri.shape[0], dtype='i4')										#Flag triangle is inside wall

		if(False):
			print("\n\nAfter removal of triangles external to the wall")
			Root.eirene_gen.Ax.triplot(RKnotAll, ZKnotAll, VTri, 'c-')
			Root.eirene_gen.Ax.plot(EWall.RWbord, EWall.ZWbord, color="blue", linestyle="-", linewidth=1, marker=".", markersize = 6)
			Rc = (RKnotAll[VTri[:,0]] + RKnotAll[VTri[:,1]] + RKnotAll[VTri[:,2]])/3.
			Zc = (ZKnotAll[VTri[:,0]] + ZKnotAll[VTri[:,1]] + ZKnotAll[VTri[:,2]])/3.
			Root.eirene_gen.Ax.plot(RKnotExtBandE[iSeqKnotsW], ZKnotExtBandE[iSeqKnotsW], color="green", linestyle="-", linewidth=1, marker="s", markersize = 4)
			Root.eirene_gen.Ax.plot(Rknots2D, Zknots2D, color="blue", linestyle="none", linewidth=1, marker="s", markersize = 4)
			Root.eirene_gen.Ax.plot(Rc, Zc, 'r.')
			Root.eirene_gen.Fig.canvas.draw()
			return

#		Fix holes on the external wall
#		=================	

		IsNotWcon = np.where(IsKnotWCon == 0)[0]
		if(len(IsNotWcon) > 0):
#			print("IsKnotWCon=",IsKnotWCon)
			for iNotWCon in range(len(IsNotWcon)):											#Loop on non-connected wall segments
				nKnot	   = iSeqKnotsW[IsNotWcon[iNotWCon]]
				nKnotLast  = iSeqKnotsW[(IsNotWcon[iNotWCon]+1)%niSeqKnotsW]
				RWs = np.array([RKnotAll[nKnot]])
				ZWs = np.array([ZKnotAll[nKnot]])
				RWe = np.array([RKnotAll[nKnotLast]])
				ZWe = np.array([ZKnotAll[nKnotLast]])
#				print("\nKnot={:d}, R={:f}, Z={:f},nKnotLast={:d}".format(nKnot,RKnotAll[nKnot],ZKnotAll[nKnot],nKnotLast))

				nKnotPrev  = iSeqKnotsW[(IsNotWcon[iNotWCon]-1)]
				ThetaWPrev = atan2(ZKnotAll[nKnotPrev]-ZKnotAll[nKnot],RKnotAll[nKnotPrev]-RKnotAll[nKnot])
				if(ThetaWPrev < 0.): ThetaWPrev += TWOPI
				
				nKnots = np.array([nKnot], dtype = 'i4')
				iMax   = 0
				while (nKnot != nKnotLast):													#Find knots from one point to the next
					iMax		 += 1
					if(iMax>100):
						print("\nATTENTION SOMETHING IS WRONG IN THE CODE!\n")
						return

					iTriKnot, iVertKnot = np.where(VTri == nKnot)
					nKnotsV		= np.empty(2*len(iTriKnot), dtype='i4')
					ThetaKnotV	= np.empty(2*len(iTriKnot), dtype='f8')
					nTriV		= np.empty(2*len(iTriKnot), dtype='i4')
					RVs			= np.empty(2*len(iTriKnot), dtype='f8') 
					ZVs			= np.empty(2*len(iTriKnot), dtype='f8') + ZKnotAll[nKnot]
					RVe			= np.empty(2*len(iTriKnot), dtype='f8')
					ZVe			= np.empty(2*len(iTriKnot), dtype='f8')

					RVs[:]		= RKnotAll[nKnot]
					ZVs[:]		= ZKnotAll[nKnot]
					iSeg		= 0
					for kTri in range(len(iVertKnot)):
						for kVert in range(3):
							if(VTri[iTriKnot[kTri],kVert] != nKnot):
								nKnotsV[iSeg] = VTri[iTriKnot[kTri],kVert]
								nTriV[iSeg]	  = iTriKnot[kTri]
								RVe[iSeg]	  = RKnotAll[nKnotsV[iSeg]]
								ZVe[iSeg]	  = ZKnotAll[nKnotsV[iSeg]]
								if((nKnot == nKnots[0]) or (nKnotsV[iSeg] == nKnots[0]) or (nKnotsV[iSeg] == nKnotLast)):				#To void intersection with first and last node
									RVs[iSeg]	= -RVs[iSeg]
									RVe[iSeg]	= -RVe[iSeg]

								ThetaKnotV[iSeg] = atan2(ZKnotAll[nKnotsV[iSeg]]-ZKnotAll[nKnot],RKnotAll[nKnotsV[iSeg]]-RKnotAll[nKnot])
								if(ThetaKnotV[iSeg] < 0.): 			ThetaKnotV[iSeg] += TWOPI
								if(ThetaKnotV[iSeg] <= ThetaWPrev):	ThetaKnotV[iSeg] += TWOPI
#								print("\tnKnotsV={:d}, R={:f}, Z={:f}".format(nKnotsV[iSeg],RKnotAll[nKnotsV[iSeg]],ZKnotAll[nKnotsV[iSeg]]))
								iSeg += 1

#					search intersection with wall segment

					Ris, Zis, is1, isV = intersect_segments(RWs, ZWs, RWe, ZWe, RVs, ZVs, RVe, ZVe)
					if(len(isV) > 0):
#						print("\tIntersection")
#						for ii in isV:
#							print("\tnKnotsV={:d}, R={:f}, Z={:f}".format(nKnotsV[ii],RKnotAll[nKnotsV[ii]],ZKnotAll[nKnotsV[ii]]))						
						IsTriOn[nTriV[isV]] = 0								            #remove triangles with intersection
						VTri[nTriV[isV],:]	= -1										#Ready to remove
						notInters		= np.ones(2*len(iTriKnot), dtype='i4')
						notInters[isV]  = 0
						IsIn = np.where(notInters == 1)[0]
						nKnotsV 	= nKnotsV[IsIn]
						ThetaKnotV	= ThetaKnotV[IsIn]

					nKnot	= nKnotsV[np.argmin(ThetaKnotV)]	
					ThetaWPrev = atan2(ZKnotAll[nKnots[-1]]-ZKnotAll[nKnot],RKnotAll[nKnots[-1]]-RKnotAll[nKnot])
					if(ThetaWPrev < 0.): ThetaWPrev += TWOPI

					nKnots	= np.append(nKnots,nKnot)
#					print("\tnKnot={:d}, R={:f}, Z={:f}".format(nKnot,RKnotAll[nKnot],ZKnotAll[nKnot]))
				
				DelTri		= Delaunay(np.array([RKnotAll[nKnots],ZKnotAll[nKnots]]).T)
				VTriNew		= nKnots[DelTri.simplices]										#Update knots number
				VTri		= np.append(VTri, VTriNew, axis=0)
				IsTriOn		= np.append(IsTriOn, np.ones(VTriNew.shape[0], dtype='i4'))

		if(False):
			print("\n\nAfter fixing hole in externa wall")
			Root.eirene_gen.Ax.plot(EWall.RWbord, EWall.ZWbord, color="blue", linestyle="-", linewidth=1, marker=".", markersize = 6)
			Root.eirene_gen.Ax.triplot(RKnotAll, ZKnotAll, VTri, 'b-')
#			Root.eirene_gen.Ax.triplot(RKnotAll, ZKnotAll, VTri, 'b-')
			Rc = (RKnotAll[VTri[:,0]] + RKnotAll[VTri[:,1]] + RKnotAll[VTri[:,2]])/3.
			Zc = (ZKnotAll[VTri[:,0]] + ZKnotAll[VTri[:,1]] + ZKnotAll[VTri[:,2]])/3.
			Root.eirene_gen.Ax.plot(Rc, Zc, 'r.')
			Root.eirene_gen.Ax.plot(RKnotExtBandE[iSeqKnotsW], ZKnotExtBandE[iSeqKnotsW], color="green", linestyle="-", linewidth=1, marker="s", markersize = 4)
			Root.eirene_gen.Fig.canvas.draw()
			return


#		Find triangles inside internall wall
#		===================

		if(len(AllIntEWalls) > 0):
			VTriExt = np.copy(VTri)
			nKnotIntBStart = nKnotIntBStart + len(EWall.RKnotsE) + len(RKnotExtB)
			for iIntEWalls in range(len(AllIntEWalls)):
				iSeqKnotsW	= np.arange(nKnotIntBStart[iIntEWalls],nKnotIntBStart[iIntEWalls+1])
				niSeqKnotsW = len(iSeqKnotsW)
				IsKnotWCon	= np.zeros(len(iSeqKnotsW), dtype='i4')								#Flag wall knots connected
				for kKnot in range(len(iSeqKnotsW)):
					nKnotW   	= iSeqKnotsW[kKnot]
					nKnotT   	= iSeqKnotsW[kKnot-1]
					ThetaWPrev	= atan2(ZKnotAll[nKnotT]-ZKnotAll[nKnotW],RKnotAll[nKnotT]-RKnotAll[nKnotW])
					if(ThetaWPrev < 0.): ThetaWPrev += TWOPI

					nKnotW   	= iSeqKnotsW[kKnot]
					nKnotT   	= iSeqKnotsW[(kKnot+1)%niSeqKnotsW]
					ThetaWNext	= atan2(ZKnotAll[nKnotT]-ZKnotAll[nKnotW],RKnotAll[nKnotT]-RKnotAll[nKnotW])
					if(ThetaWNext < 0.): 		ThetaWNext += TWOPI
					if(ThetaWNext < ThetaWPrev):	ThetaWNext += TWOPI
					DThetaW	= ThetaWNext - ThetaWPrev

					iTriKnot, iVertKnot = np.where(VTri == nKnotW)
					VTriExt[iTriKnot,iVertKnot] = -(nKnotW + 1)
					for kTri in range(len(iVertKnot)):												#Search segments of triangles starting from nknotW out of wall
						for kVert in range(3):
							if((kVert != iVertKnot[kTri]) and (VTriExt[iTriKnot[kTri],kVert] >= 0)):
								nKnotT = VTriExt[iTriKnot[kTri],kVert]								#Number of knot on Triangle
								if(nKnotT == iSeqKnotsW[kKnot-1]):									#connected to previous wall knot
									IsKnotWCon[kKnot-1] = 1											#this should happen only for kknot=0
								if(nKnotT == iSeqKnotsW[(kKnot+1)%niSeqKnotsW]):					#connected to next wall knot
									IsKnotWCon[kKnot] = 1
								else:
									ThetaKnotT	= atan2(ZKnotAll[nKnotT]-ZKnotAll[nKnotW],RKnotAll[nKnotT]-RKnotAll[nKnotW])
									if(ThetaKnotT < 0.): 		ThetaKnotT += TWOPI
									if(ThetaKnotT < ThetaWPrev):	ThetaKnotT += TWOPI
									if(ThetaKnotT-ThetaWPrev >= DThetaW-1e-10):						#Triangle inside the wall  (-1e-10 to account alligned points)
										IsTriOn[iTriKnot[kTri]]   = 0								#set ready to remove triangle
										VTriExt[iTriKnot[kTri],:] = -1								#set ready to remove triangle

		IsIn   = np.where(IsTriOn == 1)[0]
		IsOut = np.where(IsTriOn == 0)[0]
		VTri = VTri[IsIn,:]

		if(False):
			print("\n\nAfter removal of triangles inside internal wall")
			Root.eirene_gen.Ax.plot(EWall.RWbord, EWall.ZWbord, color="blue", linestyle="-", linewidth=1, marker=".", markersize = 6)
			Root.eirene_gen.Ax.triplot(RKnotAll, ZKnotAll, VTri, 'b-')
#			Root.eirene_gen.Ax.triplot(RKnotAll, ZKnotAll, VTri, 'b-')
			Rc = (RKnotAll[VTri[:,0]] + RKnotAll[VTri[:,1]] + RKnotAll[VTri[:,2]])/3.
			Zc = (ZKnotAll[VTri[:,0]] + ZKnotAll[VTri[:,1]] + ZKnotAll[VTri[:,2]])/3.
			Root.eirene_gen.Ax.plot(Rc, Zc, 'r.')
			Root.eirene_gen.Ax.plot(RKnotExtBandE[iSeqKnotsW], ZKnotExtBandE[iSeqKnotsW], color="green", linestyle="-", linewidth=1, marker="s", markersize = 4)
			Root.eirene_gen.Fig.canvas.draw()
			return

		TKnots 	= np.copy(VTri)																	#select triangles in domain and get knots number
		DirectD	= triangles_direct_or_indirect(TKnots, RKnotAll, ZKnotAll)						#get point orientation:  clockwise/counterclockwise

		if(DirectEP != DirectEP):																#reverse point orientation
			TKnots[:,1] = VTri[:,2]
			TKnots[:,2] = VTri[:,1]

#		Now update EIRENE knots and triangles
#		=======================

#		Extend EIRENE knots
#		============

		nOldKnots = len(RKnotsP)
		RKnotsP  = np.concatenate((RKnotsP, RKnotExtB, RKnotIntB, RKnotIn))
		ZKnotsP  = np.concatenate((ZKnotsP, ZKnotExtB, ZKnotIntB, ZKnotIn))
		nKnotAll  = nOldKnots + len(RKnotAll) - len(EWall.RKnotsE)

		nNewKnots = np.arange(nOldKnots - len(EWall.RKnotsE), nKnotAll, dtype = 'i4')
		nNewKnots[:len(EWall.RKnotsE)] = EWall.nKnotsE											#set plasma wall knots to their number 

		TKnots = nNewKnots[TKnots]																#Update knots number

#		Extend Triangles

		nOldTriangles = Triangles.k.shape[0]
		nNewTriangles = len(IsIn) + nOldTriangles
		Triangles = triangles_reshape(Triangles, nTriangles=nNewTriangles)

		Triangles.Area[nOldTriangles:]  = 0							#Area toward wall

		Triangles.k[nOldTriangles:]		= -1
		Triangles.i[nOldTriangles:]		= -1
		Triangles.j[nOldTriangles:]		= -1

		Triangles.np_k[nOldTriangles:]	= -1
		Triangles.np_i[nOldTriangles:]	= -1
		Triangles.np_j[nOldTriangles:]	= -1

		Triangles.p1[nOldTriangles:]	= TKnots[:,0]
		Triangles.p2[nOldTriangles:]	= TKnots[:,1]
		Triangles.p3[nOldTriangles:]	= TKnots[:,2]

		Triangles.step[nOldTriangles:]	= 1

		Triangles.PlasmaVacuum[nOldTriangles:] = 1


	print("generate_vacuum_triangles: completed")

	return Triangles, RKnotsP, ZKnotsP
	

def get_wall_knots_sequence(RWall, ZWall, RKnots, ZKnots):

	iFreeKnots = np.arange(len(RKnots), dtype='i4')
	iSeq = np.array([], dtype='i4')
	for i in range(len(RWall)):
		if(len(iFreeKnots) > 0):
			d =get_dist_from_segment(RWall[i-1], ZWall[i-1], RWall[i], ZWall[i], RKnots[iFreeKnots], ZKnots[iFreeKnots])
			OnOut = np.where(d < 1e-5, 1, 0)
			iOn   = np.where(OnOut == 1)[0]
			iOut  = np.where(OnOut == 0)[0]
			if(len(iOn) == 1):
				iSeq 		= np.append(iSeq, iFreeKnots[iOn[0]])
				iFreeKnots	= iFreeKnots[iOut]
			elif(len(iOn) > 1):
				iOn			= iFreeKnots[iOn]
				d1			= (RKnots[iOn]-RWall[i-1])**2+(ZKnots[iOn]-ZWall[i-1])**2
				iSort		= np.argsort(d1)
				iSeq		= np.append(iSeq, iOn[iSort])
				iFreeKnots	= iFreeKnots[iOut]

	if(len(iFreeKnots) > 0):
		print("\tERROR: At the end of wall len(iFreeKnots)=",len(iFreeKnots))
		exit()
		
	return iSeq

