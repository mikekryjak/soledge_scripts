import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot 		as pyp
from eirene.polygon_routines 	import polygon_inout_intersect
#from eirene.polygon_routines 	import polygon_true_intersect, polygon_inout_intersect
from routines.input_routines 	import input_int
from routines.globals			import DEBUG

# Function definition is here
#=======================================================================================================================================

#def generate_knots_on_wall(Zones, RKnots, ZKnots, WallPath, Wall, PlotCheck=0):
def generate_knots_on_wall(Root, Config, Zones, RKnots, ZKnots):

	if(DEBUG > 0):	print("generate_knots_on_wall")

#	on s'occupe des noeuds a l'intersection des mailles soledge et de la paroi


	iPwalls = Config.iPwalls
	Walls	= Config.Walls
	nZones  = len(Zones)
	nCross  = 0
	for k in range(nZones):
		for iPwall in range(Zones[k].IsCrossed.shape[3]):
			IsCrossed = np.sum(Zones[k].IsCrossed[:,:,:,iPwall], axis=2)
			ii, jj	  = np.where((Zones[k].IsAligned == 0) & (IsCrossed >= 1) & (Zones[k].InPlasma == 1))
			nCross += len(ii)

	nMaxKnotsWall = nCross*4										#max four intersections on each quadrangol
	RWallKnots	= np.empty(nMaxKnotsWall, dtype='f8')
	ZWallKnots	= np.empty(nMaxKnotsWall, dtype='f8')

#	CutSegs 	= np.zeros(2, dtype='i4')
	nKnotsWall	= 0
	for k in range(nZones):
		Zones[k].KnotE = -np.ones(Zones[k].KnotA.shape, dtype='i4')
		Zones[k].KnotF = -np.ones(Zones[k].KnotA.shape, dtype='i4')
		Zones[k].KnotG = -np.ones(Zones[k].KnotA.shape, dtype='i4')
		Zones[k].KnotH = -np.ones(Zones[k].KnotA.shape, dtype='i4') 
		KnotsM		   = [Zones[k].KnotE,     Zones[k].KnotF,     Zones[k].KnotG,     Zones[k].KnotH]
		InPlasmasC	   = [Zones[k].InPlasmaA, Zones[k].InPlasmaB, Zones[k].InPlasmaC, Zones[k].InPlasmaD, Zones[k].InPlasmaA]

		for iPwall in range(Zones[k].IsCrossed.shape[3]):
			IsCrossed 		  = np.sum(Zones[k].IsCrossed[:,:,:,iPwall], axis=2)
			ii, jj = np.where((Zones[k].IsAligned == 0) & (IsCrossed == 2) & (Zones[k].InPlasma == 1))
			for l in range(len(ii)):
				i = ii[l]
				j = jj[l]
				iKnotA = Zones[k].KnotA[i,j]
				iKnotB = Zones[k].KnotB[i,j]
				iKnotC = Zones[k].KnotC[i,j]
				iKnotD = Zones[k].KnotD[i,j]
				iKnots = [iKnotA, iKnotB, iKnotC, iKnotD, iKnotA]
				
				for ls in range(4):
	#				%segment is

					if(Zones[k].IsCrossed[i,j,ls,iPwall] == 1):
						if(InPlasmasC[ls][i,j] == 2):
							RWallKnots[nKnotsWall]	= RKnots[iKnots[ls]]
							ZWallKnots[nKnotsWall]	= ZKnots[iKnots[ls]]
						elif(InPlasmasC[ls+1][i,j] == 2):
							RWallKnots[nKnotsWall]	= RKnots[iKnots[ls+1]]
							ZWallKnots[nKnotsWall]	= ZKnots[iKnots[ls+1]]
						else:
							r_w, z_w, ind  = polygon_inout_intersect(RKnots[iKnots[ls]],   ZKnots[iKnots[ls]], \
																RKnots[iKnots[ls+1]], ZKnots[iKnots[ls+1]], Walls[iPwalls[iPwall]].WallPath)
	#						r_w, z_w, ind  = polygon_true_intersect(Root, RKnots[iKnots[ls]], ZKnots[iKnots[ls]], \
	#													 	RKnots[iKnots[ls+1]], ZKnots[iKnots[ls+1]], RWall, ZWall)
							if(ind == 0):														#segment crossing wall
								RWallKnots[nKnotsWall]	= r_w
								ZWallKnots[nKnotsWall]	= z_w
							else:
								print("\t\tCrazy!! I should not be here with ind=", ind)

						KnotsM[ls][i,j] 		= nKnotsWall
						nKnotsWall				= nKnotsWall+1
						

	if(DEBUG > 1):	print("\tgenerated nodes in wall =",nKnotsWall)
	if(DEBUG > 1):	print("\tremoving nodes with a distance lower than 1e-6")

	RWallKnots = np.copy(RWallKnots[:nKnotsWall])
	ZWallKnots = np.copy(ZWallKnots[:nKnotsWall])

#	remove nodes with a distance lower than 1e-6

	ToKnotsWall		= np.empty(nKnotsWall, dtype='i4')			#Array with true nodes number after removing near nodes
	npt		   	 	= 0
	ToKnotsWall[0]	= npt
	for n in range(1,nKnotsWall):
		npt		   	 = npt+1
		ToKnotsWall[n] = npt
		ii 		 	 = np.where(np.sqrt((RWallKnots[n]-RWallKnots[0:n])**2+(ZWallKnots[n]-ZWallKnots[0:n])**2) < 1e-6); ii=ii[0]
		if len(ii) > 0:
			ToKnotsWall[n] = ToKnotsWall[ii[0]]
			npt			 = npt-1
	
	if((npt+1 < nKnotsWall) and (DEBUG > 0)): print("\tremoved ", nKnotsWall-npt-1, "nodes of ",nKnotsWall," nodes")

#	Prepare (R,Z) and nodes array for wall nodes

	nKnotsWall = np.max(ToKnotsWall)+1
	RKnotsWall = np.empty(nKnotsWall, dtype='f8')
	ZKnotsWall = np.empty(nKnotsWall, dtype='f8')

#	Define new nodes array with node number  
	
	for nKnot in range (nKnotsWall):
		i = np.where(ToKnotsWall == nKnot); i=i[0][0]						#first element of touple & first value
		RKnotsWall[nKnot]		= RWallKnots[i]
		ZKnotsWall[nKnot]		= ZWallKnots[i]

	"""
	if PlotCheck == 1:
		print("\tPlot points")

		fig  = pyp.figure()
		ax   = fig.add_subplot(1,1,1)
		ax.set_title("Wall KnotsM")
		ax.set_xlabel("X (m)")
		ax.set_ylabel("Z (m)")

		ax.plot(RWallKnots,ZWallKnots,'ro')
		ax.plot(RKnotsWall,ZKnotsWall,'mv')
		ax.plot(Wall[:,0],Wall[:,1],'k-')
		fig.show(True)
	"""
	if(DEBUG > 0):	print("generate_knots_on_wall: Completed")

	return ToKnotsWall, RKnotsWall, ZKnotsWall

#	This routine is not used!!!!!!
#	================

def multiple_crossing_not_used(Root, Zones, RKnots, ZKnots, WallPath, RWallKnots, ZWallKnots, CutSegs, k, i, j, nKnotsWall, Message=0):

	iKnodeA = Zones[k].KnotA[i,j]
	iKnodeB = Zones[k].KnotB[i,j]
	iKnodeC = Zones[k].KnotC[i,j]
	iKnodeD = Zones[k].KnotD[i,j]
	RepeatCut = False
	
	for nSeg in range(2):
		CutSeg = CutSegs[nSeg]
		if(CutSeg == 1):
			if(Zones[k].InPlasmaA[i,j] + Zones[k].InPlasmaB[i,j] == 1):
#				r_w, z_w			= polygon_intersect(RKnots[iKnodeA], ZKnots[iKnodeA], RKnots[iKnodeB], ZKnots[iKnodeB], WallPath)
				r_w, z_w			= polygon_intersect(Root, RKnots[iKnodeA], ZKnots[iKnodeA], RKnots[iKnodeB], ZKnots[iKnodeB], RWall, ZWall)
				RWallKnots[nKnotsWall]	= r_w
				ZWallKnots[nKnotsWall]	= z_w
				KnotE				= Zones[k].KnotE[i,j]
				Zones[k].KnotE[i,j] = nKnotsWall								
				nKnotsWall 			= nKnotsWall+1
				if( Message > 1): print("\tcrossed segment AB")
			else:
				RepeatCut = True
				if(Message > 0): print("\tAttention segment AB is not crossed!")

#		segment2
		elif(CutSeg == 2):
			if(Zones[k].InPlasmaB[i,j] + Zones[k].InPlasmaC[i,j] == 1):
#				r_w, z_w  		= polygon_intersect(RKnots[iKnodeB], ZKnots[iKnodeB], RKnots[iKnodeC], ZKnots[iKnodeC], WallPath)
				r_w, z_w  		= polygon_intersect(Root, RKnots[iKnodeB], ZKnots[iKnodeB], RKnots[iKnodeC], ZKnots[iKnodeC], RWall, ZWall)
				RWallKnots[nKnotsWall]	= r_w
				ZWallKnots[nKnotsWall]	= z_w
				KnotF				= Zones[k].KnotF[i,j]
				Zones[k].KnotF[i,j] = nKnotsWall									
				nKnotsWall 			= nKnotsWall+1
				if(Message > 1): print("\tcrossed segment BC")
			else:
				RepeatCut = True
				if(Message > 0): print("\tAttention segment BC is not crossed!")

#		segment3
		elif(CutSeg == 3):
			if(Zones[k].InPlasmaC[i,j] + Zones[k].InPlasmaD[i,j] == 1):
#				r_w, z_w   		= polygon_intersect(RKnots[iKnodeC], ZKnots[iKnodeC], RKnots[iKnodeD], ZKnots[iKnodeD], WallPath)
				r_w, z_w   		= polygon_intersect(Root, RKnots[iKnodeC], ZKnots[iKnodeC], RKnots[iKnodeD], ZKnots[iKnodeD], RWall, ZWall)
				RWallKnots[nKnotsWall]	= r_w
				ZWallKnots[nKnotsWall]	= z_w
				KnotG				= Zones[k].KnotG[i,j]
				Zones[k].KnotG[i,j] = nKnotsWall
				nKnotsWall 			= nKnotsWall+1
				if(Message > 1): print("\tcrossed segment CD")
			else:
				RepeatCut = True
				if(Message > 0): print("\tAttention segment CD is not crossed!")
								
#		segment4
		else:
			if(Zones[k].InPlasmaD[i,j] + Zones[k].InPlasmaA[i,j] == 1):
#				r_w, z_w 		 	= polygon_intersect(RKnots[iKnodeD], ZKnots[iKnodeD], RKnots[iKnodeA], ZKnots[iKnodeA], WallPath)
				r_w, z_w 		 	= polygon_intersect(Root, RKnots[iKnodeD], ZKnots[iKnodeD], RKnots[iKnodeA], ZKnots[iKnodeA], RWall, ZWall)
				RWallKnots[nKnotsWall]	= r_w
				ZWallKnots[nKnotsWall]	= z_w
				KnotH				= Zones[k].KnotH[i,j]
				Zones[k].KnotH[i,j] = nKnotsWall								
				nKnotsWall 			= nKnotsWall+1
				if(Message > 1): print("\tcrossed segment DA")
			else:
				RepeatCut = True
				if(Message > 0): print("\tAttention segment DA is not crossed!")

		if RepeatCut:
			if(Message > 0): print("\tATTENTION REPEAT CUT")
			if nSeg == 1:
				CutSeg = CutSegs[0]
				if CutSeg == 1:		Zones[k].KnotE[i,j] = KnotE
				elif CutSeg == 2:	Zones[k].KnotF[i,j] = KnotF
				elif CutSeg == 3:	Zones[k].KnotG[i,j] = KnotG
				else:				Zones[k].KnotH[i,j] = KnotH
				nKnotsWall = nKnotsWall-1
			break

	return RepeatCut, nKnotsWall
