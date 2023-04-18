import numpy					as np
from eirene.triangles_routines	import triangles_copy
from interfaces.plot_and_ask	import plot_and_ask
from routines.utils_walls		import set_XYlines_walls
from routines.globals			import DEBUG, MIN_AREA_TRIANGLES, MIN_ORDER_TRIANGLES, BC_WALL, BC_CORE, DRAWING_WALL

# This function replace MATLab function remove_double_wall_traingles (calling later remove_unused_triangles)
# Function definition is here
#=======================================================================================================================================

def check_triangles_orientation_size(Root, RKnots, ZKnots, Triangles):

	if(DEBUG > 0):	print("check_triangles_orientation_size")

	iWall	= np.where(Triangles.Area == -1)[0]
	Triangles.Area[iWall] = 0

	QuadR	= np.array([RKnots[Triangles.p1[iWall]], RKnots[Triangles.p2[iWall]], RKnots[Triangles.p3[iWall]], RKnots[Triangles.p1[iWall]]])
	QuadZ	= np.array([ZKnots[Triangles.p1[iWall]], ZKnots[Triangles.p2[iWall]], ZKnots[Triangles.p3][iWall], ZKnots[Triangles.p1[iWall]]])

	AreaT   = 0.5*np.abs(np.sum((QuadZ[:1,:] + QuadZ[:-1,:])*(QuadR[:1,:] - QuadR[:-1,:]), axis = 0))
	iSmallT	= np.where(AreaT <  MIN_AREA_TRIANGLES)[0]
	QuadR	= 0; QuadZ = 0; AreaT = 0

	AR = RKnots[Triangles.p1[iWall]]*(ZKnots[Triangles.p2[iWall]] - ZKnots[Triangles.p3[iWall]]) + \
		 RKnots[Triangles.p2[iWall]]*(ZKnots[Triangles.p3[iWall]] - ZKnots[Triangles.p1[iWall]]) + \
		 RKnots[Triangles.p3[iWall]]*(ZKnots[Triangles.p1[iWall]] - ZKnots[Triangles.p2][iWall])

	iWrongO = np.where(AR < MIN_ORDER_TRIANGLES)[0]
	AR = 0

	iSmallT = iWall[iSmallT]
	iWrongO = iWall[iWrongO]

	Triangles.step = np.ones(Triangles.Area.size, dtype = 'i4')


	if((len(iWrongO) > 0) or (len(iSmallT) > 0)):
		OldTriangles			= triangles_copy(Triangles)						#save old triangles
		Triangles.step[iWrongO] = 0
		Triangles.step[iSmallT] = 0
	else:
		return Triangles

	iRemove = np.where(Triangles.step == 0)[0]
	iKeep   = np.where(Triangles.step == 1)[0]

	AddRemoved = np.empty(0, dtype='i4')
	kLoop = 0
	while (len(iRemove) > 0):

#	Fix Boundary condiction

		for kRemove in iRemove:
			iTriAtSide = np.where(Triangles.neigh1[iKeep] == kRemove)[0]
			if(len(iTriAtSide) > 0):
				iTriAtSide = iKeep[iTriAtSide] 
				Triangles.neigh1[iTriAtSide]		= -1
				Triangles.typeneigh1[iTriAtSide]	= -1
				Triangles.BC1[iTriAtSide]			= BC_WALL

			iTriAtSide = np.where(Triangles.neigh2[iKeep] == kRemove)[0]
			if(len(iTriAtSide) > 0):
				iTriAtSide = iKeep[iTriAtSide] 
				Triangles.neigh2[iTriAtSide]		= -1
				Triangles.typeneigh2[iTriAtSide]	= -1
				Triangles.BC2[iTriAtSide]			= BC_WALL


			iTriAtSide = np.where(Triangles.neigh3[iKeep] == kRemove)[0]
			if(len(iTriAtSide) > 0):
				iTriAtSide = iKeep[iTriAtSide] 
				Triangles.neigh3[iTriAtSide]		= -1
				Triangles.typeneigh3[iTriAtSide]	= -1
				Triangles.BC3[iTriAtSide]			= BC_WALL

#		Look to triangles with two sides as wall

		iKeep = np.where(Triangles.step == 1)[0]

		iWall = np.zeros(len(iKeep), dtype = 'i4')
		iWall += np.where((Triangles.BC1[iKeep] == BC_WALL) & (Triangles.neigh1[iKeep] == -1), 1, 0)
		iWall += np.where((Triangles.BC2[iKeep] == BC_WALL) & (Triangles.neigh2[iKeep] == -1), 1, 0)
		iWall += np.where((Triangles.BC3[iKeep] == BC_WALL) & (Triangles.neigh3[iKeep] == -1), 1, 0)

		iMulti = np.where(iWall > 1)[0]
		if(len(iMulti) > 0):
			iRemove = iKeep[iMulti]
			Triangles.step[iRemove] = 0
			iKeep   = np.where(Triangles.step == 1)[0]
			AddRemoved = np.append(AddRemoved, iRemove)
		else:
			iRemove = []

	if((len(iWrongO) > 0) or (len(iSmallT) > 0)):

		if(len(iWrongO) > 0):	 print("\tWarning!!!: Found small orientation for {:d} triangles".format(len(iWrongO)))
		if(len(iSmallT) > 0):	 print("\tWarning!!!: Found small size for {:d} triangles".format(len(iSmallT)))
		if(len(AddRemoved) > 0): print("\tRemoved also {:d} triangles to avoid two sides on wall".format(len(AddRemoved)))

		Xarrs,  Yarrs, lines = set_XYlines_walls(Root.Config, wall_line="c--")

		if(len(iWrongO) > 0): 
			Xarrs.append(np.array([RKnots, ZKnots]))
			Yarrs.append(np.array([Triangles.p1[iWrongO], Triangles.p2[iWrongO], Triangles.p3[iWrongO]]).T)
			lines.append("r-")

		if(len(iSmallT) > 0): 
			Xarrs.append(np.array([RKnots, ZKnots]))
			Yarrs.append(np.array([Triangles.p1[iSmallT], Triangles.p2[iSmallT], Triangles.p3[iSmallT]]).T)
			lines.append("b-")

		if(len(AddRemoved) > 0): 
			Xarrs.append(np.array([RKnots, ZKnots]))
			Yarrs.append(np.array([Triangles.p1[AddRemoved], Triangles.p2[AddRemoved], Triangles.p3[AddRemoved]]).T)
			lines.append("g-")

		LinesData = [Xarrs, Yarrs, lines]
		choice = plot_and_ask(Root, LinesData=LinesData, title="Wrong=red, small=blue, add removed=green")
		if(choice[0] == 0):
			print("keep triangles!")
			if(DEBUG > 0):	print("check_triangles_orientation_size: Completed")
			return OldTriangles

		else:
			if(DEBUG > 0):	print("check_triangles_orientation_size: Completed")
			return Triangles

	else:
		if(DEBUG > 0):	print("check_triangles_orientation_size: Completed")
		return OldTriangles

