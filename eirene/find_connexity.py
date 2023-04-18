from math						import floor
import numpy					as np
from interfaces.plot_and_ask	import plot_and_ask
from eirene.triangles_routines	import triangles_copy
from routines.utils_walls		import set_XYlines_walls
from routines.globals 			import *

# Function definition is here
#=======================================================================================================================================

def find_connexity(Root, Triangles, RKnots, ZKnots):

	if(DEBUG > 0):	print("find_connexity")

#	Set as connex by default triangles in contact with the core

	connexDomain = np.where((Triangles.BC1 == BC_CORE) | (Triangles.BC2 == BC_CORE) | (Triangles.BC3 == BC_CORE), 1, 0)

#	 Propagate connexity in while loop

	TotFound   = 1
	keep_going = True
	while(keep_going):
		iTriConnex = np.where(connexDomain == 1)[0]

		iBC1 = np.where((Triangles.BC1[iTriConnex] == BC_TRIANGLE) & (connexDomain[Triangles.neigh1[iTriConnex]]==0))[0]
		if(len(iBC1) > 0): connexDomain[Triangles.neigh1[iTriConnex[iBC1]]] = 1

		iBC2 = np.where((Triangles.BC2[iTriConnex] == BC_TRIANGLE) & (connexDomain[Triangles.neigh2[iTriConnex]]==0))[0]
		if(len(iBC2) > 0): connexDomain[Triangles.neigh2[iTriConnex[iBC2]]] = 1

		iBC3 = np.where((Triangles.BC3[iTriConnex] == BC_TRIANGLE) & (connexDomain[Triangles.neigh3[iTriConnex]]==0))[0]
		if(len(iBC3) > 0): connexDomain[Triangles.neigh3[iTriConnex[iBC3]]] = 1

		if((len(iBC1) + len(iBC2) + len(iBC3)) == 0): keep_going = False

	Triangles.step = np.copy(connexDomain)				#save in legacy mesher variable

	nTriangles = len(Triangles)
	nConnected = np.sum(Triangles.step)
	if((DEBUG > 1) and (nTriangles > nConnected > 0)):
		print("\tFound {:d} non connected triangles on {:d} triangles".format(nTriangles-nConnected, nTriangles))

		iNoCon = np.where(Triangles.step == 0)[0]

		Xarrs, Yarrs, lines = set_XYlines_walls(Root.Config,wall_line="c--")
		Xarrs.append(np.array([RKnots, ZKnots]))
		Yarrs.append(np.array([Triangles.p1[iNoCon], Triangles.p2[iNoCon], Triangles.p3[iNoCon]]).T)
		lines.append("r-")
		LinesData = [Xarrs, Yarrs, lines]
		choice = plot_and_ask(Root, LinesData=LinesData, title="Non connected triangles (red)")

	if(DEBUG > 0):	print("find_connexity: Completed")

	return