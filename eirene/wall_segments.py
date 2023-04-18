import numpy as np
from routines.globals			import BC_CORE

def wall_segments(WallTriangles, RKnots, ZKnots):

	print("wall_segments")

	triNum		= WallTriangles.ntri
	triFace		= WallTriangles.side

	triKnots = np.array([WallTriangles.p1,WallTriangles.p2,WallTriangles.p3]).T

	SideIndex1 = np.array([0,1,2])
	SideIndex2 = np.array([1,2,0])

	ii = np.arange(len(triNum))
	P1 = triKnots[ii,triFace]
	P2 = triKnots[ii,SideIndex2[triFace]]

	TriSequences = [np.array([0],dtype='i4')]
	AllSequences = [np.array([0],dtype='i4')]
	Ntri = len(triNum)
	for itri in range(1,Ntri):
		ind = np.where(P1 == P2[TriSequences[-1][-1]])[0][0]
		if(TriSequences[-1][0] != ind):									#check  if end of closed wall
			TriSequences[-1] = np.append(TriSequences[-1], ind)
			AllSequences	 = np.append(AllSequences, ind)
		else:
			for i in range(1,Ntri):
				if(len(np.where(AllSequences == i)[0]) == 0):
					TriSequences.append([np.array([i],dtype='i4')])
					AllSequences	 = np.append(AllSequences, i)
					break

	del AllSequences, P1, P2

#	Then find wall segments segments are ordered ==> to go from WallTriangles to segment TrigSeq is necessary
#	================================================================

	TriR	= np.array([RKnots[WallTriangles.p1],RKnots[WallTriangles.p2],RKnots[WallTriangles.p3]]).T
	TriZ	= np.array([ZKnots[WallTriangles.p1],ZKnots[WallTriangles.p2],ZKnots[WallTriangles.p3]]).T

	TriFaces = []
	R12		 = []
	Z12		 = []
	for iWall in range(len(TriSequences)):
		nTri		= TriSequences[iWall]
		nFace		= triFace[TriSequences[iWall]]

		TriFaces.append(nFace)
		R12.append(np.array([TriR[nTri,nFace], TriR[nTri,SideIndex2[nFace]]]).T)
		Z12.append(np.array([TriZ[nTri,nFace], TriZ[nTri,SideIndex2[nFace]]]).T)

	print("wall_segments: completed")

	return TriSequences, TriFaces, R12, Z12
