import numpy				as np
from routines.globals		import DEBUG
#=========================================================
# This routine find triangles seguence at the wall
#=========================================================

#def find_sequence(WallTriangles, RKnots, ZKnots, Wall=[0], PlotCheck=0):
def find_sequence(WallTriangles, RKnots, ZKnots):

	if(DEBUG > 0):	print("find_sequence")

	nWallTriangles = len(WallTriangles)

	RKnotsTri = np.array([RKnots[WallTriangles.p1],RKnots[WallTriangles.p2],RKnots[WallTriangles.p3]])*100.
	ZKnotsTri = np.array([ZKnots[WallTriangles.p1],ZKnots[WallTriangles.p2],ZKnots[WallTriangles.p3]])*100.

	TriSeq = np.arange(nWallTriangles, dtype='i4')

	pR1 = np.empty(nWallTriangles, dtype='f8')
	pZ1 = np.empty(nWallTriangles, dtype='f8')
	pR2 = np.empty(nWallTriangles, dtype='f8')
	pZ2 = np.empty(nWallTriangles, dtype='f8')

	SideIndex1 = np.array([0,1,2])
	SideIndex2 = np.array([1,2,0])

	pR1[TriSeq] = RKnotsTri[SideIndex1[WallTriangles.side[TriSeq]], TriSeq]
	pZ1[TriSeq] = ZKnotsTri[SideIndex1[WallTriangles.side[TriSeq]], TriSeq]

	pR2[TriSeq] = RKnotsTri[SideIndex2[WallTriangles.side[TriSeq]], TriSeq]
	pZ2[TriSeq] = ZKnotsTri[SideIndex2[WallTriangles.side[TriSeq]], TriSeq]
	
	dSqMax = 1e-16
	iSeqStart = [0]
	for iSeq  in range(nWallTriangles-2):
		dSq = (pR1[TriSeq[iSeq+1:nWallTriangles]] - pR2[TriSeq[iSeq]])**2 + (pZ1[TriSeq[iSeq+1:nWallTriangles]] - pZ2[TriSeq[iSeq]])**2
		iMin = np.argmin(d)
		iTri 					= TriSeq[iSeq+1]
		TriSeq[iSeq+1] 			= TriSeq[iMin + iSeq + 1]
		TriSeq[iMin + iSeq + 1]	= iTri

		if(dSq[iMin] > dSqMax):
			if(len(iSeqStart) == 0): 	iSeqEnd = [iSeq]
			else:						iSeqEnd.append(iSeq)
			iSeqStart.append(iSeq+1)

	if(len(iSeqStart) == 0): 	iSeqEnd = [nWallTriangles-1]
	else:						iSeqEnd.append(nWallTriangles-1)

	TriSequences = []
	for iWall in range(len(iSeqStart)):
		TriSequences.append(TriSeq[iSeqStart[iWall]:iSeqEnd[iWall]+1])

	if(DEBUG > 0):	print("find_sequence: Completed")
	
	return TriSequences
