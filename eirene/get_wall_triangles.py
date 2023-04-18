
# Function definition is here

import numpy						as np

#=========================================================
# This routine get wall triangles sequence
#=========================================================

def get_wall_triangles(Eirene, nTri=None, Side=None):

	print("get_wall_triangles")

	RKnots = Eirene.RKnots
	ZKnots = Eirene.ZKnots
	if(nTri is None): nTri	 = Eirene.WallTriangles.ntri

	WallTriKnots = np.array([Eirene.Triangles.p1[nTri], Eirene.Triangles.p2[nTri], Eirene.Triangles.p3[nTri] ]).T
	if(Side is not None):
		iKnots	= WallTriKnots[np.arange(len(Side)),Side]
		RWalls	= RKnots[iKnots]
		ZWalls	= ZKnots[iKnots]

		print("get_wall_triangles.1: Completed")

		return RWalls, ZWalls, iKnots

	WallTypeFace = np.array([Eirene.Triangles.BC1[nTri], Eirene.Triangles.BC2[nTri], Eirene.Triangles.BC3[nTri]]).T

	nWallTriangles = nTri.shape[0]
	TriSeq	= np.arange(nWallTriangles, dtype='i4')
	Side	= np.arange(nWallTriangles, dtype='i4')
	iKnots	= np.arange(nWallTriangles, dtype='i4')
	pR1		= np.empty(nWallTriangles, dtype='f8')
	pZ1		= np.empty(nWallTriangles, dtype='f8')
	pR2		= np.empty(nWallTriangles, dtype='f8')
	pZ2		= np.empty(nWallTriangles, dtype='f8')

	ind2 = [1,2,0]
	for iSide in range(3):
#		ii = np.where(WallTypeFace[:,iSide] == 1); ii = ii[0]									#TypeFace = 1 mean wall boundary
		ii = np.where(WallTypeFace[:,iSide] > 1); ii = ii[0]						#TypeFace > 1 mean wall/pump boundary
		if(len(ii) > 0):
			Side[ii]	= iSide
			iKnots[ii]	= WallTriKnots[ii,iSide]
			pR1[ii]		= RKnots[iKnots[ii]]
			pZ1[ii]		= ZKnots[iKnots[ii]]

			pR2[ii] = RKnots[WallTriKnots[ii,ind2[iSide]]]
			pZ2[ii] = ZKnots[WallTriKnots[ii,ind2[iSide]]]

	for iSeq  in range(nWallTriangles-2):
		iMin = np.argmin((pR1[TriSeq[iSeq+1:nWallTriangles]] - pR2[TriSeq[iSeq]])**2 + (pZ1[TriSeq[iSeq+1:nWallTriangles]] - pZ2[TriSeq[iSeq]])**2)
		iTri 					= TriSeq[iSeq+1]
		TriSeq[iSeq+1] 			= TriSeq[iMin + iSeq + 1]
		TriSeq[iMin + iSeq + 1]	= iTri

	RWalls	= np.copy(pR1[TriSeq])													#R & Z coordiante triangles
	ZWalls	= np.copy(pZ1[TriSeq])
	nTri	= np.copy(nTri[TriSeq])													# indexes triangles
	Side	= np.copy(Side[TriSeq])													# indexes wall triangles iSide
	iKnots	= np.copy(iKnots[TriSeq])												# indexes of knots

	print("get_wall_triangles.2: Completed")
	return RWalls, ZWalls, iKnots, nTri, Side

