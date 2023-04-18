import types
import numpy						as np
import numpy.matlib					as mat
from math							import sqrt
from matplotlib.path				import Path
from routines.intersect_segments	import intersect_segments, intersect_lines_to_segments, intersect_line_to_segment
from routines.intersect_contour		import intersect_2contours
from routines.globals				import DEBUG

def get_chords_triangle_path(Eirene, Diag):

	Diag.path = get_lines_triangles_path(Eirene, Diag.lines_start, Diag.lines_end)

	return

#===================================================================

def get_lines_triangles_path(Eirene, LinesSt, LinesEnd):

#	Extend lines to be sure to cross the wall
	
	RWmin = np.min(Eirene.Wall.R12[0][:,0])
	RWmax = np.max(Eirene.Wall.R12[0][:,0])
	ZWmin = np.min(Eirene.Wall.Z12[0][:,0])
	ZWmax = np.max(Eirene.Wall.Z12[0][:,0])
	Rmin = RWmin - 0.1*(RWmax-RWmin)
	Rmax = RWmax + 0.1*(RWmax-RWmin)
	Zmin = ZWmin - 0.1*(ZWmax-ZWmin)
	Zmax = ZWmax + 0.1*(ZWmax-ZWmin)

	RWs = np.array([Rmin, Rmax, Rmax, Rmin])
	RWe = np.array([Rmax, Rmax, Rmin, Rmin])
	ZWs = np.array([Zmin, Zmin, Zmax, Zmax])
	ZWe = np.array([Zmin, Zmax, Zmax, Zmin])

	ELinesSt = np.empty_like(LinesSt)
	ELinesEnd = np.empty_like(LinesEnd)
	for i in range(LinesSt.shape[1]):
		xis, zis, is1, is2 = intersect_lines_to_segments([LinesSt[0,i]], [LinesSt[2,i]], [LinesEnd[0,i]], [LinesEnd[2,i]], RWs, ZWs, RWe, ZWe)
		ELinesSt[0,i]  = xis[0]
		ELinesEnd[0,i] = xis[1]
		ELinesSt[1,i]  = LinesSt[1,i]
		ELinesEnd[1,i] = LinesEnd[1,i]
		ELinesSt[2,i]  = zis[0]
		ELinesEnd[2,i] = zis[1]

#	find all intersections between lines and wall

	xis, zis, ich, ims = intersect_lines_to_segments(ELinesSt[0,:], ELinesSt[2,:], ELinesEnd[0,:], ELinesEnd[2,:], 
													Eirene.Wall.R12[0][:,0], Eirene.Wall.Z12[0][:,0], 
													Eirene.Wall.R12[0][:,1], Eirene.Wall.Z12[0][:,1], remove_same = False)

	RKnots		 = Eirene.RKnots
	ZKnots		 = Eirene.ZKnots
	TriKnots	 = Eirene.TriKnots
	TriNeigh	 = Eirene.TriNeigh
	TriTypeneigh = Eirene.TriTypeneigh

#	Find intersections between lines and triangles starting from wall

	ChanRi	 = []
	ChanZi	 = []
	ChanTi	 = []
	ChanLeni = []
	for nChan in range(LinesSt.shape[1]):
		Ri	 = []
		Zi	 = []
		Ti	 = []
		Leni = []

		ii = np.where(ich == nChan)[0]
		if(len(ii) == 0):
			print("nChan=",nChan)
			print("LinesSt[:,nChan]  =",LinesSt[:,nChan])
			print("LinesEnd[:,nChan] =",LinesEnd[:,nChan])
			print("ELinesSt[:,nChan] =",ELinesSt[:,nChan])
			print("ELinesEnd[:,nChan]=",ELinesEnd[:,nChan])

		LenLine	  = sqrt((LinesEnd[0,nChan] - LinesSt[0,nChan])**2 +(LinesEnd[2,nChan] - LinesSt[2,nChan])**2)
		dWalli	= ((xis[ii] - LinesSt[0,nChan])*(LinesEnd[0,nChan] - LinesSt[0,nChan]) + \
				   (zis[ii] - LinesSt[2,nChan])*(LinesEnd[2,nChan] - LinesSt[2,nChan]))/LenLine

		ii		 = ii[np.argsort(dWalli)]															#sort on distance from starting of line
		k  		 = 0

		InPlasma = True
		while (InPlasma):
			StartTriangle = Eirene.WallTriangles.ntri[Eirene.Wall.TriSequences[Eirene.Wall.iTriSeqExtPlasma][ims[ii[k]]]]
			SubRi, SubZi, SubTi, SubSi, SubKi, SubWKi = get_triangles_on_line(RKnots, ZKnots, TriKnots, TriNeigh, TriTypeneigh, StartTriangle, 
																				ELinesSt[0,nChan], ELinesSt[2,nChan], ELinesEnd[0,nChan], ELinesEnd[2,nChan])
			Ri.append(SubRi)
			Zi.append(SubZi)
			Ti.append(SubTi)

			Leni.append(np.sqrt((SubRi[1:]-SubRi[:-1])**2 + (SubZi[1:]-SubZi[:-1])**2))
			if(k == 0):
				k += 1
				dNext = sqrt((xis[ii[k]] - SubRi[-1])**2 + (zis[ii[k]] -  SubRi[-1])**2)
				if(dNext < 1e-3): InPlasma = False
			else:
				InPlasma = False


		ChanRi.append(Ri)
		ChanZi.append(Zi)
		ChanTi.append(Ti)
		ChanLeni.append(Leni)

	path = types.SimpleNamespace()
	path.Leni = ChanLeni
	path.Ri	  = ChanRi
	path.Zi	  = ChanZi
	path.Ti	  = ChanTi

	return 	path		


def get_triangles_on_line(RKnots, ZKnots, TriKnots, TriNeigh, TriTypeneigh, OldTriangle, RLs, ZLs, RLe, ZLe):

	OtherSides  = np.array([[1,2],[0,2],[0,1]])
	Vertex	    = np.array([[0,1],[1,2],[2,0]])

	OldSide =np.where(TriNeigh[OldTriangle,:] < 0)[0]
	if(len(OldSide) != 1):
		print("get_triangles_on_line: ERROR in starting triangle")
		print("\tOldTriangle            =",OldTriangle)
		print("\tTriNeigh[OldTriangle,:]=",TriNeigh[OldTriangle,:])
		exit()
	else:
		OldSide = OldSide[0]
		
	SRi, SZi, Found = intersect_line_to_segment(RKnots[TriKnots[OldTriangle,Vertex[OldSide,0]]], ZKnots[TriKnots[OldTriangle,Vertex[OldSide,0]]],
												RKnots[TriKnots[OldTriangle,Vertex[OldSide,1]]], ZKnots[TriKnots[OldTriangle,Vertex[OldSide,1]]],
												RLs, ZLs, RLe, ZLe)

	Ri = np.array([SRi])					#Array intersection on triangles
	Zi = np.array([SZi])
	Ti = np.array([OldTriangle])			#Array of triangles
	Si = np.array([OldSide])				#Array sides of triangles

	ContinueSearch = True
	while (ContinueSearch):
		for iSide in OtherSides[OldSide,:]:
			SRi,SZi,Found = intersect_line_to_segment(RKnots[TriKnots[OldTriangle,Vertex[iSide,0]]], ZKnots[TriKnots[OldTriangle,Vertex[iSide,0]]],
													  RKnots[TriKnots[OldTriangle,Vertex[iSide,1]]], ZKnots[TriKnots[OldTriangle,Vertex[iSide,1]]],
													  RLs, ZLs, RLe, ZLe)
			
			if(Found):
				Ri = np.append(Ri, SRi)
				Zi = np.append(Zi, SZi)
				if(TriNeigh[OldTriangle,iSide] >= 0):
					OldSide     = TriTypeneigh[OldTriangle,iSide]		#Side of neighbour triangle
					OldTriangle	= TriNeigh[OldTriangle,iSide]			#Number of neighbour triangle
				else:													#found boundary
					OldSide 		= iSide
					ContinueSearch	= False

				Si = np.append(Si, OldSide)
				Ti = np.append(Ti, OldTriangle)

				break

		if(not Found):
			print("\tERROR: intersection not found before exit from triangles region")
			print("\tRLs, ZLs, RLe, ZLe                  =",RLs, ZLs, RLe, ZLe)
			for iSide in OtherSides[OldSide,:]:
				print("\tSide                                =",iSide+1)
				print("\tRKnots1, ZKnots1,RKnots2, ZKnots2   =",RKnots[TriKnots[OldTriangle,Vertex[iSide,0]]], ZKnots[TriKnots[OldTriangle,Vertex[iSide,0]]],
													  RKnots[TriKnots[OldTriangle,Vertex[iSide,1]]], ZKnots[TriKnots[OldTriangle,Vertex[iSide,1]]])

			print("\tLast Ri =",Ri[-min(len(Ri),4):])
			print("\tLast Zi =",Zi[-min(len(Zi),4):])
			ContinueSearch	= False
#			exit()

	Ki		 = np.empty((len(Ri),2), dtype = 'i4')
	WKi		 = np.empty((len(Ri),2), dtype = 'f8')
	Ki[:,0]	 = TriKnots[Ti,Vertex[Si,0]]
	Ki[:,1]  = TriKnots[Ti,Vertex[Si,1]]
	WKi[:,1] = np.sqrt((Ri - RKnots[Ki[:,0]])**2			 + (Zi - ZKnots[Ki[:,0]])**2) / \
					  np.sqrt((RKnots[Ki[:,1]]-RKnots[Ki[:,0]])**2 + (ZKnots[Ki[:,1]]-ZKnots[Ki[:,0]])**2)	
	WKi[:,0] = 1. - WKi[:,1]

	return Ri, Zi, Ti, Si, Ki, WKi
