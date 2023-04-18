import types
import numpy						as np
import numpy.matlib					as mat
from math							import sqrt
from matplotlib.path				import Path
from routines.intersect_segments	import intersect_segments, intersect_lines_to_segments
from routines.intersect_contour		import intersect_2contours
from routines.globals				import DEBUG

def get_chords_mesh_path(Config, Diag):

	Diag.path = get_lines_mesh_path(Config, Diag.lines_start, Diag.lines_end)

	return

#===================================================================

def get_lines_mesh_path(Config, LinesSt, LinesEnd):

	Zones  = Config.Zones
	nZones = len(Zones)

#	Extend lines to be sure to cross the wall
	
	RWmin = np.min(Config.Walls[0].Rwall)
	RWmax = np.max(Config.Walls[0].Rwall)
	ZWmin = np.min(Config.Walls[0].Zwall)
	ZWmax = np.max(Config.Walls[0].Zwall)
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
		ELinesSt[2,i]  = zis[0]
		ELinesEnd[2,i] = zis[1]

#	Define path to find point outside wall

	WallPath = Config.Walls[0].WallPath

	RsW	= Config.Walls[0].Rwall[:-1]
	ZsW = Config.Walls[0].Zwall[:-1]
	ReW	= Config.Walls[0].Rwall[1:]
	ZeW	= Config.Walls[0].Zwall[1:]

	isRZ   = [0,0,1,1,0]
	jsRZ   = [0,1,1,0,0]
	iSides = [0,1,2,3,0]

#	Find intersections between lines and mesh

	Xis = []
	Zis = []
	ICh = []
	Icm = []
	Jcm = []
	Ism = []
	InWall = []
	for k in range(nZones):
		nCells = Zones[k].Chi.size

		nSegs = 4*nCells
		Rs = np.empty(nSegs, dtype='f8')
		Zs = np.empty(nSegs, dtype='f8')
		Re = np.empty(nSegs, dtype='f8')
		Ze = np.empty(nSegs, dtype='f8')
		ism = np.empty(nSegs, dtype='i4')
		icm = np.empty(nSegs, dtype='i4')
		jcm = np.empty(nSegs, dtype='i4')
		Nx	   = Zones[k].Chi.shape[0]
		Nz	   = Zones[k].Chi.shape[1]
		iBase  = 0
		for iSide in range(4):
			Rs[iBase:iBase+nCells] = Zones[k].gridR[isRZ[iSide]:Nx+isRZ[iSide], jsRZ[iSide]:Nz+jsRZ[iSide]].reshape(-1)
			Zs[iBase:iBase+nCells] = Zones[k].gridZ[isRZ[iSide]:Nx+isRZ[iSide], jsRZ[iSide]:Nz+jsRZ[iSide]].reshape(-1)

			Re[iBase:iBase+nCells] = Zones[k].gridR[isRZ[iSide+1]:Nx+isRZ[iSide+1], jsRZ[iSide+1]:Nz+jsRZ[iSide+1]].reshape(-1)
			Ze[iBase:iBase+nCells] = Zones[k].gridZ[isRZ[iSide+1]:Nx+isRZ[iSide+1], jsRZ[iSide+1]:Nz+jsRZ[iSide+1]].reshape(-1)
				
			icm[iBase:iBase+nCells] = mat.repmat(np.arange(Nx), Nz, 1).T.reshape(-1)
			jcm[iBase:iBase+nCells] = mat.repmat(np.arange(Nz), Nx, 1).reshape(-1)
			ism[iBase:iBase+nCells] = iSide
			iBase += nCells

		xis, zis, ich, ims = intersect_segments(ELinesSt[0,:], ELinesSt[2,:], ELinesEnd[0,:], ELinesEnd[2,:], Rs, Zs, Re, Ze, remove_same = False)
		Xis.append(xis)
		Zis.append(zis)
		ICh.append(ich)
		Icm.append(icm[ims])
		Jcm.append(jcm[ims])
		Ism.append(ism[ims])
		InWall.append(np.where(WallPath.contains_points(np.array([xis,zis]).T), 1, -1))
		
#	Assign intersection to mesh

	LenCells = []
	kCells   = []
	iCells   = []
	jCells   = []
	XiCells	 = []
	ZiCells	 = []
	for nChan in range(LinesSt.shape[1]):
		kCells.append(np.empty(0, dtype = 'i4'))
		iCells.append(np.empty(0, dtype = 'i4'))
		jCells.append(np.empty(0, dtype = 'i4'))
		XiCells1 = np.empty(0, dtype = 'f8')
		XiCells2 = np.empty(0, dtype = 'f8')
		ZiCells1 = np.empty(0, dtype = 'f8')
		ZiCells2 = np.empty(0, dtype = 'f8')

		XCross   = np.empty(0, dtype = 'f8')
		ZCross   = np.empty(0, dtype = 'f8')
		for k in range(nZones):
			ii = np.where(ICh[k] == nChan)[0]
			if(len(ii) > 0):
				Zones[k].xci = np.zeros((Zones[k].Nx, Zones[k].Nz, 4), dtype='f8')
				Zones[k].zci = np.zeros((Zones[k].Nx, Zones[k].Nz, 4), dtype='f8')
				Zones[k].fci = np.zeros((Zones[k].Nx, Zones[k].Nz, 4), dtype='i4')

				Zones[k].xci[Icm[k][ii],Jcm[k][ii],Ism[k][ii]] = Xis[k][ii]
				Zones[k].zci[Icm[k][ii],Jcm[k][ii],Ism[k][ii]] = Zis[k][ii]
				Zones[k].fci[Icm[k][ii],Jcm[k][ii],Ism[k][ii]] = InWall[k][ii]
					
				SumFci = np.sum(np.abs(Zones[k].fci), axis = 2)

				iCross, jCross = np.where(SumFci == 1)					#remove single crossing
				if(len(iCross) > 0):
					Zones[k].xci[iCross, jCross, :] = 0.
					Zones[k].zci[iCross, jCross, :] = 0.
					Zones[k].fci[iCross, jCross, :] = 0
				
				iCross, jCross = np.where(SumFci > 2)					#remove corner crossing
				if(len(iCross) > 0):
					for ii in range(len(iCross)):
						i = iCross[ii]
						j = jCross[ii]
						for iSide in range(4):
							if((Zones[k].fci[i, j, iSide] != 0) and (Zones[k].fci[i, j, iSides[iSide]] != 0)):
								dist = sqrt((Zones[k].xci[i,j,iSide] - Zones[k].xci[i,j,iSides[iSide]])**2 + \
											(Zones[k].zci[i,j,iSide] - Zones[k].zci[i,j,iSides[iSide]])**2)
								if(dist < 2.e-5):
									Zones[k].xci[iCross, jCross, iSide] = 0.
									Zones[k].zci[iCross, jCross, iSide] = 0.
									Zones[k].fci[iCross, jCross, iSide] = 0

				SumFci = np.sum(np.abs(Zones[k].fci), axis = 2)						# check again
				iCross, jCross = np.where((SumFci != 0) & (SumFci != 2))							
				if(len(iCross) > 0):
					print("\tUnable to find only two intersection on some quadrangles")
					exit()

				iCross, jCross = np.where(np.sum(Zones[k].fci, axis = 2) == -2)		# find cells with both points outside wall
				if(len(iCross) > 0):
					Zones[k].xci[iCross, jCross, :] = 0.
					Zones[k].zci[iCross, jCross, :] = 0.
					Zones[k].fci[iCross, jCross, :] = 0
					iOut = 0
					
				iCross, jCross, iSideOut = np.where(Zones[k].fci == -1)				# find cells with one point outside wall
				if(len(iCross) > 0):
					for ii in range(len(iCross)):
						i = iCross[ii]
						j = jCross[ii]
						iSideIn = np.where(Zones[k].fci[i,j,:] == 1)				# find point inside
						xs = np.array([Zones[k].xci[i,j,iSideOut[ii]], Zones[k].xci[i,j,iSideIn]])
						zs = np.array([Zones[k].zci[i,j,iSideOut[ii]], Zones[k].zci[i,j,iSideIn]])
						
						xis, zis, is1, is2 = intersect_2contours(xs, zs, Config.Walls[0].Rwall, Config.Walls[0].Zwall)

						XCross = np.append(XCross, xis[0])
						ZCross = np.append(ZCross, zis[0])

						Zones[k].xci[i,j,iSideOut] = xis[0]
						Zones[k].zci[i,j,iSideOut] = zis[0]
						Zones[k].fci[i,j,iSideOut] = 1

				iCross, jCross = np.where(np.sum(Zones[k].fci, axis = 2) == 2)		
				nSegs = len(iCross)
				i2Sides = np.empty((nSegs, 2), dtype = 'i4')
				i2IndS  = np.zeros(nSegs, dtype = 'i4')
				for iSide in range(4):
					iiS = np.where(Zones[k].fci[iCross, jCross, iSide] == 1)[0]
					i2Sides[iiS, i2IndS[iiS]] = iSide
					i2IndS[iiS]			   += 1

				
				XiCells1 = np.append(XiCells1, Zones[k].xci[iCross,jCross,i2Sides[:,0]])
				XiCells2 = np.append(XiCells2, Zones[k].xci[iCross,jCross,i2Sides[:,1]])
				ZiCells1 = np.append(ZiCells1, Zones[k].zci[iCross,jCross,i2Sides[:,0]])
				ZiCells2 = np.append(ZiCells2, Zones[k].zci[iCross,jCross,i2Sides[:,1]])

				kCells[-1]   = np.append(kCells[-1], np.ones(nSegs, dtype = 'i4')*k)
				iCells[-1]   = np.append(iCells[-1], iCross)
				jCells[-1]   = np.append(jCells[-1], jCross)

		if(len(XiCells1) > 0):
			StartIsIn = WallPath.contains_point(np.array([LinesSt[0,nChan],LinesSt[2,nChan]]))
			XiMid	  = 0.5*(XiCells1 + XiCells2)
			ZiMid	  = 0.5*(ZiCells1 + ZiCells2)
			LenLine	  = sqrt((LinesEnd[0,nChan] - LinesSt[0,nChan])**2 + (LinesEnd[2,nChan] - LinesSt[2,nChan])**2)
			dPi		  = ((XiMid - LinesSt[0,nChan])*(LinesEnd[0,nChan] - LinesSt[0,nChan]) + \
						 (ZiMid - LinesSt[2,nChan])*(LinesEnd[2,nChan] - LinesSt[2,nChan]))/LenLine

			idSort = np.argsort(dPi)
			dPi	   = dPi[idSort]
			XiCells1 = XiCells1[idSort]
			XiCells2 = XiCells2[idSort]
			ZiCells1 = ZiCells1[idSort]
			ZiCells2 = ZiCells2[idSort]
			kCells[-1] = kCells[-1][idSort]
			iCells[-1] = iCells[-1][idSort]
			jCells[-1] = jCells[-1][idSort]
			if(len(XCross) > 2):
#				Remove segmentes after first wall crossing from line start


				dPCross	  = ((XCross - LinesSt[0,nChan])*(LinesEnd[0,nChan] - LinesSt[0,nChan]) + \
							 (ZCross - LinesSt[2,nChan])*(LinesEnd[2,nChan] - LinesSt[2,nChan]))/LenLine
				if(StartIsIn):
					dMin  = np.max(dPCross < 0., dPCross, -1e30)
					dMax  = np.min(dPCross > 0., dPCross,  1e30)
				else:
					dPCross	  = np.sort(dPCross)
					dMin  = dPCross[0]
					dMax  = dPCross[1]

				iInWall  = np.where((dPi > dMin) & (dPi < dMax))[0]
				iOutWall = np.where((dPi < dMin) | (dPi > dMax))[0]
				if(len(iInWall) < len(XiCells1)):
					XiCells1 = XiCells1[iInWall]
					XiCells2 = XiCells2[iInWall]
					ZiCells1 = ZiCells1[iInWall]
					ZiCells2 = ZiCells2[iInWall]
					kCells[-1] = kCells[-1][iInWall]
					iCells[-1] = iCells[-1][iInWall]
					jCells[-1] = jCells[-1][iInWall]

			dYFact = sqrt((LinesEnd[2,nChan] - LinesSt[2,nChan])**2 + LenLine**2)/LenLine
			LenCells.append(np.sqrt((XiCells2-XiCells1)**2 + (ZiCells2-ZiCells1)**2)*dYFact)
			XiCells.append(np.array([XiCells1, XiCells2]).T)	
			ZiCells.append(np.array([ZiCells1, ZiCells2]).T)

			dPi = ((XiCells[-1] - LinesSt[0,nChan])*(LinesEnd[0,nChan] - LinesSt[0,nChan]) + \
				   (ZiCells[-1] - LinesSt[2,nChan])*(LinesEnd[2,nChan] - LinesSt[2,nChan]))/LenLine

			idSort = np.argsort(dPi, axis=1)
			XiCells[-1] = XiCells[-1][:,idSort]
			ZiCells[-1] = ZiCells[-1][:,idSort]
		else:
			LenCells.append(np.empty(0, dtype = 'f8'))
			XiCells.append(np.empty((0,2), dtype = 'f8'))
			ZiCells.append(np.empty((0,2), dtype = 'f8'))

	path = types.SimpleNamespace()
	path.Leni = LenCells
	path.Ki	  = kCells
	path.Ii	  = iCells
	path.Ji	  = jCells
	path.Xi	  = XiCells
	path.Zi	  = ZiCells
	return 	path		
