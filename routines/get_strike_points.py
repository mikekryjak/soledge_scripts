import types
import numpy								as np
import scipy.interpolate					as interp
from routines.intersect_contour				import intersect_contour
from routines.contour_better				import contour_better
from routines.utils_walls					import intersect_walls
from routines.globals						import DEBUG
from routines.utils_walls					import walls_define_path

from mesh.get_rz_core_sep					import get_rz_core_sep

#==============================================================================
# return strike points position
#==============================================================================

def get_strike_points(Config):

	Rcore, Zcore, CoreMegazone = get_rz_core_sep(Config, core_and_sep = False)

	R0 		= 0.5*(Rcore.min() + Rcore.max())

	iMax = np.argmax(Rcore)
	Z0 = Zcore[iMax]

	X_points	= Config.X_points
	nX_points	= len(X_points)

	if(hasattr(Config,"MagPMegazones")):
#		Find strike point where mesh is alligned with the wall

		XPointsPsi = np.empty(len(X_points), dtype='f8')
		for i in range(len(X_points)): XPointsPsi[i] = X_points[i].psi

		MagPMegazones = Config.MagPMegazones
		iXPoints	= np.empty(len(MagPMegazones), dtype = 'i4')
		StrikesR	= np.empty(len(MagPMegazones), dtype='f8')
		StrikesZ	= np.empty(len(MagPMegazones), dtype='f8')

		nStrikes = 0
		f_psi   = interp.RectBivariateSpline(Config.r2D[0,:], Config.z2D[:,0], Config.flux2D.T)								#RectBivariateSpline wants [x,y] order
		for k in range(len(MagPMegazones)):
			if(MagPMegazones[k].isaligned):
				if((MagPMegazones[k].subrefpoints[0].R[-1] == MagPMegazones[k].subrefpoints[1].R[0]) and (MagPMegazones[k].subrefpoints[0].Z[-1] == MagPMegazones[k].subrefpoints[1].Z[0])):
					StrikesR[nStrikes] = MagPMegazones[k].subrefpoints[0].R[-1]
					StrikesZ[nStrikes] = MagPMegazones[k].subrefpoints[0].Z[-1]		
				else:
					StrikesR[nStrikes] = MagPMegazones[k].subrefpoints[0].R[0]
					StrikesZ[nStrikes] = MagPMegazones[k].subrefpoints[0].Z[0]		

				StrikePsi = f_psi.ev(StrikesR[nStrikes], StrikesZ[nStrikes])
				iXPoints[nStrikes] = np.argmin(np.abs(StrikePsi-XPointsPsi))
				nStrikes += 1

		iXPoints	= iXPoints[:nStrikes] 
		StrikesR	= StrikesR[:nStrikes]
		StrikesZ	= StrikesZ[:nStrikes]

		nInXPoints = int(nStrikes/2)
		if(nInXPoints*2 != nStrikes):
			print("\t\tERROR: Odd number of alligned strike points")
			exit()

		iXSort = np.argsort(iXPoints)
		InXPoints	= iXPoints[iXSort]

		InXPoints	= np.empty(nInXPoints, dtype='i4')
		RStrikes	= np.empty((nInXPoints,2), dtype='f8')
		ZStrikes	= np.empty((nInXPoints,2), dtype='f8')
		for i in range(nInXPoints):
			if(iXPoints[iXSort[2*i]] != iXPoints[iXSort[2*i+1]]):
				print("\t\tERROR: Only one strike point for an X-point")
				exit()

			InXPoints[i]	= iXPoints[iXSort[2*i]]
			for j in range(2):
				RStrikes[i,j]	= StrikesR[iXSort[2*i+j]]
				ZStrikes[i,j]	= StrikesZ[iXSort[2*i+j]]
	else:

		InXPoints	= np.array([0])
		if(nX_points > 1):
			walls_define_path(Config)				#Find x points inside wall

			PXPoints  = np.array([[X_points[k].R for k in range(nX_points)], [X_points[k].Z for k in range(nX_points)]]).T
			InXPoints = np.where(Config.WallPath.contains_points(PXPoints))[0]

			nInXPoints = len(InXPoints)

			RStrikes	= np.empty((nInXPoints,2), dtype='f8')
			ZStrikes	= np.empty((nInXPoints,2), dtype='f8')

#		Find strike points on wall

		if(nInXPoints > 1):

#			Find sequence separatrixes from internal to external
			cEq 	   		= types.SimpleNamespace()
			cEq.arc 		= [types.SimpleNamespace()]
			cEq.arc[0].x	= np.array([R0, R0+100.])
			cEq.arc[0].y	= np.array([Z0, Z0])

			dSep		= np.zeros(nInXPoints, dtype='f8')
			Xint		= []
			for k in range(nInXPoints):
				cSep = contour_better(Config.r2D, Config.z2D, Config.flux2D, [Config.X_points[InXPoints[k]].psi])
				X	 = intersect_contour(cEq,cSep)
				dSep[k]	= X[0].x - R0

				Xint.append(intersect_walls(Config,cSep))					#interserction with wall

			iSep	 	= np.argsort(dSep)
			dSep	 	= dSep[iSep]
			dSep	   -= dSep[0]
			InXPoints	= InXPoints[iSep]
			Xint		= Xint[iSep]
			RStrikes	= np.empty((nInXPoints,4), dtype='f8')
			ZStrikes	= np.empty((nInXPoints,4), dtype='f8')
			nMaxInt		= 2
			for k in range(nInXPoints):
				if(len(Xint[k]) > nMaxInt):
					Rint = np.array([Xint[k][i].x for i in range(len(Xint[k]))])
					Zint = np.array([Xint[k][i].y for i in range(len(Xint[k]))])
					ClosestToXpts = np.argsort((Rint - Config.X_points[InXPoints[k]].R)**2 + (Zint - Config.X_points[InXPoints[k]].Z)**2)	#find closest to X-points
					RStrikes[k,:nMaxInt] = Rint[ClosestToXpts[:nMaxInt]]
					ZStrikes[k,:nMaxInt] = Rint[ClosestToXpts[:nMaxInt]]			
				else:
					for i in range(nMaxInt):
						RStrikes[k,i] = Xint[k][i].x
						ZStrikes[k,i] = Xint[k][i].y

	#		Find strike points far from previous

			for k in range(1,nInXPoints):
				dMin = np.empty(4, dtype='f8')
				for i in range(4):
					dMin[i] = min(sqrt((RStrikes[k,i] - RStrikes[k-1,0])**2 + (ZStrikes[k,i] - ZStrikes[k-1,0])**2),
								  sqrt((RStrikes[k,i] - RStrikes[k-1,1])**2 + (ZStrikes[k,i] - ZStrikes[k-1,1])**2))
				iMax = np.argsort(dMin)[::-1]
				RStrikes[k,:] = RStrikes[k,iMax]
				ZStrikes[k,:] = ZStrikes[k,iMax]

			RStrikes = RStrikes[:,:2]									#Keep first two intersections with wall
			ZStrikes = ZStrikes[:,:2]
		else:
			RStrikes	= np.empty((1,2), dtype='f8')
			ZStrikes	= np.empty((1,2), dtype='f8')
			cSep = contour_better(Config.r2D, Config.z2D, Config.flux2D, [Config.X_points[InXPoints[0]].psi])
			X	 = intersect_walls(Config,cSep)							#interserction with walls

			if(len(X) > 2):
				Rint = np.array([X[i].x for i in range(len(X))])
				Zint = np.array([X[i].y for i in range(len(X))])
				ClosestToXpts = np.argsort((Rint - Config.X_points[InXPoints[k]].R)**2 + (Zint - Config.X_points[InXPoints[k]].Z)**2)	#find closest to X-points
				RStrikes[0,:] = Rint[ClosestToXpts[:2]]
				ZStrikes[0,:] = Rint[ClosestToXpts[:2]]			
			else:
				for i in range(2):
					RStrikes[0,i] = X[i].x
					ZStrikes[0,i] = X[i].y

	return RStrikes, ZStrikes, InXPoints



