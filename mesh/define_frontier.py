import types
import numpy					as np
from scipy						import interpolate
from routines.intersect_contour	import intersect_contour

def define_frontier(Config, list_res):
			

	X_points = Config.X_points
	P1		= np.empty((0,3), dtype='f8')
	c1		= types.SimpleNamespace()
	c1.arc	= [types.SimpleNamespace()]
	c2		= types.SimpleNamespace()
	c2.arc	= [types.SimpleNamespace()]
	for k in range(len(X_points)):
		for n in range(4):
			c1.num		= 1
			c1.arc[0].x = X_points[k].cut[n].psiR
			c1.arc[0].y = X_points[k].cut[n].psiZ
			c2.num		= 1
			c2.arc[0].x = list_res[:,0]
			c2.arc[0].y = list_res[:,1]
			inter = intersect_contour(c1,c2)
			if(len(inter) == 1):
#				check if psi surface already tested
				inSurf = False
				m = P1.shape[0]
				for k1 in range(m):
					if(P1[k1,2] == X_points[k].cut[n].psilim): inSurf = True

#				if no doublon
				if(not inSurf): P1 = np.append(P1,[[inter[0].x ,inter[0].y, X_points[k].cut[n].psilim]], axis=0)

	Frontier = types.SimpleNamespace()

	f_flux2D		= interpolate.RectBivariateSpline(Config.r2D[0,:], Config.z2D[:,0], Config.flux2D.T) #RectBivariateSpline want [x,y] order
	psipt			= f_flux2D.ev(list_res[:,0], list_res[:,1])

	b				= np.argsort(psipt)
	Frontier.R	 	= list_res[b,0]
	Frontier.Z	 	= list_res[b,1]
	Frontier.P1	 	= P1
	Frontier.sel	= False
	Frontier.psimin	= np.min(P1[:,2])
	Frontier.psimax = np.max(P1[:,2])

	Frontier.psiA	= float(f_flux2D.ev(Frontier.R[0],  Frontier.Z[0]))								#float to convert from numpy.ndarray class
	Frontier.psiB	= float(f_flux2D.ev(Frontier.R[-1], Frontier.Z[-1]))							#float to convert from numpy.ndarray class
	f_flux2D		= 0
	
	return Frontier
