import types
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot			as pyp
import numpy 						as np
from routines.show_no_intersection	import show_no_intersection
from routines.contour_better		import contour_better
from routines.intersect_contour		import intersect_contour
from routines.globals				import DEBUG

def getArcLim(Root, Config):

	if(DEBUG > 0): print("getArcLim")
	
	X_points  = Config.X_points
	nX_points = len(X_points)
	for k in range(nX_points):
		for n in range(4):
			c1			= types.SimpleNamespace()
			c1.arc		= [types.SimpleNamespace()]
			c1.num		= 1
			c1.arc[0].x	= X_points[k].cut[n].R
			c1.arc[0].y	= X_points[k].cut[n].Z
			c2			= contour_better(Config.r2D, Config.z2D, Config.flux2D, [X_points[k].cut[n].psilim])
			X			= intersect_contour(c1,c2)

			if(len(X) == 0):																	#If no intersection try to extend cut
				c1.arc[0].x = np.append(c1.arc[0].x, 2*c1.arc[0].x[-1]-c1.arc[0].x[-2])
				c1.arc[0].y = np.append(c1.arc[0].y, 2*c1.arc[0].y[-1]-c1.arc[0].y[-2])
				X  = intersect_contour(c1,c2)

			if(len(X) == 0):
				print("\tnX_points, k_point, n_cut=",nX_points, k+1, n+1)
				print("\tX_points[k].cut[n].psilim=",X_points[k].cut[n].psilim)
				show_no_intersection(Root, Config, c1, c2)
				return False
			else:
				X_points[k].cut[n].psiR = c2.arc[X[0].arc2].x
				X_points[k].cut[n].psiZ = c2.arc[X[0].arc2].y

	if(DEBUG > 0): print("getArcLim: Completed")
	return True
