import types
import numpy							as np
from tkinter							import messagebox
import scipy.interpolate				as interp
from routines.contour_better			import contour_better
from routines.intersect_contour			import intersect_contour
from routines.globals					import CIRCLE_RADIUS, MAX_N_XPOINTS, DEBUG, ASK_XPOINTS
from interfaces.ask_for_Xpoint			import ask_for_Xpoint

def find_Xpoints(Root, Config, Fig=0):

	if(DEBUG > 0): print("find_Xpoints")

	theta = np.linspace(0.,2.*np.pi,100)

#	 compute first derivatives
	dp_dR	= (Config.flux2D[1:-1,2:] - Config.flux2D[1:-1,:-2])/(Config.r2D[1:-1,2:] - Config.r2D[1:-1,:-2])
	dp_dZ	= (Config.flux2D[2:,1:-1] - Config.flux2D[:-2,1:-1])/(Config.z2D[2:,1:-1] - Config.z2D[:-2,1:-1])

	rred	= Config.r2D[1:-1,1:-1]
	zred	= Config.z2D[1:-1,1:-1]
	d2p_dR2	=(dp_dR[1:-1,2:]-dp_dR[1:-1,:-2])/(rred[1:-1,2:]-rred[1:-1,:-2])
	d2p_dZ2	=(dp_dZ[2:,1:-1]-dp_dZ[:-2,1:-1])/(zred[2:,1:-1]-zred[:-2,1:-1])

#	find contour lines where the first derivatives are zero
	cx	= contour_better(Config.r2D[1:-1,1:-1], Config.z2D[1:-1,1:-1], dp_dR, [0.])
	cy  = contour_better(Config.r2D[1:-1,1:-1], Config.z2D[1:-1,1:-1], dp_dZ, [0.])

	f_psi	  = interp.RectBivariateSpline(Config.r2D[0,:], Config.z2D[:,0], Config.flux2D.T) #RectBivariateSpline wants [x,y] order
	f_d2p_dR2 = interp.RectBivariateSpline(rred[1,1:-1], zred[1:-1,1], d2p_dR2.T)
	f_d2p_dZ2 = interp.RectBivariateSpline(rred[1,1:-1], zred[1:-1,1], d2p_dZ2.T)

	c2			= types.SimpleNamespace()
	c2.arc		= [types.SimpleNamespace()]
	c2.num		= 1

	X_points  = []
	x_cxcy = intersect_contour(cx,cy)

	k			= 0
	XSearch		= True
	while (XSearch): 
		
		if(CIRCLE_RADIUS > 0.):
			r = CIRCLE_RADIUS
		else:
			r = -CIRCLE_RADIUS*x_cxcy[k].x
		
		psiguess	= float(f_psi.ev(x_cxcy[k].x, x_cxcy[k].y))								#float to convert from numpy.ndarray class 
		c1			= contour_better(Config.r2D, Config.z2D, Config.flux2D, [psiguess])
		c2.arc[0].x = x_cxcy[k].x+r*np.cos(theta)
		c2.arc[0].y = x_cxcy[k].y+r*np.sin(theta)
		I			= intersect_contour(c1,c2)
		if(len(I) == 4):
			d2p_dR2_ = float(f_d2p_dR2.ev(x_cxcy[k].x, x_cxcy[k].y))						#float to convert from numpy.ndarray class
			d2p_dZ2_ = float(f_d2p_dZ2.ev(x_cxcy[k].x, x_cxcy[k].y))
		
#			threshold for merdiques cases
			th = np.abs(np.max(Config.flux2D) - np.min(Config.flux2D))/np.abs(np.max(Config.r2D) - np.min(Config.r2D))**2*0.6
			if((np.abs(d2p_dR2_) < th) or (np.abs(d2p_dZ2_) < th)):			#manual
				
				if(ASK_XPOINTS == 1):
					choiche = ask_for_Xpoint(Root, Config, psiguess, x_cxcy[k])
				else:
					choiche = 'X'

				if(choiche == 'X'):
					X_points.append(types.SimpleNamespace())
					X_points[-1].R		= x_cxcy[k].x
					X_points[-1].Z		= x_cxcy[k].y
					X_points[-1].psi	= float(psiguess)
					if(Fig != 0): Fig.plot(x_cxcy[k].x, x_cxcy[k].y, 'bo')
				if(choiche == 'O'):
					print("Point rejected")
			else:															#auto
	
				if(d2p_dR2_*d2p_dZ2_ < 0):
					X_points.append(types.SimpleNamespace())
					X_points[-1].R		= x_cxcy[k].x
					X_points[-1].Z		= x_cxcy[k].y
					X_points[-1].psi	= psiguess
					if(Fig != 0): Fig.plot(x_cxcy[k].x, x_cxcy[k].y, 'bo')
		k += 1
		if(k >= len(x_cxcy)):
			XSearch = False
		elif(len(X_points) >= MAX_N_XPOINTS):
			XSearch = False
			messagebox.showwarning("Grid-Gen", "Reached max number of X point = {:d}".format(MAX_N_XPOINTS))


	psi = np.array([X_points[k].psi for k in range(len(X_points))])
	psi	= np.sort(psi)
	for k in range(len(X_points)):
		c = np.where(psi == X_points[k].psi)[0][0]
		X_points[k].index = c
#		X_points[k].cut[1:4].psilim  = np.zeros(4, dtype='f8')     #Definito in surround_Xpoints.py

	Config.X_points = X_points

	if(DEBUG > 0): print("find_Xpoints: Completed")
