import types
import numpy						as np
import scipy.interpolate			as interp
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot 			as pyp
from routines.intersect_contour		import intersect_contour
from routines.contour_better		import contour_better
from routines.globals				import CIRCLE_RADIUS, DEBUG
from interfaces.ask_for_Xpoint		import ask_for_Xpoint

def surround_Xpoints(Root, Config, CheckPlot=False, RadArroundXp=CIRCLE_RADIUS):

	if(DEBUG > 0): print("surround_Xpoints")

#	r = 2*(Config.r2D[0,1]-Config.r2D[0,0])
	theta = np.linspace(0.,2.*np.pi,100)

	X_points = Config.X_points
	
	f_psi	= interp.RectBivariateSpline(Config.r2D[0,:], Config.z2D[:,0], Config.flux2D.T)										#RectBivariateSpline wants [x,y] order

	for k in range(len(X_points)):
		if(RadArroundXp > 0.):
			r = RadArroundXp
		else:
			r = -RadArroundXp*X_points[k].R

		c1			= contour_better(Config.r2D, Config.z2D, Config.flux2D, [X_points[k].psi])
		c2			= types.SimpleNamespace()
		c2.arc		= [types.SimpleNamespace()]
		c2.num		= 1
		c2.arc[0].x = X_points[k].R+r*np.cos(theta)
		c2.arc[0].y = X_points[k].Z+r*np.sin(theta)
		I			= intersect_contour(c1,c2)
		arg			= np.zeros(4)
		if(len(I)==0):
#			print("r=",r)
			point   = types.SimpleNamespace()
			point.x = X_points[k].R 
			point.y = X_points[k].Z
			choiche = ask_for_Xpoint(Root, Config, X_points[k].psi, point)
		for n in range(4):
			cplx	= (I[n].x-X_points[k].R)+complex(0,1)*(I[n].y-X_points[k].Z)
			arg[n]	= np.angle(cplx)


		if(DEBUG > 1):
			print("\n\tX_points = ",k+1," R=",X_points[k].R," Z=",X_points[k].Z)
			print("\tangles   = ",arg*180./np.pi)

		b = np.argsort(arg)
		X_points[k].branch = []
		for n in range(4):
			X_points[k].branch.append(types.SimpleNamespace())
			X_points[k].branch[n].R		= I[b[n]].x
			X_points[k].branch[n].Z		= I[b[n]].y
			X_points[k].branch[n].theta = arg[b[n]]
			X_points[k].branch[n].arc	= I[b[n]].arc1
			if(DEBUG > 1):
				print("\tX_points[{:d}].branch[{:d}].theta = {:f} ".format(k+1,n+1,X_points[k].branch[n].theta*180./np.pi))

		B = np.roll(np.arange(4),-1)
		X_points[k].cut = []
		for n in range(4):
			X_points[k].cut.append(types.SimpleNamespace())
			X_points[k].cut[n].theta	= (X_points[k].branch[n].theta + X_points[k].branch[B[n]].theta)/2
			if(np.abs(X_points[k].cut[n].theta - X_points[k].branch[n].theta) > np.pi/2):
				X_points[k].cut[n].theta = X_points[k].cut[n].theta + np.pi
																						# choose the small angle
			X_points[k].cut[n].Rs		= X_points[k].R+r*np.cos(X_points[k].cut[n].theta)
			X_points[k].cut[n].Zs		= X_points[k].Z+r*np.sin(X_points[k].cut[n].theta)
			X_points[k].cut[n].psi		= float(f_psi.ev(X_points[k].cut[n].Rs, X_points[k].cut[n].Zs))								#float to convert from numpy.ndarray class
			X_points[k].cut[n].psilim   = 0.
			if(DEBUG > 1):
				print("\tX_points[{:d}].cut[{:d}].theta={:f} ".format(k+1,n+1,X_points[k].cut[n].theta*180./np.pi))
				print("\tX_points[{:d}].cut[{:d}].Rs-R={:f}, X_points[{:d}].cut[{:d}].Zs-Z = {:f}".format(k+1,n+1,X_points[k].cut[n].Rs-X_points[k].R,k+1,n+1,X_points[k].cut[n].Zs-X_points[k].Z))


		if(X_points[k].cut[0].psi > X_points[k].cut[1].psi):							# reorder
			if(DEBUG > 1):
				print("\tReorder cuts")

			R_		= X_points[k].branch[3].R
			Z_		= X_points[k].branch[3].Z
			theta_	= X_points[k].branch[3].theta
			for i in range(2,-1,-1):
				X_points[k].branch[i+1].R		= X_points[k].branch[i].R
				X_points[k].branch[i+1].Z		= X_points[k].branch[i].Z
				X_points[k].branch[i+1].theta	= X_points[k].branch[i].theta

			X_points[k].branch[0].R		= R_
			X_points[k].branch[0].Z		= Z_
			X_points[k].branch[0].theta	= theta_
			R_							= X_points[k].cut[3].Rs
			Z_							= X_points[k].cut[3].Zs
			psi_						= X_points[k].cut[3].psi
			theta_						= X_points[k].cut[3].theta
			for i in range(2,-1,-1):
				X_points[k].cut[i+1].Rs		= X_points[k].cut[i].Rs
				X_points[k].cut[i+1].Zs		= X_points[k].cut[i].Zs
				X_points[k].cut[i+1].theta	= X_points[k].cut[i].theta
				X_points[k].cut[i+1].psi	= X_points[k].cut[i].psi

			X_points[k].cut[0].Rs		= R_
			X_points[k].cut[0].Zs		= Z_
			X_points[k].cut[0].theta	= theta_
			X_points[k].cut[0].psi		= psi_

		if(CheckPlot):
			for n in range(len(c1.arc)): pyp.plot(c1.arc[n].x, c1.arc[n].y,'b-')
			pyp.plot(c2.arc[0].x, c2.arc[0].y,'r-')
	
			pyp.text(X_points[k].branch[0].R,	X_points[k].branch[0].Z,'1')
			pyp.text(X_points[k].branch[1].R,	X_points[k].branch[1].Z,'2')
			pyp.text(X_points[k].branch[2].R,	X_points[k].branch[2].Z,'3')
			pyp.text(X_points[k].branch[3].R,	X_points[k].branch[3].Z,'4')
			pyp.text(X_points[k].cut[0].Rs,		X_points[k].cut[0].Zs,	'A')
			pyp.text(X_points[k].cut[1].Rs,		X_points[k].cut[1].Zs,	'B')
			pyp.text(X_points[k].cut[2].Rs,		X_points[k].cut[2].Zs,	'C')
			pyp.text(X_points[k].cut[3].Rs,		X_points[k].cut[3].Zs,	'D')

			pyp.show(block=False)

	if(DEBUG > 0): print("surround_Xpoints: Completed")
