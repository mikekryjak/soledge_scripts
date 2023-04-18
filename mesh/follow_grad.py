import types
import numpy						as np
import scipy.interpolate			as interp
from routines.globals				import STEP_FOLLOW_GRAD
from routines.contour_better		import contour_better
from routines.intersect_contour		import intersect_contour

def follow_grad(Config, Rstart, Zstart, psi_objective, step=STEP_FOLLOW_GRAD):

#	compute first derivatives

	dp_dR	= (Config.flux2D[1:-1,2:]-Config.flux2D[1:-1,:-2])/(Config.r2D[1:-1,2:]-Config.r2D[1:-1,:-2])
	dp_dZ	= (Config.flux2D[2:,1:-1]-Config.flux2D[:-2,1:-1])/(Config.z2D[2:,1:-1]-Config.z2D[:-2,1:-1])

	rred	= np.copy(Config.r2D[1:-1,1:-1])
	zred	= np.copy(Config.z2D[1:-1,1:-1])

	f_psi	= interp.RectBivariateSpline(Config.r2D[0,:], Config.z2D[:,0], Config.flux2D.T)	#RectBivariateSpline want [x,y] order
	f_dp_dR	= interp.RectBivariateSpline(rred[0,:], zred[:,0], dp_dR.T)
	f_dp_dZ	= interp.RectBivariateSpline(rred[0,:], zred[:,0], dp_dZ.T)

	psistart = float(f_psi.ev(Rstart, Zstart))								#float to convert from numpy.ndarray class
	if(psistart > psi_objective):
		signe = 1
	else:
		signe=-1

	R	= np.array([Rstart])
	Z	= np.array([Zstart])
	psi = psistart


	c1		= types.SimpleNamespace()
	c1.arc	= [types.SimpleNamespace()]
	c2 = contour_better(Config.r2D, Config.z2D, Config.flux2D, [psi_objective])

	cutContour = 0
	while(cutContour == 0):
		dpdR_	= float(f_dp_dR.ev(R[-1], Z[-1]))								#float to convert from numpy.ndarray class
		dpdZ_	= float(f_dp_dZ.ev(R[-1], Z[-1]))								#float to convert from numpy.ndarray class

		Rnew	= R[-1] - signe*dpdR_/np.sqrt(dpdR_**2+dpdZ_**2)*step
		Znew	= Z[-1] - signe*dpdZ_/np.sqrt(dpdR_**2+dpdZ_**2)*step

		R		= np.append(R, Rnew)
		Z		= np.append(Z, Znew)

		psiold	= psi
		psi		= float(f_psi.ev(R[-1], Z[-1]))									#float to convert from numpy.ndarray class
		if(signe*(psi-psi_objective) <= 0.):
#			Checks whether we have actually crossed the iso-contour we are  interested in
#			It might not be the case due to resolution/interpolation errors
#			But it has to be true for the rest of the algorithm so we need to  push the descent until this is true
			c1.arc[0].x = R[-2:]
			c1.arc[0].y = Z[-2:]
			X = intersect_contour(c1,c2)
			if (len(X) > 0): cutContour = 1

		if(signe*(psi-psiold) > 0):
			print("follow_grad: Attention!!")
			print("\tUnable to reach psi=",psi_objective)
			print("\tstarting from psi  =",psistart)
			print("\treached psi=",psi," after ",len(R)-1," steps")
			print("\tTry to increas STEP_FOLLOW_GRAD in globals routine or smoothing")
			break

																				#change sign in descent
	return R, Z
