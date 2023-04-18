import numpy 					as np
from interfaces.progressbar 	import ProgressBar
from mesh.follow_grad			import follow_grad
from routines.globals			import DEBUG, STEP_FOLLOW_GRAD

def plot_ortho_lines(Root, Config, step=STEP_FOLLOW_GRAD):

	if(DEBUG > 0): print("plot_ortho_lines")

	prog_bar = ProgressBar(Root, Label="Please wait: Xpoint number 1", Value=0.) 
	
	X_points	= Config.X_points
	nX_points	= len(X_points)
	Nsteps		= 4*nX_points
	for k in range(nX_points):

		for n in range(4):
																			# start point
			R_ 		= X_points[k].cut[n].Rs
			Z_		= X_points[k].cut[n].Zs
			psi_	= X_points[k].cut[n].psi
			psiX	= X_points[k].psi
			if(psi_ < psiX):												# psi descendant --> psimin = psicore
				Rv	= np.append(X_points[k].R, R_)
				Zv	= np.append(X_points[k].Z, Z_)
			else:															# psi montant --> psimax = psiout
				Rv	= np.append(X_points[k].R, R_)
				Zv	= np.append(X_points[k].Z, Z_)


			X_points[k].cut[n].R		= np.copy(Rv)
			X_points[k].cut[n].Z		= np.copy(Zv)
			X_points[k].cut[n].type		= 1

		psilim		= np.empty(4)
		psilim[0]	= X_points[k].cut[0].psilim
		psilim[1]	= Config.psiout
		psilim[2]	= X_points[k].cut[2].psilim
		psilim[3]	= Config.psiout
		for n in range(4):
			if(DEBUG > 1): 
				print("\tnX_points, k_point, n_cut=",nX_points, k+1, n+1)
				print("\tpsilim=",psilim[n])

			R, Z =	follow_grad(Config, X_points[k].cut[n].R[-1], X_points[k].cut[n].Z[-1], psilim[n], step=step)
			X_points[k].cut[n].R	= np.append(X_points[k].cut[n].R, R[1:])
			X_points[k].cut[n].Z	= np.append(X_points[k].cut[n].Z, Z[1:])
			Nstep					= (k-1)*4+n

		prog_bar.Update(Nstep/Nsteps, Label="Please wait: Xpoint number {:d}".format(k))

	prog_bar = 0

	if(DEBUG > 0): print("plot_ortho_lines:Completed")

	return True
