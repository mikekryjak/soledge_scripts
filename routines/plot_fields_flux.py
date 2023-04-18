
# Function definition is here

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot 	as pyp
import numpy 				as np
import scipy.interpolate	as interpolate
from routines.utils_walls	import plot2d_walls
from routines.globals		import DEBUG

#=========================================================
# This routine plot magnetic fields as provided
#=========================================================
#

def plot_fields_flux(Config, plot_type=0, one_plot=1):

	if((plot_type == 0) and (not hasattr(Config,'r2'))):		return					#fields structure not defined
	if((plot_type == 1) and (not hasattr(Config,'r2D'))):	return					#fields structure not defined
	

	ax	   = []
	fig    = pyp.figure()
	fig.patch.set_facecolor('white')
	for i in range(1,5): ax.append(fig.add_subplot(2,2,i))

	im = []
	for i in range(len(ax)):
		ax[i].set_aspect(1.)
		ax[i].set_xlabel("$R\ (m)$")
		ax[i].set_ylabel("$Z\ (m)$")
#		ax[i].set_title("Fields & flux")
		ax[i].autoscale(enable=True, axis='both', tight=True)
	
	if(plot_type == 0):
		if(DEBUG > 0):print("\tplotting loaded magnetic field")

		if(DEBUG > 1):
			print("Min(Br2)  (T)  ={:8.3f}, Max(Br2)   (T)  ={:8.3f}".format(np.min(Config.Br2),np.max(Config.Br2)))
			print("Min(Bz2)  (T)  ={:8.3f}, Max(Bz2)   (T)  ={:8.3f}".format(np.min(Config.Bz2),np.max(Config.Bz2)))
			print("Min(Bphi2)(T)  ={:8.3f}, Max(Bphi2) (T)  ={:8.3f}".format(np.min(Config.Bphi2),np.max(Config.Bphi2)))
			print("Min(flux2)(Wb) ={:8.3f}, Max(flux2) (Wb) ={:8.3f}".format(np.min(Config.flux2),np.max(Config.flux2)))

		im.append(ax[0].pcolormesh(Config.r2, Config.z2, Config.Br2,   shading='gouraud'))
		ax[0].contour(Config.r2, Config.z2, Config.Br2, colors="r", linestyles='-', linewidths=4, levels=[0.])

		im.append(ax[1].pcolormesh(Config.r2, Config.z2, Config.Bz2,   shading='gouraud'))
		ax[1].contour(Config.r2, Config.z2, Config.Bz2, colors="r", linestyles='-', linewidths=4, levels=[0.])

		im.append(ax[2].pcolormesh(Config.r2, Config.z2, Config.Bphi2, shading='gouraud'))
		im.append(ax[3].pcolormesh(Config.r2, Config.z2, Config.flux2, shading='gouraud'))
		for k in range(len(Config.flux2_x_psi)):
				ax[3].contour(Config.r2, Config.z2, Config.flux2, colors="r", linestyles='-', linewidths=4, levels=[Config.flux2_x_psi[k]])

	elif(plot_type == 1):
		if(DEBUG > 0):print("\tplotting plot_type magnetic field")

		nExt =  len(np.where(Config.r2D[0,:] < Config.r2[0,0])[0])

		if(nExt > 0):
			iS = nExt
			iEr = -nExt; iEz = -nExt
		else:
			iS = 0
			iEr = Config.r2D.shape[0]; iEz = Config.r2D.shape[1]
			
		if(DEBUG > 1):
			print("Min(Br2D)  (T)  ={:8.3f}, Max(Br2D)  (T)  ={:8.3f}".format(np.min(Config.Br2D[iS:iEr,iS:iEz]),np.max(Config.Br2D[iS:iEr,iS:iEz])))
			print("Min(Bz2D)  (T)  ={:8.3f}, Max(Bz2D)  (T)  ={:8.3f}".format(np.min(Config.Bz2D[iS:iEr,iS:iEz]),np.max(Config.Bz2D[iS:iEr,iS:iEz])))
			print("Min(Bphi2D)(T)  ={:8.3f}, Max(Bphi2D)(T)  ={:8.3f}".format(np.min(Config.Bphi2D[iS:iEr,iS:iEz]),np.max(Config.Bphi2D[iS:iEr,iS:iEz])))
			print("Min(flux2D)(Wb) ={:8.3f}, Max(flux2D)(Wb) ={:8.3f}".format(np.min(Config.flux2D[iS:iEr,iS:iEz]),np.max(Config.flux2D[iS:iEr,iS:iEz])))

		im.append(ax[0].pcolormesh(Config.r2D[iS:iEr,iS:iEz], Config.z2D[iS:iEr,iS:iEz], Config.Br2D[iS:iEr,iS:iEz],   shading='gouraud'))
		ax[0].contour(Config.r2D[iS:iEr,iS:iEz], Config.z2D[iS:iEr,iS:iEz], Config.Br2D[iS:iEr,iS:iEz], colors="r", linestyles='-', linewidths=4, levels=[0.])

		im.append(ax[1].pcolormesh(Config.r2D[iS:iEr,iS:iEz], Config.z2D[iS:iEr,iS:iEz], Config.Bz2D[iS:iEr,iS:iEz],   shading='gouraud'))
		ax[1].contour(Config.r2D[iS:iEr,iS:iEz], Config.z2D[iS:iEr,iS:iEz], Config.Bz2D[iS:iEr,iS:iEz], colors="r", linestyles='-', linewidths=4, levels=[0.])
		im.append(ax[2].pcolormesh(Config.r2D[iS:iEr,iS:iEz], Config.z2D[iS:iEr,iS:iEz], Config.Bphi2D[iS:iEr,iS:iEz], shading='gouraud'))
		im.append(ax[3].pcolormesh(Config.r2D[iS:iEr,iS:iEz], Config.z2D[iS:iEr,iS:iEz], Config.flux2D[iS:iEr,iS:iEz], shading='gouraud'))
		for k in range(len(Config.X_points)):
			ax[3].contour(Config.r2D[iS:iEr,iS:iEz], Config.z2D[iS:iEr,iS:iEz], Config.flux2D[iS:iEr,iS:iEz], colors="r", linestyles='-', linewidths=4, levels=[Config.X_points[k].psi])

	elif(plot_type == 2):
		if(DEBUG > 0):print("\tplotting plot_type magnetic difference")

		nExt =  len(np.where(Config.r2D[0,:] < Config.r2[0,0])[0])

		if(nExt > 0):
			iS = nExt
			iEr = -nExt; iEz = -nExt
		else:
			iS = 0
			iEr = Config.r2D.shape[0]; iEz = Config.r2D.shape[1]
			
		rz2 = np.array([Config.r2.reshape(-1),Config.z2.reshape(-1)]).T

		Bdiff = Config.Br2D[iS:iEr,iS:iEz] - interpolate.griddata(rz2, Config.Br2.reshape(-1), (Config.r2D[iS:iEr,iS:iEz], Config.z2D[iS:iEr,iS:iEz]), method='linear')
		im.append(ax[0].pcolormesh(Config.r2D[iS:iEr,iS:iEz], Config.z2D[iS:iEr,iS:iEz], Bdiff,   shading='gouraud'))

		Bdiff = Config.Bz2D[iS:iEr,iS:iEz] - interpolate.griddata(rz2, Config.Bz2.reshape(-1), (Config.r2D[iS:iEr,iS:iEz], Config.z2D[iS:iEr,iS:iEz]), method='linear')
		im.append(ax[1].pcolormesh(Config.r2D[iS:iEr,iS:iEz], Config.z2D[iS:iEr,iS:iEz], Bdiff,   shading='gouraud'))

		Bdiff = Config.Bphi2D[iS:iEr,iS:iEz] - interpolate.griddata(rz2, Config.Bphi2.reshape(-1), (Config.r2D[iS:iEr,iS:iEz], Config.z2D[iS:iEr,iS:iEz]), method='linear')
		im.append(ax[2].pcolormesh(Config.r2D[iS:iEr,iS:iEz], Config.z2D[iS:iEr,iS:iEz], Bdiff, shading='gouraud'))

		Bdiff = Config.flux2D[iS:iEr,iS:iEz] - interpolate.griddata(rz2, Config.flux2.reshape(-1), (Config.r2D[iS:iEr,iS:iEz], Config.z2D[iS:iEr,iS:iEz]), method='linear')
		im.append(ax[3].pcolormesh(Config.r2D[iS:iEr,iS:iEz], Config.z2D[iS:iEr,iS:iEz], Bdiff, shading='gouraud'))

		for k in range(len(Config.flux2_x_psi)):
				ax[3].contour(Config.r2, Config.z2, Config.flux2, colors="r", linestyles='-', linewidths=4, levels=[Config.flux2_x_psi[k]])

		for k in range(len(Config.X_points)):
				ax[3].contour(Config.r2D[iS:iEr,iS:iEz], Config.z2D[iS:iEr,iS:iEz], Config.flux2D[iS:iEr,iS:iEz], colors="b", linestyles='-', linewidths=4, levels=[Config.X_points[k].psi])

	if(hasattr(Config,'Walls')):												# I can plot wall also
		for i in range(len(ax)):
			plot2d_walls(ax[i], Config.Walls)											# plot boundary

	cb = []
	for i in range(len(im)):
		cb.append(fig.colorbar(im[i],   ax=ax[i]))

	labels = ["$B_R\ (T)$", "$B_Z\ (T)$", "$B_{\phi}\ (T)$", "$Flux\ (Wb)$"]
	for i in range(len(cb)): cb[i].set_label(labels[i])

	pyp.show(block=False)


