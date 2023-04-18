from tkinter 						import messagebox
from tkinter.messagebox 			import Message
import numpy						as np
from routines.RectBivariateSpline	import RectBivariateSpline
from routines.utils_walls			import get_min_max_walls
from routines.globals				import DEBUG

def smooth_fields(Config):

	if(DEBUG > 0): print("smooth_fields")

	nr							= Config.nr
	nz							= Config.nz
	smooth						= Config.Smooth2D
	OutLength2					= Config.OutLength2

	if(nr == 0): nr = Config.r2.shape[1]
	if(nz == 0): nz = Config.r2.shape[0]

	if((nr == Config.r2.shape[1]) and (nz == Config.r2.shape[0]) and 
		(smooth == Config.Smooth2) and (OutLength2 == Config.OutLength2)):				# No modification 
		Config.r2D   	= np.copy(Config.r2)
		Config.z2D		= np.copy(Config.z2)
		Config.Br2D		= np.copy(Config.Br2)
		Config.Bphi2D	= np.copy(Config.Bphi2)
		Config.Bz2D		= np.copy(Config.Bz2)
		Config.flux2D	= np.copy(Config.flux2)
		return

	
	if((smooth == Config.Smooth2) and (OutLength2 != Config.OutLength2)):				# check if only mesh size reduction
		r_min, r_max, z_min, z_max = get_min_max_walls(Config,plasma_wall=True, eirene_wall=True, extra_wall=False)
		r_min 	   -= OutLength2
		r_max 	   += OutLength2
		z_min 	   -= OutLength2
		z_max 	   += OutLength2
		i_min 		= np.where(Config.r2[0,:] >= r_min)[0][0]
		i_max 		= np.where(Config.r2[0,:] <= r_max)[0][-1]
		j_min 		= np.where(Config.z2[:,0] >= z_min)[0][0]
		j_max 		= np.where(Config.z2[:,0] <= z_max)[0][-1]
		if((nr == i_max - i_min +1) and (nz == j_max - j_min +1)):
			Config.r2D   	= Config.r2[j_min:j_max+1, i_min:i_max+1]
			Config.z2D		= Config.z2[j_min:j_max+1, i_min:i_max+1]
			Config.Br2D		= Config.Br2[j_min:j_max+1, i_min:i_max+1]
			Config.Bphi2D	= Config.Bphi2[j_min:j_max+1, i_min:i_max+1]
			Config.Bz2D		= Config.Bz2[j_min:j_max+1, i_min:i_max+1]
			Config.flux2D	= Config.flux2[j_min:j_max+1, i_min:i_max+1]
			return


	r_min 	= Config.r2.min()
	r_max 	= Config.r2.max()
	z_min 	= Config.z2.min()
	z_max 	= Config.z2.max()
	if(OutLength2 != Config.OutLength2):												# Mesh Size reduction
		rw_min, rw_max, zw_min, zw_max = get_min_max_walls(Config,plasma_wall=True, eirene_wall=True, extra_wall=False)
		r_min 		= max(rw_min - OutLength2, r_min)
		r_max 		= min(rw_max + OutLength2, r_max)
		z_min 		= max(zw_min - OutLength2, z_min)
		z_max 		= min(zw_max + OutLength2, z_max)

	r1D	= np.linspace(r_min, r_max, nr)
	z1D	= np.linspace(z_min, z_max, nz)
	
	Config.r2D, Config.z2D		= np.meshgrid(r1D, z1D)
	
	r2D_mat1		= Config.r2D.reshape(Config.r2D.size)
	z2D_mat1		= Config.z2D.reshape(Config.z2D.size)

	Config.flux2D	= np.reshape(RectBivariateSpline(Config.r2[0,:], Config.z2[:,0], Config.flux2.T, r2D_mat1, z2D_mat1, s=smooth), Config.r2D.shape)
	Config.Bphi2D	= np.reshape(RectBivariateSpline(Config.r2[0,:], Config.z2[:,0], Config.Bphi2.T, r2D_mat1, z2D_mat1, s=smooth), Config.r2D.shape)

#	Check for normalize flux

	NormFlux = False
	for k in range(len(Config.flux2_x_psi)):
		if(abs(Config.flux2_x_psi[k] - 1.0) < 0.1):
			NormFlux = True
			break

	Br2D = np.reshape(RectBivariateSpline(Config.r2[0,:], Config.z2[:,0], Config.Br2.T,   r2D_mat1, z2D_mat1, s=smooth), Config.r2D.shape)		#	Python indexes matrix Bz[z,r]
	Bz2D = np.reshape(RectBivariateSpline(Config.r2[0,:], Config.z2[:,0], Config.Bz2.T,   r2D_mat1, z2D_mat1, s=smooth), Config.r2D.shape)

	if(NormFlux):								# for normalized flux use smoothed Br and Bz
		Config.Br2D = Br2D
		Config.Bz2D = Bz2D

	else:											# for not normalized flux use flux to compute Br and Bz 


		Config.Br2D			 = np.empty_like(Config.r2D)
		Config.Bz2D			 = np.empty_like(Config.r2D)

		Dr	 = r1D[2] - r1D[0]
		Dz	 = z1D[2] - z1D[0]
		Config.Bz2D[:,1:-1]	 = -(Config.flux2D[:,2:] - Config.flux2D[:,:-2])/Dr/Config.r2D[:,1:-1]
		Config.Br2D[1:-1,:]	 =  (Config.flux2D[2:,:] - Config.flux2D[:-2,:])/Dz/Config.r2D[1:-1,:]

		Config.Bz2D[:,0]	 = -(Config.flux2D[:,1]  - Config.flux2D[:,0])/(0.5*Dr)/r1D[0]
		Config.Bz2D[:,-1]	 = -(Config.flux2D[:,-1] - Config.flux2D[:,-2])/(0.5*Dr)/r1D[-1]
		Config.Br2D[0,:]	 =  (Config.flux2D[1,:]  - Config.flux2D[0,:])/(0.5*Dz)/r1D
		Config.Br2D[-1,:]	 =  (Config.flux2D[-1,:] - Config.flux2D[-2,:])/(0.5*Dz)/r1D

		r1D = 0; z1D = 0

#		Correct sign if it was reversed relatively to flux evaluation

		if(np.sum(np.abs(Config.Br2D-Br2D)) > np.sum(np.abs(Config.Br2D+Br2D))):
			print("\tReverse sign to Br:")
			print("\t\tnp.abs(Config.Br2D-Br2D) =",np.sum(np.abs(Config.Br2D-Br2D)))
			print("\t\tnp.abs(Config.Br2D+Br2D) =",np.sum(np.abs(Config.Br2D+Br2D)))
			Config.Br2D = -Config.Br2D
		if(np.sum(np.abs(Config.Bz2D-Bz2D)) > np.sum(np.abs(Config.Bz2D+Bz2D))):
			print("\tReverse sign to Bz:")
			print("\t\tnp.abs(Config.Bz2D-Bz2D) =",np.sum(np.abs(Config.Bz2D-Bz2D)))
			print("\t\tnp.abs(Config.Bz2D+Bz2D) =",np.sum(np.abs(Config.Bz2D+Bz2D)))
			Config.Bz2D = -Config.Bz2D

	DSmooth = np.array([[(Config.Br2D.min()-Config.Br2.min())/(Config.Br2.max()-Config.Br2.min()),			(Config.Br2D.max()-Config.Br2.max())/(Config.Br2.max()-Config.Br2.min())],
						[(Config.Bz2D.min()-Config.Bz2.min())/(Config.Bz2.max()-Config.Bz2.min()),			(Config.Bz2D.max()-Config.Bz2.max())/(Config.Bz2.max()-Config.Bz2.min())],
						[(Config.Bphi2D.min()-Config.Bphi2.min())/(Config.Bphi2.max()-Config.Bphi2.min()),	(Config.Bphi2D.max()-Config.Bphi2.max())/(Config.Bphi2.max()-Config.Bphi2.min())],
						[(Config.flux2D.min()-Config.flux2.min())/(Config.flux2.max()-Config.flux2.min()),	(Config.flux2D.max()-Config.flux2.max())/(Config.flux2.max()-Config.flux2.min())]])*100

		
	if((DEBUG > 1) or ((smooth > 0.) and (np.max(np.abs(DSmooth)) > 2*smooth))):
		print("\tMin(Br2D)   (T)  ={:8.3f}, Max(Br2D)   (T)  ={:8.3f}".format(Config.Br2D.min(),  Config.Br2D.max()))
		print("\tMin(Bz2D)   (T)  ={:8.3f}, Max(Bz2d)   (T)  ={:8.3f}".format(Config.Bz2D.min(),  Config.Bz2D.max()))
		print("\tMin(Bphi2D) (T)  ={:8.3f}, Max(Bphi2D) (T)  ={:8.3f}".format(Config.Bphi2D.min(),Config.Bphi2D.max()))
		print("\tMin(flux2D) (Wb) ={:8.3f}, Max(flux2D) (Wb) ={:8.3f}".format(Config.flux2D.min(),     Config.flux2D.max()))
		print("\tDMin(Br2D)  (%)  ={:8.1f}, DMax(Br2D)  (%)  ={:8.1f}".format(DSmooth[0,0], DSmooth[0,1]))
		print("\tDMin(Bz2D)  (%)  ={:8.1f}, DMax(Bz2d)  (%)  ={:8.1f}".format(DSmooth[1,0], DSmooth[1,1]))
		print("\tDMin(Bphi2D)(%)  ={:8.1f}, DMax(Bphi2D)(%)  ={:8.1f}".format(DSmooth[2,0], DSmooth[2,1]))
		print("\tDMin(flux2D)(%)  ={:8.1f}, DMax(flux2D)(%)  ={:8.1f}".format(DSmooth[3,0], DSmooth[3,1]))

	if(smooth > 0.):
		if((np.abs(np.sum(np.abs(Br2D))/Br2D.size - np.sum(np.abs(Config.Br2))/Config.Br2.size) > np.sum(np.abs(Config.Br2))/Config.Br2.size) or
		   (np.abs(np.sum(np.abs(Bz2D))/Bz2D.size - np.sum(np.abs(Config.Bz2))/Config.Bz2.size) > np.sum(np.abs(Config.Br2))/Config.Bz2.size)):
			messagebox.showerror("Smooth fields", "POSSIBLE ERROR IN SMOTTHING\nCHECK ORIGINAL AND SMOOTHED FIELDS")

		elif(np.max(np.abs(DSmooth)) > 2*smooth):
			if(np.max(np.abs(DSmooth)) < 999.):
				messagebox.showwarning("Smooth fields", "Max diff. in field/flux smoothing = {:3.1f} (%)".format(np.max(np.abs(DSmooth))))
			else:
				messagebox.showwarning("Smooth fields", "Max diff. in field/flux smoothing = {:.3e} (%)".format(np.max(np.abs(DSmooth))))

	if(DEBUG > 0): print("smooth_fields: Completed")
