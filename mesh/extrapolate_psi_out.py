import numpy 						as np
from routines.RectBivariateSpline	import RectBivariateSpline
from routines.globals				import DEBUG

def extrapolate_psi_out(Config, nExt=100):

	if(DEBUG > 0):print("extrapolate_psi_out")

#	nExt		num. of points to extrapolate on all sides.

	extrapol_val = Config.extrapol_val
	raise_value	 = Config.raise_value
	raise_power	 = Config.raise_power
	smooth		 = Config.Smooth2DExt
	fluxmin		 = Config.flux2D.min()
	fluxmax		 = Config.flux2D.max()

	Rmin = Config.r2D.min()
	Zmin = Config.z2D.min()
	Rmax = Config.r2D.max()
	Zmax = Config.z2D.max()

	dpsiN = (fluxmax-fluxmin)*raise_value[0]
	dpsiS = (fluxmax-fluxmin)*raise_value[1]
	dpsiE = (fluxmax-fluxmin)*raise_value[2]
	dpsiW = (fluxmax-fluxmin)*raise_value[3]

	nExtZ = max(nExt, int(extrapol_val*Config.r2D.shape[0]+0.5))
	nExtR = max(nExt, int(extrapol_val*Config.r2D.shape[1]+0.5))
	shape	= [Config.r2D.shape[0]+2*nExtZ, Config.r2D.shape[1]+2*nExtR]
	r2D		= np.zeros(shape,	dtype=Config.r2D.dtype)
	z2D		= np.zeros(shape,	dtype=Config.z2D.dtype)
	flux2D	= np.zeros(shape,	dtype=Config.flux2D.dtype)
	Br2D	= np.zeros(shape,	dtype=Config.Br2D.dtype)
	Bz2D	= np.zeros(shape,	dtype=Config.Bz2D.dtype)
	Bphi2D	= np.zeros(shape,	dtype=Config.Bphi2D.dtype)

	r2D[nExtZ:-nExtZ,	nExtR:-nExtR]	= Config.r2D
	z2D[nExtZ:-nExtZ,	nExtR:-nExtR]	= Config.z2D
	flux2D[nExtZ:-nExtZ,nExtR:-nExtR]	= Config.flux2D
	Br2D[nExtZ:-nExtZ,	nExtR:-nExtR]	= Config.Br2D
	Bz2D[nExtZ:-nExtZ,	nExtR:-nExtR]	= Config.Bz2D
	Bphi2D[nExtZ:-nExtZ,nExtR:-nExtR]	= Config.Bphi2D

#	50/200 ==> +25% de flux

#	extrapolate radial

	k = Config.r2D.shape[1]+nExtR
	for i in range(1,nExtR+1):
		r2D[nExtZ:-nExtZ,k]	 	=	Rmax+i/nExtR*((Rmax-Rmin)*extrapol_val)
		z2D[nExtZ:-nExtZ,k]	 	=	z2D[nExtZ:-nExtZ,k-1]
		flux2D[nExtZ:-nExtZ,k]	=	2*flux2D[nExtZ:-nExtZ,k-1] - flux2D[nExtZ:-nExtZ,k-2]+dpsiW*(i/nExtR)**raise_power
		Br2D[nExtZ:-nExtZ,k]	=	2*Br2D[nExtZ:-nExtZ,k-1]	 - Br2D[nExtZ:-nExtZ,k-2]
		Bz2D[nExtZ:-nExtZ,k]	=	2*Bz2D[nExtZ:-nExtZ,k-1]	 - Bz2D[nExtZ:-nExtZ,k-2]
		Bphi2D[nExtZ:-nExtZ,k]	=	2*Bphi2D[nExtZ:-nExtZ,k-1] - Bphi2D[nExtZ:-nExtZ,k-2]
		k +=1

	k = nExtR - 1
	for i in range(1,nExtR+1):
		r2D[nExtZ:-nExtZ,k]		=	Rmin-i/nExtR*((Rmax-Rmin)*extrapol_val)
		z2D[nExtZ:-nExtZ,k]	 	=	z2D[nExtZ:-nExtZ,k+1]
		flux2D[nExtZ:-nExtZ,k]	=	2*flux2D[nExtZ:-nExtZ,k+1] - flux2D[nExtZ:-nExtZ,k+2]+dpsiE*(i/nExtR)**raise_power
		Br2D[nExtZ:-nExtZ,k]	=	2*Br2D[nExtZ:-nExtZ,k+1]	- Br2D[nExtZ:-nExtZ,k+2]
		Bz2D[nExtZ:-nExtZ,k]	=	2*Bz2D[nExtZ:-nExtZ,k+1]	- Bz2D[nExtZ:-nExtZ,k+2]
		Bphi2D[nExtZ:-nExtZ,k]	=	2*Bphi2D[nExtZ:-nExtZ,k+1] - Bphi2D[nExtZ:-nExtZ,k+2]
		k -= 1

#	extrapolate vertical

	k = Config.r2D.shape[0]+nExtZ
	for i in range(1,nExtZ+1):
		r2D[k,:]	 = r2D[k-1,:]
		z2D[k,:]	 = Zmax			   + i/nExtZ*((Zmax-Zmin)*extrapol_val)
		flux2D[k,:]  = 2*flux2D[k-1,:] - flux2D[k-2,:]+dpsiN*(i/nExtZ)**raise_power
		Br2D[k,:]	 = 2*Br2D[k-1,:]   - Br2D[k-2,:]
		Bz2D[k,:]	 = 2*Bz2D[k-1,:]   - Bz2D[k-2,:]
		Bphi2D[k,:]  = 2*Bphi2D[k-1,:] - Bphi2D[k-2,:]
		k +=1

	k = nExtZ - 1
	for i in range(1,nExtZ+1):
		r2D[k,:]	 = r2D[k+1,:]
		z2D[k,:]	 = Zmin				- i/nExtZ*((Zmax-Zmin)*extrapol_val)
		flux2D[k,:]  = 2*flux2D[k+1,:]	- flux2D[k+2,:]+dpsiS*(i/nExtZ)**raise_power
		Br2D[k,:]	 = 2*Br2D[k+1,:]    - Br2D[k+2,:]
		Bz2D[k,:]	 = 2*Bz2D[k+1,:]    - Bz2D[k+2,:]
		Bphi2D[k,:]  = 2*Bphi2D[k+1,:]	- Bphi2D[k+2,:]
		k -= 1

	if(smooth != 0.):
		r2D_mat1		= r2D.reshape(r2D.size)
		z2D_mat1		= z2D.reshape(z2D.size)
		flux2D	= np.reshape(RectBivariateSpline(r2D[0,:], z2D[:,0], flux2D.T, r2D_mat1, z2D_mat1, s=smooth), r2D.shape)

	Config.r2D		= r2D
	Config.z2D		= z2D
	Config.flux2D	= flux2D
	Config.Br2D		= Br2D
	Config.Bz2D		= Bz2D
	Config.Bphi2D	= Bphi2D


	if(DEBUG > 1):
		print("Min(Br2)   (T)  ={:8.3f}, Max(Br2)    (T)  ={:8.3f}".format(np.min(Config.Br2),np.max(Config.Br2)))
		print("Min(Bz2)   (T)  ={:8.3f}, Max(Bz2)    (T)  ={:8.3f}".format(np.min(Config.Bz2),np.max(Config.Bz2)))
		print("Min(Bphi2) (T)  ={:8.3f}, Max(Bphi2)  (T)  ={:8.3f}".format(np.min(Config.Bphi2),np.max(Config.Bphi2)))
		print("Min(flux2) (Wb) ={:8.3f}, Max(flux2)  (Wb) ={:8.3f}".format(np.min(Config.flux2),np.max(Config.flux2)))

		print("Min(Br2D)  (T)  ={:8.3f}, Max(Br2D)  (T)  ={:8.3f}".format(np.min(Config.Br2D),np.max(Config.Br2D)))
		print("Min(Bz2D)  (T)  ={:8.3f}, Max(Bz2D)  (T)  ={:8.3f}".format(np.min(Config.Bz2D),np.max(Config.Bz2D)))
		print("Min(Bphi2D)(T)  ={:8.3f}, Max(Bphi2D)(T)  ={:8.3f}".format(np.min(Config.Bphi2D),np.max(Config.Bphi2D)))
		print("Min(flux2D)(Wb) ={:8.3f}, Max(flux2D)(Wb) ={:8.3f}".format(np.min(Config.flux2D),np.max(Config.flux2D)))

	if(DEBUG > 0):print("extrapolate_psi_out: Completed")
