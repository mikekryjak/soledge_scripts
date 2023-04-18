
import numpy as np

def broadcast_mesh(Zones):

	for k in range(len(Zones)):
		Zones[k].xmin	= Zones[k].xb[0, 0]
		Zones[k].xmax	= Zones[k].xb[-1,0]
		Zones[k].zmin	= Zones[k].zb[0, 0]
		Zones[k].zmax	= Zones[k].zb[0,-1]
		Zones[k].xm		= Zones[k].xb[:-1, :-1]
		Zones[k].xp		= Zones[k].xb[1:,  :-1]
		Zones[k].zm		= Zones[k].zb[:-1, :-1]
		Zones[k].zp		= Zones[k].zb[:-1, 1:]

	for k in range(len(Zones)):
		
		Nx = Zones[k].x.shape[0]
		Nz = Zones[k].x.shape[1]
		x  = np.empty((Nx+2, Nz+2), dtype='f8')
		z  = np.empty((Nx+2, Nz+2), dtype='f8')
		xm = np.empty((Nx+2, Nz+2), dtype='f8')
		xp = np.empty((Nx+2, Nz+2), dtype='f8')
		zm = np.empty((Nx+2, Nz+2), dtype='f8')
		zp = np.empty((Nx+2, Nz+2), dtype='f8')

		x[1:-1,1:-1]	= Zones[k].x
		z[1:-1,1:-1]	= Zones[k].z
		xm[1:-1,1:-1]	= Zones[k].xm
		xp[1:-1,1:-1]	= Zones[k].xp
		zm[1:-1,1:-1]	= Zones[k].zm
		zp[1:-1,1:-1]	= Zones[k].zp


#		North boundary
		North = Zones[k].Neighbour.north
		if(North < -1):
			x[-1,:]		= Zones[k].xmax
			z[-1,:1:-1] = Zones[k].z[-1,:]
			zm[-1,1:-1] = Zones[k].zm[-1,:]
			zp[-1,1:-1] = Zones[k].zp[-1,:]
		else:
			x[-1,1:-1]	= Zones[k].xmax + Zones[North].x[0,:] - Zones[North].xmax
			z[-1,:1:-1] = Zones[North].z[0,:]
			zm[-1,1:-1] = Zones[North].zm[-1,:]
			zp[-1,1:-1] = Zones[North].zp[-1,:]		
#			corners
			x[-1,0]		=  Zones[k].xmax + Zones[North].x[0,0] - Zones[North].xmin
			x[-1,-1]	=  Zones[k].xmax + Zones[North].x[0,-1] - Zones[North].xmin

#		South boundary
		South = Zones[k].Neighbour.south
		if(South < -1):
			x[0,:]		= Zones[k].xmin
			z[0,1:-1]	= Zones[k].z[0,:]
			zm[0,1:-1]	= Zones[k].zm[0,:]
			zp[0,1:-1]	= Zones[k].zp[0,:]
		else:
			x[0,1:-1]	= Zones[k].xmin + Zones[South].x[-1,:] -  Zones[South].xmax
			z[0,1:-1]	= Zones[South].z[-1,:]
			zm[0,1:-1]	= Zones[South].zm[-1,:]
			zp[0,1:-1]	= Zones[South].zp[-1,:]
#			corners
			x[0,0]		= Zones[k].xmin + Zones[South].x[-1,0] - Zones[South].xmax
			x[0,-1]		= Zones[k].xmin + Zones[South].x[-1,0] - Zones[South].xmax

#		East boundary
		East = Zones[k].Neighbour.east
		if(East < 1):
			zp[1:-1,-1]	= 2.*Zones[k].zp[:,-1]-Zones[k].zp[:,-2]
			z[1:-1,-1]	= (3.*Zones[k].zp[:,-1]-Zones[k].zp[:,-2])/2.
			zm[1:-1,-1] = Zones[k].zp[:,-1]
			x[1:-1,-1]  = Zones[k].x[:,-1]
			xm[1:-1,-1] = Zones[k].xm[:,-1]
			xp[1:-1,-1]	= Zones[k].xp[:,-1]
		else:
			zp[1:-1,-1]	= Zones[k].zmax + Zones[East].zp[:,0] - Zones[East].zmin
			z[1:-1,-1]	= Zones[k].zmax + Zones[East].z[:,0]  - Zones[East].zmin
			zm[1:-1,-1]	= Zones[k].zmax + Zones[East].zm[:,0] - Zones[East].zmin
			x[1:-1,-1]	= Zones[East].x[:,0]
			xm[1:-1,-1] = Zones[East].xm[:,0]
			xp[1:-1,-1] = Zones[East].xp[:,0]

#		West boundary
		West = Zones[k].Neighbour.west
		if(West < -1):
			zp[1:-1,0]	= Zones[k].zm[:,0]
			z[1:-1,0]	= (3.*Zones[k].zm[:,0]-Zones[k].zm[:,1])/2.
			zm[1:-1,0]	= 2.*Zones[k].zm[:,0]-Zones[k].zm[:,1]
			x[1:-1,0]	= Zones[k].x[:,0]
			xm[1:-1,0]	= Zones[k].xm[:,0]
			xp[1:-1,0]	= Zones[k].xp[:,0]
		else:
			zp[1:-1,0]	= Zones[k].zmin + Zones[West].zp[:,-1] - Zones[West].zmax
			z[1:-1,0]	= Zones[k].zmin + Zones[West].z[:,-1]  - Zones[West].zmax
			zm[1:-1,0]	= Zones[k].zmin + Zones[West].zm[:,-1] - Zones[West].zmax
			x[1:-1,0]	= Zones[West].x[:,-1]
			xm[1:-1,0]	= Zones[West].xm[:,-1]
			xp[1:-1,0]	= Zones[West].xp[:,-1]


		Zones[k].bc_x  = x
		Zones[k].bc_z  = z
		Zones[k].bc_xm = xm
		Zones[k].bc_xm = xm
		Zones[k].bc_zp = zp		
		Zones[k].bc_zm = zm		

	
		Rcenter				= np.zeros((Nx+2,Nz+2), dtype='f8')
		Rcenter[1:-1,1:-1]	= Zones[k].gridRc

		Rcenter[0,1:-1]		= Zones[k].SouthP.R
		Rcenter[-1,1:-1]	= Zones[k].NorthP.R
		Rcenter[1:-1,0]		= Zones[k].WestP.R
		Rcenter[1:-1,-1]	= Zones[k].EastP.R
	
		Zcenter				= np.zeros((Nx+2,Nz+2), dtype='f8')
		Zcenter[1:-1,1:-1]	= Zones[k].gridZc
		Zcenter[0,1:-1]		= Zones[k].SouthP.Z
		Zcenter[-1,1:-1]	= Zones[k].NorthP.Z
		Zcenter[1:-1,0]		= Zones[k].WestP.Z
		Zcenter[1:-1,-1]	= Zones[k].EastP.Z

		Zones[k].bc_gridRc	= Rcenter
		Zones[k].bc_gridZc	= Zcenter

		Br					= np.zeros((Nx+2,Nz+2), dtype='f8')
		Br[1:-1,1:-1]		= Zones[k].Br
		Br[0,1:-1]			= Zones[k].SouthP.Br
		Br[-1,1:-1]			= Zones[k].NorthP.Br
		Br[1:-1,0]			= Zones[k].WestP.Br
		Br[1:-1,-1]			= Zones[k].EastP.Br
	
		Bz					= np.zeros((Nx+2,Nz+2), dtype='f8')
		Bz[1:-1,1:-1]		= Zones[k].Bz
		Bz[0,1:-1]			= Zones[k].SouthP.Bz
		Bz[-1,1:-1]			= Zones[k].NorthP.Bz
		Bz[1:-1,0]			= Zones[k].WestP.Bz
		Bz[1:-1,-1]			= Zones[k].EastP.Bz
	
		Bphi				= np.zeros((Nx+2,Nz+2), dtype='f8')
		Bphi[1:-1,1:-1]		= Zones[k].Bphi
		Bphi[0,1:-1]		= Zones[k].SouthP.Bphi
		Bphi[-1,1:-1]		= Zones[k].NorthP.Bphi
		Bphi[1:-1,0]		= Zones[k].WestP.Bphi
		Bphi[1:-1,-1]		= Zones[k].EastP.Bphi

		Zones[k].bc_Br		= Br
		Zones[k].bc_Bz		= Bz
		Zones[k].bc_Bphi	= Bphi



