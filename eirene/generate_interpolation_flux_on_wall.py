import os
import numpy				as np
import numpy.matlib			as mat
from routines.globals		import *

# Function definition is here
#=======================================================================================================================================

def generate_interpolation_flux_on_wall(Zones):

	if(DEBUG > 0):	print("generate_interpolation_flux_on_wall")

	nZones = len(Zones)
	for k in range(nZones):
		Nx = Zones[k].Nx
		Nz = Zones[k].Nz
		Zones[k].Chi2 = np.zeros((Nx+2,Nz+2), dtype='f8')
		Zones[k].R2   = np.zeros((Nx+2,Nz+2), dtype='f8')
		Zones[k].Z2   = np.zeros((Nx+2,Nz+2), dtype='f8')

		Zones[k].Chi2	[1:-1, 1:-1]	= np.copy(Zones[k].Chi)
		Zones[k].R2		[1:-1, 1:-1] 	= np.copy(Zones[k].gridRc)
		Zones[k].Z2		[1:-1, 1:-1]	= np.copy(Zones[k].gridZc)

		if Zones[k].Neighbour.north >= 0:														#North
			kNe  = Zones[k].Neighbour.north
			Zones[k].Chi2	[-1, 1:-1] = Zones[kNe].Chi		[0, :]
			Zones[k].R2		[-1, 1:-1] = Zones[kNe].gridRc	[0, :]
			Zones[k].Z2		[-1, 1:-1] = Zones[kNe].gridZc	[0, :]
		else:
			Zones[k].Chi2	[-1, 1:-1] = Zones[k].Chi2	[-2, 1:-1]
			Zones[k].R2		[-1, 1:-1] = Zones[k].R2	[-2, 1:-1]
			Zones[k].Z2		[-1, 1:-1] = Zones[k].Z2	[-2, 1:-1]

		if Zones[k].Neighbour.south >= 0:														#South
			kNe  = Zones[k].Neighbour.south
			NxNe = Zones[kNe].Nx - 1
			Zones[k].Chi2	[0, 1:-1] = Zones[kNe].Chi		[NxNe, :]
			Zones[k].R2		[0, 1:-1] = Zones[kNe].gridRc	[NxNe, :]
			Zones[k].Z2		[0, 1:-1] = Zones[kNe].gridZc	[NxNe, :]
		else:
			Zones[k].Chi2	[0, 1:-1] = Zones[k].Chi2	[1, 1:-1]
			Zones[k].R2		[0, 1:-1] = Zones[k].R2		[1, 1:-1]
			Zones[k].Z2		[0, 1:-1] = Zones[k].Z2		[1, 1:-1]

		if Zones[k].Neighbour.east >= 0:														#East
			kNe  = Zones[k].Neighbour.east
			Zones[k].Chi2	[1:-1, -1] = Zones[kNe].Chi		[:, 0]
			Zones[k].R2		[1:-1, -1] = Zones[kNe].gridRc	[:, 0]
			Zones[k].Z2		[1:-1, -1] = Zones[kNe].gridZc	[:, 0]
		else:
			Zones[k].Chi2	[1:-1, -1] = Zones[k].Chi2	[1:-1, -2]
			Zones[k].R2		[1:-1, -1] = Zones[k].R2	[1:-1, -2]
			Zones[k].Z2		[1:-1, -1] = Zones[k].Z2	[1:-1, -2]

		if Zones[k].Neighbour.west >= 0:														#West
			kNe  = Zones[k].Neighbour.west
			NzNe = Zones[kNe].Nz - 1
			Zones[k].Chi2 	[1:-1, 0] = Zones[kNe].Chi		[:, NzNe]
			Zones[k].R2		[1:-1, 0] = Zones[kNe].gridRc	[:, NzNe]
			Zones[k].Z2		[1:-1, 0] = Zones[kNe].gridZc	[:, NzNe]
		else:
			Zones[k].Chi2	[1:-1, 0] = Zones[k].Chi2	[1:-1, 1]
			Zones[k].R2		[1:-1, 0] = Zones[k].R2		[1:-1, 1]
			Zones[k].Z2		[1:-1, 0] = Zones[k].Z2		[1:-1, 1]

#	remove point isolated (only one point along theta) masked points

	for k in range(nZones):
		ii,jj = np.where((Zones[k].Chi2[1:-1,0:-2] == 0) & (Zones[k].Chi2[1:-1,1:-1] == 1) & (Zones[k].Chi2[1:-1,2:] == 0))		#isolated masked inside zone 
		if(len(ii) > 0): 
			if(DEBUG > 1):	
				print("\tRemoved isolated chi=1 at Zone=",k+1," i=",ii+1," j=",jj+1)
				print("\t\tR=",Zones[k].gridRc[ii,jj]," Z=",Zones[k].gridZc[ii,jj])		
				print("\t\tfrom=",Zones[k].Chi2[ii+1,jj+1])
			Zones[k].Chi2[ii+1,jj+1] = 0
			if(DEBUG > 1):	
				print("\t\tto  =",Zones[k].Chi2[ii+1,jj+1])

		ii = np.where((Zones[k].Chi2[1:-1,0] == 1) & (Zones[k].Chi2[1:-1,1] == 0))[0]			#isolated masked point point at zone gegin
		if(len(ii) > 0):
			kWest	= Zones[k].Neighbour.west
			ii_west = np.where(Zones[kWest].Chi2[ii+1,-3] == 0)[0]
			if(len(ii_west) > 0):
				if(DEBUG > 1):	
					print("\tRemoved isolated chi=1 at Zone=",kWest+1," i=",ii[ii_west]+1," j=end")
					print("\t\tR=",Zones[kWest].gridRc[ii[ii_west],-1]," Z=",Zones[kWest].gridZc[ii[ii_west],-1])
					print("\t\tfrom=",Zones[k].Chi2[ii[ii_west]+1,0])
				Zones[k].Chi2[ii[ii_west]+1,0] = 0
				Zones[kWest].Chi2[ii[ii_west]+1,-2] = 0
				if(DEBUG > 1):	
					print("\t\tto  =",Zones[k].Chi2[ii[ii_west]+1,0])


		ii = np.where((Zones[k].Chi2[1:-1,-1] == 1) & (Zones[k].Chi2[1:-1,-2] == 0))[0]		#isolated masked point point at zone end
		if(len(ii) > 0):
			kEast  = Zones[k].Neighbour.east
			ii_est = np.where(Zones[kEast].Chi2[ii+1,2] == 0)[0]
			if(len(ii_est) > 0):
				if(DEBUG > 1):	
					print("\tRemoved isolated chi=1 at Zone=",kEast+1," i=",ii[ii_est]+1," j=1]")
					print("\t\tR=",Zones[kEast].gridRc[ii[ii_est],0]," Z=",Zones[kEast].gridZc[ii[ii_est],0])
					print("\t\tfrom=",Zones[k].Chi2[ii[ii_est]+1,-1])
				Zones[k].Chi2[ii[ii_est]+1,-1] = 0
				Zones[kEast].Chi2[ii[ii_est]+1,1] = 0
				if(DEBUG > 1):	
					print("\t\tto  =",Zones[k].Chi2[ii[ii_est]+1,-1])


	for k in range(nZones):
		Nx = Zones[k].Nx
		Nz = Zones[k].Nz

		iVec = np.arange(Nx)
		jVec = np.arange(Nz-1)
		kMat = np.zeros((len(iVec), len(jVec)), dtype='i4') + k
		iMat = np.transpose(mat.repmat(iVec, len(jVec), 1))
		jMat = mat.repmat(jVec, len(iVec), 1)

		Zones[k].point_east_k = np.empty((Nx, Nz), dtype = 'i4')
		Zones[k].point_east_i = np.empty((Nx, Nz), dtype = 'i4')
		Zones[k].point_east_j = np.empty((Nx, Nz), dtype = 'i4')
										
		Zones[k].point_east_k[:, :-1] = kMat
		Zones[k].point_east_i[:, :-1] = iMat
		Zones[k].point_east_j[:, :-1] = jMat+1

		if Zones[k].Neighbour.east >= 0:														#East
			Zones[k].point_east_k[:, -1] = Zones[k].Neighbour.east
			Zones[k].point_east_i[:, -1] = iVec
			Zones[k].point_east_j[:, -1] = np.zeros(Nx)
		else:
			Zones[k].point_east_k[:, -1] = -np.ones((Nx), dtype='i4')
			Zones[k].point_east_i[:, -1] = -np.ones((Nx), dtype='i4')
			Zones[k].point_east_j[:, -1] = -np.ones((Nx), dtype='i4')

		iVec = np.arange(Nx)																#West
		jVec = np.arange(1,Nz)
		kMat = np.zeros((len(iVec), len(jVec)), dtype='i4') + k
		iMat = np.transpose(mat.repmat(iVec, len(jVec), 1))
		jMat = mat.repmat(jVec, len(iVec), 1)

		Zones[k].point_west_k = np.empty((Nx, Nz), dtype = 'i4')
		Zones[k].point_west_i = np.empty((Nx, Nz), dtype = 'i4')
		Zones[k].point_west_j = np.empty((Nx, Nz), dtype = 'i4')
										
		Zones[k].point_west_k[:, 1:] = kMat
		Zones[k].point_west_i[:, 1:] = iMat
		Zones[k].point_west_j[:, 1:] = jMat-1

		if Zones[k].Neighbour.west >= 0:
			Zones[k].point_west_k[:, 0] = Zones[k].Neighbour.west
			Zones[k].point_west_i[:, 0] = iVec
			Zones[k].point_west_j[:, 0] = Zones[Zones[k].Neighbour.west].Nz-1
		else:
			Zones[k].point_west_k[:, 0] = -np.ones((Nx), dtype='i4')
			Zones[k].point_west_i[:, 0] = -np.ones((Nx), dtype='i4')
			Zones[k].point_west_j[:, 0] = -np.ones((Nx), dtype='i4')


		iVec = np.arange(Nx-1)																#North
		jVec = np.arange(Nz)
		kMat = np.zeros((len(iVec), len(jVec)), dtype='i4') + k
		iMat = np.transpose(mat.repmat(iVec, len(jVec), 1))
		jMat = mat.repmat(jVec, len(iVec), 1)

		Zones[k].point_north_k = np.empty((Nx, Nz), dtype = 'i4')
		Zones[k].point_north_i = np.empty((Nx, Nz), dtype = 'i4')
		Zones[k].point_north_j = np.empty((Nx, Nz), dtype = 'i4')

		Zones[k].point_north_k[:-1, :] = kMat
		Zones[k].point_north_i[:-1, :] = iMat+1
		Zones[k].point_north_j[:-1, :] = jMat

		if Zones[k].Neighbour.north >= 0:
			Zones[k].point_north_k[-1, :] = Zones[k].Neighbour.north
			Zones[k].point_north_i[-1, :] = np.zeros(Nz)
			Zones[k].point_north_j[-1, :] = jVec
		else:
			Zones[k].point_north_k[-1, :] = -np.ones((Nz), dtype='i4')
			Zones[k].point_north_i[-1, :] = -np.ones((Nz), dtype='i4')
			Zones[k].point_north_j[-1, :] = -np.ones((Nz), dtype='i4')


		iVec = np.arange(1,Nx)																#South
		jVec = np.arange(Nz)
		kMat = np.zeros((len(iVec), len(jVec)), dtype='i4') + k
		iMat = np.transpose(mat.repmat(iVec, len(jVec), 1))
		jMat = mat.repmat(jVec, len(iVec), 1)

		Zones[k].point_south_k = np.empty((Nx, Nz), dtype = 'i4')
		Zones[k].point_south_i = np.empty((Nx, Nz), dtype = 'i4')
		Zones[k].point_south_j = np.empty((Nx, Nz), dtype = 'i4')

		Zones[k].point_south_k[1:, :] = kMat
		Zones[k].point_south_i[1:, :] = iMat-1
		Zones[k].point_south_j[1:, :] = jMat

		if Zones[k].Neighbour.south >= 0:
			Zones[k].point_south_k[0, :] = Zones[k].Neighbour.south
			Zones[k].point_south_i[0, :] = Zones[Zones[k].Neighbour.south].Nx-1
			Zones[k].point_south_j[0, :] = jVec
		else:
			Zones[k].point_south_k[0, :] = -np.ones((Nz), dtype='i4')
			Zones[k].point_south_i[0, :] = -np.ones((Nz), dtype='i4')
			Zones[k].point_south_j[0, :] = -np.ones((Nz), dtype='i4')

	if(DEBUG > 0):	print("generate_interpolation_flux_on_wall: Completed")
