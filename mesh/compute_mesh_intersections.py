#import matplotlib.pyplot 				as pyp

import numpy							as np
import scipy.interpolate				as interp
from math								import sqrt

from routines.intersect_contour			import intersect_2contours
from routines.contour_better			import contour_better
from routines.find_closest_segment		import find_closest_segment
from routines.get_in_out				import get_in_out
from routines.globals					import DEBUG, LENGTH_TYPES
from mesh.get_rz_core_sep				import get_rz_core_sep
from mesh.mag_to_zones_pos				import mag_to_zones_pos


def compute_mesh_intersections(Config, Cut, also_pos=False, use_mag_zones = True):

	if(DEBUG > 0): print("compute_mesh_intersections")

	MagZones	 = Config.MagZones
	MagMegazones = Config.MagMegazones

#	CutCore = False
	nLen	= 0
	for iMagZone in Cut.nz:
		nLen	+= MagZones[iMagZone].gridRc.shape[0]
#		if(MagMegazones[MagZones[iMagZone].mz].isperiodic): CutCore = True


	Lengths = np.empty((nLen, len(LENGTH_TYPES)), dtype='f8')						#Lengths array [Lenght, X, psi/psi_sep, rho_pol=sqrt(psi/psi_sep), area/area_sep]
	if(also_pos):
		INtCEll = np.empty((nLen,3), dtype = 'i4')
		IntRZ   = np.empty((nLen,2),  dtype = 'f8')

#	if(CutCore):
	RZcore, RZsep, CoreMegazone, SepMegazone, jSep = get_rz_core_sep(Config)
	RZcore = 0
	iSep = MagMegazones[SepMegazone].list[0]

	pts_in, pts_out = get_in_out(RZsep[:,0], RZsep[:,1], Config.r2D.reshape(Config.r2D.size), Config.z2D.reshape(Config.z2D.size))
	PsiAxis = np.min(Config.flux2D.reshape(Config.flux2D.size)[pts_in])
	RZsep   = 0; pts_out = 0; pts_in  = 0

	f_psi   = interp.RectBivariateSpline(Config.r2D[0,:], Config.z2D[:,0], Config.flux2D.T)								#RectBivariateSpline wants [x,y] order
	PsiSep, AreaSep = compute_psi_area(Config, f_psi, iSep, jSep)

	nLen = 0
	for i in range(Cut.nz.size):
		iMagZone	= Cut.nz[i]		

		nGrid = MagZones[iMagZone].gridR.shape[0]
		rz = np.empty((nGrid, 2), dtype='f8')
		rz[0, 0]  = Cut.r12[i,0]					
		rz[0, 1]  = Cut.z12[i,0]					
		rz[-1, 0] = Cut.r12[i,1]					
		rz[-1, 1] = Cut.z12[i,1]	
						
		for iGrid in range(1,nGrid-1):																#grid contours
			gridR, gridZ = extend_RZ_east_west(MagZones, iMagZone, iGrid)
			Ri, Zi, is1, is2  = intersect_2contours(Cut.r12[i,:], Cut.z12[i,:], gridR, gridZ)
			rz[iGrid, 0] = Ri[0]					
			rz[iGrid, 1] = Zi[0]
								
		dGrid = np.sqrt((rz[:,0] - Cut.r12[i,0])**2 + (rz[:,1] - Cut.z12[i,0])**2) + Cut.d12[i,0]
		dGrid = 0.5*(dGrid[:-1] + dGrid[1:])
		rz	  = 0.5*(rz[:-1,:] + rz[1:,:])

		if(Cut.dir[i] != 0):																		#reverse
			dGrid = dGrid[::-1]
			rz	  = rz[::-1,:]

		Nx		= MagZones[iMagZone].gridRc.shape[0]
		Lengths[nLen:nLen+Nx, 0] = dGrid															#0: Lenght
		Lengths[nLen:nLen+Nx, 1] = MagZones[iMagZone].x[:,0]										#1: X mesh

		if(also_pos):
			IntRZ[nLen:nLen+Nx, 0] = rz[:,0]
			IntRZ[nLen:nLen+Nx, 1] = rz[:,1]

		i1 	 = nLen
		for kk in range(Nx):
			if(Cut.dir[i] == 0): k = kk
			else:				 k = Nx - 1 - kk
			
			Psi, Area = compute_psi_area(Config, f_psi, iMagZone, k)

			Lengths[i1,2] = np.abs((Psi - PsiAxis)/(PsiSep - PsiAxis))								#2: Normalized poloidal flux does not need sqrt
			Lengths[i1,3] = np.sqrt(Lengths[i1,2])													#3: rho of poloidal flux needs sqrt

			if(MagMegazones[MagZones[iMagZone].mz].isperiodic): 									#inside separatrix
				Lengths[i1,4] = sqrt(Area/AreaSep)													#4: Rho area
			else:
				Lengths[i1,4] = np.nan																#4: Rho area

			if(also_pos):
				INtCEll[i1,0] = iMagZone
				INtCEll[i1,1] = k
				d, INtCEll[i1,2] = find_closest_segment(rz[kk,0], rz[kk,1], MagZones[iMagZone].gridRc[k,:], MagZones[iMagZone].gridZc[k,:])

			i1 += 1

		nLen += Nx

	if(DEBUG > 0): print("compute_mesh_intersections:Completed")
		

#	print("Lengths=",Lengths)
	if(also_pos):
		if(not use_mag_zones): INtCEll = mag_to_zones_pos(Config, INtCEll)

		return Lengths, IntRZ, INtCEll
	else:
		return Lengths

def compute_psi_area(Config, f_psi, iMagZone, kx):

	MagZones		= Config.MagZones

	Psi			= np.mean(f_psi.ev(MagZones[iMagZone].gridRc[kx,:], MagZones[iMagZone].gridZc[kx,:]))
	C_gridRc	= contour_better(Config.r2D, Config.z2D, Config.flux2D, [Psi])

	nC_gridRc = len(C_gridRc.arc)
	if(nC_gridRc > 1):
		ky = int(MagZones[iMagZone].gridZc.shape[1]/2)
		d  = np.empty(nC_gridRc, dtype = 'f8')
		for iC in range(nC_gridRc):
			d[iC] = np.min(np.sqrt((C_gridRc.arc[iC].x - MagZones[iMagZone].gridRc[kx,ky])**2 + (C_gridRc.arc[iC].y - MagZones[iMagZone].gridZc[kx,ky])**2))
		iC = np.argmin(d)
	else:
		iC = 0
	Area = abs(np.trapz(C_gridRc.arc[iC].y, x=C_gridRc.arc[iC].x))

	return Psi, Area

def extend_RZc_east_west(MagZones, iMagZone, kx):

	Rc = MagZones[iMagZone].gridRc[kx,:]
	Zc = MagZones[iMagZone].gridZc[kx,:]

	if((MagZones[iMagZone].Neighbour.east >= 0) and  (MagZones[iMagZone].Neighbour.east != iMagZone)):
		Rc = np.append(Rc,MagZones[MagZones[iMagZone].Neighbour.east].gridRc[kx,0])
		Zc = np.append(Zc,MagZones[MagZones[iMagZone].Neighbour.east].gridZc[kx,0])

	if((MagZones[iMagZone].Neighbour.west >= 0) and  (MagZones[iMagZone].Neighbour.west != iMagZone)):
		Rc = np.append(MagZones[MagZones[iMagZone].Neighbour.west].gridRc[kx,-1], Rc)
		Zc = np.append(MagZones[MagZones[iMagZone].Neighbour.west].gridZc[kx,-1], Zc)

	return Rc, Zc



def extend_RZ_east_west(MagZones, iMagZone, kx):

	R = MagZones[iMagZone].gridR[kx,:]
	Z = MagZones[iMagZone].gridZ[kx,:]

	if((MagZones[iMagZone].Neighbour.east >= 0) and  (MagZones[iMagZone].Neighbour.east != iMagZone)):
		R = np.append(R,MagZones[MagZones[iMagZone].Neighbour.east].gridR[kx,0])
		Z = np.append(Z,MagZones[MagZones[iMagZone].Neighbour.east].gridZ[kx,0])

	if((MagZones[iMagZone].Neighbour.west >= 0) and  (MagZones[iMagZone].Neighbour.west != iMagZone)):
		R = np.append(MagZones[MagZones[iMagZone].Neighbour.west].gridR[kx,-1], R)
		Z = np.append(MagZones[MagZones[iMagZone].Neighbour.west].gridZ[kx,-1], Z)

	return R, Z
