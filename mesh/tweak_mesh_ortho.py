import numpy			as np
from routines.globals	import DEBUG

def tweak_mesh_ortho(Config, dmin_ad, drmin_ad):

	if(DEBUG > 0): print("tweak_mesh_ortho")

	MagZones		= Config.MagZones
	MagMegazones	= Config.MagMegazones
	MagPMegazones	= Config.MagPMegazones

#	poloidal tweak only orthogonal mesh

	for k in range(len(MagMegazones)):
		IsEmpty	= True
		for k1 in range(len(MagMegazones[k].list)):
			nz = MagMegazones[k].list[k1]
			if(IsEmpty): 
				R			= np.copy(MagZones[nz].gridR[:,:-1])
				Z			= np.copy(MagZones[nz].gridZ[:,:-1])
				meshortho	= np.ones((MagZones[nz].gridZ.shape[1] - 1), dtype='i4')*MagZones[nz].meshortho
				IsEmpty	= False
			else:
				R		  = np.append(R, MagZones[nz].gridR[:,:-1], axis=1)
				Z  		  = np.append(Z, MagZones[nz].gridZ[:,:-1], axis=1)
				meshortho = np.append(meshortho, np.ones((MagZones[nz].gridZ.shape[1] - 1), dtype='i4')*MagZones[nz].meshortho)

		R = np.append(R, MagZones[nz].gridR[:,-1].reshape(R.shape[0],1), axis=1)
		Z = np.append(Z, MagZones[nz].gridZ[:,-1].reshape(Z.shape[0],1), axis=1)
		m = R.shape[0]
		p = R.shape[1]
		for i in range(m):
			R_ = R[i,:]
			Z_ = Z[i,:]
			dist = np.append(0., np.cumsum(np.sqrt((R_[1:]-R_[:-1])**2+(Z_[1:]-Z_[:-1])**2)))

			dmin	= dmin_ad*1.e-3/dist[-1]
			dist	= dist/dist[-1]
			distnew = np.copy(dist)
			for j in range(1,len(distnew)):
				dd = distnew[j]-distnew[j-1]
				if((dd < dmin) and (meshortho[j-1] == 1)): distnew[j]=distnew[j-1]+dmin

			distnew = distnew/distnew[-1]

			R_new = np.interp(distnew, dist, R_)
			Z_new = np.interp(distnew, dist, Z_)
			R[i,:] = R_new
			Z[i,:] = Z_new

		for k1 in range(len(MagMegazones[k].list)):
			nz				= MagMegazones[k].list[k1]
			m				= MagZones[nz].gridR.shape[0]
			p				= MagZones[nz].gridR.shape[1]
			MagZones[nz].gridR = R[:m,:p]
			MagZones[nz].gridZ = Z[:m,:p]
			R = R[:,p-1:]
			Z = Z[:,p-1:]

	if(DEBUG > 0): print("tweak_mesh_ortho: completed")
