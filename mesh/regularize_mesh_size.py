import numpy				as np
from routines.globals		import DEBUG

def regularize_mesh_size(Config):

	if(DEBUG > 0): print("regularize_mesh_size")

	MagZones		= Config.MagZones
	MagMegazones	= Config.MagMegazones

	for k in range(len(MagMegazones)):
		IsEmpty	= True
		for k1 in range(len(MagMegazones[k].list)):
			nz = MagMegazones[k].list[k1]
			if(IsEmpty): 
				R		= np.copy(MagZones[nz].gridR[:,:-1])
				Z		= np.copy(MagZones[nz].gridZ[:,:-1])
				IsEmpty	= False
			else:
				R  = np.append(R, MagZones[nz].gridR[:,:-1], axis=1)			#Same size on axis=0
				Z  = np.append(Z, MagZones[nz].gridZ[:,:-1], axis=1)
				
		R = np.append(R, MagZones[nz].gridR[:,-1].reshape(R.shape[0],1),axis=1)
		Z = np.append(Z, MagZones[nz].gridZ[:,-1].reshape(Z.shape[0],1),axis=1)
		m = R.shape[0]
		p = R.shape[1]
		for i in range(m):
			R_	 	= R[i,:]
			Z_	 	= Z[i,:]
			dist 	= np.append(0., np.cumsum(np.sqrt((R_[1:]-R_[:-1])**2+(Z_[1:]-Z_[:-1])**2)))
			dist 	= dist/dist[-1]
			distnew = np.copy(dist)
			for j in range(1,len(distnew)-1):
				dd  = distnew[j]-distnew[j-1]
				dd2 = distnew[j+1]-distnew[j-1]
				rat = dd/dd2
				if(rat < 0.4):
					ratn = 0.4
				elif(rat > 0.6):
					ratn = 0.6
				else:
					ratn = rat

				distnew[j] = distnew[j-1]+ratn*dd2

			distnew = distnew/distnew[-1]

			R_new	= np.interp(distnew, dist, R_)
			Z_new	= np.interp(distnew, dist, Z_)
			R[i,:]	= R_new
			Z[i,:]	= Z_new

		for k1 in range(len(MagMegazones[k].list)):
			nz	= MagMegazones[k].list[k1]
			m	= MagZones[nz].gridR.shape[0]
			p	= MagZones[nz].gridR.shape[1]
			MagZones[nz].gridR	= R[:m,:p]
			MagZones[nz].gridZ = Z[:m,:p]
			R	= R[:,p-1:]
			Z	= Z[:,p-1:]

	if(DEBUG > 0): print("regularize_mesh_size: completed")
