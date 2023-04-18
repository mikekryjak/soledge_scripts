from routines.globals	import DEBUG

def conform_mesh(Config):

	if(DEBUG > 0): print("conform_mesh")


	X_points		= Config.X_points
	MagZones		= Config.MagZones
	MagMegazones	= Config.MagMegazones
	MagPMegazones	= Config.MagPMegazones

	for k in range(len(MagPMegazones)):
		for k1 in range(len(MagPMegazones[k].list)-1):
			nz1 = MagPMegazones[k].list[k1]
			nz2 = MagPMegazones[k].list[k1+1]
			MagZones[nz2].gridR[0,:] = MagZones[nz1].gridR[-1,:]
			MagZones[nz2].gridZ[0,:] = MagZones[nz1].gridZ[-1,:]

	for k in range(len(MagMegazones)):
		for k1 in range(len(MagMegazones[k].list)-1):
			nz1 = MagMegazones[k].list[k1]
			nz2 = MagMegazones[k].list[k1+1]
			MagZones[nz2].gridR[:,0] = MagZones[nz1].gridR[:,-1]
			MagZones[nz2].gridZ[:,0] = MagZones[nz1].gridZ[:,-1]


	for k in range(len(MagZones)):
		if((int(MagZones[k].pA.coord[2])  ==  -1) and (int(MagZones[k].pA.coord[3]) == 1)):
			MagZones[k].gridR[0,0] = X_points[int(MagZones[k].pA.coord[1])].R
			MagZones[k].gridZ[0,0] = X_points[int(MagZones[k].pA.coord[1])].Z

		if((int(MagZones[k].pB.coord[2]) == -1) and (int(MagZones[k].pB.coord[3]) == 1)):
			MagZones[k].gridR[-1,0] = X_points[int(MagZones[k].pB.coord[1])].R
			MagZones[k].gridZ[-1,0] = X_points[int(MagZones[k].pB.coord[1])].Z

		if((int(MagZones[k].pC.coord[2]) == -1) and (int(MagZones[k].pC.coord[3]) == 1)):
			MagZones[k].gridR[-1,-1] = X_points[int(MagZones[k].pC.coord[1])].R
			MagZones[k].gridZ[-1,-1] = X_points[int(MagZones[k].pC.coord[1])].Z

		if((int(MagZones[k].pD.coord[2]) == -1) and (int(MagZones[k].pD.coord[3]) == 1)):
			MagZones[k].gridR[0,-1] = X_points[int(MagZones[k].pD.coord[1])].R
			MagZones[k].gridZ[0,-1] = X_points[int(MagZones[k].pD.coord[1])].Z


	if(DEBUG > 0): print("conform_mesh: completed")
