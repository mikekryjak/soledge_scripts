import numpy				as np
from routines.globals		import DEBUG
#=========================================================
# This routine find triangles seguence at the wall
#=========================================================

def assign_non_penalized_cell(Zones, Triangles):

	if(DEBUG > 0):	print("assign_non_penalized_cell")

	Triangles.np_k = np.copy(Triangles.k)
	Triangles.np_i = np.copy(Triangles.i)
	Triangles.np_j = np.copy(Triangles.j)

	nstepmax = 2
	for kZone in range(len(Zones)):
		ik = np.where(Triangles.np_k == kZone)[0]
		if(len(ik) > 0):
			ik1 = np.where(Zones[kZone].Chi[Triangles.np_i[ik],Triagles,np_j[ik]] == 1)[0]
			if(len(ik1) > 0): ik = ik[ik1]
			else:			  ik = []

#		move in the poloidal direction to find a non penalized point

		for n in range(len(ik)):
			found		= False
			moveleft	= True
			moveright	= True
			step		= 0
			jlocl		= Triangles.np_j[ik[n]]
			jlocr		= Triangles.np_j[ik[n]]
			klocl		= kZone
			klocr		= kZone
			while((not found) and (step < nstepmax)):
				step +=1
				if(moveleft):
					if(jlocl > 0):												#not on the edge
						jlocl -= 1
					else:
						klocl = Zones[klocl].Neighbour.west
						if(klocl > -1):	jlocl = Zones[klocl].Nz-1
						else:			moveleft = False

					if(moveleft):												#test new point
						if(Zones[klocl].Chi[iloc,jlocl] == 0):
							found = True
							Triangles.np_k[ik[n]] = klocl
							Triangles.np_j[ik[n]] = jlocl

				if(moveright):
					if(jlocr < Zones[klocr].Nz-1):								#not on the edge
						jlocr += 1
					else:
						klocr=Zones[klocr].Neighbour.east
						if(klocr > -1):	jlocr = 0
						else: 			moveright = False

					if(moveright):											#test new point
						if(Zones[klocr].Chi[iloc,jlocr] == 0):
							found				  = True
							Triangles.np_k[ik[n]] = klocr
							Triangles.np_j[ik[n]] = jlocr


			if(not found):														#test radial direction
				moveleft	= True
				moveright	= True
				step		= 0
				ilocl		= Triangles.np_i[ik[n]]
				ilocr		= Triangles.np_i[ik[n]]
				klocl		= kZone
				klocr		= kZone
				while((not found) and (step < nstepmax)):
					step += 1
					if(moveleft):
						if(ilocl > 0):											#not on the edge
							ilocl -= 1
						else:
							klocl = Zones[kloc].Neighbour.south
							if(klocl > -1): ilocl 	 = Zones[kloc].Nx-1
							else:			moveleft = False

						if(moveleft):											#test new point
							if(Zones[kloc].Chi[ilocl,jloc] == 0):
								found				  = True
								Triangles.np_k[ik[n]] = klocl
								Triangles.np_i[ik[n]] = ilocl

					if(moveright):
						if(ilocr < Zones[klocr].Nx-1):							#not on the edge
							ilocr += 1
						else:
							klocr = Zones[klocr].Neighbour.north
							if(klocr > -1):	ilocr=1
							else:			moveright = False

						if(moveright):										#test new point
							if(Zones[klocr].Chi[ilocr,jloc] == 0):
								found				  =true
								Triangles.np_k[ik[n]] = klocr
								Triangles.np_i[ik[n]] = ilocr

#	second pass
#	=======

	for kZone in range(len(Zones)):
		ik = np.where(Triangles.np_k == kZone)[0]
		if(len(ik) > 0):
			ik1 = np.where(Zones[kZone].Chi[Triangles.np_i[ik],Triagles,np_j[ik]] == 1)[0]
			if(len(ik1) > 0): ik = ik[ik1]
			else:			  ik = []

		for n in range(len(ik)):
			found		= False

			itri = Triangles.neigh1[ik[n]]							#use surrounding triangles
			if(itri > -1):
				kloc2 = Triangles.np_k[itri]
				iloc2 = Triangles.np_i[itri]
				jloc2 = Triangles.np_j[itri]
				if(Zones[kloc2].Chi[iloc2,jloc2] == 0):
					Triangles.np_k[ik[n]] = kloc2
					Triangles.np_i[ik[n]] = iloc2
					Triangles.np_j[ik[n]] = jloc2
					found				  = True

			if(not found):
				itri = Triangles.neigh2[ik[n]]
				if(itri > -1):
					kloc2 = Triangles.np_k[itri]
					iloc2 = Triangles.np_i[itri]
					jloc2 = Triangles.np_j[itri]
					if(Zones[kloc2].Chi[iloc2,jloc2] == 0):
						Triangles.np_k[ik[n]] = kloc2
						Triangles.np_i[ik[n]] = iloc2
						Triangles.np_j[ik[n]] = jloc2

			if(not found):
				itri = Triangles.neigh3[ik[n]]
				if(itri > -1):
					kloc2 = Triangles.np_k[itri]
					iloc2 = Triangles.np_i[itri]
					jloc2 = Triangles.np_j[itri]
					if(Zones[kloc2].Chi[iloc2,jloc2] == 0):
						Triangles.np_k[ik[n]] = kloc2
						Triangles.np_i[ik[n]] = iloc2
						Triangles.np_j[ik[n]] = jloc2

			if(not found):
				print("Bad chi=1 removing with triangles ",ik[n])


	if(DEBUG > 0):	print("assign_non_penalized_cell: Completed")
