import types
import numpy as np
from routines.globals			import DEBUG

def generate_knots(Zones, X_points):

	nZones		= len(Zones)
	nX_points	= len(X_points)

	for k in range(nZones):
		Nx				= Zones[k].Nx
		Nz				= Zones[k].Nz
		Zones[k].nKnots	= np.zeros((Nx+1,Nz+1), dtype='i4')
		Zones[k].Rknots	= np.zeros((Nx+1,Nz+1), dtype='f8')
		Zones[k].Zknots	= np.zeros((Nx+1,Nz+1), dtype='f8')
		Zones[k].V		= -np.ones((4), dtype='i4')						# tables to define if X point in zone corners
		Zones[k].V2		= -np.ones((4), dtype='i4')						# tables to define if normal point in zone corners

	tol = 1.e-5
	for nx in range(nX_points):
		for k in range(nZones):
			d = np.sqrt((Zones[k].gridR[0,0]-X_points[nx].R)**2+(Zones[k].gridZ[0,0]-X_points[nx].Z)**2)
			if(d < tol): Zones[k].V[0] = nx

			d = np.sqrt((Zones[k].gridR[-1,0]-X_points[nx].R)**2+(Zones[k].gridZ[-1,0]-X_points[nx].Z)**2)
			if(d < tol): Zones[k].V[1] = nx

			d = np.sqrt((Zones[k].gridR[-1,-1]-X_points[nx].R)**2+(Zones[k].gridZ[-1,-1]-X_points[nx].Z)**2)
			if(d < tol): Zones[k].V[2] = nx

			d = np.sqrt((Zones[k].gridR[0,-1]-X_points[nx].R)**2+(Zones[k].gridZ[0,-1]-X_points[nx].Z)**2)
			if(d < tol): Zones[k].V[3] = nx


#	count normal points

	RN = np.empty(0, dtype='f8')
	ZN = np.empty(0, dtype='f8')
	for k in range(nZones):
		if((Zones[k].Neighbour.north > 0) and (Zones[k].Neighbour.east > 0)):
			RN = np.append(RN, Zones[k].gridR[-1,-1])
			ZN = np.append(ZN, Zones[k].gridZ[-1,-1])

#	remove doublon
	RN_	= np.array([RN[0]])
	ZN_	= np.array([ZN[0]])
	for k in range(1,len(RN)):
		d = np.sqrt((RN[k]-RN_)**2+(ZN[k]-ZN_)**2)
		if(d.min() > tol):
			RN_ = np.append(RN_, RN[k])
			ZN_ = np.append(ZN_, ZN[k])

#	remove Xpoints
	RN = np.empty(0, dtype='f8')
	ZN = np.empty(0, dtype='f8')
	for k in range(len(RN_)):
		d = 1.e10
		for nx in range(nX_points):
			d = min(d, np.sqrt((RN_[k]-X_points[nx].R)**2+(ZN_[k]-X_points[nx].Z)**2))

		if(d > tol):
			RN = np.append(RN, RN_[k])
			ZN = np.append(ZN, ZN_[k])

	nN = len(RN)
	for nn in range(nN):
		for k in range(nZones):
			d = np.sqrt((Zones[k].gridR[0,0]-RN[nn])**2+(Zones[k].gridZ[0,0]-ZN[nn])**2)
			if(d < tol): Zones[k].V2[0] = nn

			d = np.sqrt((Zones[k].gridR[-1,0]-RN[nn])**2+(Zones[k].gridZ[-1,0]-ZN[nn])**2)
			if(d < tol): Zones[k].V2[1] = nn

			d = np.sqrt((Zones[k].gridR[-1,-1]-RN[nn])**2+(Zones[k].gridZ[-1,-1]-ZN[nn])**2)
			if(d < tol): Zones[k].V2[2] = nn

			d = np.sqrt((Zones[k].gridR[0,-1]-RN[nn])**2+(Zones[k].gridZ[0,-1]-ZN[nn])**2)
			if(d < tol): Zones[k].V2[3] = nn

	Nnknot	= np.arange(nN)	+ nX_points		 						# indexes 0...nN-1 (Python style)

#	X		= np.arange(nX_points, dtype='f8')
	Xnknot	= np.arange(nX_points, dtype='i4')
	RX		= np.array([X_points[k].R for k in range(nX_points)])
	ZX		= np.array([X_points[k].Z for k in range(nX_points)])

	nKnots	= nX_points + nN										# nX_points xpoints and nN normal intersection

	for k in range(nZones):
		Nx = Zones[k].Nx
		Nz = Zones[k].Nz
		Zones[k].RKnots = np.empty((Nx+1,Nz+1), dtype='f8')
		Zones[k].ZKnots = np.empty((Nx+1,Nz+1), dtype='f8')
		Zones[k].nKnots = np.ones((Nx+1,Nz+1), dtype='i4')*100000

		for i in range(1,Nx):
			for j in range(1,Nz):
				Zones[k].RKnots[i,j] = Zones[k].gridR[i,j]
				Zones[k].ZKnots[i,j] = Zones[k].gridZ[i,j]
				Zones[k].nKnots[i,j] = nKnots
				nKnots			    += 1

#		treating north

		if(Zones[k].MagNeighbour.north == 0):
			if(Zones[k].Neighbour.north < k): 
				iz = Zones[k].Neighbour.north
				for j in range(1,Nz):
					Zones[k].RKnots[Nx,j] = Zones[iz].RKnots[0,j]
					Zones[k].ZKnots[Nx,j] = Zones[iz].ZKnots[0,j]
					Zones[k].nKnots[Nx,j] = Zones[iz].nKnots[0,j]

			else:
				for j in range(1,Nz):
					Zones[k].RKnots[Nx,j] = Zones[k].gridR[Nx,j]
					Zones[k].ZKnots[Nx,j] = Zones[k].gridZ[Nx,j]
					Zones[k].nKnots[Nx,j] = nKnots
					nKnots				 += 1

		else:
			for j in range(1,Nz):
				Zones[k].RKnots[Nx,j] = Zones[k].gridR[Nx,j]
				Zones[k].ZKnots[Nx,j] = Zones[k].gridZ[Nx,j]
				Zones[k].nKnots[Nx,j] = nKnots
				nKnots				 += 1

#		treating south

		if(Zones[k].MagNeighbour.south == 0):
			if(Zones[k].Neighbour.south < k):
				iz = Zones[k].Neighbour.south
				imx = Zones[iz].Nx
				for j in range(1,Nz):
					Zones[k].RKnots[0,j] = Zones[iz].RKnots[imx,j]
					Zones[k].ZKnots[0,j] = Zones[iz].ZKnots[imx,j]
					Zones[k].nKnots[0,j] = Zones[iz].nKnots[imx,j]
			else:
				for j in range(1,Nz):
					Zones[k].RKnots[0,j] = Zones[k].gridR[0,j]
					Zones[k].ZKnots[0,j] = Zones[k].gridZ[0,j]
					Zones[k].nKnots[0,j] = nKnots
					nKnots			   += 1
		else:
			for j in range(1,Nz):
				Zones[k].RKnots[0,j] = Zones[k].gridR[0,j]
				Zones[k].ZKnots[0,j] = Zones[k].gridZ[0,j]
				Zones[k].nKnots[0,j] = nKnots
				nKnots			    += 1
		
#		treating east

		if(Zones[k].MagNeighbour.east == 0):
			if(Zones[k].Neighbour.east < k):
				iz = Zones[k].Neighbour.east
				for i in range(1,Nx):
					Zones[k].RKnots[i,Nz] = Zones[iz].RKnots[i,0]
					Zones[k].ZKnots[i,Nz] = Zones[iz].ZKnots[i,0]
					Zones[k].nKnots[i,Nz] = Zones[iz].nKnots[i,0]

			else:
				for i in range(1,Nx):
					Zones[k].RKnots[i,Nz] = Zones[k].gridR[i,Nz]
					Zones[k].ZKnots[i,Nz] = Zones[k].gridZ[i,Nz]
					Zones[k].nKnots[i,Nz] = nKnots
					nKnots				+= 1

		else:
			for i in range(1,Nx):
				Zones[k].RKnots[i,Nz] = Zones[k].gridR[i,Nz]
				Zones[k].ZKnots[i,Nz] = Zones[k].gridZ[i,Nz]
				Zones[k].nKnots[i,Nz] = nKnots
				nKnots				+= 1

#		treating west

		if(Zones[k].MagNeighbour.west == 0):
			if(Zones[k].Neighbour.west <= k):
				iz = Zones[k].Neighbour.west
				imz = Zones[iz].Nz
				for i in range(1,Nx):
					Zones[k].RKnots[i,0] = Zones[iz].RKnots[i,imz]
					Zones[k].ZKnots[i,0] = Zones[iz].ZKnots[i,imz]
					Zones[k].nKnots[i,0] = Zones[iz].nKnots[i,imz]

			else:
				for i in range(1,Nx):
					Zones[k].RKnots[i,0] = Zones[k].gridR[i,0]
					Zones[k].ZKnots[i,0] = Zones[k].gridZ[i,0]
					Zones[k].nKnots[i,0] = nKnots
					nKnots			   += 1

		else:
			for i in range(1,Nx):
				Zones[k].RKnots[i,0] = Zones[k].gridR[i,0]
				Zones[k].ZKnots[i,0] = Zones[k].gridZ[i,0]
				Zones[k].nKnots[i,0] = nKnots
				nKnots			   += 1

		
#		treating corner SW

		if (Zones[k].V2[0] == -1) and (Zones[k].V[0] == -1):
			nex = 0
			if(Zones[k].MagNeighbour.west == 0):										#something west
				if(Zones[k].Neighbour.west < k):										#corner SW = corner Neigh SE
					nex = 1
					iz  = Zones[k].Neighbour.west
					imz = Zones[iz].Nz
					Zones[k].RKnots[0,0] = Zones[iz].RKnots[0,imz]
					Zones[k].ZKnots[0,0] = Zones[iz].ZKnots[0,imz]
					Zones[k].nKnots[0,0] = Zones[iz].nKnots[0,imz]
					
			if(Zones[k].MagNeighbour.south == 0):											#something south
				if(Zones[k].Neighbour.south < k):
					nex = 1
					iz  = Zones[k].Neighbour.south
					imx = Zones[iz].Nx
					Zones[k].RKnots[0,0] = Zones[iz].RKnots[imx,0]
					Zones[k].ZKnots[0,0] = Zones[iz].ZKnots[imx,0]
					Zones[k].nKnots[0,0] = Zones[iz].nKnots[imx,0]
			if(nex == 0):																	#create a new knot
				Zones[k].RKnots[0,0] = Zones[k].gridR[0,0]
				Zones[k].ZKnots[0,0] = Zones[k].gridZ[0,0]
				Zones[k].nKnots[0,0] = nKnots
				nKnots			   += 1
		else:
			if(Zones[k].V2[0] >= 0):
				iz = Zones[k].V2[0]
				Zones[k].RKnots[0,0] = RN[iz]
				Zones[k].ZKnots[0,0] = ZN[iz]
				Zones[k].nKnots[0,0] = Nnknot[iz]

			if(Zones[k].V[0] >= 0):
				iz = Zones[k].V[0]
				Zones[k].RKnots[0,0] = RX[iz]
				Zones[k].ZKnots[0,0] = ZX[iz]
				Zones[k].nKnots[0,0] = Xnknot[iz]

#		treating corner SE

		if (Zones[k].V2[3] == -1) and (Zones[k].V[3] == -1):
			nex = 0
			if(Zones[k].MagNeighbour.east == 0):											#something east
				if(Zones[k].Neighbour.east <= k):
					nex = 1
					iz = Zones[k].Neighbour.east
					Zones[k].RKnots[0,Nz] = Zones[iz].RKnots[0,0]
					Zones[k].ZKnots[0,Nz] = Zones[iz].ZKnots[0,0]
					Zones[k].nKnots[0,Nz] = Zones[iz].nKnots[0,0]
					
			if(Zones[k].MagNeighbour.south == 0):											#something south
				if(Zones[k].Neighbour.south < k):
					nex = 1
					iz  = Zones[k].Neighbour.south
					imx = Zones[iz].Nx
					imz = Zones[iz].Nz
					Zones[k].RKnots[0,Nz] = Zones[iz].RKnots[imx,imz]
					Zones[k].ZKnots[0,Nz] = Zones[iz].ZKnots[imx,imz]
					Zones[k].nKnots[0,Nz] = Zones[iz].nKnots[imx,imz]
				
			if(nex == 0):																	#create a new knot
				Zones[k].RKnots[0,Nz] = Zones[k].gridR[0,Nz]
				Zones[k].ZKnots[0,Nz] = Zones[k].gridZ[0,Nz]
				Zones[k].nKnots[0,Nz] = nKnots
				nKnots				+= 1
		else:
			if(Zones[k].V2[3] >= 0):
				iz = Zones[k].V2[3]
				Zones[k].RKnots[0,Nz] = RN[iz]
				Zones[k].ZKnots[0,Nz] = ZN[iz]
				Zones[k].nKnots[0,Nz] = Nnknot[iz]

			if(Zones[k].V[3] >= 0):
				iz = Zones[k].V[3]
				Zones[k].RKnots[0,Nz] = RX[iz]
				Zones[k].ZKnots[0,Nz] = ZX[iz]
				Zones[k].nKnots[0,Nz] = Xnknot[iz]
		
#		treating corner NE

		if (Zones[k].V2[2] == -1) and (Zones[k].V[2] == -1):
			nex = 0
			if(Zones[k].MagNeighbour.east == 0):											#something east
				if(Zones[k].Neighbour.east < k):
					nex = 1
					iz	= Zones[k].Neighbour.east
					imx = Zones[iz].Nx
					Zones[k].RKnots[Nx,Nz] = Zones[iz].RKnots[imx,0]
					Zones[k].ZKnots[Nx,Nz] = Zones[iz].ZKnots[imx,0]
					Zones[k].nKnots[Nx,Nz] = Zones[iz].nKnots[imx,0]
					
			if(Zones[k].MagNeighbour.north == 0):											#something north
				if Zones[k].Neighbour.north < k:
					nex = 1
					iz 	= Zones[k].Neighbour.north
					imz = Zones[iz].Nz
					Zones[k].RKnots[Nx,Nz] = Zones[iz].RKnots[0, imz]
					Zones[k].ZKnots[Nx,Nz] = Zones[iz].ZKnots[0, imz]
					Zones[k].nKnots[Nx,Nz] = Zones[iz].nKnots[0, imz]
					
			if(nex == 0):																	#create a new knot
				Zones[k].RKnots[Nx,Nz] = Zones[k].gridR[Nx,Nz]
				Zones[k].ZKnots[Nx,Nz] = Zones[k].gridZ[Nx,Nz]
				Zones[k].nKnots[Nx,Nz] = nKnots
				nKnots				 += 1

		else:
			if(Zones[k].V2[2] >= 0):
				iz = Zones[k].V2[2]
				Zones[k].RKnots[Nx,Nz] = RN[iz]
				Zones[k].ZKnots[Nx,Nz] = ZN[iz]
				Zones[k].nKnots[Nx,Nz] = Nnknot[iz]

			if(Zones[k].V[2] >= 0):
				iz = Zones[k].V[2]
				Zones[k].RKnots[Nx,Nz] = RX[iz]
				Zones[k].ZKnots[Nx,Nz] = ZX[iz]
				Zones[k].nKnots[Nx,Nz] = Xnknot[iz]
		
#		treating corner NW

		if (Zones[k].V2[1] == -1) and (Zones[k].V[1] == -1):
			nex = 0
			if(Zones[k].MagNeighbour.west == 0):										#something west
				if(Zones[k].Neighbour.west <= k):
					nex = 1
					iz = Zones[k].Neighbour.west
					imx = Zones[iz].Nx
					imz = Zones[iz].Nz
					Zones[k].RKnots[Nx,0] = Zones[iz].RKnots[imx,imz]
					Zones[k].ZKnots[Nx,0] = Zones[iz].ZKnots[imx,imz]
					Zones[k].nKnots[Nx,0] = Zones[iz].nKnots[imx,imz]
					
			if(Zones[k].MagNeighbour.north == 0):										#something north
				if(Zones[k].Neighbour.north < k):
					nex = 1
					iz 	= Zones[k].Neighbour.north
					imz = Zones[iz].Nz
					Zones[k].RKnots[Nx,0] = Zones[iz].RKnots[0, 0]
					Zones[k].ZKnots[Nx,0] = Zones[iz].ZKnots[0, 0]
					Zones[k].nKnots[Nx,0] = Zones[iz].nKnots[0, 0]
					
			if(nex == 0):																	#create a new knot
				Zones[k].RKnots[Nx,0] = Zones[k].gridR[Nx,0]
				Zones[k].ZKnots[Nx,0] = Zones[k].gridZ[Nx,0]
				Zones[k].nKnots[Nx,0] = nKnots
				nKnots				+= 1
		else:
			if(Zones[k].V2[1] >= 0):
				iz = Zones[k].V2[1]
				Zones[k].RKnots[Nx,0] = RN[iz]
				Zones[k].ZKnots[Nx,0] = ZN[iz]
				Zones[k].nKnots[Nx,0] = Nnknot[iz]

			if(Zones[k].V[1] >= 0):
				iz = Zones[k].V[1]
				Zones[k].RKnots[Nx,0] = RX[iz]
				Zones[k].ZKnots[Nx,0] = ZX[iz]
				Zones[k].nKnots[Nx,0] = Xnknot[iz]

#	save the eirene nodes in the corners of the soledge mesh

	for k in range(nZones):
		Nx = Zones[k].Nx
		Nz = Zones[k].Nz
		Zones[k].KnotA = Zones[k].nKnots[:Nx,	:Nz]
		Zones[k].KnotB = Zones[k].nKnots[1:Nx+1, :Nz]
		Zones[k].KnotC = Zones[k].nKnots[1:Nx+1, 1:Nz+1]
		Zones[k].KnotD = Zones[k].nKnots[:Nx,	1:Nz+1]


#		remaining node (normal intersection of 4 zones)

	if(DEBUG > 1): print("\tstored node number and (R,Z) coordinates for all nodes")
	
	nKnots = np.array([], dtype='i4')
	RKnots = np.array([], dtype='f8')
	ZKnots = np.array([], dtype='f8')
	
	for k in range(nZones):
		ii, jj = np.where(Zones[k].nKnots > -1)
		nKnots = np.append(nKnots, Zones[k].nKnots[ii,jj])
		RKnots = np.append(RKnots, Zones[k].RKnots[ii,jj])
		ZKnots = np.append(ZKnots, Zones[k].ZKnots[ii,jj])

	iSort	= np.argsort(nKnots)
	nKnots	= nKnots[iSort]
	RKnots	= RKnots[iSort]
	ZKnots	= ZKnots[iSort]

	iDiff	= nKnots[1:len(nKnots)] - nKnots[:len(nKnots)-1]
	iDiff	= np.append(np.array([1]), iDiff)
	ll		= np.where(iDiff != 0); ll=ll[0]					#index python style
	
	return nKnots[ll], RKnots[ll], ZKnots[ll]