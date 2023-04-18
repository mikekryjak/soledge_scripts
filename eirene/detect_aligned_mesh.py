import numpy as np
from routines.globals			import DEBUG
#=======================================================================================================================================
# Function definition is here

def detect_aligned_mesh(Zones, Walls, CheckAlligned=False):

	if(DEBUG > 0):	print("detect_aligned_mesh")


	nZones	= len(Zones)
	if(CheckAlligned):

		nsw		= Wall.shape[0]

		eps=1.e-3
		for k in range(nZones):
			Nx = Zones[k].Nx
			Nz = Zones[k].Nz
			for i in range(Nx):
				for j in range(Nz):
					if	(Zones[k].IsCrossed[i,j] == 1) or (Zones[k].IsCrossed[np.max(1,i-1),j] == 1) or \
						(Zones[k].IsCrossed[np.min(i+1,Nx-1),j] == 1) or (Zones[k].IsCrossed[i,np.max(1,j-1)]==1) or \
						(Zones[k].IsCrossed[i,np.min(j+1,Nz-1)] == 1):
						for n in range(nsw-1):
							if (Wall[n,0]-Wall[n+1,0]) != 0:
#								segment Wall
								a = (Wall[n,1]-Wall[n+1,1])/(Wall[n,0]-Wall[n+1,0])
								b = Wall[n,1]-a*Wall[n,0]

#								segment AB
								if (Zones[k].Rr[i+1,j]-Zones[k].Rr[i,j]) != 0.:
									a1 = (Zones[k].Zr[i+1,j] - Zones[k].Zr[i,j])/(Zones[k].Rr[i+1,j] - Zones[k].Rr[i,j])
									b1 = Zones[k].Zr[i,j] - a1*Zones[k].Rr[i,j]
									if np.sqrt([a-a1]^2/a^2+[b-b1]^2/b^2) < eps: Zones.IsAligned[i,j] = 1

#								segment BC
								if (Zones[k].Rr[i+1,j]-Zones[k].Rr[i+1,j+1]) != 0.:
									a1 = (Zones[k].Zr[i+1,j] - Zones[k].Zr[i+1,j+1])/(Zones[k].Rr[i+1,j] - Zones[k].Rr[i+1,j+1])
									b1 = Zones[k].Zr[i+1,j]-a1*Zones[k].Rr[i+1,j]
									if np.sqrt([a-a1]^2/a^2+[b-b1]^2/b^2) < eps: Zones[k].IsAligned[i,j] = 1

#								segment CD
								if (Zones[k].Rr[i+1,j+1]-Zones[k].Rr[i,j+1]) != 0.:
									a1 = (Zones[k].Zr[i+1,j+1] - Zones[k].Zr[i,j+1])/(Zones[k].Rr[i+1,j+1] - Zones[k].Rr[i,j+1])
									b1 = Zones[k].Zr[i+1,j+1] - a1*Zones[k].Rr[i+1,j+1]
									if  np.sqrt((a-a1)^2/a^2+(b-b1)^2/b^2) < eps: Zones[k].IsAligned[i,j]=1

#								segment AD
								if (Zones[k].Rr[i,j+1]-Zones[k].Rr[i,j]) != 0.:
									a1 = (Zones[k].Zr[i,j+1] - Zones[k].Zr[i,j])/(Zones[k].Rr[i,j+1] - Zones[k].Rr[i,j])
									b1 = Zones[k].Zr[i,j] - a1*Zones[k].Rr[i,j]
									if np.sqrt((a-a1)^2/a^2+(b-b1)^2/b^2) < eps: Zones[k].IsAligned[i,j] = 1
							else:
#								segment AB
								if (Zones[k].Rr[i+1,j] - Zones[k].Rr[i,j]) == 0.:
									if np.sqrt((Zones[k].Rr[i,j] - Wall[n,0])^2) < eps: Zones[k].IsAligned[i,j] = 1

#								segment BC
								if (Zones[k].Rr[i+1,j] - Zones[k].Rr[i+1,j+1]) == 0.:
									if np.sqrt((Zones[k].Rr[i+1,j] - Wall[n,0])^2) < eps:  Zones[k].IsAligned[i,j] = 1

#								segment CD
								if (Zones[k].Rr[i+1,j+1] - Zones[k].Rr[i,j+1]) == 0.:
									if np.sqrt((Zones[k].Rr[i+1,j+1] - Wall[n,0])^2) < eps: Zones[k].IsAligned[i,j] = 1

#								segment AD
								if (Zones[k].Rr[i,j] - Zones[k].Rr[i,j+1]) == 0.:
									if np.sqrt((Zones[k].Rr[i,j] - Wall[n,0])^2) < eps: 	Zones[k].IsAligned[i,j] = 1

#	Allignement checking non necessary
	else:
		if(DEBUG > 1):	print("\tnot checking alligned mesh")
		for k in range(nZones):
			Zones[k].IsAligned = np.zeros(Zones[k].gridRc.shape, dtype='i4')

	if(DEBUG > 0):	print("detect_aligned_mesh: Completed")
