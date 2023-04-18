import numpy as np
from routines.globals			import DEBUG
# Function definition is here
#=======================================================================================================================================

def finalization_remove(Zones, ToKnotsPlasma, RKnotsPlasma, ZKnotsPlasma, ToKnotsWall, RKnotsWall, ZKnotsWall):

	if(DEBUG > 0):	print("finalization_remove")

	nZones		 = len(Zones)
	nKnotsPlasma = RKnotsPlasma.shape[0]
	nKnotsWall   = RKnotsWall.shape[0]
	
#	Convert from old plasma node number to the new node number (after removal of nodes outside wall)

	for k in range(nZones):
		Zones[k].KnotA = ToKnotsPlasma[Zones[k].KnotA]
		Zones[k].KnotB = ToKnotsPlasma[Zones[k].KnotB]
		Zones[k].KnotC = ToKnotsPlasma[Zones[k].KnotC]
		Zones[k].KnotD = ToKnotsPlasma[Zones[k].KnotD]
	
#	Convert old wall nodes number to new node number (after removal of to close nodes)
	
		ii,jj = np.where(Zones[k].KnotE >= 0)
		if len(ii) > 0: Zones[k].KnotE[ii,jj] = ToKnotsWall[Zones[k].KnotE[ii,jj]]+ nKnotsPlasma
				
		ii,jj = np.where(Zones[k].KnotF >= 0)
		if len(ii) > 0: Zones[k].KnotF[ii,jj] = ToKnotsWall[Zones[k].KnotF[ii,jj]]+ nKnotsPlasma
				
		ii,jj = np.where(Zones[k].KnotG >= 0)
		if len(ii) > 0: Zones[k].KnotG[ii,jj] = ToKnotsWall[Zones[k].KnotG[ii,jj]]+ nKnotsPlasma
			
		ii,jj = np.where(Zones[k].KnotH >= 0)
		if len(ii) > 0: Zones[k].KnotH[ii,jj] = ToKnotsWall[Zones[k].KnotH[ii,jj]]+ nKnotsPlasma

#	Merge (R,Z) knots arrays inside wall and at the wall

	nKnots = nKnotsPlasma + nKnotsWall
	
	RKnots = np.empty(nKnots, dtype='f8')
	ZKnots = np.empty(nKnots, dtype='f8')

	RKnots[0:nKnotsPlasma] = RKnotsPlasma; RKnots[nKnotsPlasma:nKnots] = RKnotsWall
	ZKnots[0:nKnotsPlasma] = ZKnotsPlasma; ZKnots[nKnotsPlasma:nKnots] = ZKnotsWall
	
	Knots = np.arange(nKnots, dtype='i4')

#	Fix quadrangles with a single corner on wall

	NextSide= [1,2,3,0]
	for k in range(nZones):
		KnotsC	   = [Zones[k].KnotB, Zones[k].KnotC, Zones[k].KnotD, Zones[k].KnotA]
		KnotsM	   = [Zones[k].KnotE, Zones[k].KnotF, Zones[k].KnotG, Zones[k].KnotH, Zones[k].KnotE]
		InPlasmasC = [Zones[k].InPlasmaB, Zones[k].InPlasmaC, Zones[k].InPlasmaD, Zones[k].InPlasmaA]
		for iPwall in range(Zones[k].IsCrossed.shape[3]):
			for ls in range(4):
				iOnWall, jOnWall  = np.where(InPlasmasC[ls] == 2)
				for kOn in range(len(iOnWall)):
					i = iOnWall[kOn]
					j = jOnWall[kOn]
					if((Zones[k].IsCrossed[i,j,ls,iPwall] == 1) and (Zones[k].IsCrossed[i,j,NextSide[ls],iPwall] == 1)):
						InPlasmasC[ls] = 1									#set node inside plasma						
						Zones[k].IsCrossed[i,j,ls,iPwall] = 0				#remove crossing (but not used after this routine)
						Zones[k].IsCrossed[i,j,NextSide[ls],iPwall] = 0		# ""
						KnotsC[ls][i,j]		= KnotsM[ls][i,j]				#set node at quadrange corner
						KnotsM[ls][i,j]		= -1							#remove nodes in quadrangla sides
						KnotsM[ls+1][i,j]	= -1
						if(DEBUG > 2):
							Corner  = ["B","C","D","A"]
							print("\tIsolated point on wall at: k,i,j",k+1,i+1,j+1) 
							print("\tCorner:           ",Corner[ls])
							print("\tR isolated point: ",RKnots[KnotsC[ls][i,j]])
							print("\tZ isolated point: ",ZKnots[KnotsC[ls][i,j]])

	if(DEBUG > 0):	print("finalization_remove: Completed")

	return Knots, RKnots, ZKnots