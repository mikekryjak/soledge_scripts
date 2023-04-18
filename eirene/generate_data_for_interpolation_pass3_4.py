import numpy 				as np
from routines.globals		import DEBUG
# Function definition is here
#=======================================================================================================================================

def generate_data_for_interpolation_pass3_4(Zones, RKnots, ZKnots, KnotsInterp):

	if(DEBUG > 0):	print("generate_data_for_interpolation_pass3_4")
	nZones	= len(Zones)
	NxMax	= 0
	NzMax	= 0
	for k in range(nZones):
		NxMax	= max(NxMax, Zones[k].gridRc.shape[0])
		NzMax	= max(NzMax, Zones[k].gridRc.shape[1])

	nKnots = len(RKnots)
	
#	treating pass 3 points
#	For each node look for closest node
	
	nn = np.where(KnotsInterp.nsol == 2); nn=nn[0]

	if(DEBUG > 1):	print("\t3 points conditions: processing",len(nn)," points")

	for n in nn:
		d		 = np.zeros(nKnots, dtype='f8') + 10000.
		nAssp	 = np.where(KnotsInterp.assp == 1); nAssp = nAssp[0]
		d[nAssp] = np.sqrt((RKnots[n]-RKnots[nAssp])**2+(ZKnots[n]-ZKnots[nAssp])**2) 
		d[n] 	 = 10000.

		KnotsInterp.eir[n, 0] = np.argmin(d)
		KnotsInterp.neir[n]   = 1
		
#	For each node look for closest soledge node

	nn = np.where(KnotsInterp.nsol == 1); nn=nn[0]
	if(DEBUG > 1):	print("\tSearch closest SOLEDGE node: processing",len(nn)," points")
	for n in nn:
		d		 = np.zeros(nKnots, dtype='f8') + 10000.
		nAssp	 = np.where(KnotsInterp.assp == 1); nAssp = nAssp[0]
		d[nAssp] = np.sqrt((RKnots[n]-RKnots[nAssp])**2+(ZKnots[n]-ZKnots[nAssp])**2) 
		d[n] 	 = 10000.
		
		p = np.argmin(d); m = d[p]
		d[p]  = 10000.	
		p2    = np.argmin(d); m2 = d[p2]
		dMins = np.array([m, m2, 0.])						#distance of two closest eirene points 

		d2	   = np.zeros((nZones,NxMax,NzMax), dtype='f8') + 10000.
		for k in range(nZones):
			ii, jj = np.where(Zones[k].Chi == 0)
			d2[k,ii,jj] = np.sqrt((RKnots[n]-Zones[k].gridRc[ii, jj])**2+(ZKnots[n]-Zones[k].gridZc[ii, jj])**2)

		kSol = KnotsInterp.sol[n, 0, 0]
		iSol = KnotsInterp.sol[n, 0, 1]
		jSol = KnotsInterp.sol[n, 0, 2]
		d2[kSol,iSol,jSol] = 10000.
		
		k1,i1,j1  = np.unravel_index(np.argmin(d2), d2.shape)
		dMins[2]  = d2[k1, i1, j1]

#		taking the two closest points
		B = np.argsort(dMins)
		if B[0] == 0:
			KnotsInterp.eir[n, KnotsInterp.neir[n]]    = p
			KnotsInterp.neir[n]						  += 1
		elif B[0] == 1:
#			impossible
			KnotsInterp.eir[n, KnotsInterp.neir[n]]    = p2
			KnotsInterp.neir[n]						  += 1
		elif B[0] == 2:
			KnotsInterp.sol[n, KnotsInterp.nsol[n], :] = np.array([k1, i1, j1])
			KnotsInterp.nsol[n]						  += 1

		if B[1] == 0:
			KnotsInterp.eir[n, KnotsInterp.neir[n]]    = p
			KnotsInterp.neir[n]						  += 1
		elif B[1] == 1:
			KnotsInterp.eir[n, KnotsInterp.neir[n]]    = p2
			KnotsInterp.neir[n]						  += 1
		elif B[1] == 2:
			KnotsInterp.sol[n, KnotsInterp.nsol[n], :] = np.array([k1, i1, j1])
			KnotsInterp.nsol[n]						  += 1

#	treating pass 4 points

	nn = np.where(KnotsInterp.assp == 4); nn=nn[0]
	if(DEBUG > 1):	print("\t4 points conditions: processing",len(nn)," points")
	for n in nn:
		d = np.zeros(nKnots, dtype='f8') + 10000.
		d		 = np.zeros(nKnots, dtype='f8') + 10000.
		nAssp	 = np.where(KnotsInterp.assp == 1); nAssp = nAssp[0]
		d[nAssp] = np.sqrt((RKnots[n]-RKnots[nAssp])**2+(ZKnots[n]-ZKnots[nAssp])**2) 

		d[n] 	 = 10000.
		p = np.argmin(d); m = d[p]

		d[p]  = 10000.
		p2    = np.argmin(d); m2 = d[p2]

		dMins = np.array([m, m2, 0., 0.])
	
#		looking for closest soledge points
		d2	   = np.zeros((nZones,NxMax,NzMax)) + 10000.
		for k in range(nZones):
			ii, jj = np.where(Zones[k].Chi == 0)
			d2[k, ii,jj] = np.sqrt((RKnots[n]-Zones[k].gridRc[ii, jj])**2+(ZKnots[n]-Zones[k].gridZc[ii, jj])**2)

		k1,i1,j1  = np.unravel_index(np.argmin(d2),d2.shape)
		dMins[2]  = d2[k1,i1,j1]

		d2[k1,i1,j1] = 10000.
		k2,i2,j2  = np.unravel_index(np.argmin(d2),d2.shape)

		dMins[3]  = d2[k2,i2,j2]

#		taking the three closest points
		B = np.argsort(dMins)
		if B[0] == 0:
			KnotsInterp.eir[n, KnotsInterp.neir[n]]    = p
			KnotsInterp.neir[n]						  += 1
		elif B[0] == 1:
#			impossible
			KnotsInterp.eir[n, KnotsInterp.neir[n]]	   = p2
			KnotsInterp.neir[n]						  += 1
		elif B[0] == 2:
			KnotsInterp.sol[n, KnotsInterp.nsol[n],:] =  np.array([k1, i1, j1])
			KnotsInterp.nsol[n]						 += 1
		elif B[0] == 3:
#			impossible
			KnotsInterp.sol[n, KnotsInterp.nsol[n], :] = np.array([k2, i2, j2])
			KnotsInterp.nsol[n]						  += 1

		if B[1] == 0:
			KnotsInterp.eir[n, KnotsInterp.neir[n]]    = p
			KnotsInterp.neir[n]						  += 1
		elif B[1] == 1:
			KnotsInterp.eir[n, KnotsInterp.neir[n]]	   = p2
			KnotsInterp.neir[n]						  += 1
		elif B[1] == 2:
			KnotsInterp.sol[n, KnotsInterp.nsol[n], :] = np.array([k1, i1, j1])
			KnotsInterp.nsol[n]						  += 1
		elif B[1] == 3:
			KnotsInterp.sol[n, KnotsInterp.nsol[n], :] = np.array([k2, i2, j2])
			KnotsInterp.nsol[n]						  += 1

		if B[2] == 0:
			KnotsInterp.eir[n, KnotsInterp.neir[n]]    = p
			KnotsInterp.neir[n]						  += 1
		elif B[2] == 1:
			KnotsInterp.eir[n, KnotsInterp.neir[n]]	   = p2
			KnotsInterp.neir[n]						  += 1
		elif B[2] == 2:
			KnotsInterp.sol[n, KnotsInterp.nsol[n], :] = np.array([k1, i1, j1])
			KnotsInterp.nsol[n]						  += 1
		elif B[2] == 3:
			KnotsInterp.sol[n, KnotsInterp.nsol[n], :] = np.array([k2, i2, j2])
			KnotsInterp.nsol[n]						  += 1

	if(DEBUG > 1):	print("\tmax_nsol=",np.max(KnotsInterp.nsol)," max_neir=",np.max(KnotsInterp.neir))
	if(DEBUG > 0):	print("generate_data_for_interpolation_pass3_4: Completed")

