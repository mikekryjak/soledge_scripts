import numpy as np

def find_center_position(Zones):

	Rce		= 0.
	Zce		= 0.
	weight	= 0.      
	for k in range(len(Zones)):
		if(Zones[k].Neighbour.south==-2):
			w		= np.abs(Zones[k].zb[0,1:] - Zones[k].zb[0,:-1])
			weight += np.sum(w)
			Rce	   += np.sum(Zones[k].gridRc[0,:]*w)
			Zce	   += np.sum(Zones[k].gridZc[0,:]*w)

	Rce /= weight
	Zce /= weight

	return Rce, Zce
