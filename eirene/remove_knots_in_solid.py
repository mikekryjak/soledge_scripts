import numpy as np
from routines.globals			import DEBUG

# Function definition is here
#=======================================================================================================================================

def remove_knots_in_solid(Zones, RKnots, ZKnots, InPWallsKnots):

	if(DEBUG > 0):	print("remove_knots_in_solid")

	nZones			= len(Zones)
	nKnots			= len(RKnots)

	KnotsPlasma = np.sort(np.where(InPWallsKnots == 1)[0])

	RKnotsPlasma = np.copy(RKnots[KnotsPlasma])
	ZKnotsPlasma = np.copy(ZKnots[KnotsPlasma])
		
	ToKnotsPlasma				= np.zeros(nKnots, dtype='i4')				#Array to convert from old knots number to new KnotsPlasma number
	ToKnotsPlasma[KnotsPlasma]	= 1
	KnotsNone					= np.where(ToKnotsPlasma == 0)[0]
	ToKnotsPlasma 				= np.cumsum(ToKnotsPlasma)
	ToKnotsPlasma[KnotsNone] 	= 0
	ToKnotsPlasma			 	= np.int32(ToKnotsPlasma) - 1  				#index python style from 0

#	*** with respect MatLab version InPlasma was computed in routine find_inplasma.py ***
		
	if(DEBUG > 0):	print("\tnodes in plasma =",len(RKnotsPlasma))
	if(DEBUG > 1):	print("remove_knots_in_solid: Completed")

	return ToKnotsPlasma, RKnotsPlasma, ZKnotsPlasma
