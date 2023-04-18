from routines.globals			import DEBUG

#=======================================================================================================================================
# Function definition is here

#=========================================================
# This routine check direct or indirect orientation
#=========================================================
#
def zones_direct_or_indirect(Zones, RKnots, ZKnots):

	if(DEBUG > 0):	print("zones_direct_or_indirect")

	vector1x = RKnots[Zones[0].KnotC[0,0]] - RKnots[Zones[0].KnotA[0,0]]
	vector1y = ZKnots[Zones[0].KnotC[0,0]] - ZKnots[Zones[0].KnotA[0,0]]

	vector2x = RKnots[Zones[0].KnotB[0,0]] - RKnots[Zones[0].KnotA[0,0]]
	vector2y = ZKnots[Zones[0].KnotB[0,0]] - ZKnots[Zones[0].KnotA[0,0]]
	vecprod  = vector1x*vector2y-vector1y*vector2x

	if(vecprod >= 0):
		if(DEBUG > 1):	print("\tindirect orientation")
		direct=0
	else:
		if(DEBUG > 1):	print("\tdirect orientation")
		direct=1

	if(DEBUG > 0):	print("zones_direct_or_indirect: Completed")

	return direct

def triangles_direct_or_indirect(TriKnots, RKnots, ZKnots):

	if(DEBUG > 0):	print("triangles_direct_or_indirect")

	vector1x = RKnots[TriKnots[0,1]] - RKnots[TriKnots[0,0]]
	vector1y = ZKnots[TriKnots[0,1]] - ZKnots[TriKnots[0,0]]

	vector2x = RKnots[TriKnots[0,2]] - RKnots[TriKnots[0,0]]
	vector2y = ZKnots[TriKnots[0,2]] - ZKnots[TriKnots[0,0]]
	vecprod  = vector1x*vector2y-vector1y*vector2x

	if(vecprod >= 0):
		if(DEBUG > 1):	print("\tindirect orientation")
		direct=0
	else:
		if(DEBUG > 1):	print("\tdirect orientation")
		direct=1

	if(DEBUG > 0):	print("direct_or_indirect: Completed")

	return direct