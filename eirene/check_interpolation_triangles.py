import numpy as np
from routines.globals			import DEBUG

# Function definition is here
#=======================================================================================================================================

#def check_interpolation_triangles(Zones, RKnots, ZKnots, KnotsInterp, Wall, PlotCheck=0):
def check_interpolation_triangles(Zones, RKnots, ZKnots, KnotsInterp):

	if(DEBUG > 0):	print("check_interpolation_triangles")

#	checking if points are not aligned

	nn=np.where(KnotsInterp.nsol <= 3); nn = nn[0]
	if(DEBUG > 1):	print("\tchecking interpolation for ",len(nn)," triangles")

	for n in nn:
		if KnotsInterp.nsol[n] == 3:
			k	= KnotsInterp.sol[n, 0, 0]
			i	= KnotsInterp.sol[n, 0, 1]
			j	= KnotsInterp.sol[n, 0, 2]
			p1r	= Zones[k].gridRc[i,j]
			p1z	= Zones[k].gridZc[i,j]

			k	= KnotsInterp.sol[n, 1, 0]
			i	= KnotsInterp.sol[n, 1, 1]
			j	= KnotsInterp.sol[n, 1, 2]
			p2r	= Zones[k].gridRc[i,j]
			p2z	= Zones[k].gridZc[i,j]

			k	= KnotsInterp.sol[n, 2, 0]
			i	= KnotsInterp.sol[n, 2, 1]
			j	= KnotsInterp.sol[n, 2, 2]
			p3r	= Zones[k].gridRc[i,j]
			p3z	= Zones[k].gridZc[i,j]
		elif KnotsInterp.nsol[n] == 2:
			k	= KnotsInterp.sol[n, 0, 0]
			i	= KnotsInterp.sol[n, 0, 1]
			j	= KnotsInterp.sol[n, 0, 2]
			p1r	= Zones[k].gridRc[i,j]
			p1z	= Zones[k].gridZc[i,j]

			k	= KnotsInterp.sol[n, 1, 0]
			i	= KnotsInterp.sol[n, 1, 1]
			j	= KnotsInterp.sol[n, 1, 2]
			p2r	= Zones[k].gridRc[i,j]
			p2z	= Zones[k].gridZc[i,j]

			p3r	= RKnots[KnotsInterp.eir[n, 0]]
			p3z	= ZKnots[KnotsInterp.eir[n, 0]]
		elif KnotsInterp.nsol[n] == 1:
			k	= KnotsInterp.sol[n, 0, 0]
			i	= KnotsInterp.sol[n, 0, 1]
			j	= KnotsInterp.sol[n, 0, 2]
			p1r	= Zones[k].gridRc[i,j]
			p1z	= Zones[k].gridZc[i,j]

			p2r	= RKnots[KnotsInterp.eir[n, 0]]
			p2z	= ZKnots[KnotsInterp.eir[n, 0]]
			p3r	= RKnots[KnotsInterp.eir[n, 1]]
			p3z = ZKnots[KnotsInterp.eir[n, 1]]
		else:
			if KnotsInterp.neir[n] != 0:
				p1r	= RKnots[KnotsInterp.eir[n, 0]]
				p1z	= ZKnots[KnotsInterp.eir[n, 0]]
				p2r	= RKnots[KnotsInterp.eir[n, 1]]
				p2z	= ZKnots[KnotsInterp.eir[n, 1]]
				p3r	= RKnots[KnotsInterp.eir[n, 2]]
				p3z	= ZKnots[KnotsInterp.eir[n, 2]]

		v1r	= p2r - p1r
		v1z	= p2z - p1z
		v2r	= p3r - p1r
		v2z	= p3z - p1z
		if (v1z != 0.) and (v2r != 0.):
			ratio = (v2z*v1r)/(v2r*v1z)
		elif (v1r != 0.) and (v2z != 0.):
			ratio = (v2r*v1z)/(v2z*v1r)
		elif (v1r != v2r) and (v1z != v2z):
			ratio = 2.
		else:
			ratio = 1.

		if((np.abs(ratio-1) <= 1e-3) and (DEBUG > 0)):
			print("\tpoints almost aligned for triangle ",n, " nsol=",KnotsInterp.nsol[n])
	"""
	if PlotCheck > 0:
		print("\tplot nodes")
		print("\t\tred:        nodes on one soledge quadrangle")
		print("\t\tblue:       nodes on two soledge quadrangles")
		print("\t\tyellow:     nodes on three soledge quadrangles")
		print("\t\tgreen:      nodes on four soledge quadrangles")
		print("\t\tmagenta-x:  nodes without soledge quadrangle")

		fig = pyp.figure()
		ax  = fig.add_subplot(111)
		ax.set_aspect(1.)
		ax.set_xlabel("R (m)")
		ax.set_ylabel("Z (m)")
		ax.set_title("Knots")

		nn=np.where(KnotsInterp.nsol == 1); nn = nn[0]
		if len(nn) > 0:	ax.plot(RKnots[nn], ZKnots[nn],'r.')

		nn=np.where(KnotsInterp.nsol == 2); nn = nn[0]
		if len(nn) > 0:	ax.plot(RKnots[nn], ZKnots[nn],'b.')

		nn=np.where(KnotsInterp.nsol == 3); nn = nn[0]
		if len(nn) > 0:	ax.plot(RKnots[nn], ZKnots[nn],'y.')

		nn=np.where(KnotsInterp.nsol == 4); nn = nn[0]
		if len(nn) > 0:	ax.plot(RKnots[nn], ZKnots[nn],'g.')

		nn=np.where((KnotsInterp.nsol == 0) & (KnotsInterp.neir != 0)); nn = nn[0]
		if len(nn) > 0:	ax.plot(RKnots[nn], ZKnots[nn],'mx')

		ax.plot(Wall[:,0],Wall[:,1],'b-')

		pyp.show()
	"""
	if(DEBUG > 0):	print("check_interpolation_triangles: Completed")
