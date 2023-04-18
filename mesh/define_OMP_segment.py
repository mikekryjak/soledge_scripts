import types
import numpy						as np
from routines.find_core				import find_core
from routines.contour_better		import contour_better
from routines.intersect_contour		import intersect_contour
from mesh.follow_grad				import follow_grad
from mesh.get_magzone_of_point		import get_magzone_of_point
from routines.globals			import DEBUG

def define_OMP_segment(Config):

	if(DEBUG > 0): print("define_OMP_segment")

	if(hasattr(Config, "OMPseg")): return True

	MagZones		= Config.MagZones
	X_points		= Config.X_points
	MagMegazones	= Config.MagMegazones
	MagPMegazones	= Config.MagPMegazones
	Config.OMPseg	= types.SimpleNamespace()
	OMPseg			= Config.OMPseg

	Rcore, Zcore, Bcore = find_core(Config)						# Start by locating the core

#	Now find periodic Megazone

	iMzCore = -1
	for iMz in range(len(MagMegazones)):
		if(MagMegazones[iMz].isperiodic): 
			iMzCore = iMz
			break

	if(iMzCore == -1):
		print("\tERROR: Unable to find core")
		return false

#	Now find the intersection between the horizontal line from the core and north side of core Zone

	ccore 	   		= types.SimpleNamespace()
	ccore.arc 		= [types.SimpleNamespace()]
	ccore.arc[0].x	= np.array([Rcore, np.max(Config.r2D)])
	ccore.arc[0].y	= np.array([Zcore, Zcore])

	c2 	   		= types.SimpleNamespace()
	c2.arc 		= [types.SimpleNamespace()]
	for iz in MagMegazones[iMzCore].list:
		c2.arc[0].x = MagZones[iz].north.R
		c2.arc[0].y = MagZones[iz].north.Z
		X = intersect_contour(ccore, c2)
		if(len(X) >  0):
			if(len(X) >  1): ROMP = np.min(np.array([x.x for x in X]))
			else: 				ROMP = X[0].x
			izoneOMP = iz
			break

#	Find the Pmegazone this point belongs to, as well as the order list of zones and megazones

	pmz = MagZones[izoneOMP].pmz
	OMPseg.zonelist = np.empty(0, dtype='i4')
	OMPseg.mzlist	= np.empty(0, dtype='i4')
	for iz in range(len(MagPMegazones[pmz].list)):
		izone = MagPMegazones[pmz].list[iz]
		OMPseg.zonelist = np.append(OMPseg.zonelist, izone)
		OMPseg.mzlist	= np.append(OMPseg.mzlist, MagZones[izone].mz)

#	Extract psimin and psimax of Pmegazone

	psimin = MagZones[MagPMegazones[pmz].list[ 0]].pA.coord[0]
	psimax = MagZones[MagPMegazones[pmz].list[-1]].pB.coord[0]

#	Follow the gradient of psi between psimin and psimax to define the segment
	R1, Z1 = follow_grad(Config, ROMP, Zcore, psimin)
	R2, Z2 = follow_grad(Config, ROMP, Zcore, psimax)
	OMPseg.R = np.append(R1[:0:-1], R2)
	OMPseg.Z = np.append(Z1[:0:-1], Z2)

#	Find intersections of the segment with zone sides on its way
	c2 	   		= types.SimpleNamespace()
	c2.arc 		= [types.SimpleNamespace()]
	c1 	   		= types.SimpleNamespace()
	c1.arc 		= [types.SimpleNamespace()]
	c1.arc[0].x = OMPseg.R
	c1.arc[0].y = OMPseg.Z
	NzonesOMP = len(OMPseg.zonelist)
	OMPseg.Intersec_R = np.empty((NzonesOMP+1), dtype = 'f8')
	OMPseg.Intersec_Z = np.empty((NzonesOMP+1), dtype = 'f8')
	for iz in range(NzonesOMP):
		izone = OMPseg.zonelist[iz]
		c2.arc[0].x = MagZones[izone].south.R
		c2.arc[0].y = MagZones[izone].south.Z
		X = intersect_contour(c1,c2)
		if(len(X) > 0):
			if(len(X) > 1): 	iRmin = np.argmin(np.array([x.x for x in X]))
			else:					iRmin = 0
			OMPseg.Intersec_R[iz] = X[iRmin].x
			OMPseg.Intersec_Z[iz] = X[iRmin].y
		else:
			print("\tERRROR define_OMP_segment: not found intersection with zone=",izone+1)

	izone = OMPseg.zonelist[-1]
	c2.arc[0].x = MagZones[izone].north.R
	c2.arc[0].y = MagZones[izone].north.Z
	X = intersect_contour(c1,c2)
	if(len(X) > 1): iRmin = np.argmin(np.array([x.x for x in X]))
	else:			iRmin = 0
	OMPseg.Intersec_R[-1] = X[iRmin].x
	OMPseg.Intersec_Z[-1] = X[iRmin].y

	OMPseg.R[0]  = OMPseg.Intersec_R[0]					# Correct segment so that its end point are exactly at the intersections with extreme zone sides
	OMPseg.Z[0]	 = OMPseg.Intersec_Z[0]
	OMPseg.R[-1] = OMPseg.Intersec_R[-1]
	OMPseg.Z[-1] = OMPseg.Intersec_Z[-1]

	OMPseg.ismeshed = True									# Check if at this stage the special OMP segment should be considered as meshed
	for imz in range(len(OMPseg.mzlist)):
		if(not MagMegazones[OMPseg.mzlist[imz]].ismeshed): OMPseg.ismeshed = False

	if(DEBUG > 0): print("\ndefine_OMP_segment: completed\n")

	return True