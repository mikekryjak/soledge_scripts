import numpy			as np
import itertools 		as it
from routines.globals				import DEBUG
from eirene.triangles_routines		import KnotsInterp_new

# Function definition is here
#=======================================================================================================================================

def generate_data_for_interpolation_pass1_2(Zones, Knots,  RKnots, ZKnots, Triangles, TriKnots):

	if(DEBUG > 0):	print("generate_data_for_interpolation_pass1_2")

	nKnots = len(Knots)

	KnotsInterp = KnotsInterp_new(nKnots)

#	For all nodes in soledge get mesh indexes

	nZones = len(Zones)
	NxMax  = 0
	NzMax  = 0
	for k in range(nZones):
		NxMax  = max(NxMax, Zones[k].Nx)
		NzMax  = max(NzMax, Zones[k].Nz)
		
#	Get mesh of triangles (indexes k,i,j at all knots)
	
	lnZones = "{:d}".format(int(np.log10(nZones)+1))
	lNxMax  = "{:d}".format(int(np.log10(NxMax)+1))
	lNzMax  = "{:d}".format(int(np.log10(NzMax)+1))
	MeshFormat = "{:0"+lnZones+"d},{:0"+lNxMax+"d},{:0"+lNzMax+"d}"

	for nKnot in range(nKnots):
		if((int(nKnot/10000)*10000 == nKnot) and (DEBUG > 0)):
			print("\tprocessing node n.={:6d} of nodes={:6d}".format(nKnot,nKnots))
#		
		kTri, iTri = np.where(TriKnots == nKnot)

		if(len(kTri) == 0):
			print("\tError no triangle with node = ",nKnot)
			print("\tat position ({:.4f},{:.4f}".format(RKnots[nKnot], ZKnots[nKnot]))
			continue

#		Old Method
		kijMat	   = np.array([Triangles.k[kTri], Triangles.i[kTri], Triangles.j[kTri]])
		MeshString = np.unique(np.sort(np.array(list(it.starmap(MeshFormat.format,kijMat.T)))))
		if(len(MeshString) > 0):
			SolMesh    = np.loadtxt(MeshString,dtype='i4',delimiter=',')
			if len(MeshString) == 1: SolMesh=np.reshape(SolMesh, (1,3))

			iChi = []
			for k in range(SolMesh.shape[0]):
				if(Zones[SolMesh[k, 0]].Chi[SolMesh[k, 1],SolMesh[k, 2]] == 0): iChi.append(k)	
		else:
			iChi = []
			print("\tnKnots,nKnot=",nKnots,nKnot)
			print("\tkijMat=",kijMat)
			
		
			"""
#		New Method
		OnlyOne = True
		if(len(kTri) > 1):
			Triangles_Knot = np.empty((len(kTri)), dtype=[('k','i4'), ('i','i4'), ('j', 'i4')])
			Triangles_Knot.k = Triangles.k[kTri]
			Triangles_Knot.i = Triangles.i[kTri]
			Triangles_Knot.j = Triangles.j[kTri]
			Triangle_Knot    = np.sort(Triangles_Knot, order=['k', 'i', 'j'])
			sum_diff		 = Triangles_Knot.k[1:] - Triangles_Knot.k[0:-1] + \
							   Triangles_Knot.i[1:] - Triangles_Knot.i[0:-1] + \
							   Triangles_Knot.j[1:] - Triangles_Knot.j[0:-1]
			iUnique   		 = np.where(sum_diff > 0); iUnique = iUnique[0]
			if(len(iUnique) > 0):
				iUnique = np.append(0, iUnique+1)
				SolMesh = np.array([Triangles.k[iUnique], Triangles.i[iUnique], Triangles.j[iUnique]]).T
				OnlyOne = False
		if(OnlyOne):
			SolMesh = np.array([[Triangles.k[0], Triangles.i[0], Triangles.j[0]]])
		"""
				
#		selection of the soledge quadrangle in the plasma
#		iChi	   = np.where(Zones[SolMesh[:, 0]].Chi[SolMesh[:, 1],SolMesh[:, 2]] == 0);	iChi = iChi[0]

		if(len(iChi) > 0):
			SolMesh2   = np.copy(SolMesh[iChi, :])
			nSolMesh2  = SolMesh2.shape[0]
		else:
			nSolMesh2  = 0

		if(nSolMesh2 >= 4):													#The knot is connected with more than 4 soledge quadrangle in the plasma => first pass candidate
			KnotsInterp.assp[nKnot]					= 1
			KnotsInterp.nsol[nKnot]					= nSolMesh2
			KnotsInterp.sol[nKnot, 0:nSolMesh2, :]	= np.copy(SolMesh2)

		elif(nSolMesh2 == 3):												# The knot is connected with 3 soledge quadrangle in the plasma => second pass candidate
			KnotsInterp.assp[nKnot]					= 2
			KnotsInterp.nsol[nKnot]					= nSolMesh2
			KnotsInterp.sol[nKnot, 0:nSolMesh2, :]	= np.copy(SolMesh2)
		elif(nSolMesh2 >= 1):
			KnotsInterp.assp[nKnot]					= 3
			KnotsInterp.nsol[nKnot]					= nSolMesh2
			KnotsInterp.sol[nKnot, 0:nSolMesh2,:]	= np.copy(SolMesh2)
			KnotsInterp.neir[nKnot]					= 0
		else:
			KnotsInterp.assp[nKnot]					= 4						#no soledge point (corner in penalized quadrangle typically or center node)
			KnotsInterp.nsol[nKnot]					= 0

	if(DEBUG > 0):	print("generate_data_for_interpolation_pass1_2: Completed")

	return KnotsInterp
