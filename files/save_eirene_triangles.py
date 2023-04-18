# Function definition is here

import numpy as np
import os
import h5py
from routines.h5_routines 	import h5_write
from routines.globals		import DEBUG, FILE_TYPE_SOLEDGE2D, FILE_TYPE_SOLEDGE3X
#=========================================================
# This routine to write data in H5DF file
#=========================================================

def save_eirene_triangles(Dir, Eirene, FileType=FILE_TYPE_SOLEDGE2D):

	if(DEBUG > 0):	print("save_eirene_triangles")

	Version = 2

	try:
		os.mkdir(Dir)
	except OSError:
		try:
			os.remove(Dir + "triangles.h5")
		except OSError:
			pass

	Triangles		= Eirene.Triangles
	WallTriangles	= Eirene.WallTriangles
	KnotsInterp		= Eirene.KnotsInterp
	RKnots			= Eirene.RKnots
	ZKnots			= Eirene.ZKnots

	nTriangles = len(Triangles)
	nKnots     = len(RKnots)

	if(DEBUG > 1):	print("\tsaving to: "+Dir + "triangles.h5")

	if_triangles = h5py.File(Dir + "triangles.h5", "w")

	h5_write(if_triangles, 'Version',	    Version)

	h5_write(if_triangles, 'ToTheCenter',	Eirene.ToTheCenter)

	h5_write(if_triangles,"Ntriangles", np.array([nTriangles],dtype='i4'), order='F')
	h5_write(if_triangles,"Nknots", np.array([nKnots],dtype='i4'), order='F')

	GrouName   = "knots"
	KnotsGroup = if_triangles.create_group(GrouName)

	RKnots = 100.*RKnots																#change to cm for eirene
	ZKnots = 100.*ZKnots
	h5_write(if_triangles,GrouName+"/R", RKnots, order='F')
	h5_write(if_triangles,GrouName+"/Z", ZKnots, order='F')

	h5_write(if_triangles,GrouName+"/pass", KnotsInterp.assp, order='F')
	h5_write(if_triangles,GrouName+"/nsol", KnotsInterp.nsol, order='F')
	h5_write(if_triangles,GrouName+"/neir", KnotsInterp.neir, order='F')

	sol = np.empty(KnotsInterp.sol.shape, dtype='i4')
	sol[:, :,0] = KnotsInterp.sol[:, :, 1] + 1									#python to Matlab/Fortan indexes
	sol[:, :,1] = KnotsInterp.sol[:, :, 2] + 1									#nearest soledge nodes
	sol[:, :,2] = KnotsInterp.sol[:, :, 0] + 1
	eir			= KnotsInterp.eir + 1											#python to Matlab/Fortan indexes
	h5_write(if_triangles,GrouName+"/sol",  sol, order='F')
	h5_write(if_triangles,GrouName+"/eir",  eir, order='F')	


	GrouName   = "wall"
	WallGroup = if_triangles.create_group(GrouName)

	nWallTriangles = len(WallTriangles)
	h5_write(if_triangles,GrouName+"/Ntriangles", np.array([nWallTriangles], dtype='i4'), order='F')

	if(FileType == FILE_TYPE_SOLEDGE2D):
		code_name = "_s2d_"
		WallTri	     = np.empty((nWallTriangles,4),'i4')
		WallTri[:,0] = WallTriangles.ntri+1													#python to Matlab/Fortan indexes
		WallTri[:,1] = WallTriangles.k+1													#indexes soledge mesh (NO USED IN CODE)
		WallTri[:,2] = WallTriangles.i+1													# (NO USED IN CODE)
		WallTri[:,3] = WallTriangles.j+1													# (NO USED IN CODE)
		h5_write(if_triangles,GrouName+"/back_interp", WallTri, order='F')
		h5_write(if_triangles,GrouName+"/side", WallTriangles.side, order='F')			#NO MORE USED

	else:
		code_name = "_cell_"
		h5_write(if_triangles,GrouName+"/triNum", WallTriangles.ntri+1, order='F')		#python to Matlab/Fortan indexes
		h5_write(if_triangles,GrouName+"/triFace", WallTriangles.side+1, order='F')
		nSeqs = len(Eirene.TriSequences)
		h5_write(if_triangles,GrouName+"/nseq", np.array([nSeqs], dtype='i4'), order='F')
		for iSeq in range(nSeqs):
			TriSeq = Eirene.TriSequences[iSeq]
			SeqName  = "wall/triSequence{:d}".format(iSeq+1)
			SeqGroup = if_triangles.create_group(SeqName)
			h5_write(if_triangles,SeqName+"/tri", TriSeq+1, order='F')
			h5_write(if_triangles,SeqName+"/typemat", Eirene.Wall.TypeMat[TriSeq]+1, order='F')
			h5_write(if_triangles,SeqName+"/ispump", Eirene.Wall.IsPump[TriSeq], order='F')

	h5_write(if_triangles,GrouName+"/north"+code_name+"to_use", WallTriangles.list_n + 1, order='F')	#python to Matlab/Fortan indexes
	h5_write(if_triangles,GrouName+"/north"+code_name+"weight", WallTriangles.weight_n,   order='F')

	h5_write(if_triangles,GrouName+"/south"+code_name+"to_use", WallTriangles.list_s + 1, order='F')	#python to Matlab/Fortan indexes
	h5_write(if_triangles,GrouName+"/south"+code_name+"weight", WallTriangles.weight_s,   order='F')

	h5_write(if_triangles,GrouName+"/east"+code_name+"to_use", WallTriangles.list_e + 1, order='F')	#python to Matlab/Fortan indexes
	h5_write(if_triangles,GrouName+"/east"+code_name+"weight", WallTriangles.weight_e,   order='F')

	h5_write(if_triangles,GrouName+"/west"+code_name+"to_use", WallTriangles.list_w + 1, order='F')	#python to Matlab/Fortan indexes
	h5_write(if_triangles,GrouName+"/west"+code_name+"weight", WallTriangles.weight_w,   order='F')

	GrouName   = "triangles"
	TrianglesGroup = if_triangles.create_group(GrouName)

	TriKnots = np.array([Triangles.p1+1, Triangles.p2+1, Triangles.p3+1], dtype='i4').T											#python to Matlab/Fortan indexes
	h5_write(if_triangles,GrouName+"/tri_knots",   TriKnots, order='F')

	TypeFace = np.array([Triangles.BC1, Triangles.BC2, Triangles.BC3], dtype='i4').T
	h5_write(if_triangles,GrouName+"/type_face", TypeFace, order='F')		#Type of boundary condiction (core, wall, puff, pump)

	NeighTri = np.array([Triangles.neigh1+1, Triangles.neigh2+1, Triangles.neigh3+1], dtype='i4').T								#python to Matlab/Fortan indexes
	h5_write(if_triangles,GrouName+"/neigh_tri", 	 NeighTri, order='F')	#number of neighbour triangle on sides

	NeighFace = np.array([Triangles.typeneigh1+1, Triangles.typeneigh2+1, Triangles.typeneigh3+1], dtype='i4').T					#python to Matlab/Fortan indexes
	h5_write(if_triangles,GrouName+"/neigh_face",  NeighFace, order='F')	#number of side on neighbour triangles


	if(FileType == FILE_TYPE_SOLEDGE2D): AllTri = np.array([Triangles.k+1,Triangles.i+1,Triangles.j+1], dtype='i4').T				#python to Matlab/Fortan indexes
	else:								 AllTri = np.array([Triangles.np_k+1,Triangles.np_i+1,Triangles.np_j+1], dtype='i4').T		#python to Matlab/Fortan indexes
	h5_write(if_triangles,GrouName+"/back_interp", AllTri, order='F')

	h5_write(if_triangles,GrouName+"/plasma_vacuum", Triangles.PlasmaVacuum, order='F')

#	Data stored to easily reload configuration
#	########################

	GrouName   = "wall_config"
	ConfigGroup = if_triangles.create_group(GrouName)

	nMaterials = len(Eirene.Wall.Material)
	h5_write(if_triangles,GrouName+"/nMaterials", nMaterials)
	for i in range(nMaterials):
		h5_write(if_triangles,GrouName+"/Material_Name{:d}".format(i+1), Eirene.Wall.Material[i].Name)
		
	nPumps = len(Eirene.Wall.Pump)
	h5_write(if_triangles,GrouName+"/nPumps", nPumps)
	for i in range(nPumps):
		h5_write(if_triangles,GrouName+"/Pump_Name{:d}".format(i+1), Eirene.Wall.Pump[i].Name)
		
	nPuffs = len(Eirene.Wall.Puff)
	h5_write(if_triangles,GrouName+"/nPuffs", nPuffs)
	if(nPuffs > 0):
		for i in range(nPuffs):
			h5_write(if_triangles,GrouName+"/Puff_Name{:d}".format(i+1), Eirene.Wall.Puff[i].Name)

		PuffR = np.array([Eirene.Wall.Puff[k].R for k in range(nPuffs)])
		PuffZ = np.array([Eirene.Wall.Puff[k].Z for k in range(nPuffs)])
	
		h5_write(if_triangles,GrouName+"/PuffR", PuffR)
		h5_write(if_triangles,GrouName+"/PuffZ", PuffZ)

	if(FileType == FILE_TYPE_SOLEDGE2D):
		nWallSeq = len(Eirene.Wall.TypeMat)
		h5_write(if_triangles,GrouName+"/nWallSeq", np.array([nWallSeq], dtype='i4'), order='F')
		for iSeq in range(nWallSeq):
			SeqName  = "wall_config/wall_seq_{:d}".format(iSeq+1)
			SeqGroup = if_triangles.create_group(SeqName)
			h5_write(if_triangles,SeqName+"/TypeMat", Eirene.Wall.TypeMat[iSeq], order='F')
			h5_write(if_triangles,SeqName+"/IsPump", Eirene.Wall.IsPump[iSeq], order='F')

#		write Eirene parameters for wall 

	if(len(Eirene.Wall.EWalls) > 0):
		WallsExtName = "walls_eirene_external"
		WallExtGroup = if_triangles.create_group(WallsExtName)
		nExtWalls = len(Eirene.Wall.EWalls)
		h5_write(if_triangles,WallsExtName+"/nExtWalls", nExtWalls)
		EWalls = Eirene.Wall.EWalls
		for iExtWall in range(nExtWalls):
			WallName = WallsExtName + "/wall_{:d}".format(iExtWall+1)
			if_triangles.create_group(WallName)

			EWall = EWalls[iExtWall]
			h5_write(if_triangles, WallName+"/iWall",	EWall.iWall+1)							#Index of wall in Config

			ESubZones  = EWall.ESubZones
			nESubZones = len(ESubZones)
			h5_write(if_triangles,WallName+"/nESubZones", nESubZones)
			for iSub in range(nESubZones):
				SubName = WallName+"/ESubZone_{:d}".format(iSub+1)
				SubGroup = if_triangles.create_group(SubName)
				h5_write(if_triangles,SubName+"/dr_knots", ESubZones[iSub].dr_knots)
				h5_write(if_triangles,SubName+"/dz_knots", ESubZones[iSub].dz_knots)
				if(len(ESubZones) > 1):
					h5_write(if_triangles,SubName+"/Rbord", ESubZones[iSub].Rbord, order='F')
					h5_write(if_triangles,SubName+"/Zbord", ESubZones[iSub].Zbord, order='F')

	if(len(Eirene.Surfaces) > 0):
		GrouName    = "surfaces"
		ConfigGroup = if_triangles.create_group(GrouName)
		nSurfaces = len(Eirene.Surfaces)
		h5_write(if_triangles,GrouName+"/nSurfaces", np.array([nSurfaces], dtype='i4'), order='F')
		for iSurf in range(nSurfaces):
			SurfName  = "surfaces/surface_{:d}".format(iSurf+1)
			SurfGroup = if_triangles.create_group(SurfName)
			h5_write(if_triangles,SurfName+"/Name",  Eirene.Surfaces[iSurf].Name)
			h5_write(if_triangles,SurfName+"/Rcoef", Eirene.Surfaces[iSurf].Rcoef)
			h5_write(if_triangles,SurfName+"/R",      Eirene.Surfaces[iSurf].R)
			h5_write(if_triangles,SurfName+"/Z",      Eirene.Surfaces[iSurf].Z)

	if_triangles.close()

	if(DEBUG > 0):	print("save_eirene_triangles: Completed")
