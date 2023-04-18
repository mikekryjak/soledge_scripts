# Function definition is here

import numpy as np
import os
import types
import h5py
from routines.h5_routines			import h5_read
from routines.globals				import DEBUG, FILE_TYPE_SOLEDGE2D, FILE_TYPE_SOLEDGE3X
from files.load_neighbors			import load_neighbors
from eirene.wall_segments			import wall_segments
from eirene.set_wall_puff			import set_wall_puff

#=========================================================
# This routine to write data in H5DF file
#=========================================================

def load_eirene_triangles(TriFile, FileType=FILE_TYPE_SOLEDGE2D):

	if(DEBUG > 0): print("load_eirene_triangles: reading file ",TriFile)
	
	Eirene = types.SimpleNamespace()

	if_triangles = h5py.File(TriFile, "r")

	Version = 2
	try:
		FileVersion	= h5_read(if_triangles, 'Version',	keep_array=False)
	except:
		FileVersion = 0

	if(FileVersion != Version):
		print("load_eirene_triangles: file version = ",FileVersion," is older than current version = ",Version)

	
	GrouName   = "knots"

	RKnots = h5_read(if_triangles,GrouName+"/R", order='F')/100.
	ZKnots = h5_read(if_triangles,GrouName+"/Z", order='F')/100.

	KnotsInterp = np.zeros((len(RKnots)), dtype=[('assp','i4'), ('nsol','i4'), ('neir', 'i4'), ('sol', 'i4', (8,3)), ('eir', 'i4', (3))])
	KnotsInterp = KnotsInterp.view(np.recarray)

	KnotsInterp.assp = h5_read(if_triangles,GrouName+"/pass", order='F')
	KnotsInterp.nsol = h5_read(if_triangles,GrouName+"/nsol", order='F')
	KnotsInterp.neir = h5_read(if_triangles,GrouName+"/neir", order='F')

	sol = h5_read(if_triangles,GrouName+"/sol",  order='F')
	eir = h5_read(if_triangles,GrouName+"/eir",  order='F')	
	KnotsInterp.sol[:, :, 1] = sol[:, :,0] - 1									#Matlab/Fortan to python indexes
	KnotsInterp.sol[:, :, 2] = sol[:, :,1] - 1									#nearest soledge nodes
	KnotsInterp.sol[:, :, 0] = sol[:, :,2] - 1
	KnotsInterp.eir			 = eir -  1											#Matlab/Fortan to python indexes


	GrouName   = "wall"
	
	nWallTriangles = h5_read(if_triangles,GrouName+"/Ntriangles", keep_array = False)
	WallTriangles  = np.empty((nWallTriangles), dtype=[('k','i4'), ('i','i4'), ('j','i4'), ('ntri','i4'), 
												('p1','i4'), ('p2','i4'), ('p3','i4'), ('side','i4'), 
												('list_e', 'i4',(3)), ('weight_e','f8'),
												('list_w', 'i4',(3)), ('weight_w','f8'),
												('list_n', 'i4',(3)), ('weight_n','f8'),
												('list_s', 'i4',(3)), ('weight_s','f8')])
	WallTriangles  = WallTriangles.view(np.recarray)

	if(FileType == FILE_TYPE_SOLEDGE2D):
		code_name = "_s2d_"
		WallTri = h5_read(if_triangles,GrouName+"/back_interp", order='F')
		WallTriangles.ntri = WallTri[:,0] - 1											#Matlab/Fortan to python indexes
		WallTriangles.k	   = WallTri[:,1] - 1											#indexes soledge mesh
		WallTriangles.i	   = WallTri[:,2] - 1
		WallTriangles.j	   = WallTri[:,3] - 1
		WallTriangles.side = h5_read(if_triangles,GrouName+"/side", order='F')
	else:
		code_name = "_cell_"
		WallTriangles.ntri = h5_read(if_triangles,GrouName+"/triNum", order='F') - 1	#python to Matlab/Fortan indexes
		WallTriangles.side = h5_read(if_triangles,GrouName+"/triFace",order='F') - 1
		nSeqs = len(Eirene.TriSequences)
		nSeqs = h5_read(if_triangles,GrouName+"/nseq", keep_array = False)
		Eirene.TriSequences = []
		for iSeq in range(nSeqs):
			SeqName  = "wall/triSequence{:d}".format(iSeq+1)
			TriSeq = h5_read(if_triangles,SeqName+"/tri", order='F') - 1
			Eirene.Wall.TypeMat[TriSeq] = h5_read(if_triangles,SeqName+"/typemat", order='F') - 1
			Eirene.Wall.IsPump[TriSeq]  = h5_read(if_triangles,SeqName+"/ispump", order='F')
			Eirene.TriSequences.append(TriSeq)

	WallTriangles.list_n   = h5_read(if_triangles,GrouName+"/north"+code_name+"to_use", order='F') - 1	#Matlab/Fortan to python indexes
	WallTriangles.weight_n = h5_read(if_triangles,GrouName+"/north"+code_name+"weight", order='F')

	WallTriangles.list_s   = h5_read(if_triangles,GrouName+"/south"+code_name+"to_use", order='F') - 1	#Matlab/Fortan to python indexes
	WallTriangles.weight_s = h5_read(if_triangles,GrouName+"/south"+code_name+"weight", order='F')

	WallTriangles.list_e   = h5_read(if_triangles,GrouName+"/east"+code_name+"to_use", order='F') - 1		#Matlab/Fortan to python indexes
	WallTriangles.weight_e = h5_read(if_triangles,GrouName+"/east"+code_name+"weight", order='F')

	WallTriangles.list_w	= h5_read(if_triangles,GrouName+"/west"+code_name+"to_use",order='F') - 1		#Matlab/Fortan to python indexes
	WallTriangles.weight_w	= h5_read(if_triangles,GrouName+"/west"+code_name+"weight",order='F')

	GrouName   = "triangles"
	nTriangles = h5_read(if_triangles,"Ntriangles", keep_array = False) 

	Triangles = np.empty((nTriangles), dtype=[('k','i4'), ('i','i4'), ('j','i4'),('p1','i4'), ('p2','i4'), ('p3','i4'), 
												('neigh1','i4'), ('typeneigh1','i4'), ('BC1','i4'),
												('neigh2','i4'), ('typeneigh2','i4'), ('BC2','i4'),
												('neigh3','i4'), ('typeneigh3','i4'), ('BC3','i4')])
	Triangles  = Triangles.view(np.recarray)
	
	TriKnots = h5_read(if_triangles,GrouName+"/tri_knots",   order='F')
	Triangles.p1 = TriKnots[:,0] - 1											#Matlab/Fortan to python indexes
	Triangles.p2 = TriKnots[:,1] - 1											#indexes triangles knots
	Triangles.p3 = TriKnots[:,2] - 1

	TypeFace = h5_read(if_triangles,GrouName+"/type_face",   order='F')
	Triangles.BC1 = TypeFace[:,0] 												#Matlab/Fortan to python indexes
	Triangles.BC2 = TypeFace[:,1] 
	Triangles.BC3 = TypeFace[:,2] 

	NeighTri  = h5_read(if_triangles,GrouName+"/neigh_tri", 	  order='F')
	Triangles.neigh1 = NeighTri[:,0] - 1										#Matlab/Fortan to python indexes
	Triangles.neigh2 = NeighTri[:,1] - 1
	Triangles.neigh3 = NeighTri[:,2] - 1
	
	NeighFace = h5_read(if_triangles,GrouName+"/neigh_face",  order='F')
	Triangles.typeneigh1 = NeighFace[:,0] - 1									#Matlab/Fortan to python indexes
	Triangles.typeneigh2 = NeighFace[:,1] - 1
	Triangles.typeneigh3 = NeighFace[:,2] - 1

	AllTri  = h5_read(if_triangles,GrouName+"/back_interp", order='F')
	if(FileType == FILE_TYPE_SOLEDGE2D):
		Triangles.k = AllTri[:,0] - 1												#Matlab/Fortan to python indexes
		Triangles.i = AllTri[:,1] - 1												#indexes soledge mesh
		Triangles.j = AllTri[:,2] - 1
	else:
		Triangles.np_k = AllTri[:,0] - 1											#Matlab/Fortan to python indexes
		Triangles.np_i = AllTri[:,1] - 1											#indexes soledge mesh
		Triangles.np_j = AllTri[:,2] - 1

	WallTriangles.p1 = Triangles.p1[WallTriangles.ntri]
	WallTriangles.p2 = Triangles.p2[WallTriangles.ntri]
	WallTriangles.p3 = Triangles.p3[WallTriangles.ntri]

	try:
		Triangles.PlasmaVacuum  = h5_read(if_triangles,GrouName+"/plasma_vacuum", order='F')
	except:
		Triangles.PlasmaVacuum = np.zeros(len(Triangles.p1), dtype='i4')

	Eirene.Eirene_OK = True

#	Define Eirene structure

	Eirene.Triangles		= Triangles
	Eirene.WallTriangles	= WallTriangles
	Eirene.KnotsInterp		= KnotsInterp
	Eirene.RKnots			= RKnots
	Eirene.ZKnots			= ZKnots
	Eirene.Wall				= types.SimpleNamespace()
	Eirene.Wall.TriSequences, Eirene.Wall.TriFaces, Eirene.Wall.R12, Eirene.Wall.Z12 = wall_segments(WallTriangles, RKnots, ZKnots)
	
#	Data stored to easily reload configuration
#	##############################################

	try:
		Eirene.ToTheCenter = h5_read(if_triangles, "ToTheCenter", 	keep_array=False)
	except:
		Eirene.ToTheCenter = 0
	
	GrouName   = "wall_config"

	Wall	   	  = Eirene.Wall
	Wall.Material = []
	Wall.Pump	  = []
	Wall.Puff	  = []
	Wall.TypeMat  = []
	Wall.IsPump   = []
	try:	
		nMaterials =  h5_read(if_triangles,GrouName+"/nMaterials", keep_array = False)
		for i in range(nMaterials):
			Wall.Material.append(types.SimpleNamespace())
			Wall.Material[-1].Name = h5_read(if_triangles,GrouName+"/Material_Name{:d}".format(i+1), keep_array = False)
			
		nPumps = h5_read(if_triangles,GrouName+"/nPumps", keep_array = False)
		for i in range(nPumps):
			Wall.Pump.append(types.SimpleNamespace())
			Wall.Pump[-1].Name = h5_read(if_triangles,GrouName+"/Pump_Name{:d}".format(i+1), keep_array = False)
		

		nPuffs = h5_read(if_triangles,GrouName+"/nPuffs", keep_array = False)
		PuffR  = h5_read(if_triangles,GrouName+"/PuffR")
		PuffZ  = h5_read(if_triangles,GrouName+"/PuffZ")

		for i in range(nPuffs):
			Wall.Puff.append(types.SimpleNamespace())
			try:	
				Wall.Puff[-1].Name = h5_read(if_triangles,GrouName+"/Puff_Name{:d}".format(i+1), keep_array = False)
			except:
				Wall.Puff[-1].Name = "{:d}".format(i+1)
			Wall.Puff[-1].R = PuffR[i]
			Wall.Puff[-1].Z = PuffZ[i]

		for i in range(len(Wall.Puff)):
			set_wall_puff(Eirene, i, Wall.Puff[i].R, Wall.Puff[i].Z)

		if(FileType == FILE_TYPE_SOLEDGE2D):
			try:
				nWallSeq = h5_read(if_triangles,GrouName+"/nWallSeq", keep_array = False)
				Eirene.Wall.TypeMat = []
				Eirene.Wall.IsPump = []
				for iSeq in range(nWallSeq):
					SeqName  = "wall_config/wall_seq_{:d}".format(iSeq+1)
					Eirene.Wall.TypeMat.append(h5_read(if_triangles,SeqName+"/TypeMat", order='F'))
					Eirene.Wall.IsPump.append(h5_read(if_triangles,SeqName+"/IsPump", order='F'))
			except:
				Eirene.Wall.TypeMat = [h5_read(if_triangles, GrouName+"/TypeMat")]
				Eirene.Wall.IsPump	= [h5_read(if_triangles, GrouName+"/IsPump")]
	except:
		Wall.Material.append(types.SimpleNamespace())
		Wall.Material[-1].Name 	= "W"

		Wall.TypeMat = []
		Wall.IsPump  = []
		for iWall in range(len(Eirene.Wall.TriSequences)):
			Wall.TypeMat.append(np.zeros((len(Eirene.Wall.TriSequences[iWall])), dtype='i4'))
			Wall.IsPump.append( np.zeros((len(Eirene.Wall.TriSequences[iWall])), dtype='bool'))

	Wall.EWalls = []
	try:
		nExtWalls	= h5_read(if_triangles, "/walls_eirene_external/nExtWalls", keep_array=False)
		for iExtWall in range(nExtWalls):
			WallName = "/walls_eirene_external/wall_{:d}".format(iExtWall+1)
			Wall.EWalls.append(types.SimpleNamespace())
			Wall.EWalls[-1].iWall	= h5_read(if_triangles, WallName+"/iWall", keep_array=False) - 1
			nESubZones				= h5_read(if_triangles, WallName+"/nESubZones", keep_array=False)
			if(nESubZones > 0):
				Wall.EWalls[-1].ESubZones = []
				for iSub in range(nESubZones):
					SubName = WallName+"/ESubZone_{:d}".format(iSub+1)
					Wall.EWalls[-1].ESubZones.append(types.SimpleNamespace())
					Wall.EWalls[-1].ESubZones[-1].dr_knots = h5_read(if_triangles, SubName+"/dr_knots", keep_array=False)
					Wall.EWalls[-1].ESubZones[-1].dz_knots = h5_read(if_triangles, SubName+"/dz_knots", keep_array=False)
					if(nESubZones > 1):
						Wall.EWalls[-1].ESubZones[-1].Rbord = h5_write(if_triangles,SubName+"/Rbord",  order='F')
						Wall.EWalls[-1].ESubZones[-1].Zbord = h5_write(if_triangles,SubName+"/Zbord",  order='F')
	except:
		pass

	Eirene.Surfaces = []
	try:	
		GrouName    = "surfaces"
		nSurfaces = h5_read(if_triangles,GrouName+"/nSurfaces", keep_array = False)
		for iSurf in range(nSurfaces):
			SurfName  = "surfaces/surface_{:d}".format(iSurf+1)
			Eirene.Surfaces.append(types.SimpleNamespace())
			Eirene.Surfaces[-1].Name	= h5_read(if_triangles,SurfName+"/Name", keep_array = False)
			Eirene.Surfaces[-1].Rcoef 	= h5_read(if_triangles,SurfName+"/Rcoef", keep_array = False)
			Eirene.Surfaces[-1].R		= h5_read(if_triangles,SurfName+"/R")
			Eirene.Surfaces[-1].Z		= h5_read(if_triangles,SurfName+"/Z")
	except:
		pass

	if_triangles.close()

	if(DEBUG > 0):	print("load_eirene_triangles: Completed")
	
	return Eirene
