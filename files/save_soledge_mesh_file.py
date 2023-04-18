#!/usr/bin/python

# Function definition is here

import numpy as np
import os
import h5py
from routines.globals	  import EW_BORDER_NEIGHBOUR, CUSTOM_PLOTS_TYPES, FILE_TYPE_SOLEDGE2D, FILE_TYPE_SOLEDGE3X
from routines.h5_routines import h5_write

#=========================================================
# This routine to write data to the H5DF mesh file
#=========================================================

def save_soledge_mesh_file(MeshFile, Config, FileType=FILE_TYPE_SOLEDGE2D): 

	if(FileType == FILE_TYPE_SOLEDGE2D):	Version = 3
	else:									Version = 1000
	
	print("save_soledge_mesh_file: Saving to ",MeshFile)
	
	if_mesh = h5py.File(MeshFile, "w")

	MagZones	  = Config.MagZones
	MagMegazones  = Config.MagMegazones

	h5_write(if_mesh, 'Version',	Version)
	
	mesher_fields = "/mesher/Fields"
	if(Config.in_equ_OK):
		if_mesh.create_group("mesher")
	
		if_mesh.create_group(mesher_fields)
		
		h5_write(if_mesh, mesher_fields+ '/r2',				Config.r2,			order='F')
		h5_write(if_mesh, mesher_fields+ '/z2',				Config.z2,			order='F')
		h5_write(if_mesh, mesher_fields+ '/Br2',				Config.Br2,			order='F')
		h5_write(if_mesh, mesher_fields+ '/Bz2',				Config.Bz2, 		order='F')
		h5_write(if_mesh, mesher_fields+ '/Bphi2',			Config.Bphi2,		order='F')
		h5_write(if_mesh, mesher_fields+ '/psi2',				Config.flux2,		order='F')
		h5_write(if_mesh, mesher_fields+ '/psi2_x_psi',		Config.flux2_x_psi,	order='F')

		h5_write(if_mesh, mesher_fields+ '/OutLength2',		Config.OutLength2)
		h5_write(if_mesh, mesher_fields+ '/Smooth2',			Config.Smooth2)
		h5_write(if_mesh, mesher_fields+ '/Smooth2D',		Config.Smooth2D)
		h5_write(if_mesh, mesher_fields+ '/Nr',				Config.nr)
		h5_write(if_mesh, mesher_fields+ '/Nz',				Config.nz)

	if(Config.extrapol_OK):

		h5_write(if_mesh, mesher_fields+ '/ExtrapolValue',	Config.extrapol_val)
		h5_write(if_mesh, mesher_fields+ '/Smooth2DExt',		Config.Smooth2DExt)
		h5_write(if_mesh, mesher_fields+ '/RaiseValue',		Config.raise_value)
		h5_write(if_mesh, mesher_fields+ '/RaisePower',		Config.raise_power)
		
		mesher_x_points = "/mesher/XPoints"
		if_mesh.create_group(mesher_x_points)

		h5_write(if_mesh, mesher_x_points+ '/RadArroundXp',	Config.RadArroundXp)
		h5_write(if_mesh, mesher_x_points+ '/DescentStep',		Config.DescentStep)

		nX_points = len(Config.X_points)
		h5_write(if_mesh, mesher_x_points +"/nXPoints", nX_points)
		for k in range(nX_points):
			x_point = mesher_x_points + "/XPoint{:d}".format(k+1)
			if_mesh.create_group(x_point)
			h5_write(if_mesh,x_point + "/R",		Config.X_points[k].R)										#From find_Xpoints (in post_loaded_fields)
			h5_write(if_mesh,x_point + "/Z",		Config.X_points[k].Z)
			h5_write(if_mesh,x_point + "/psi",		Config.X_points[k].psi)
			h5_write(if_mesh,x_point + "/index",	Config.X_points[k].index)
			
			for n in range(4):
				cut = x_point+"/cut{:d}".format(n+1)
				if_mesh.create_group(cut)
				h5_write(if_mesh,cut+"/Rs", 	Config.X_points[k].cut[n].Rs)									#From surround_Xpoints (in post_loaded_fields)
				h5_write(if_mesh,cut+"/Zs", 	Config.X_points[k].cut[n].Zs)
				h5_write(if_mesh,cut+"/psi", 	Config.X_points[k].cut[n].psi)
				
			for n in range(4):
				branch = x_point+"/branch{:d}".format(n+1)
				if_mesh.create_group(branch)
				h5_write(if_mesh,branch+"/R", 		Config.X_points[k].branch[n].R)
				h5_write(if_mesh,branch+"/Z", 		Config.X_points[k].branch[n].Z)
				h5_write(if_mesh,branch+"/theta", 	Config.X_points[k].branch[n].theta)
				h5_write(if_mesh,branch+"/arc", 	Config.X_points[k].branch[n].arc)

#	X points have been found
#-----------------------------

	if(Config.xPoints_OK):
		h5_write(if_mesh, mesher_fields+ '/PsiCore',			Config.psicore)
		h5_write(if_mesh, mesher_fields+ '/PsiEdge',			Config.psiout)

		mesher_x_points = "/mesher/XPoints"

		for k in range(nX_points):
			x_point = mesher_x_points + "/XPoint{:d}".format(k+1)
			for n in range(4):
				cut = x_point+"/cut{:d}".format(n+1)
				h5_write(if_mesh,cut+"/R", 		Config.X_points[k].cut[n].R,		order='F')					#From plot_ortho_lines (in cmd_split_zones)
				h5_write(if_mesh,cut+"/Z", 		Config.X_points[k].cut[n].Z,		order='F')
				h5_write(if_mesh,cut+"/type", 	Config.X_points[k].cut[n].type)
				h5_write(if_mesh,cut+"/psiR", 	Config.X_points[k].cut[n].psiR,		order='F')					#From getArcLim (in cmd_split_zones)
				h5_write(if_mesh,cut+"/psiZ", 	Config.X_points[k].cut[n].psiZ,		order='F')
				h5_write(if_mesh,cut+"/psilim", Config.X_points[k].cut[n].psilim)

#	frontiers have been defined
#-----------------------------

	if(Config.Frontiers_OK):
		mesher_frontiers = "/mesher/Frontiers"
		if_mesh.create_group(mesher_frontiers)
		nFrontiers =len(Config.Frontiers)
		h5_write(if_mesh,mesher_frontiers +'/nFrontiers',	nFrontiers)
		for k in range(nFrontiers):
			frontier = mesher_frontiers + "/Frontier{:d}".format(k+1)
			if_mesh.create_group(frontier)
			h5_write(if_mesh, frontier+"/R",	Config.Frontiers[k].R,	order='F')								#From define_frontier
			h5_write(if_mesh, frontier+"/Z",	Config.Frontiers[k].Z,	order='F')
			h5_write(if_mesh, frontier+"/P1",	Config.Frontiers[k].P1,	order='F')
			h5_write(if_mesh, frontier+"/psiA",	Config.Frontiers[k].psiA)
			h5_write(if_mesh, frontier+"/psiB",	Config.Frontiers[k].psiB)

#	domains have been defined
#-----------------------------
	if(Config.MagZones_OK):
		mesher_mag_zones = "/mesher/MagZones"
		if_mesh.create_group(mesher_mag_zones)
		nMagZones = len(MagZones)
		h5_write(if_mesh,mesher_mag_zones +'/nMagZones', nMagZones)
		for k in range(nMagZones):
			mag_zone = mesher_mag_zones + "/MagZone{:d}".format(k+1)
			if_mesh.create_group(mag_zone)
			h5_write(if_mesh, mag_zone+"/east_R",			MagZones[k].east.R,		order='F')
			h5_write(if_mesh, mag_zone+"/east_Z",			MagZones[k].east.Z,		order='F')
			h5_write(if_mesh, mag_zone+"/west_R",			MagZones[k].west.R,		order='F')
			h5_write(if_mesh, mag_zone+"/west_Z",			MagZones[k].west.Z,		order='F')
			h5_write(if_mesh, mag_zone+"/north_R",		MagZones[k].north.R,	order='F')
			h5_write(if_mesh, mag_zone+"/north_Z",		MagZones[k].north.Z,	order='F')
			h5_write(if_mesh, mag_zone+"/south_R",		MagZones[k].south.R,	order='F')
			h5_write(if_mesh, mag_zone+"/south_Z",		MagZones[k].south.Z,	order='F')
			h5_write(if_mesh, mag_zone+"/coord",			MagZones[k].coord,		order='F')
			h5_write(if_mesh, mag_zone+"/pA_coord",		MagZones[k].pA.coord,	order='F')
			h5_write(if_mesh, mag_zone+"/pB_coord",		MagZones[k].pB.coord,	order='F')
			h5_write(if_mesh, mag_zone+"/pC_coord",		MagZones[k].pC.coord,	order='F')
			h5_write(if_mesh, mag_zone+"/pD_coord",		MagZones[k].pD.coord,	order='F')
			
			h5_write(if_mesh, mag_zone+"/mz",				MagZones[k].mz)
			h5_write(if_mesh, mag_zone+"/pmz",				MagZones[k].pmz)
			h5_write(if_mesh, mag_zone+"/Xtype_east",		MagZones[k].Xtype_east)
			h5_write(if_mesh, mag_zone+"/Xtype_west",		MagZones[k].Xtype_west)	
			
			h5_write(if_mesh, mag_zone+"/Neighbour_east",	MagZones[k].Neighbour.east)
			h5_write(if_mesh, mag_zone+"/Neighbour_west",	MagZones[k].Neighbour.west)
			h5_write(if_mesh, mag_zone+"/Neighbour_north",	MagZones[k].Neighbour.north)
			h5_write(if_mesh, mag_zone+"/Neighbour_south",	MagZones[k].Neighbour.south)
			
		mesher_pmegazones = "/mesher/MagPMegazones"
		if_mesh.create_group(mesher_pmegazones)
		MagPMegazones = Config.MagPMegazones
		nMagPMegazones =len(MagPMegazones)
		h5_write(if_mesh,mesher_pmegazones +'/nMagPMegazones', nMagPMegazones)
		for k in range(nMagPMegazones):
			pmegazone = mesher_pmegazones + "/MagPMegazone{:d}".format(k+1)
			if_mesh.create_group(pmegazone)
			h5_write(if_mesh, pmegazone+"/list",			MagPMegazones[k].list,	order='F')

		mesher_megazones = "/mesher/MagMegazones"
		if_mesh.create_group(mesher_megazones)
		
		nMagMegazones = len(MagMegazones)
		h5_write(if_mesh,mesher_megazones +'/nMagMegazones', nMagMegazones)
		for k in range(nMagMegazones):
			megazone = mesher_megazones + "/MagMegazone{:d}".format(k+1)
			if_mesh.create_group(megazone)
			h5_write(if_mesh, megazone+"/list",			MagMegazones[k].list,	order='F')
			h5_write(if_mesh, megazone+"/isperiodic",		MagMegazones[k].isperiodic)

#	Segments have been defined
#-----------------------------

	if(Config.Segments_OK):	
		mesher_pmegazones = "/mesher/MagPMegazones"
		nMagPMegazones =len(MagPMegazones)
		for k in range(nMagPMegazones):
			pmegazone = mesher_pmegazones + "/MagPMegazone{:d}".format(k+1)
			h5_write(if_mesh, pmegazone+"/isaligned", MagPMegazones[k].isaligned)
			
			pmegazone_refpoints = pmegazone + "/refpoints"
			if_mesh.create_group(pmegazone_refpoints)
			if(MagPMegazones[k].isaligned):
				h5_write(if_mesh, pmegazone+"/iAlignWall", MagPMegazones[k].iAlignWall)

				h5_write(if_mesh, pmegazone_refpoints+"/R_sub1",			MagPMegazones[k].subrefpoints[0].R,	order='F')
				h5_write(if_mesh, pmegazone_refpoints+"/Z_sub1",			MagPMegazones[k].subrefpoints[0].Z,	order='F')
				h5_write(if_mesh, pmegazone_refpoints+"/R_sub2",			MagPMegazones[k].subrefpoints[1].R,	order='F')
				h5_write(if_mesh, pmegazone_refpoints+"/Z_sub2",			MagPMegazones[k].subrefpoints[1].Z,	order='F')
				h5_write(if_mesh, pmegazone_refpoints+"/align_psimin",	MagPMegazones[k].align_psimin)
				h5_write(if_mesh, pmegazone_refpoints+"/align_psimax",	MagPMegazones[k].align_psimax)
				
			h5_write(if_mesh, pmegazone_refpoints+"/R",		MagPMegazones[k].refpoints.R,	order='F')
			h5_write(if_mesh, pmegazone_refpoints+"/Z",		MagPMegazones[k].refpoints.Z,	order='F')
			h5_write(if_mesh, pmegazone_refpoints+"/nz",		MagPMegazones[k].refpoints.nz)
			h5_write(if_mesh, pmegazone_refpoints+"/nzB",		MagPMegazones[k].refpoints.nzB)

			if(MagPMegazones[k].isaligned):
				for n in range(2):
					sub = "_sub{:d}".format(n+1)
					h5_write(if_mesh, pmegazone_refpoints+"/nPoints"+sub,		MagPMegazones[k].subrefpoints[n].nPoints)							#Config parameters
					h5_write(if_mesh, pmegazone_refpoints+"/RefineType"+sub,	MagPMegazones[k].subrefpoints[n].RefineType)
					h5_write(if_mesh, pmegazone_refpoints+"/RefineSide"+sub,	MagPMegazones[k].subrefpoints[n].RefineSide)
					h5_write(if_mesh, pmegazone_refpoints+"/AdjustMode"+sub,	MagPMegazones[k].subrefpoints[n].AdjustMode)
					h5_write(if_mesh, pmegazone_refpoints+"/ParamL"+sub,		MagPMegazones[k].subrefpoints[n].ParamL)
					h5_write(if_mesh, pmegazone_refpoints+"/ParamR"+sub,		MagPMegazones[k].subrefpoints[n].ParamR)
					if(MagPMegazones[k].subrefpoints[n].RefineType == 2):
						h5_write(if_mesh, pmegazone_refpoints+"/xBezier"+sub,	MagPMegazones[k].subrefpoints[n].xBezier,	order='F')
						h5_write(if_mesh, pmegazone_refpoints+"/yBezier"+sub,	MagPMegazones[k].subrefpoints[n].yBezier,	order='F')
			else:
				h5_write(if_mesh, pmegazone_refpoints+"/nPoints",		MagPMegazones[k].refpoints.nPoints)							#Config parameters
				h5_write(if_mesh, pmegazone_refpoints+"/RefineType",	MagPMegazones[k].refpoints.RefineType)
				h5_write(if_mesh, pmegazone_refpoints+"/RefineSide",	MagPMegazones[k].refpoints.RefineSide)
				h5_write(if_mesh, pmegazone_refpoints+"/AdjustMode",	MagPMegazones[k].refpoints.AdjustMode)
				h5_write(if_mesh, pmegazone_refpoints+"/ParamL",		MagPMegazones[k].refpoints.ParamL)
				h5_write(if_mesh, pmegazone_refpoints+"/ParamR",		MagPMegazones[k].refpoints.ParamR)
				if(MagPMegazones[k].refpoints.RefineType == 2):
					h5_write(if_mesh, pmegazone_refpoints+"/xBezier",	MagPMegazones[k].refpoints.xBezier,	order='F')
					h5_write(if_mesh, pmegazone_refpoints+"/yBezier",	MagPMegazones[k].refpoints.yBezier,	order='F')

			h5_write(if_mesh, pmegazone_refpoints+"/ForceOrtho",	MagPMegazones[k].ForceOrtho)
		
		mesher_megazones = "/mesher/MagMegazones"
		MagMegazones  = Config.MagMegazones
		nMagMegazones =len(MagMegazones)
		for k in range(nMagMegazones):
			megazone_refpoints = mesher_megazones + "/MagMegazone{:d}/refpoints".format(k+1)
			if_mesh.create_group(megazone_refpoints)
			h5_write(if_mesh, megazone_refpoints+"/R",		MagMegazones[k].refpoints.R,	order='F')
			h5_write(if_mesh, megazone_refpoints+"/Z",		MagMegazones[k].refpoints.Z,	order='F')
			h5_write(if_mesh, megazone_refpoints+"/psi",	MagMegazones[k].refpoints.psi,	order='F')

			h5_write(if_mesh, megazone_refpoints+"/nPoints",			MagMegazones[k].refpoints.nPoints)							#Config parameters
			h5_write(if_mesh, megazone_refpoints+"/RefineType",		MagMegazones[k].refpoints.RefineType)
			h5_write(if_mesh, megazone_refpoints+"/RefineSide",		MagMegazones[k].refpoints.RefineSide)
			h5_write(if_mesh, megazone_refpoints+"/AdjustMode",		MagMegazones[k].refpoints.AdjustMode)
			h5_write(if_mesh, megazone_refpoints+"/ParamL",			MagMegazones[k].refpoints.ParamL)
			h5_write(if_mesh, megazone_refpoints+"/ParamR",			MagMegazones[k].refpoints.ParamR)
			if(MagMegazones[k].refpoints.RefineType == 2):
				h5_write(if_mesh, megazone_refpoints+"/xBezier",	MagMegazones[k].refpoints.xBezier,	order='F')
				h5_write(if_mesh, megazone_refpoints+"/yBezier",	MagMegazones[k].refpoints.yBezier,	order='F')

#		write Outer Mid Plane segment

		OMPseg  = Config.OMPseg
		mesher_OMP_segment = "/mesher/OMP_segment"
		if_mesh.create_group(mesher_OMP_segment)
		h5_write(if_mesh, mesher_OMP_segment+"/zonelist",	OMPseg.zonelist,	order='F')
		h5_write(if_mesh, mesher_OMP_segment+"/mzlist",		OMPseg.mzlist,		order='F')
		h5_write(if_mesh, mesher_OMP_segment+"/Intersec_R",	OMPseg.Intersec_R,	order='F')
		h5_write(if_mesh, mesher_OMP_segment+"/Intersec_Z",	OMPseg.Intersec_Z,	order='F')
		h5_write(if_mesh, mesher_OMP_segment+"/R",				OMPseg.R,	order='F')
		h5_write(if_mesh, mesher_OMP_segment+"/Z",				OMPseg.Z,	order='F')

		if(hasattr(OMPseg, "refpoints")):
			OMP_segment_refpoints = mesher_OMP_segment+"/refpoints"
			if_mesh.create_group(OMP_segment_refpoints)
			h5_write(if_mesh, OMP_segment_refpoints+"/nPoints", 		OMPseg.refpoints.nPoints)							#Mesher parameters
			h5_write(if_mesh, OMP_segment_refpoints+"/RefineType",	OMPseg.refpoints.RefineType)
			h5_write(if_mesh, OMP_segment_refpoints+"/RefineSide",	OMPseg.refpoints.RefineSide)
			h5_write(if_mesh, OMP_segment_refpoints+"/AdjustMode",	OMPseg.refpoints.AdjustMode)
			h5_write(if_mesh, OMP_segment_refpoints+"/ParamL", 		OMPseg.refpoints.ParamL)
			h5_write(if_mesh, OMP_segment_refpoints+"/ParamR", 		OMPseg.refpoints.ParamR)
			if(OMPseg.refpoints.RefineType == 2):
				h5_write(if_mesh, OMP_segment_refpoints+"/xBezier",	OMPseg.refpoints.xBezier, order='F')
				h5_write(if_mesh, OMP_segment_refpoints+"/yBezier",	OMPseg.refpoints.yBezier, order='F')

		mesher_mag_zones = "/mesher/MagZones"
		for k in range(nMagZones):
			mag_zone = mesher_mag_zones + "/MagZone{:d}".format(k+1)
			h5_write(if_mesh, mag_zone+"/northaligned", MagZones[k].northaligned)
			h5_write(if_mesh, mag_zone+"/southaligned", MagZones[k].southaligned)
			if(MagZones[k].northaligned or MagZones[k].southaligned): h5_write(if_mesh, mag_zone+"/iAlignWall", MagZones[k].iAlignWall)
			
			if(MagZones[k].northaligned):
				for n in range(2):
					sub = "_sub{:d}".format(n+1)
					h5_write(if_mesh, mag_zone+"/North_R"+sub,	MagZones[k].subNorth[n].R,	order='F')
					h5_write(if_mesh, mag_zone+"/North_Z"+sub,	MagZones[k].subNorth[n].Z,	order='F')

			if(MagZones[k].southaligned):
				for n in range(2):
					sub = "_sub{:d}".format(n+1)
					h5_write(if_mesh, mag_zone+"/South_R"+sub,	MagZones[k].subSouth[n].R,	order='F')
					h5_write(if_mesh, mag_zone+"/South_Z"+sub,	MagZones[k].subSouth[n].Z,	order='F')

	if(Config.MagGrid_OK):
		for k in range(nMagZones):
			mag_zone = mesher_mag_zones + "/MagZone{:d}".format(k+1)
			h5_write(if_mesh,	mag_zone+"/meshortho",	MagZones[k].meshortho)
			h5_write(if_mesh,	mag_zone+"/Rcorner",		MagZones[k].gridR, order='F')				#Soledge mesh corners
			h5_write(if_mesh,	mag_zone+"/Zcorner",		MagZones[k].gridZ, order='F')

		Optimization = Config.Optimization
		mesher_optim = "/mesher/Optimization"
		if_mesh.create_group(mesher_optim)
		h5_write(if_mesh, mesher_optim+"/Ortho",		Optimization.orthoptim)
		h5_write(if_mesh, mesher_optim+"/dMin",		Optimization.dmin_ad)
		h5_write(if_mesh, mesher_optim+"/Range", 		Optimization.range_ad)
		h5_write(if_mesh, mesher_optim+"/Length", 		Optimization.length_ad)
		h5_write(if_mesh, mesher_optim+"/Alpha", 		Optimization.alpha_ad)
		h5_write(if_mesh, mesher_optim+"/Target", 		Optimization.target_ad)
		h5_write(if_mesh, mesher_optim+"/drMin", 		Optimization.drmin_ad)

	if(Config.wall_OK):
		if_mesh.create_group("walls")
		nWalls = int(len(Config.Walls))
		h5_write(if_mesh, "/walls/Nwalls",	 nWalls)							#Number of walls
		for k in range(nWalls):
			wall = "walls/wall{:d}".format(k+1)
			if_mesh.create_group(wall)
			h5_write(if_mesh, wall+ "/R", Config.Walls[k].Rwall, order='F')
			h5_write(if_mesh, wall+ "/Z", Config.Walls[k].Zwall, order='F')
			h5_write(if_mesh, wall+ "/type", Config.Walls[k].Type)
			h5_write(if_mesh, wall+ "/line_type", Config.Walls[k].LineType, order='F')

	if(Config.equ_OK):
		if_mesh.create_group("config")
	
		h5_write(if_mesh,'/config/r',		Config.r2D,		order='F')
		h5_write(if_mesh,'/config/z',		Config.z2D,		order='F')
		h5_write(if_mesh,'/config/psi',	Config.flux2D,	order='F')
		if(Config.new_equ_OK):
			h5_write(if_mesh,'/config/Br',	Config.Br2D,	order='F')
			h5_write(if_mesh,'/config/Bz',	Config.Bz2D, 	order='F')
			h5_write(if_mesh,'/config/Bphi',	Config.Bphi2D,	order='F')

	if(len(Config.X_points) > 0):
		nX_points = len(Config.X_points)
		h5_write(if_mesh,'/config/nsep', nX_points)
	
		for k in range(len(Config.X_points)):
			h5_write(if_mesh,"/config/psisep{:d}".format(k+1),	Config.X_points[k].psi)
	
		if(hasattr(Config.X_points[0],'psicore')):
			h5_write(if_mesh,'/config/psicore',	Config.X_points[0].psicore)
		elif(hasattr(Config,'psicore')):
			h5_write(if_mesh,'/config/psicore',	Config.psicore)

	if(Config.Mesh_OK):
		if(FileType == FILE_TYPE_SOLEDGE2D):	write_mesh_S2D(if_mesh, Config)
		else:									write_mesh_S3X(if_mesh, Config)

#		write Zones data to recover splitting configuration

		mesher_zones = "/mesher/Zones"
		if_mesh.create_group(mesher_zones)

		Zones  = Config.Zones
		nZones = len(Zones)
		for k in range (nZones):
			zone = mesher_zones + "/Zone{:d}".format(k+1)
			if_mesh.create_group(zone)
			h5_write(if_mesh,zone+"/MagZone",	Zones[k].magz + 1, order='F')	#Python index

		if(len(Config.Split) > 0): h5_write(if_mesh,"mesher/Split",	Config.Split, order='F')

#	Transport profiles have been defined

	if(Config.transp_prof_OK):
		if_mesh.create_group("/mesher/TranspCuts")
		transp_cuts= "/mesher/TranspCuts"

		TranspCuts		= Config.TranspCuts
		nTranspCuts		= len(TranspCuts)
		nTranspProfiles = len(TranspCuts[0].Flux.Profiles)
		h5_write(if_mesh, transp_cuts+"/nCuts", nTranspCuts)
		h5_write(if_mesh, transp_cuts+"/nProfiles", nTranspProfiles)
		
		for i in range(nTranspProfiles):
			h5_write(if_mesh, transp_cuts+"/Prof_Name{:d}".format(i+1),  TranspCuts[0].Flux.Profiles[i].Name,  order='F')
			
		for i in range(nTranspCuts):
			cut = transp_cuts+"/cut{:d}".format(i+1)
			Cgroup = if_mesh.create_group(cut)

			h5_write(if_mesh, cut+"/BallooningMode",  TranspCuts[i].BallooningMode)
			h5_write(if_mesh, cut+"/ZoneMode",  		  TranspCuts[i].ZoneMode)
			
			h5_write(if_mesh, cut+"/Flux_nz",  TranspCuts[i].Flux.nz,  order='F')
			h5_write(if_mesh, cut+"/Flux_dir", TranspCuts[i].Flux.dir, order='F')
			h5_write(if_mesh, cut+"/Flux_d12", TranspCuts[i].Flux.d12, order='F')
			h5_write(if_mesh, cut+"/Flux_r12", TranspCuts[i].Flux.r12, order='F')
			h5_write(if_mesh, cut+"/Flux_z12", TranspCuts[i].Flux.z12, order='F')
			for k in range(len(TranspCuts[i].Flux.nz)):
				h5_write(if_mesh, cut+"/Flux_Zones{:d}".format(k+1), TranspCuts[i].Flux.MagZones[k], order='F')
						
			for k in range(len(TranspCuts[i].Flux.Profiles)):
				h5_write(if_mesh, cut+"/Flux_"+ TranspCuts[i].Flux.Profiles[k].Name+"_xValues", TranspCuts[i].Flux.Profiles[k].xValues, order='F')
				h5_write(if_mesh, cut+"/Flux_"+ TranspCuts[i].Flux.Profiles[k].Name+"_Values",	TranspCuts[i].Flux.Profiles[k].Values,	order='F')
	
			h5_write(if_mesh, cut+"/Theta_nz", TranspCuts[i].Theta.nz, 	order='F')
			h5_write(if_mesh, cut+"/Theta_d12", TranspCuts[i].Theta.d12,	order='F')	
						
			for k in range(len(TranspCuts[i].Theta.Profiles)):
				h5_write(if_mesh, cut+"/Theta_"+ TranspCuts[i].Theta.Profiles[k].Name+"_xValues", TranspCuts[i].Theta.Profiles[k].xValues,	order='F')
				h5_write(if_mesh, cut+"/Theta_"+ TranspCuts[i].Theta.Profiles[k].Name+"_Values", TranspCuts[i].Theta.Profiles[k].Values,	order='F')

	if(Config.transp_values_OK):
		TranspCuts	= Config.TranspCuts
		nZones		= len(Zones)
		for i in range (nZones):
			zone = "zone{:d}".format(i+1)
			for k in range(len(TranspCuts[0].Flux.Profiles)):
				h5_write(if_mesh,zone+"/ballooning"+ TranspCuts[0].Flux.Profiles[k].Name, Zones[i].Ballooning[k], order='F')			#Chi e values at soledge center

#	FeedbackTransp has been defined

	if(Config.feedback_transp_OK):
		mesher_feedback_transp	= "/mesher/FeedbackTransp"
		if_mesh.create_group(mesher_feedback_transp)

		FeedbackTransp = Config.FeedbackTransp

		h5_write(if_mesh, mesher_feedback_transp +"/Dmin",		FeedbackTransp.Data.Dmin)
		h5_write(if_mesh, mesher_feedback_transp +"/Dmax",		FeedbackTransp.Data.Dmax)
		h5_write(if_mesh, mesher_feedback_transp +"/Keep",		FeedbackTransp.Data.Keep)
		h5_write(if_mesh, mesher_feedback_transp +"/Gain",		FeedbackTransp.Data.Gain)
		h5_write(if_mesh, mesher_feedback_transp +"/GainG",		FeedbackTransp.Data.GainG)
		h5_write(if_mesh, mesher_feedback_transp +"/Lambda",		FeedbackTransp.Data.Lambda)
		h5_write(if_mesh, mesher_feedback_transp +"/UsePattern",	FeedbackTransp.Data.UsePattern)

		h5_write(if_mesh, mesher_feedback_transp +"/Nmin",		FeedbackTransp.Data.Nmin)
		h5_write(if_mesh, mesher_feedback_transp +"/Tmin",		FeedbackTransp.Data.Tmin)
		h5_write(if_mesh, mesher_feedback_transp +"/Nref",		FeedbackTransp.Data.NRef)
		h5_write(if_mesh, mesher_feedback_transp +"/Tref",		FeedbackTransp.Data.TRef)

		h5_write(if_mesh, mesher_feedback_transp +"/InterpMode",		FeedbackTransp.Data.InterpMode)

		nFeedbackProfiles = len(FeedbackTransp.Cut.Flux.Profiles)			
		h5_write(if_mesh, mesher_feedback_transp +"/nProfiles", nFeedbackProfiles)
			
		for k in range(nFeedbackProfiles):
			h5_write(if_mesh, mesher_feedback_transp +"/Prof_Name{:d}".format(k+1), FeedbackTransp.Cut.Flux.Profiles[k].Name, order='F')
		
		cut = mesher_feedback_transp +"/cut"
		Cgroup = if_mesh.create_group(cut)
		
		h5_write(if_mesh, cut+"/Flux_nz",  FeedbackTransp.Cut.Flux.nz,  order='F')
		h5_write(if_mesh, cut+"/Flux_dir", FeedbackTransp.Cut.Flux.dir, order='F')
		h5_write(if_mesh, cut+"/Flux_d12", FeedbackTransp.Cut.Flux.d12, order='F')
		h5_write(if_mesh, cut+"/Flux_r12", FeedbackTransp.Cut.Flux.r12, order='F')
		h5_write(if_mesh, cut+"/Flux_z12", FeedbackTransp.Cut.Flux.z12, order='F')
		
#		for k in range(len(FeedbackTransp.Cut.Flux.nz)):
#			h5_write(if_mesh, cut+"/Flux_Zones{:d}".format(k+1),FeedbackTransp.Cut.Flux.Zones[k])
						
		for k in range(nFeedbackProfiles):
			h5_write(if_mesh, cut+"/Flux_"+ FeedbackTransp.Cut.Flux.Profiles[k].Name+"_xValues",	FeedbackTransp.Cut.Flux.Profiles[k].xValues, order='F')
			h5_write(if_mesh, cut+"/Flux_"+ FeedbackTransp.Cut.Flux.Profiles[k].Name+"_Values",		FeedbackTransp.Cut.Flux.Profiles[k].Values, order='F')

#	FeedbackPuffing has been defined

	if(Config.feedback_puffing_OK):
		mesher_feedback_puffing		= "/mesher/FeedbackPuffing"
		if_mesh.create_group(mesher_feedback_puffing)

		h5_write(if_mesh, mesher_feedback_puffing +"/AutoTarget",	FeedbackPuffing.AutoTarget)
		h5_write(if_mesh, mesher_feedback_puffing +"/nTarget",		FeedbackPuffing.nTarget)
		h5_write(if_mesh, mesher_feedback_puffing +"/iTarget",		FeedbackPuffing.iTarget+1)
		h5_write(if_mesh, mesher_feedback_puffing +"/jTarget",		FeedbackPuffing.jTarget+1)
		h5_write(if_mesh, mesher_feedback_puffing +"/kTarget",		FeedbackPuffing.kTarget+1)
		h5_write(if_mesh, mesher_feedback_puffing +"/Gain",			FeedbackPuffing.Gain)
		h5_write(if_mesh, mesher_feedback_puffing +"/Tau_i",			FeedbackPuffing.Tau_i)
		h5_write(if_mesh, mesher_feedback_puffing +"/Tau_d",			FeedbackPuffing.Tau_d)
		h5_write(if_mesh, mesher_feedback_puffing +"/MaxPuff",		FeedbackPuffing.MaxPuff)
		h5_write(if_mesh, mesher_feedback_puffing +"/MinPuff",		FeedbackPuffing.MinPuff)


#	CustomPlots has been defined

	if(Config.CustomPlots_OK):
		custom_plots = "/CustomPlots"
		if_mesh.create_group(custom_plots)

		CustomPlots = Config.CustomPlots
		for iType in range(len(CustomPlots)):
			if(len(CustomPlots[iType]) > 0):
				custom_plots_type = custom_plots+"/{:s}".format(CUSTOM_PLOTS_TYPES[iType])
				if_mesh.create_group(custom_plots_type)		

				h5_write(if_mesh, custom_plots_type +"/nPlots", len(CustomPlots[iType]))
				for iPlot in range(len(CustomPlots[iType])):
					custom_plots_number = custom_plots_type+"/plot_{:d}".format(iPlot+1)
					if_mesh.create_group(custom_plots_number)
			
					h5_write(if_mesh, custom_plots_number +"/Name",		CustomPlots[iType][iPlot].Name)
					h5_write(if_mesh, custom_plots_number +"/kMagZones",	CustomPlots[iType][iPlot].kMagZones+1)
					h5_write(if_mesh, custom_plots_number +"/MagCoords",	CustomPlots[iType][iPlot].MagCoords+1)
					h5_write(if_mesh, custom_plots_number +"/nSubPlots", len(CustomPlots[iType][iPlot].SubPlots))

					if(len(CustomPlots[iType][iPlot].SubPlots) > 0):
						h5_write(if_mesh, custom_plots_number +"/nRows",	CustomPlots[iType][iPlot].nRows)
						h5_write(if_mesh, custom_plots_number +"/nCols",	CustomPlots[iType][iPlot].nCols)
						h5_write(if_mesh, custom_plots_number +"/SameX",	CustomPlots[iType][iPlot].SameX)
						for iSub in range(len(CustomPlots[iType][iPlot].SubPlots)):
							custom_sub_plot = custom_plots_number+"/sub_plot_{:d}".format(iSub+1)
							if_mesh.create_group(custom_sub_plot)

							h5_write(if_mesh, custom_sub_plot +"/Ions",		CustomPlots[iType][iPlot].SubPlots[iSub].Ions)
							h5_write(if_mesh, custom_sub_plot +"/ParLabels", CustomPlots[iType][iPlot].SubPlots[iSub].ParLabels)
							h5_write(if_mesh, custom_sub_plot +"/xLabel",		CustomPlots[iType][iPlot].SubPlots[iSub].xLabel+" ")
							h5_write(if_mesh, custom_sub_plot +"/yLabel",		CustomPlots[iType][iPlot].SubPlots[iSub].yLabel+" ")
							h5_write(if_mesh, custom_sub_plot +"/xMin",		CustomPlots[iType][iPlot].SubPlots[iSub].xMin+" ")
							h5_write(if_mesh, custom_sub_plot +"/xMax",		CustomPlots[iType][iPlot].SubPlots[iSub].xMax+" ")
							h5_write(if_mesh, custom_sub_plot +"/yMin",		CustomPlots[iType][iPlot].SubPlots[iSub].yMin+" ")
							h5_write(if_mesh, custom_sub_plot +"/yMax",		CustomPlots[iType][iPlot].SubPlots[iSub].yMax+" ")
							h5_write(if_mesh, custom_sub_plot +"/xLog",		CustomPlots[iType][iPlot].SubPlots[iSub].xLog)
							h5_write(if_mesh, custom_sub_plot +"/yLog",		CustomPlots[iType][iPlot].SubPlots[iSub].yLog)
							h5_write(if_mesh, custom_sub_plot +"/nParams", len(CustomPlots[iType][iPlot].SubPlots[iSub].Params))

							for iPar in range(len(CustomPlots[iType][iPlot].SubPlots[iSub].Params)):
								custom_sub_plots_param = custom_sub_plot+"/parameter_{:d}".format(iPar+1)
								if_mesh.create_group(custom_sub_plots_param)

								h5_write(if_mesh, custom_sub_plots_param +"/Name", CustomPlots[iType][iPlot].SubPlots[iSub].Params[iPar].Name)
								h5_write(if_mesh, custom_sub_plots_param +"/LineType", CustomPlots[iType][iPlot].SubPlots[iSub].Params[iPar].LineType)
								h5_write(if_mesh, custom_sub_plots_param +"/LineColor", CustomPlots[iType][iPlot].SubPlots[iSub].Params[iPar].LineColor)

	if_mesh.close()

	return

#	CODE SPECIFIC MESH DATA: SOLEDGE2D
#	========================

def write_mesh_S2D(if_mesh, Config): 

	Zones  = Config.Zones
	nZones = len(Zones)
	h5_write(if_mesh,"NZones", nZones)
	
	for k in range (nZones):
		Nx = Zones[k].Nx
		Nz = Zones[k].Nz
		zone = "zone{:d}".format(k+1)
		Zgroup = if_mesh.create_group(zone)

		h5_write(if_mesh,zone+"/Nx", Zones[k].Nx)
		h5_write(if_mesh,zone+"/Nz", Zones[k].Nz)
		h5_write(if_mesh,zone+"/x", Zones[k].x, order='F')
		h5_write(if_mesh,zone+"/z", Zones[k].z, order='F')
	
		h5_write(if_mesh,zone+"/xm", Zones[k].xb[:-1, :-1], order='F')
		h5_write(if_mesh,zone+"/xp", Zones[k].xb[1:,  :-1], order='F')
		h5_write(if_mesh,zone+"/zm", Zones[k].zb[:-1, :-1], order='F')
		h5_write(if_mesh,zone+"/zp", Zones[k].zb[:-1, 1:],  order='F')
	
		h5_write(if_mesh,zone+"/xmin", Zones[k].xb[0, 0])
		h5_write(if_mesh,zone+"/xmax", Zones[k].xb[-1,0])
		h5_write(if_mesh,zone+"/zmin", Zones[k].zb[0, 0])
		h5_write(if_mesh,zone+"/zmax", Zones[k].zb[0,-1])
		if(Config.use_penalization):
			h5_write(if_mesh,zone+"/chi", Zones[k].Chi, order='F')											#		optional : only if penalisation
		else:
			h5_write(if_mesh,zone+"/chi", np.zeros(Zones[k].Chi.shape,dtype='f8'), order='F')				#		optional : only if penalisation
				
		Ngroup = Zgroup.create_group("Neighbors")
		Neighbors = zone+"/Neighbors"
	
		h5_write(if_mesh,Neighbors+"/North", Zones[k].Neighbour.north + 1, order='F')			#Python index
		h5_write(if_mesh,Neighbors+"/South", Zones[k].Neighbour.south + 1, order='F')			#Python index
		if(Zones[k].Neighbour.east == EW_BORDER_NEIGHBOUR):											#Python index
			if(Config.use_penalization):				
				h5_write(if_mesh,Neighbors+"/East",   Zones[k].Neighbour.east + 1, order='F')		#Python index
			else:
				h5_write(if_mesh,Neighbors+"/East",   np.array([-3],dtype='i4'), order='F')		#Python index
		else:
			h5_write(if_mesh,Neighbors+"/East",   Zones[k].Neighbour.east + 1, order='F')
	
		if(Zones[k].Neighbour.west == EW_BORDER_NEIGHBOUR):															#Python index
			if(Config.use_penalization):				
				h5_write(if_mesh,Neighbors+"/West",   Zones[k].Neighbour.west + 1, order='F')
			else:
				h5_write(if_mesh,Neighbors+"/West",   np.array([-3],dtype='i4'), order='F')
		else:
			h5_write(if_mesh,Neighbors+"/West",   Zones[k].Neighbour.west + 1, order='F')
						
	
		MNgroup = Zgroup.create_group("MagNeighbors")
		MagNeighbors = zone+"/MagNeighbors"
	
		h5_write(if_mesh,MagNeighbors+"/North", Zones[k].MagNeighbour.north, order='F')
		h5_write(if_mesh,MagNeighbors+"/South", Zones[k].MagNeighbour.south, order='F')
		h5_write(if_mesh,MagNeighbors+"/East",  Zones[k].MagNeighbour.east,  order='F')
		h5_write(if_mesh,MagNeighbors+"/West",  Zones[k].MagNeighbour.west,  order='F')
	
		Br					= np.zeros((Nx+2,Nz+2), dtype='f8')
		Br[1:-1,1:-1]		= Zones[k].Br
		Br[0,1:-1]			= Zones[k].SouthP.Br
		Br[-1,1:-1]			= Zones[k].NorthP.Br
		Br[1:-1,0]			= Zones[k].WestP.Br
		Br[1:-1,-1]			= Zones[k].EastP.Br
	
		Bz					= np.zeros((Nx+2,Nz+2), dtype='f8')
		Bz[1:-1,1:-1]		= Zones[k].Bz
		Bz[0,1:-1]			= Zones[k].SouthP.Bz
		Bz[-1,1:-1]			= Zones[k].NorthP.Bz
		Bz[1:-1,0]			= Zones[k].WestP.Bz
		Bz[1:-1,-1]			= Zones[k].EastP.Bz
	
		Bphi				= np.zeros((Nx+2,Nz+2), dtype='f8')
		Bphi[1:-1,1:-1]		= Zones[k].Bphi
		Bphi[0,1:-1]		= Zones[k].SouthP.Bphi
		Bphi[-1,1:-1]		= Zones[k].NorthP.Bphi
		Bphi[1:-1,0]		= Zones[k].WestP.Bphi
		Bphi[1:-1,-1]		= Zones[k].EastP.Bphi
	
		Rcenter				= np.zeros((Nx+2,Nz+2), dtype='f8')
		Rcenter[1:-1,1:-1]	= Zones[k].gridRc
		Rcenter[0,1:-1]		= Zones[k].SouthP.R
		Rcenter[-1,1:-1]	= Zones[k].NorthP.R
		Rcenter[1:-1,0]		= Zones[k].WestP.R
		Rcenter[1:-1,-1]	= Zones[k].EastP.R
	
		Zcenter				= np.zeros((Nx+2,Nz+2), dtype='f8')
		Zcenter[1:-1,1:-1]	= Zones[k].gridZc
		Zcenter[0,1:-1]		= Zones[k].SouthP.Z
		Zcenter[-1,1:-1]	= Zones[k].NorthP.Z
		Zcenter[1:-1,0]		= Zones[k].WestP.Z
		Zcenter[1:-1,-1]	= Zones[k].EastP.Z
	
		h5_write(if_mesh,zone+"/Br",	 Br,   order='F')							#values at soledge center
		h5_write(if_mesh,zone+"/Bz",	 Bz,   order='F')
		h5_write(if_mesh,zone+"/Bphi",  Bphi, order='F')
		h5_write(if_mesh,zone+"/Rgeom", Rcenter, order='F')
		h5_write(if_mesh,zone+"/Zgeom", Zcenter, order='F')
	
		h5_write(if_mesh,zone+"/Rcorner", Zones[k].gridR, order='F')				#Soledge mesh corners
		h5_write(if_mesh,zone+"/Zcorner", Zones[k].gridZ, order='F')
		
	Megazones  = Config.Megazones	
	nMegazones = len(Megazones)
	h5_write(if_mesh,"NMegazones", nMegazones)
	
	for k in range(nMegazones):
		Mzone  = "megazone{:d}".format(k+1)
		Mgroup = if_mesh.create_group(Mzone)
	
		size	= len(Megazones[k].list)
		h5_write(if_mesh,Mzone+"/size",	 size)
		h5_write(if_mesh,Mzone+"/configuration", 	Megazones[k].list + 1, order='F')	#Python index
		h5_write(if_mesh,Mzone+"/isperiodic",			Megazones[k].isperiodic)

	return


#	CODE SPECIFIC MESH DATA: SOLEDGE3X
#	========================

def write_mesh_S3X(if_mesh, Config): 

	Nphi	= Config.Soledge3x.Nphi
	Rcore	= Config.Soledge3x.R0
	Bcore	= Config.Soledge3x.Bcore
	deltaR	= Config.Soledge3x.deltaR
	h5_write(if_mesh,zone+"/R0", Rcore)
	h5_write(if_mesh,zone+"/a0", deltaR)
	h5_write(if_mesh,zone+"/B0", Bcore)

	Zones  = Config.Zones
	nZones = len(Zones)
	h5_write(if_mesh,"NZones", nZones)
	
	phinodes = np.linspace(0,2.*np.pi,Nphi+1)

	for k in range (nZones):
		Nx = Zones[k].Nx
		Nz = Zones[k].Nz
		zone = "zone{:d}".format(k+1)
		Zgroup = if_mesh.create_group(zone)

		h5_write(if_mesh,zone+"/Nx", Zones[k].Nx)
		h5_write(if_mesh,zone+"/Nz", Zones[k].Nz)
		h5_write(if_mesh,zone+"/Nphi", Nphi)
		h5_write(if_mesh,zone+"/phicorners", phinodes)

		h5_write(if_mesh,zone+"/Bphi", (Zones[k].Bphi/Bcore)[:,:,np.newaxis])						#Auto transpose writing in C row-colums order
		h5_write(if_mesh,zone+"/Br", (Zones[k].Br/Bcore)[:,:,np.newaxis])
		h5_write(if_mesh,zone+"/Bz", (Zones[k].Bz/Bcore)[:,:,np.newaxis])

		h5_write(if_mesh,zone+"/Rcorners", ((Zones[k].gridR-Rcore)/deltaR)[:,:,np.newaxis])
		h5_write(if_mesh,zone+"/Zcorners", (Zones[k].gridZ/deltaR)[:,:,np.newaxis])
		h5_write(if_mesh,zone+"/psi", Zones[k].x[:,0], order='F')

	Neighbortab = np.zeros((nZones,6), dtype='i8')
	for k in range (nZones):
		if(Zones[k].Neighbour.south > -1): Neighbortab[k,0] = Zones[k].Neighbour.south+1				#Python index
		else:
			if(Zones[k].Neighbour.south == -2): Neighbortab[k,0] = -2
			else:								Neighbortab[k,0] = -3

		if(Zones[k].Neighbour.north > -1):		Neighbortab[k,1] = Zones[k].Neighbour.north + 1			#Python index
		else:
			if(zonesE.zone(k).Neighbour.north == -2): Neighbortab[k,1] = -2
			else: 									  Neighbortab[k,1] = -3

		if(Zones[k].Neighbour.west > -1):		Neighbortab[k,2] = Zones[k].Neighbour.west+1			#Python index
		else:									Neighbortab[k,2] = -1

		if(Zones[k].Neighbour.east > -1):		Neighbortab[k,3] = Zones[k].Neighbour.east+1			#Python index
		else:							        Neighbortab[k,3] = -1

	Neighbortab[k,5] = k
	Neighbortab[k,5] = k

	h5_write(if_mesh,zone+"/neighbors", Neighbortab[:,:,np.newaxis])									#Auto transpose writing in C row-colums order

	return