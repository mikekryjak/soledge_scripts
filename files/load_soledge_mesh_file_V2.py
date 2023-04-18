# Function definition is here

import types
import numpy as np
import os
import h5py
from mesh.recompute_megazones	import recompute_megazones
from routines.h5_routines		import h5_read
from routines.reset_config		import reset_config
from routines.globals			import DEBUG, DEF_BALLOONING_MODE, BALLOONING_NAMES, DEF_ZONE_MODE

#=========================================================
# This load mesh and config data from H5DF mesh file
#=========================================================

def load_soledge_mesh_file_V2(MeshFile): 

	if(DEBUG > 1): print("load_soledge_mesh_file: reading file ",MeshFile)
	try:	
		if_mesh = h5py.File(MeshFile, "r")
	except:
		if(DEBUG > 0): print("\tload_soledge_mesh_file: Unable to open ",MeshFile)
		return

#	Clear data

	Config			= reset_config()
	Config.MeshFile	= MeshFile
	X_points		= []
	Frontiers 		= []
	MagZones		= []
	MagMegazones	= []
	MagPMegazones	= []
	Zones			= []
	Megazones		= []
	TranspCuts		= []
	FeedbackTransp		= []
	
##	start reading

	MinVersion = 1
	try:
		FileVersion	= h5_read(if_mesh, 'Version',	keep_array=False)
	except:
		FileVersion = 0
		
	if(FileVersion < MinVersion):
		print("load_soledge_mesh_file: mesh file too old, you must try to update it with old_mesh_to_new.py")
		if_mesh.close()
		exit()
	
	try:
		mesher_fields		= "/mesher/Fields"
		Config.r2			= h5_read(if_mesh, mesher_fields+ '/r2',			order='F')
		Config.z2			= h5_read(if_mesh, mesher_fields+ '/z2',			order='F')
		Config.Br2			= h5_read(if_mesh, mesher_fields+ '/Br2',			order='F')
		Config.Bz2			= h5_read(if_mesh, mesher_fields+ '/Bz2',			order='F')
		Config.Bphi2		= h5_read(if_mesh, mesher_fields+ '/Bphi2',		order='F')
		Config.flux2		= h5_read(if_mesh, mesher_fields+ '/psi2',			order='F')
		try:
			Config.OutLength2	= h5_read(if_mesh, mesher_fields+ '/OutLength2',	keep_array=False)
		except:
			Config.OutLength2	= -1.
		Config.Smooth2		= h5_read(if_mesh, mesher_fields+ '/Smooth2',		keep_array=False)
		Config.Smooth2D		= h5_read(if_mesh, mesher_fields+ '/Smooth2D',	keep_array=False)
		Config.nr			= h5_read(if_mesh, mesher_fields+ '/Nr',			keep_array=False)
		Config.nz			= h5_read(if_mesh, mesher_fields+ '/Nz',			keep_array=False)
		Config.in_equ_OK	= True
		if(DEBUG > 1): print("\tload_soledge_mesh_file: in_equ = OK")
	except:
#	else:
		Config.in_equ_OK	= False
		if(DEBUG > 1): print("\tload_soledge_mesh_file: in_equ = NO")

	if(Config.in_equ_OK):
		try:
			mesher_fields = "/mesher/Fields"
			Config.extrapol_val	= h5_read(if_mesh, mesher_fields+ '/ExtrapolValue',	keep_array=False)
			Config.raise_value	= h5_read(if_mesh, mesher_fields+ '/RaiseValue', 	keep_array=False)

			mesher_x_points = "/mesher/XPoints"
			nX_points = h5_read(if_mesh, mesher_x_points +"/nXPoints", keep_array=False)
			X_points = []
			for k in range(nX_points):
				x_point = mesher_x_points + "/XPoint{:d}".format(k+1)
				X_points.append(types.SimpleNamespace())	
				X_points[-1].R		= h5_read(if_mesh,x_point + "/R",		keep_array=False)	
				X_points[-1].Z		= h5_read(if_mesh,x_point + "/Z",		keep_array=False)	
				X_points[-1].psi 	= h5_read(if_mesh,x_point + "/psi",		keep_array=False)	
				X_points[-1].index	= h5_read(if_mesh,x_point + "/index",	keep_array=False)			
				X_points[-1].sel		= False
	
				X_points[-1].cut = []
				for n in range(4):
					cut = x_point+"/cut{:d}".format(n+1)
					X_points[-1].cut.append(types.SimpleNamespace())
					X_points[-1].cut[-1].Rs		= h5_read(if_mesh,cut+"/Rs", 		keep_array=False)
					X_points[-1].cut[-1].Zs		= h5_read(if_mesh,cut+"/Zs", 		keep_array=False)
#					X_points[-1].cut[-1].Rs		= h5_read(if_mesh,cut+"/R", 		keep_array=False)[1]
#					X_points[-1].cut[-1].Zs		= h5_read(if_mesh,cut+"/Z", 		keep_array=False)[1]
					X_points[-1].cut[-1].psi	= h5_read(if_mesh,cut+"/psi", 		keep_array=False)
	
				X_points[-1].branch = []
				for n in range(4):
					branch = x_point+"/branch{:d}".format(n+1)
					X_points[-1].branch.append(types.SimpleNamespace())
					X_points[-1].branch[-1].R		= h5_read(if_mesh,branch+"/R",		keep_array=False)
					X_points[-1].branch[-1].Z		= h5_read(if_mesh,branch+"/Z",		keep_array=False)
					X_points[-1].branch[-1].theta	= h5_read(if_mesh,branch+"/theta",	keep_array=False)
					X_points[-1].branch[-1].arc		= h5_read(if_mesh,branch+"/arc",	keep_array=False)
			
			Config.extrapol_OK	= True
			if(DEBUG > 1): print("\tload_soledge_mesh_file: extrapol = OK")
		except:
			Config.X_points		= []
			Config.extrapol_OK	= False
			if(DEBUG > 1): print("\tload_soledge_mesh_file: extrapol = NO")

	if(Config.extrapol_OK):
		try:
			mesher_fields = "/mesher/Fields"
			Config.psicore		= h5_read(if_mesh, mesher_fields+ '/PsiCore', keep_array=False)
			Config.psiout		= h5_read(if_mesh, mesher_fields+ '/PsiEdge', keep_array=False)
			
			mesher_x_points = "/mesher/XPoints"
			for k in range(nX_points):
				x_point = mesher_x_points + "/XPoint{:d}".format(k+1)
				for n in range(4):
					cut = x_point+"/cut{:d}".format(n+1)
					X_points[k].cut[n].R		= h5_read(if_mesh,cut+"/R",			order='F')
					X_points[k].cut[n].Z		= h5_read(if_mesh,cut+"/Z",			order='F')
					X_points[k].cut[n].type		= h5_read(if_mesh,cut+"/type", 		keep_array=False)
					X_points[k].cut[n].psiR		= h5_read(if_mesh,cut+"/psiR",			order='F')
					X_points[k].cut[n].psiZ		= h5_read(if_mesh,cut+"/psiZ", 		order='F')		
					X_points[k].cut[n].psilim	= h5_read(if_mesh,cut+"/psilim",		keep_array=False)
			Config.xPoints_OK = True
			if(DEBUG > 1): print("\tload_soledge_mesh_file: xPoints = OK")
		except:
			Config.xPoints_OK	= False
			if(DEBUG > 1): print("\tload_soledge_mesh_file: xPoints = NO")

	if(not Config.extrapol_OK):
		try:
			X_points	= []
			nX_points 	= h5_read(if_mesh,"/config/nsep", keep_array=False)
			for k in range(nX_points):
				if(len(X_points) <= k): X_points.append(types.SimpleNamespace())
				X_points[k].psi = h5_read(if_mesh,"/config/psisep{:d}".format(k+1), keep_array=False)
			
			if(not Config.xPoints_OK): X_points[0].psicore = h5_read(if_mesh,"/config/psicore", keep_array=False)
		except:
			X_points	= []


	if(Config.xPoints_OK):
		try:
			Frontiers = []
			mesher_frontiers = "/mesher/Frontiers"
			nFrontiers = h5_read(if_mesh,mesher_frontiers +'/nFrontiers', keep_array=False)
			for k in range(nFrontiers):
				frontier = mesher_frontiers + "/Frontier{:d}".format(k+1)
				Frontiers.append(types.SimpleNamespace())
				Frontiers[-1].R		 = h5_read(if_mesh, frontier+"/R",    order='F')
				Frontiers[-1].Z		 = h5_read(if_mesh, frontier+"/Z",	   order='F')
				Frontiers[-1].P1	 = h5_read(if_mesh, frontier+"/P1",   order='F')
				Frontiers[-1].psiA	 = h5_read(if_mesh, frontier+"/psiA", keep_array=False	)
				Frontiers[-1].psiB	 = h5_read(if_mesh, frontier+"/psiB", keep_array=False	)
				Frontiers[-1].sel	 = False
				Frontiers[-1].psimin = np.min(Frontiers[-1].P1[:,2])
				Frontiers[-1].psimax = np.max(Frontiers[-1].P1[:,2])
				
			Config.Frontiers_OK	= True
			if(DEBUG > 1): print("\tload_soledge_mesh_file: Frontiers = OK")
		except:
			Frontiers	= []
			Config.Frontiers_OK	= False
			if(DEBUG > 1): print("\tload_soledge_mesh_file: Frontiers = NO")


	if(Config.Frontiers_OK):
		try:
			mesher_mag_zones = "/mesher/MagZones"
			nMagZones = h5_read(if_mesh, mesher_mag_zones +"/nMagZones", keep_array=False)
			MagZones = []
			for k in range(nMagZones):
				mag_zone = mesher_mag_zones + "/MagZone{:d}".format(k+1)
				MagZones.append(types.SimpleNamespace())	
				MagZones[k].east		= types.SimpleNamespace()
				MagZones[k].west		= types.SimpleNamespace()
				MagZones[k].north		= types.SimpleNamespace()
				MagZones[k].south		= types.SimpleNamespace()
				MagZones[k].pA			= types.SimpleNamespace()
				MagZones[k].pB			= types.SimpleNamespace()
				MagZones[k].pC			= types.SimpleNamespace()
				MagZones[k].pD			= types.SimpleNamespace()
				MagZones[k].Neighbour	= types.SimpleNamespace()
				
				MagZones[k].east.R		= h5_read(if_mesh, mag_zone+"/east_R",			order='F')
				MagZones[k].east.Z		= h5_read(if_mesh, mag_zone+"/east_Z",			order='F')
				MagZones[k].west.R		= h5_read(if_mesh, mag_zone+"/west_R",			order='F')
				MagZones[k].west.Z		= h5_read(if_mesh, mag_zone+"/west_Z",			order='F')
				MagZones[k].north.R		= h5_read(if_mesh, mag_zone+"/north_R",		order='F')
				MagZones[k].north.Z		= h5_read(if_mesh, mag_zone+"/north_Z",		order='F')
				MagZones[k].south.R		= h5_read(if_mesh, mag_zone+"/south_R",		order='F')
				MagZones[k].south.Z		= h5_read(if_mesh, mag_zone+"/south_Z",		order='F')
				MagZones[k].coord		= h5_read(if_mesh, mag_zone+"/coord",			order='F')
				MagZones[k].pA.coord	= h5_read(if_mesh, mag_zone+"/pA_coord",		order='F')
				MagZones[k].pB.coord	= h5_read(if_mesh, mag_zone+"/pB_coord",		order='F')
				MagZones[k].pC.coord	= h5_read(if_mesh, mag_zone+"/pC_coord",		order='F')
				MagZones[k].pD.coord	= h5_read(if_mesh, mag_zone+"/pD_coord",		order='F')
				
				MagZones[k].mz				= h5_read(if_mesh, mag_zone+"/mz", 	  			keep_array=False)
				MagZones[k].pmz				= h5_read(if_mesh, mag_zone+"/pmz", 	  			keep_array=False)
				MagZones[k].Xtype_east		= h5_read(if_mesh, mag_zone+"/Xtype_east", 	   	keep_array=False)
				MagZones[k].Xtype_west		= h5_read(if_mesh, mag_zone+"/Xtype_west",		keep_array=False)	
				
				MagZones[k].Neighbour.east	= h5_read(if_mesh, mag_zone+"/Neighbour_east",	keep_array=False)
				MagZones[k].Neighbour.west	= h5_read(if_mesh, mag_zone+"/Neighbour_west",	keep_array=False)
				MagZones[k].Neighbour.north	= h5_read(if_mesh, mag_zone+"/Neighbour_north",	keep_array=False)
				MagZones[k].Neighbour.south	= h5_read(if_mesh, mag_zone+"/Neighbour_south",	keep_array=False)
				
	
			mesher_pmegazones = "/mesher/MagPMegazones"
			nMagPMegazones = h5_read(if_mesh, mesher_pmegazones +"/nMagPMegazones", keep_array=False)
			MagPMegazones	= []
			for k in range(nMagPMegazones):
				pmegazone = mesher_pmegazones + "/MagPMegazone{:d}".format(k+1)
				MagPMegazones.append(types.SimpleNamespace())	
				MagPMegazones[k].list		= h5_read(if_mesh, pmegazone+"/list",		order='F')
						
			mesher_megazones = "/mesher/MagMegazones"
			nMagMegazones = h5_read(if_mesh, mesher_megazones +"/nMagMegazones", keep_array=False)
			MagMegazones				= []
			for k in range(nMagMegazones):
				megazone = mesher_megazones + "/MagMegazone{:d}".format(k+1)
				MagMegazones.append(types.SimpleNamespace())	
				MagMegazones[k].list		= h5_read(if_mesh, megazone+"/list",		order='F')
				MagMegazones[k].isperiodic	= h5_read(if_mesh, megazone+"/isperiodic",	keep_array=False)
	
			Config.MagZones_OK	= True
			if(DEBUG > 1): print("\tload_soledge_mesh_file: Limits = OK")
		except:
			Config.MagZones_OK	= False
			if(DEBUG > 1): print("\tload_soledge_mesh_file: Limits = NO")

#	Try if segments have been defined

	if(Config.MagZones_OK):
		try:
			mesher_pmegazones = "/mesher/MagPMegazones"
			nMagPMegazones =len(MagPMegazones)
			for k in range(nMagPMegazones):
				pmegazone = mesher_pmegazones + "/MagPMegazone{:d}".format(k+1)
				MagPMegazones[k].isaligned = h5_read(if_mesh, pmegazone+"/isaligned", 	keep_array=False)
	
				pmegazone_refpoints = pmegazone + "/refpoints"
				MagPMegazones[k].refpoints = types.SimpleNamespace()		
				if(MagPMegazones[k].isaligned):
					MagPMegazones[k].subrefpoints	   = [types.SimpleNamespace(), types.SimpleNamespace()]
					MagPMegazones[k].subrefpoints[0].R = h5_read(if_mesh, pmegazone_refpoints+"/R_sub1",	order='F')
					MagPMegazones[k].subrefpoints[0].Z = h5_read(if_mesh, pmegazone_refpoints+"/Z_sub1",	order='F')
					MagPMegazones[k].subrefpoints[1].R = h5_read(if_mesh, pmegazone_refpoints+"/R_sub2",	order='F')
					MagPMegazones[k].subrefpoints[1].Z = h5_read(if_mesh, pmegazone_refpoints+"/Z_sub2",	order='F')
					MagPMegazones[k].align_psimin	   = h5_read(if_mesh, pmegazone_refpoints+"/align_psimin", 	keep_array=False)
					MagPMegazones[k].align_psimax	   = h5_read(if_mesh, pmegazone_refpoints+"/align_psimax", 	keep_array=False)
		#				
				MagPMegazones[k].refpoints.R	= h5_read(if_mesh, pmegazone_refpoints+"/R",	order='F')
				MagPMegazones[k].refpoints.Z	= h5_read(if_mesh, pmegazone_refpoints+"/Z",	order='F')
				MagPMegazones[k].refpoints.nz	= h5_read(if_mesh, pmegazone_refpoints+"/nz", 	keep_array=False)
				MagPMegazones[k].refpoints.nzB  = h5_read(if_mesh, pmegazone_refpoints+"/nzB",	keep_array=False)
				
				MagPMegazones[k].ismeshed	  = True
	
				if(MagPMegazones[k].isaligned):
					for n in range(2):
						sub = "_sub{:d}".format(n+1)
						MagPMegazones[k].subrefpoints[n].nPoints	= h5_read(if_mesh, pmegazone_refpoints+"/nPoints"+sub, 		keep_array=False)							#Mesher parameters
						MagPMegazones[k].subrefpoints[n].RefineType	= h5_read(if_mesh, pmegazone_refpoints+"/RefineType"+sub,	keep_array=False)
						MagPMegazones[k].subrefpoints[n].RefineSide	= h5_read(if_mesh, pmegazone_refpoints+"/RefineSide"+sub, 	keep_array=False)
						MagPMegazones[k].subrefpoints[n].AdjustMode	= h5_read(if_mesh, pmegazone_refpoints+"/AdjustMode"+sub, 	keep_array=False)
						MagPMegazones[k].subrefpoints[n].ParamL		= h5_read(if_mesh, pmegazone_refpoints+"/ParamL"+sub,	 	keep_array=False)
						MagPMegazones[k].subrefpoints[n].ParamR		= h5_read(if_mesh, pmegazone_refpoints+"/ParamR"+sub, 		keep_array=False)
						if(MagPMegazones[k].subrefpoints[n].RefineType == 2):
							MagPMegazones[k].subrefpoints[n].xBezier = h5_read(if_mesh, pmegazone_refpoints+"/xBezier"+sub,	order='F')
							MagPMegazones[k].subrefpoints[n].yBezier = h5_read(if_mesh, pmegazone_refpoints+"/yBezier"+sub,	order='F')
				else:
					MagPMegazones[k].refpoints.nPoints		= h5_read(if_mesh, pmegazone_refpoints+"/nPoints", 		keep_array=False)							#Mesher parameters
					MagPMegazones[k].refpoints.RefineType	= h5_read(if_mesh, pmegazone_refpoints+"/RefineType", 	keep_array=False)
					MagPMegazones[k].refpoints.RefineSide	= h5_read(if_mesh, pmegazone_refpoints+"/RefineSide", 	keep_array=False)
					MagPMegazones[k].refpoints.AdjustMode	= h5_read(if_mesh, pmegazone_refpoints+"/AdjustMode", 	keep_array=False)
					MagPMegazones[k].refpoints.ParamL		= h5_read(if_mesh, pmegazone_refpoints+"/ParamL", 		keep_array=False)
					MagPMegazones[k].refpoints.ParamR		= h5_read(if_mesh, pmegazone_refpoints+"/ParamR", 		keep_array=False)
					if(MagPMegazones[k].refpoints.RefineType == 2):
						MagPMegazones[k].refpoints.xBezier = h5_read(if_mesh, pmegazone_refpoints+"/xBezier",	order='F')
						MagPMegazones[k].refpoints.yBezier = h5_read(if_mesh, pmegazone_refpoints+"/yBezier",	order='F')
	
				MagPMegazones[k].meshchanged = True
				
			mesher_megazones = "/mesher/MagMegazones"
			nMagMegazones =len(MagMegazones)
			for k in range(nMagMegazones):
				megazone_refpoints = mesher_megazones + "/MagMegazone{:d}/refpoints".format(k+1)
				MagMegazones[k].refpoints		= types.SimpleNamespace()		
	
				MagMegazones[k].refpoints.R		= h5_read(if_mesh, megazone_refpoints+"/R",		order='F')
				MagMegazones[k].refpoints.Z		= h5_read(if_mesh, megazone_refpoints+"/Z",		order='F')
				MagMegazones[k].refpoints.psi	= h5_read(if_mesh, megazone_refpoints+"/psi",	order='F')
				MagMegazones[k].ismeshed	  	= True
	
				MagMegazones[k].refpoints.nPoints		= h5_read(if_mesh, megazone_refpoints+"/nPoints", 		keep_array=False)							#Mesher parameters
				MagMegazones[k].refpoints.RefineType	= h5_read(if_mesh, megazone_refpoints+"/RefineType", 	keep_array=False)
				MagMegazones[k].refpoints.RefineSide	= h5_read(if_mesh, megazone_refpoints+"/RefineSide", 	keep_array=False)
				MagMegazones[k].refpoints.AdjustMode	= h5_read(if_mesh, megazone_refpoints+"/AdjustMode", 	keep_array=False)
				MagMegazones[k].refpoints.ParamL		= h5_read(if_mesh, megazone_refpoints+"/ParamL", 		keep_array=False)
				MagMegazones[k].refpoints.ParamR		= h5_read(if_mesh, megazone_refpoints+"/ParamR", 		keep_array=False)
				if(MagMegazones[k].refpoints.RefineType == 2):
					MagMegazones[k].refpoints.xBezier	= h5_read(if_mesh, megazone_refpoints+"/xBezier",	order='F')
					MagMegazones[k].refpoints.yBezier	= h5_read(if_mesh, megazone_refpoints+"/yBezier",	order='F')
	
			mesher_mag_zones = "/mesher/MagZones"
			for k in range(len(MagZones)):
				mag_zone = mesher_mag_zones + "/MagZone{:d}".format(k+1)
	
				MagZones[k].northaligned		= h5_read(if_mesh, mag_zone+"/northaligned",    	keep_array=False)
				MagZones[k].southaligned		= h5_read(if_mesh, mag_zone+"/southaligned",    	keep_array=False)
	
				if(MagZones[k].northaligned):
					MagZones[k].subNorth = [types.SimpleNamespace(), types.SimpleNamespace()]
					for n in range(2):
						sub = "_sub{:d}".format(n+1)
						MagZones[k].subNorth[n].R = h5_read(if_mesh, mag_zone+"/North_R"+sub,	order='F')
						MagZones[k].subNorth[n].Z = h5_read(if_mesh, mag_zone+"/North_Z"+sub,	order='F')
						MagZones[k].subNorth[n].ismeshed = True
						
				if(MagZones[k].southaligned):
					MagZones[k].subSouth = [types.SimpleNamespace(), types.SimpleNamespace()]
					for n in range(2):
						sub = "_sub{:d}".format(n+1)
						MagZones[k].subSouth[n].R = h5_read(if_mesh, mag_zone+"/South_R"+sub,	order='F')
						MagZones[k].subSouth[n].Z = h5_read(if_mesh, mag_zone+"/South_R"+sub,	order='F')
						MagZones[k].subSouth[n].ismeshed = True
				
				MagZones[k].east.ismeshed  		= True
				MagZones[k].west.ismeshed  		= True
				MagZones[k].north.ismeshed		= True
				MagZones[k].south.ismeshed		= True
				MagZones[k].orthomeshchanged	= True
	
			Config.Segments_OK	= True
			if(DEBUG > 1): print("\tload_soledge_mesh_file: Segments = OK")
		except:
			if(Config.MagZones_OK):
				for k in range(len(MagZones)):
					MagZones[k].northaligned	= False
					MagZones[k].southaligned	= False
					MagZones[k].east.ismeshed	= False
					MagZones[k].west.ismeshed	= False
					MagZones[k].north.ismeshed	= False
					MagZones[k].south.ismeshed	= False
					
				for k in range(len(MagMegazones)):
					MagMegazones[k].ismeshed 	= False
				
				for k in range(len(MagPMegazones)):
					MagPMegazones[k].ismeshed 		= False
					MagPMegazones[k].isaligned 		= False
					MagPMegazones[k].subrefpoints	= [types.SimpleNamespace(), types.SimpleNamespace()]
					for n in range(2):
						MagPMegazones[k].subrefpoints[n].R = np.array([])
						MagPMegazones[k].subrefpoints[n].Z = np.array([])
	
			Config.Segments_OK	= False
			if(DEBUG > 1): print("\tload_soledge_mesh_file: Segments = NO")

	if(Config.Segments_OK):		
		try:
			mesher_mag_zones = "/mesher/MagZones"
			for k in range(len(MagZones)):
				mag_zone = mesher_mag_zones + "/MagZone{:d}".format(k+1)
				MagZones[k].meshortho = h5_read(if_mesh, mag_zone+"/meshortho",    	keep_array=False)
				MagZones[k].gridR	= h5_read(if_mesh,mag_zone+"/Rcorner", order='F')				#Soledge mesh corners
				MagZones[k].gridZ	= h5_read(if_mesh,mag_zone+"/Zcorner", order='F')
				
			Config.MagGrid_OK	= True
			if(DEBUG > 1): print("\tload_soledge_mesh_file: MagGrid_OK = OK")
		except:
			Config.MagGrid_OK	= False
			if(DEBUG > 1): print("\tload_soledge_mesh_file: MagGrid_OK = NO")


#	Reading base (what is needed to soledge) mesh data	
	
	try:
		Config.Rwall 	= h5_read(if_mesh,"/wall/R", order='F')
		Config.Zwall 	= h5_read(if_mesh,"/wall/Z", order='F')
		Config.wall_OK	= True
	except:
		Config.wall_OK	= False

	try:
		Config.r2D		= h5_read(if_mesh,'/config/r',		order='F')
		Config.z2D		= h5_read(if_mesh,'/config/z',		order='F')
		Config.flux2D	= h5_read(if_mesh,'/config/psi',	order='F')
		try:
			Config.Br2D		= h5_read(if_mesh,'/config/Br',		order='F')
			Config.Bz2D		= h5_read(if_mesh,'/config/Bz',		order='F')
			Config.Bphi2D	= h5_read(if_mesh,'/config/Bphi',	order='F')
			Config.new_equ_OK	= True
			if(DEBUG > 1): print("\tload_soledge_mesh_file: new_equ = OK")
		except:
			Config.new_equ_OK	= False
			if(DEBUG > 1): print("\tload_soledge_mesh_file: new_equ = NO")
			
		Config.equ_OK	= True
		if(DEBUG > 1): print("\tload_soledge_mesh_file: equ = OK")
	except:
		Config.equ_OK	= False
		if(DEBUG > 1): print("\tload_soledge_mesh_file: equ = NO")
	
	try:
		Config.use_penalization = False

		nZones = h5_read(if_mesh,"NZones", keep_array=False)
		Zones  = []
		for k in range(nZones):
			Zones.append(types.SimpleNamespace())
	
			zone = "zone{:d}".format(k+1)
			Zones[k].Nx = h5_read(if_mesh,zone+"/Nx", keep_array=False)
			Zones[k].Nz = h5_read(if_mesh,zone+"/Nz", keep_array=False)

			Zones[k].x	= h5_read(if_mesh,zone+"/x", order='F')
			Zones[k].z	= h5_read(if_mesh,zone+"/z", order='F')

			Zones[k].xb = np.empty((Zones[k].Nx+1, Zones[k].Nz+1), dtype='f8')
			Zones[k].zb = np.empty((Zones[k].Nx+1, Zones[k].Nz+1), dtype='f8')
			
			Zones[k].xb[:-1, :-1]	= h5_read(if_mesh,zone+"/xm", order='F') 
			Zones[k].xb[1:,  :-1]	= h5_read(if_mesh,zone+"/xp", order='F')
			Zones[k].zb[:-1, :-1]	= h5_read(if_mesh,zone+"/zm", order='F') 
			Zones[k].zb[:-1, 1:]	= h5_read(if_mesh,zone+"/zp", order='F') 
			Zones[k].Chi = h5_read(if_mesh,zone+"/chi", order='F')	
			if(Zones[k].Chi.max() > 0.): Config.use_penalization = True
	
			Neighbors = zone+"/Neighbors"
	
			Zones[k].Neighbour		 = types.SimpleNamespace()
			Zones[k].Neighbour.north = h5_read(if_mesh,Neighbors+"/North", keep_array=False) - 1		#Python index
			Zones[k].Neighbour.south = h5_read(if_mesh,Neighbors+"/South", keep_array=False) - 1
			Zones[k].Neighbour.east  = h5_read(if_mesh,Neighbors+"/East", keep_array=False) - 1
			Zones[k].Neighbour.west	 = h5_read(if_mesh,Neighbors+"/West", keep_array=False) - 1 
			
			MagNeighbors = zone+"/MagNeighbors"
	
			Zones[k].MagNeighbour		= types.SimpleNamespace()
			Zones[k].MagNeighbour.north = h5_read(if_mesh,MagNeighbors+"/North", keep_array=False)
			Zones[k].MagNeighbour.south = h5_read(if_mesh,MagNeighbors+"/South", keep_array=False)
			Zones[k].MagNeighbour.east  = h5_read(if_mesh,MagNeighbors+"/East", keep_array=False)
			Zones[k].MagNeighbour.west  = h5_read(if_mesh,MagNeighbors+"/West", keep_array=False)
	
			Br		= h5_read(if_mesh,zone+"/Br",	  order='F')							#values at soledge center
			Bz		= h5_read(if_mesh,zone+"/Bz",	  order='F')
			Bphi	= h5_read(if_mesh,zone+"/Bphi",  order='F')
			Rcenter = h5_read(if_mesh,zone+"/Rgeom", order='F')
			Zcenter	= h5_read(if_mesh,zone+"/Zgeom", order='F')
	
			Zones[k].SouthP			= types.SimpleNamespace()
			Zones[k].NorthP			= types.SimpleNamespace()
			Zones[k].WestP			= types.SimpleNamespace()
			Zones[k].EastP			= types.SimpleNamespace()
			
			Zones[k].Br				= Br[1:-1,1:-1]
			Zones[k].SouthP.Br		= Br[0,1:-1] 
			Zones[k].NorthP.Br		= Br[-1,1:-1] 
			Zones[k].WestP.Br		= Br[1:-1,0] 
			Zones[k].EastP.Br		= Br[1:-1,-1] 
	
			Zones[k].Bz				= Bz[1:-1,1:-1] 
			Zones[k].SouthP.Bz		= Bz[0,1:-1] 
			Zones[k].NorthP.Bz		= Bz[-1,1:-1] 
			Zones[k].WestP.Bz		= Bz[1:-1,0] 
			Zones[k].EastP.Bz		= Bz[1:-1,-1] 
	
			Zones[k].Bphi			= Bphi[1:-1,1:-1] 
			Zones[k].SouthP.Bphi	= Bphi[0,1:-1] 
			Zones[k].NorthP.Bphi	= Bphi[-1,1:-1] 
			Zones[k].WestP.Bphi		= Bphi[1:-1,0] 
			Zones[k].EastP.Bphi		= Bphi[1:-1,-1] 
	
			Zones[k].gridRc			= Rcenter[1:-1,1:-1]
			Zones[k].SouthP.R		= Rcenter[0,1:-1] 
			Zones[k].NorthP.R		= Rcenter[-1,1:-1] 
			Zones[k].WestP.R		= Rcenter[1:-1,0] 
			Zones[k].EastP.R		= Rcenter[1:-1,-1] 
	
			Zones[k].gridZc			= Zcenter[1:-1,1:-1] 
			Zones[k].SouthP.Z		= Zcenter[0,1:-1] 
			Zones[k].NorthP.Z		= Zcenter[-1,1:-1]
			Zones[k].WestP.Z		= Zcenter[1:-1,0] 
			Zones[k].EastP.Z		= Zcenter[1:-1,-1]
	
			Zones[k].gridR		= h5_read(if_mesh,zone+"/Rcorner", order='F')				#Soledge mesh corners
			Zones[k].gridZ		= h5_read(if_mesh,zone+"/Zcorner", order='F')
	
		nMegazones = h5_read(if_mesh,"NMegazones", keep_array=False)
		Megazones  = []
		for k in range(nMegazones):
			if(len(Megazones) <= k): Megazones.append(types.SimpleNamespace())
	
			Mzone  = "megazone{:d}".format(k+1)
			Megazones[k].list		= h5_read(if_mesh,Mzone+"/configuration", order='F', keep_array=True) - 1	# Python index
			Megazones[k].isperiodic = h5_read(if_mesh,Mzone+"/isperiodic", keep_array=False)
			for iZone in Megazones[k].list: Zones[iZone].mz = k

#		read Zones data to recover splitting configuration
			
		try:
			mesher_zones = "/mesher/Zones"
			for k in range (nZones):
				zone = mesher_zones + "/Zone{:d}".format(k+1)
				Zones[k].magz  = h5_read(if_mesh, zone+"/MagZone", order='F')	- 1 			#Python index
		except:
			for k in range (nZones): Zones[k].magz = np.array([k,0,0], dtype = 'i4')

		Config.Mesh_OK			= True
		if(DEBUG > 1): print("\tload_soledge_mesh_file: Mesh = OK")
	except:
		Config.Mesh_OK			= False
		if(DEBUG > 1): print("\tload_soledge_mesh_file: Mesh = NO")

	if(Config.Mesh_OK):

#		Transport data
	
		try:
			nTranspCuts		= h5_read(if_mesh, "/mesher/TranspCuts/nCuts", 	   keep_array=False)
			nTranspProfiles = h5_read(if_mesh, "/mesher/TranspCuts/nProfiles", keep_array=False)	
			
			TranspCuts = []
			for i in range(nTranspCuts):
				cut = "/mesher/TranspCuts/cut{:d}".format(i+1)
				TranspCuts.append(types.SimpleNamespace())
	
				try:
					TranspCuts[i].BallooningMode = h5_read(if_mesh, cut+"/BallooningMode", keep_array=False)	
				except:
					TranspCuts[i].BallooningMode = DEF_BALLOONING_MODE

				try:
					TranspCuts[i].ZoneMode = h5_read(if_mesh, cut+"/ZoneMode", keep_array=False)	
				except:
					TranspCuts[i].ZoneMode  = DEF_ZONE_MODE
					
				TranspCuts[i].Flux		= types.SimpleNamespace()
				TranspCuts[i].Flux.nz	= h5_read(if_mesh, cut+"/Flux_nz",  order='F')

				TranspCuts[i].Flux.d12	= h5_read(if_mesh, cut+"/Flux_d12", order='F')
				TranspCuts[i].Flux.r12	= h5_read(if_mesh, cut+"/Flux_r12", order='F')
				TranspCuts[i].Flux.z12	= h5_read(if_mesh, cut+"/Flux_z12", order='F')
				try:
					TranspCuts[i].Flux.dir	= h5_read(if_mesh, cut+"/Flux_dir", order='F')
				except:
					TranspCuts[i].Flux.dir	= np.zeros(TranspCuts[i].Flux.d12.shape[0], dtype='f8')

				TranspCuts[i].Flux.MagZones = []
				for k in range(len(TranspCuts[i].Flux.nz)):
					TranspCuts[i].Flux.MagZones.append(h5_read(if_mesh, cut+"/Flux_Zones{:d}".format(k+1), order='F'))
				
				TranspCuts[i].Flux.Profiles = []
				for k in range(nTranspProfiles):
					TranspCuts[i].Flux.Profiles.append(types.SimpleNamespace())
					if(i == 0):
						TranspCuts[i].Flux.Profiles[k].Name = h5_read(if_mesh, "/mesher/TranspCuts/Prof_Name{:d}".format(k+1), keep_array=False)
					else:
						TranspCuts[i].Flux.Profiles[k].Name = TranspCuts[0].Flux.Profiles[k].Name
							
				for k in range(nTranspProfiles):
					TranspCuts[i].Flux.Profiles[k].xValues	= h5_read(if_mesh, cut+"/Flux_"+ TranspCuts[i].Flux.Profiles[k].Name+"_xValues", order='F')
					TranspCuts[i].Flux.Profiles[k].Values	= h5_read(if_mesh, cut+"/Flux_"+ TranspCuts[i].Flux.Profiles[k].Name+"_Values", order='F')
		
				TranspCuts[i].Theta		= types.SimpleNamespace()
				TranspCuts[i].Theta.nz	= h5_read(if_mesh, cut+"/Theta_nz", order='F')
				TranspCuts[i].Theta.d12	= h5_read(if_mesh, cut+"/Theta_d12", order='F')	
							
				TranspCuts[i].Theta.Profiles = []
				for k in range(nTranspProfiles):
					TranspCuts[i].Theta.Profiles.append(types.SimpleNamespace())
					TranspCuts[i].Theta.Profiles[k].Name = TranspCuts[i].Flux.Profiles[k].Name
					TranspCuts[i].Theta.Profiles[k].xValues = h5_read(if_mesh, cut+"/Theta_"+ TranspCuts[i].Theta.Profiles[k].Name+"_xValues", order='F')
					TranspCuts[i].Theta.Profiles[k].Values	= h5_read(if_mesh, cut+"/Theta_"+ TranspCuts[i].Theta.Profiles[k].Name+"_Values", order='F')
			Config.transp_prof_OK = True
			if(DEBUG > 1): print("\tload_soledge_mesh_file: transp_prof = OK")
		except:
			TranspCuts = []
			Config.transp_prof_OK = False
			if(DEBUG > 1): print("\tload_soledge_mesh_file: transp_prof = NO")
		
		try:
			if(not Config.transp_prof_OK):
				TranspCuts					= [types.SimpleNamespace()]
				TranspCuts[0].Flux			= types.SimpleNamespace()
				TranspCuts[0].Flux.Profiles	= []
				for i in range(len(BALLOONING_NAMES)):
					TranspCuts[0].Flux.Profiles.append(types.SimpleNamespace())
					TranspCuts[0].Flux.Profiles[-1].Name = BALLOONING_NAMES[i]

			nZones = len(Zones)
			for i in range (nZones):
				zone = "zone{:d}".format(i+1)
				Zones[i].Ballooning = []
				for k in range(len(TranspCuts[0].Flux.Profiles)):
					Zones[i].Ballooning.append(h5_read(if_mesh,zone+"/ballooning"+ TranspCuts[0].Flux.Profiles[k].Name, order='F'))
			Config.transp_values_OK = True
			if(DEBUG > 1): print("\tload_soledge_mesh_file: transp_values = OK")
		except:
			if(not Config.transp_prof_OK): TranspCuts = []
			Config.transp_values_OK = False
			if(DEBUG > 1): print("\tload_soledge_mesh_file: transp_values = NO")
		
#		FeedbackTransp data

		try:
			mesher_feedback		= "/mesher/Feedback"
	
			FeedbackTransp					= types.SimpleNamespace()
			FeedbackTransp.Data				= types.SimpleNamespace()
			FeedbackTransp.Data.Dmin		= h5_read(if_mesh, mesher_feedback +"/Dmin",		keep_array=False)
			FeedbackTransp.Data.Dmax		= h5_read(if_mesh, mesher_feedback +"/Dmax",		keep_array=False)
			FeedbackTransp.Data.Keep		= h5_read(if_mesh, mesher_feedback +"/Keep",		keep_array=False)
			FeedbackTransp.Data.Gain		= h5_read(if_mesh, mesher_feedback +"/Gain",		keep_array=False)
			FeedbackTransp.Data.GainG		= h5_read(if_mesh, mesher_feedback +"/GainG",		keep_array=False)
			FeedbackTransp.Data.Lambda		= h5_read(if_mesh, mesher_feedback +"/Lambda",		keep_array=False)
			FeedbackTransp.Data.UsePattern	= h5_read(if_mesh, mesher_feedback +"/UsePattern",	keep_array=False)
			FeedbackTransp.Data.UsePattern	= ("No", "1/Bt", "1/B").index(FeedbackTransp.Data.UsePattern)
	
			FeedbackTransp.Data.Nmin		= h5_read(if_mesh, mesher_feedback +"/Nmin",	keep_array=False)
			FeedbackTransp.Data.Tmin		= h5_read(if_mesh, mesher_feedback +"/Tmin",	keep_array=False)
			FeedbackTransp.Data.NRef		= h5_read(if_mesh, mesher_feedback +"/Nref",	keep_array=False)
			FeedbackTransp.Data.TRef		= h5_read(if_mesh, mesher_feedback +"/Tref",	keep_array=False)
	
			FeedbackTransp.Data.InterpMode	= h5_read(if_mesh, mesher_feedback +"/InterpMode",		keep_array=False)
	
			
			cut = mesher_feedback +"/cut"
			FeedbackTransp.Cut			= types.SimpleNamespace()
			FeedbackTransp.Cut.Flux		= types.SimpleNamespace()
			FeedbackTransp.Cut.Flux.nz	= h5_read(if_mesh, cut+"/Flux_nz",  order='F')
			FeedbackTransp.Cut.Flux.d12	= h5_read(if_mesh, cut+"/Flux_d12", order='F')
			FeedbackTransp.Cut.Flux.r12	= h5_read(if_mesh, cut+"/Flux_r12", order='F')
			FeedbackTransp.Cut.Flux.z12	= h5_read(if_mesh, cut+"/Flux_z12", order='F')
			try:
				FeedbackTransp.Cut.Flux.dir	= h5_read(if_mesh, cut+"/Flux_dir", order='F')
			except:
				FeedbackTransp.Cut.Flux.dir	= np.zeros(FeedbackTransp.Cut.Flux.d12.shape[0], dtype='f8')

#			FeedbackTransp.Cut.Flux.Zones = []
#			for k in range(len(FeedbackTransp.Cut.Flux.nz)):
#			FeedbackTransp.Cut.Flux.Zones.append(h5_read(if_mesh, cut+"/Flux_Zones{:d}".format(k+1), order='F'))
				
			nFeedbackProfiles = h5_read(if_mesh, mesher_feedback +"/nProfiles", keep_array=False)
				
			FeedbackTransp.Cut.Flux.Profiles = []
			for k in range(nFeedbackProfiles):
				FeedbackTransp.Cut.Flux.Profiles.append(types.SimpleNamespace())
				FeedbackTransp.Cut.Flux.Profiles[k].Name = h5_read(if_mesh, "/mesher/Feedback/Prof_Name{:d}".format(k+1), keep_array=False)
							
			for k in range(nFeedbackProfiles):
				FeedbackTransp.Cut.Flux.Profiles[k].xValues	= h5_read(if_mesh, cut+"/Flux_"+ FeedbackTransp.Cut.Flux.Profiles[k].Name+"_xValues", order='F')
				FeedbackTransp.Cut.Flux.Profiles[k].Values	= h5_read(if_mesh, cut+"/Flux_"+ FeedbackTransp.Cut.Flux.Profiles[k].Name+"_Values", order='F')
		
			Config.feedback_transp_OK = True
			if(DEBUG > 1): print("\tload_soledge_mesh_file: feedback_tranps = OK")
		except:
			FeedbackTransp		   	  = 0
			Config.feedback_transp_OK = False
			if(DEBUG > 1): print("\tload_soledge_mesh_file: feedback_transp = NO")

		
#		FeedbackPuffing data

		FeedbackPuffing		   	   = 0
		Config.feedback_puffing_OK = False
		if(DEBUG > 1): print("\tload_soledge_mesh_file: feedback_puffing = NO")


	else:
		TranspCuts					= []
		FeedbackTransp				= 0
		FeedbackPuffing				= 0
		Config.transp_prof_OK		= False
		Config.transp_values_OK 	= False
		Config.feedback_transp_OK	= False
		Config.feedback_puffing_OK	= False

	if_mesh.close()
	
	Config.X_points			= X_points
	Config.Frontiers		= Frontiers

	if(Config.MagZones_OK):
		Config.MagZones			= MagZones
		Config.MagMegazones		= MagMegazones
		Config.MagPMegazones	= MagPMegazones
	else:
		Config.MagZones			= Zones
		Config.MagMegazones		= Megazones

	Config.Zones			= Zones
	Config.Megazones		= Megazones
	Config.TranspCuts		= TranspCuts
	Config.FeedbackTransp	= FeedbackTransp
	Config.FeedbackPuffing	= FeedbackPuffing

	if(Config.Mesh_OK): recompute_megazones(Config)

	if(DEBUG > 1): print("load_soledge_mesh_file: Completed")
	
	return Config
