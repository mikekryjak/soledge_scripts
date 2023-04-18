# Function definition is here

import types
import numpy as np
import os
import h5py
from routines.h5_routines		import h5_read
from routines.reset_config		import reset_config
from routines.globals			import DEBUG, DEF_BALLOONING_MODE, BALLOONING_NAMES, DEF_ZONE_MODE
from mesh.recompute_megazones	import recompute_megazones

#=========================================================
# This load mesh and config data from H5DF mesh file
#=========================================================

def load_soledge_mesh_file_V0(MeshFile): 

	print("load_soledge_mesh_file: reading file ",MeshFile)
	try:	
		if_mesh = h5py.File(MeshFile, "r")
	except:
		if(DEBUG > 1): print("\tload_soledge_mesh_file: Unable to open ",MeshFile)
		return

#	Clear data

	Mesher			= reset_config()
	Mesher.MeshFile	= MeshFile
	X_points		= []
	Zones			= []
	Megazones		= []
	PMegazones		= []
	
##	start reading
	try:
		FileVersion	= h5_read(if_mesh, 'Version',	keep_array=False)
	except:
		FileVersion = 0

	if(FileVersion > 1):
		print("load_soledge_mesh_file: Im am sorry but mesh file is non managed from this version routine")
		if_mesh.close()
		exit()

	try:
#	if(True):
		mesher_fields = "/mesher/Fields"
		Mesher.r2			= h5_read(if_mesh, mesher_fields+ '/r2',			order='F')
		Mesher.z2			= h5_read(if_mesh, mesher_fields+ '/z2',			order='F')
		Mesher.Br2			= h5_read(if_mesh, mesher_fields+ '/Br2',			order='F')
		Mesher.Bz2			= h5_read(if_mesh, mesher_fields+ '/Bz2',			order='F')
		Mesher.Bphi2		= h5_read(if_mesh, mesher_fields+ '/Bphi2',		order='F')
		Mesher.flux2		= h5_read(if_mesh, mesher_fields+ '/psi2',			order='F')
		try:
			Mesher.OutLength2	= h5_read(if_mesh, mesher_fields+ '/OutLength2',	keep_array=False)
		except:
			Mesher.OutLength2	= -1.

		Mesher.Smooth2		= h5_read(if_mesh, mesher_fields+ '/Smooth2',		keep_array=False)
		Mesher.Smooth2D		= h5_read(if_mesh, mesher_fields+ '/Smooth2D',	keep_array=False)
		Mesher.nr			= h5_read(if_mesh, mesher_fields+ '/Nr',			keep_array=False)
		Mesher.nz			= h5_read(if_mesh, mesher_fields+ '/Nz',			keep_array=False)
		Mesher.in_equ_OK	= True
		if(DEBUG > 1): print("\tload_soledge_mesh_file: in_equ = OK")
	except:
#	else:
		Mesher.in_equ_OK	= False
		if(DEBUG > 1): print("\tload_soledge_mesh_file: in_equ = NO")

	if(Mesher.in_equ_OK):
		try:
			mesher_fields = "/mesher/Fields"
			Mesher.extrapol_val	= h5_read(if_mesh, mesher_fields+ '/ExtrapolValue',	keep_array=False)
			Mesher.raise_value	= h5_read(if_mesh, mesher_fields+ '/RaiseValue', 	keep_array=False)

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
			
			Mesher.extrapol_OK	= True
			if(DEBUG > 1): print("\tload_soledge_mesh_file: extrapol = OK")
		except:
			X_points		= []
			Mesher.extrapol_OK	= False
			if(DEBUG > 1): print("\tload_soledge_mesh_file: extrapol = NO")

	if(Mesher.extrapol_OK):
		try:
			mesher_fields = "/mesher/Fields"
			Mesher.psicore		= h5_read(if_mesh, mesher_fields+ '/PsiCore', keep_array=False)
			Mesher.psiout		= h5_read(if_mesh, mesher_fields+ '/PsiEdge', keep_array=False)
			
			mesher_x_points = "/mesher/XPoints"
			for k in range(nX_points):
				x_point = mesher_x_points + "/XPoint{:d}".format(k+1)
				for n in range(4):
					cut = x_point+"/cut{:d}".format(n+1)
					X_points[k].cut[n].R		= h5_read(if_mesh,cut+"/R",			order='F')
					X_points[k].cut[n].Z		= h5_read(if_mesh,cut+"/Z",			order='F')
					X_points[k].cut[n].type		= h5_read(if_mesh,cut+"/type", 		keep_array=False)
					X_points[k].cut[n].psiR		= h5_read(if_mesh,cut+"/psiR",		order='F')
					X_points[k].cut[n].psiZ		= h5_read(if_mesh,cut+"/psiZ", 		order='F')		
					X_points[k].cut[n].psilim	= h5_read(if_mesh,cut+"/psilim",	keep_array=False)
			Mesher.xPoints_OK = True
			if(DEBUG > 1): print("\tload_soledge_mesh_file: xPoints = OK")
		except:
			Mesher.xPoints_OK	= False
			if(DEBUG > 1): print("\tload_soledge_mesh_file: xPoints = NO")

	if(not Mesher.extrapol_OK):
		try:
			nX_points 	= h5_read(if_mesh,"/config/nsep", keep_array=False)
			for k in range(nX_points):
				if(len(X_points) <= k): X_points.append(types.SimpleNamespace())
				X_points[k].psi = h5_read(if_mesh,"/config/psisep{:d}".format(k+1), keep_array=False)
			
			if(not Mesher.xPoints_OK): X_points[0].psicore = h5_read(if_mesh,"/config/psicore", keep_array=False)
		except:
			X_points	= []


	if(Mesher.xPoints_OK):
		try:
			Frontiers = []
			mesher_frontiers = "/mesher/Frontiers"
			nFrontiers = h5_read(if_mesh,mesher_frontiers +'/nFrontiers', keep_array=False)
			for k in range(nFrontiers):
				frontier = mesher_frontiers + "/Frontier{:d}".format(k+1)
				Frontiers.append(types.SimpleNamespace())
				Frontiers[-1].R		 = h5_read(if_mesh, frontier+"/R", 	  order='F')
				Frontiers[-1].Z		 = h5_read(if_mesh, frontier+"/Z",	  order='F')
				Frontiers[-1].P1	 	 = h5_read(if_mesh, frontier+"/P1",	  order='F')
				Frontiers[-1].psiA	 = h5_read(if_mesh, frontier+"/psiA", keep_array=False	)
				Frontiers[-1].psiB	 = h5_read(if_mesh, frontier+"/psiB", keep_array=False	)
				Frontiers[-1].sel	 = False
				Frontiers[-1].psimin	 = np.min(Frontiers[-1].P1[:,2])
				Frontiers[-1].psimax  = np.max(Frontiers[-1].P1[:,2])
				
			Mesher.Frontiers_OK	= True
			if(DEBUG > 1): print("\tload_soledge_mesh_file: Frontiers = OK")
		except:
			Frontiers	= []
			Mesher.Frontiers_OK	= False
			if(DEBUG > 1): print("\tload_soledge_mesh_file: Frontiers = NO")
	else:	Frontiers	= []

	if(Mesher.Frontiers_OK):
		try:
			mesher_zones = "/mesher/Zones"
			nZones = h5_read(if_mesh, mesher_zones +"/nZones", keep_array=False)
			Zones = []
			for k in range(nZones):
				zone = mesher_zones + "/Zone{:d}".format(k+1)
				Zones.append(types.SimpleNamespace())	
				Zones[k].east		= types.SimpleNamespace()
				Zones[k].west		= types.SimpleNamespace()
				Zones[k].north		= types.SimpleNamespace()
				Zones[k].south		= types.SimpleNamespace()
				Zones[k].pA			= types.SimpleNamespace()
				Zones[k].pB			= types.SimpleNamespace()
				Zones[k].pC			= types.SimpleNamespace()
				Zones[k].pD			= types.SimpleNamespace()
				Zones[k].Neighbour	= types.SimpleNamespace()
				
				Zones[k].east.R		= h5_read(if_mesh, zone+"/east_R",			order='F')
				Zones[k].east.Z		= h5_read(if_mesh, zone+"/east_Z",			order='F')
				Zones[k].west.R		= h5_read(if_mesh, zone+"/west_R",			order='F')
				Zones[k].west.Z		= h5_read(if_mesh, zone+"/west_Z",			order='F')
				Zones[k].north.R	= h5_read(if_mesh, zone+"/north_R",			order='F')
				Zones[k].north.Z	= h5_read(if_mesh, zone+"/north_Z",			order='F')
				Zones[k].south.R	= h5_read(if_mesh, zone+"/south_R",			order='F')
				Zones[k].south.Z	= h5_read(if_mesh, zone+"/south_Z",			order='F')
				Zones[k].coord		= h5_read(if_mesh, zone+"/coord",			order='F')
				Zones[k].pA.coord	= h5_read(if_mesh, zone+"/pA_coord",		order='F')
				Zones[k].pB.coord	= h5_read(if_mesh, zone+"/pB_coord",		order='F')
				Zones[k].pC.coord	= h5_read(if_mesh, zone+"/pC_coord",		order='F')
				Zones[k].pD.coord	= h5_read(if_mesh, zone+"/pD_coord",		order='F')
				
				Zones[k].mz					= h5_read(if_mesh, zone+"/mz", 	  				keep_array=False)
				Zones[k].pmz				= h5_read(if_mesh, zone+"/pmz", 	  			keep_array=False)
				Zones[k].Xtype_east			= h5_read(if_mesh, zone+"/Xtype_east", 	   		keep_array=False)
				Zones[k].Xtype_west			= h5_read(if_mesh, zone+"/Xtype_west",			keep_array=False)	
				
				Zones[k].Neighbour.east		= h5_read(if_mesh, zone+"/Neighbour_east",		keep_array=False)
				Zones[k].Neighbour.west		= h5_read(if_mesh, zone+"/Neighbour_west",		keep_array=False)
				Zones[k].Neighbour.north	= h5_read(if_mesh, zone+"/Neighbour_north",		keep_array=False)
				Zones[k].Neighbour.south	= h5_read(if_mesh, zone+"/Neighbour_south",		keep_array=False)
				
	
				
			mesher_pmegazones = "/mesher/PMegazones"
			nPMegazones = h5_read(if_mesh, mesher_pmegazones +"/nPMegazones", keep_array=False)
			PMegazones	= []
			for k in range(nPMegazones):
				pmegazone = mesher_pmegazones + "/PMegazone{:d}".format(k+1)
				PMegazones.append(types.SimpleNamespace())	
				PMegazones[k].list		= h5_read(if_mesh, pmegazone+"/list",		order='F')
					
			mesher_megazones = "/mesher/Megazones"
			nMegazones = h5_read(if_mesh, mesher_megazones +"/nMegazones", keep_array=False)
			Megazones				= []
			for k in range(nMegazones):
				megazone = mesher_megazones + "/Megazone{:d}".format(k+1)
				Megazones.append(types.SimpleNamespace())	
				Megazones[k].list		= h5_read(if_mesh, megazone+"/list",		order='F')
				Megazones[k].isperiodic	= h5_read(if_mesh, megazone+"/isperiodic",	keep_array=False)
	
			Mesher.Zones_OK	= True
			if(DEBUG > 1): print("\tload_soledge_mesh_file: Limits = OK")
		except:
			Zones				= []
			PMegazones			= []
			Megazones			= []
			Mesher.Zones_OK	= False
			if(DEBUG > 1): print("\tload_soledge_mesh_file: Limits = NO")

#	Try if segments have been defined

	if(Mesher.Zones_OK):
		try:
			mesher_pmegazones = "/mesher/PMegazones"
			nPMegazones =len(PMegazones)
			for k in range(nPMegazones):
				pmegazone = mesher_pmegazones + "/PMegazone{:d}".format(k+1)
				PMegazones[k].isaligned = h5_read(if_mesh, pmegazone+"/isaligned", 	keep_array=False)
	
				pmegazone_refpoints = pmegazone + "/refpoints"
				PMegazones[k].refpoints = types.SimpleNamespace()		
				if(PMegazones[k].isaligned):
					PMegazones[k].subrefpoints		= [types.SimpleNamespace(), types.SimpleNamespace()]
					PMegazones[k].subrefpoints[0].R = h5_read(if_mesh, pmegazone_refpoints+"/R_sub1",	order='F')
					PMegazones[k].subrefpoints[0].Z = h5_read(if_mesh, pmegazone_refpoints+"/Z_sub1",	order='F')
					PMegazones[k].subrefpoints[1].R = h5_read(if_mesh, pmegazone_refpoints+"/R_sub2",	order='F')
					PMegazones[k].subrefpoints[1].Z = h5_read(if_mesh, pmegazone_refpoints+"/Z_sub2",	order='F')
					PMegazones[k].align_psimin		= h5_read(if_mesh, pmegazone_refpoints+"/align_psimin", 	keep_array=False)
					PMegazones[k].align_psimax		= h5_read(if_mesh, pmegazone_refpoints+"/align_psimax", 	keep_array=False)
		#				
				PMegazones[k].refpoints.R	= h5_read(if_mesh, pmegazone_refpoints+"/R",	order='F')
				PMegazones[k].refpoints.Z	= h5_read(if_mesh, pmegazone_refpoints+"/Z",	order='F')
				PMegazones[k].refpoints.nz	= h5_read(if_mesh, pmegazone_refpoints+"/nz", 	keep_array=False)
				PMegazones[k].refpoints.nzB = h5_read(if_mesh, pmegazone_refpoints+"/nzB", 	keep_array=False)
				
				PMegazones[k].ismeshed	  = True
	
				if(PMegazones[k].isaligned):
					for n in range(2):
						sub = "_sub{:d}".format(n+1)
						PMegazones[k].subrefpoints[n].nPoints		= h5_read(if_mesh, pmegazone_refpoints+"/nPoints"+sub, 		keep_array=False)							#Mesher parameters
						PMegazones[k].subrefpoints[n].RefineType	= h5_read(if_mesh, pmegazone_refpoints+"/RefineType"+sub,	keep_array=False)
						PMegazones[k].subrefpoints[n].RefineSide	= h5_read(if_mesh, pmegazone_refpoints+"/RefineSide"+sub, 	keep_array=False)
						PMegazones[k].subrefpoints[n].AdjustMode	= h5_read(if_mesh, pmegazone_refpoints+"/AdjustMode"+sub, 	keep_array=False)
						PMegazones[k].subrefpoints[n].ParamL		= h5_read(if_mesh, pmegazone_refpoints+"/ParamL"+sub,	 	keep_array=False)
						PMegazones[k].subrefpoints[n].ParamR		= h5_read(if_mesh, pmegazone_refpoints+"/ParamR"+sub, 		keep_array=False)
						if(PMegazones[k].subrefpoints[n].RefineType == 2):
							PMegazones[k].subrefpoints[n].xBezier	= h5_read(if_mesh, pmegazone_refpoints+"/xBezier"+sub,	order='F')
							PMegazones[k].subrefpoints[n].yBezier	= h5_read(if_mesh, pmegazone_refpoints+"/yBezier"+sub,	order='F')
				else:
					PMegazones[k].refpoints.nPoints		= h5_read(if_mesh, pmegazone_refpoints+"/nPoints", 		keep_array=False)							#Mesher parameters
					PMegazones[k].refpoints.RefineType	= h5_read(if_mesh, pmegazone_refpoints+"/RefineType", 	keep_array=False)
					PMegazones[k].refpoints.RefineSide	= h5_read(if_mesh, pmegazone_refpoints+"/RefineSide", 	keep_array=False)
					PMegazones[k].refpoints.AdjustMode	= h5_read(if_mesh, pmegazone_refpoints+"/AdjustMode", 	keep_array=False)
					PMegazones[k].refpoints.ParamL		= h5_read(if_mesh, pmegazone_refpoints+"/ParamL", 		keep_array=False)
					PMegazones[k].refpoints.ParamR		= h5_read(if_mesh, pmegazone_refpoints+"/ParamR", 		keep_array=False)
					if(PMegazones[k].refpoints.RefineType == 2):
						PMegazones[k].refpoints.xBezier = h5_read(if_mesh, pmegazone_refpoints+"/xBezier",	order='F')
						PMegazones[k].refpoints.yBezier	= h5_read(if_mesh, pmegazone_refpoints+"/yBezier",	order='F')
	
				PMegazones[k].meshchanged = True
				
			mesher_megazones = "/mesher/Megazones"
			nMegazones =len(Megazones)
			for k in range(nMegazones):
				megazone_refpoints = mesher_megazones + "/Megazone{:d}/refpoints".format(k+1)
				Megazones[k].refpoints		= types.SimpleNamespace()		
	
				Megazones[k].refpoints.R	= h5_read(if_mesh, megazone_refpoints+"/R",		order='F')
				Megazones[k].refpoints.Z	= h5_read(if_mesh, megazone_refpoints+"/Z",		order='F')
				Megazones[k].refpoints.psi	= h5_read(if_mesh, megazone_refpoints+"/psi",	order='F')
				Megazones[k].ismeshed	  	= True
	
				Megazones[k].refpoints.nPoints		= h5_read(if_mesh, megazone_refpoints+"/nPoints", 		keep_array=False)							#Mesher parameters
				Megazones[k].refpoints.RefineType	= h5_read(if_mesh, megazone_refpoints+"/RefineType", 	keep_array=False)
				Megazones[k].refpoints.RefineSide	= h5_read(if_mesh, megazone_refpoints+"/RefineSide", 	keep_array=False)
				Megazones[k].refpoints.AdjustMode	= h5_read(if_mesh, megazone_refpoints+"/AdjustMode", 	keep_array=False)
				Megazones[k].refpoints.ParamL		= h5_read(if_mesh, megazone_refpoints+"/ParamL", 		keep_array=False)
				Megazones[k].refpoints.ParamR		= h5_read(if_mesh, megazone_refpoints+"/ParamR", 		keep_array=False)
				if(Megazones[k].refpoints.RefineType == 2):
					Megazones[k].refpoints.xBezier	= h5_read(if_mesh, megazone_refpoints+"/xBezier",	order='F')
					Megazones[k].refpoints.yBezier	= h5_read(if_mesh, megazone_refpoints+"/yBezier",	order='F')
	
			mesher_zones = "/mesher/Zones"
			for k in range(len(Zones)):
				zone = mesher_zones + "/Zone{:d}".format(k+1)
	
				Zones[k].northaligned		= h5_read(if_mesh, zone+"/northaligned",    	keep_array=False)
				Zones[k].southaligned		= h5_read(if_mesh, zone+"/southaligned",    	keep_array=False)
	
				if(Zones[k].northaligned):
					Zones[k].subNorth = [types.SimpleNamespace(), types.SimpleNamespace()]
					for n in range(2):
						sub = "_sub{:d}".format(n+1)
						Zones[k].subNorth[n].R = h5_read(if_mesh, zone+"/North_R"+sub,	order='F')
						Zones[k].subNorth[n].Z = h5_read(if_mesh, zone+"/North_Z"+sub,	order='F')
						Zones[k].subNorth[n].ismeshed = True
						
				if(Zones[k].southaligned):
					Zones[k].subSouth = [types.SimpleNamespace(), types.SimpleNamespace()]
					for n in range(2):
						sub = "_sub{:d}".format(n+1)
						Zones[k].subSouth[n].R = h5_read(if_mesh, zone+"/South_R"+sub,	order='F')
						Zones[k].subSouth[n].Z = h5_read(if_mesh, zone+"/South_R"+sub,	order='F')
						Zones[k].subSouth[n].ismeshed = True
				
				Zones[k].east.ismeshed  	= True
				Zones[k].west.ismeshed  	= True
				Zones[k].north.ismeshed		= True
				Zones[k].south.ismeshed		= True
				Zones[k].orthomeshchanged	= True
	
			Mesher.Segments_OK	= True
			if(DEBUG > 1): print("\tload_soledge_mesh_file: Segments = OK")
	#	else:
		except:
			if(Mesher.Zones_OK):
				for k in range(len(Zones)):
					Zones[k].northaligned	= False
					Zones[k].southaligned	= False
					Zones[k].east.ismeshed  = False
					Zones[k].west.ismeshed  = False
					Zones[k].north.ismeshed = False
					Zones[k].south.ismeshed = False
					
				for k in range(len(Megazones)):
					Megazones[k].ismeshed 			= False
				
				PMegazones = Mesher.PMegazones
				for k in range(len(PMegazones)):
					PMegazones[k].ismeshed 			= False
					PMegazones[k].isaligned 		= False
					PMegazones[k].subrefpoints		= [types.SimpleNamespace(), types.SimpleNamespace()]
					for n in range(2):
						PMegazones[k].subrefpoints[n].R = np.array([])
						PMegazones[k].subrefpoints[n].Z = np.array([])
	
			Mesher.Segments_OK	= False
			if(DEBUG > 1): print("\tload_soledge_mesh_file: Segments = NO")

	if(Mesher.Segments_OK):		
		try:
			mesher_zones = "/mesher/Zones"
			for k in range(len(Zones)):
				zone = mesher_zones + "/Zone{:d}".format(k+1)
				Zones[k].meshortho = h5_read(if_mesh, zone+"/meshortho",    	keep_array=False)

			for k in range(len(Zones)):
				zone = "zone{:d}".format(k+1)
				Zones[k].gridR	= h5_read(if_mesh,zone+"/Rcorner", order='F')				#Soledge mesh corners
				Zones[k].gridZ	= h5_read(if_mesh,zone+"/Zcorner", order='F')
				
			Mesher.MagGrid_OK	= True
			if(DEBUG > 1): print("\tload_soledge_mesh_file: Grid = OK")
		except:
			Mesher.MagGrid_OK	= False
			if(DEBUG > 1): print("\tload_soledge_mesh_file: Grid = NO")


#	Reading base (what is needed to soledge) mesh data	
	
	try:
		Mesher.Rwall 	= h5_read(if_mesh,"/wall/R", order='F')
		Mesher.Zwall 	= h5_read(if_mesh,"/wall/Z", order='F')
		Mesher.wall_OK	= True
	except:
		Mesher.wall_OK	= False

	try:
		Mesher.r2D		= h5_read(if_mesh,'/config/r',		order='F')
		Mesher.z2D		= h5_read(if_mesh,'/config/z',		order='F')
		if(len(Mesher.r2D.shape) == 1): Mesher.r2D, Mesher.z2D = np.meshgrid(Mesher.r2D, Mesher.z2D)

		try:
			Mesher.flux2D	= h5_read(if_mesh,'/config/psi',	order='F')
		except:
			Mesher.flux2D	= h5_read(if_mesh,'/config/Psi',	order='F')

		try:
			Mesher.Br2D		= h5_read(if_mesh,'/config/Br',		order='F')
			Mesher.Bz2D		= h5_read(if_mesh,'/config/Bz',		order='F')
			Mesher.Bphi2D	= h5_read(if_mesh,'/config/Bphi',	order='F')
			Mesher.new_equ_OK	= True
			if(DEBUG > 1): print("\tload_soledge_mesh_file: new_equ = OK")
		except:
			if(DEBUG > 1): print("\tload_soledge_mesh_file: new_equ = NO")
			Mesher.new_equ_OK	= False
			
		Mesher.equ_OK	= True
		if(DEBUG > 1): print("\tload_soledge_mesh_file: equ = OK")
	except:
		Mesher.equ_OK	= False
		if(DEBUG > 1): print("\tload_soledge_mesh_file: equ = NO")
	
	try:
		Mesher.use_penalization = False

		nZones = h5_read(if_mesh,"NZones", keep_array=False)
		for k in range(nZones):
			if(len(Zones) <= k): Zones.append(types.SimpleNamespace())
	
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
			if(Zones[k].Chi.max() > 0.): Mesher.use_penalization = True
	
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
	
			if(not Mesher.MagGrid_OK):
				Zones[k].gridR		= h5_read(if_mesh,zone+"/Rcorner", order='F')				#Soledge mesh corners
				Zones[k].gridZ		= h5_read(if_mesh,zone+"/Zcorner", order='F')
	
		nMegazones = h5_read(if_mesh,"NMegazones", keep_array=False)
		
		for k in range(nMegazones):
			if(len(Megazones) <= k): Megazones.append(types.SimpleNamespace())
	
			Mzone  = "megazone{:d}".format(k+1)
			Megazones[k].list		= h5_read(if_mesh,Mzone+"/configuration", order='F', keep_array=True) - 1	# Python index
			Megazones[k].isperiodic = h5_read(if_mesh,Mzone+"/isperiodic", keep_array=False)


			if(len(Megazones[k].list.shape) > 1): Megazones[k].list = Megazones[k].list[0,:]						#Pro old mesh
			
			for iZone in Megazones[k].list: Zones[iZone].mz = k
			
		Mesher.Mesh_OK			= True
		if(DEBUG > 1): print("\tload_soledge_mesh_file: Mesh = OK")
	except:
		Mesher.Mesh_OK			= False
		if(DEBUG > 1): print("\tload_soledge_mesh_file: Mesh = NO")

	if(Mesher.Mesh_OK):

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

				TranspCuts[i].Flux.Zones = []
				for k in range(len(TranspCuts[i].Flux.nz)):
					TranspCuts[i].Flux.Zones.append(h5_read(if_mesh, cut+"/Flux_Zones{:d}".format(k+1), order='F'))
				
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
			Mesher.transp_prof_OK = True
			if(DEBUG > 1): print("\tload_soledge_mesh_file: transp_prof = OK")
		except:
			TranspCuts = []
			Mesher.transp_prof_OK = False
			if(DEBUG > 1): print("\tload_soledge_mesh_file: transp_prof = NO")
		
		try:
			if(not Mesher.transp_prof_OK):
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
			Mesher.transp_values_OK = True
			if(DEBUG > 1): print("\tload_soledge_mesh_file: transp_values = OK")
		except:
			TranspCuts = []
			Mesher.transp_values_OK = False
			if(DEBUG > 1): print("\tload_soledge_mesh_file: transp_values = NO")
		
#		FeedbackTransp data

		try:
			mesher_feedback		= "/mesher/FeedbackTransp"
	
			FeedbackTransp				= types.SimpleNamespace()
			FeedbackTransp.Data			= types.SimpleNamespace()
			FeedbackTransp.Data.Dmin		= h5_read(if_mesh, mesher_feedback +"/Dmin",		keep_array=False)
			FeedbackTransp.Data.Dmax		= h5_read(if_mesh, mesher_feedback +"/Dmax",		keep_array=False)
			FeedbackTransp.Data.Keep		= h5_read(if_mesh, mesher_feedback +"/Keep",		keep_array=False)
			FeedbackTransp.Data.Gain		= h5_read(if_mesh, mesher_feedback +"/Gain",		keep_array=False)
			FeedbackTransp.Data.GainG		= h5_read(if_mesh, mesher_feedback +"/GainG",		keep_array=False)
			FeedbackTransp.Data.Lambda	= h5_read(if_mesh, mesher_feedback +"/Lambda",		keep_array=False)
			FeedbackTransp.Data.UsePattern= h5_read(if_mesh, mesher_feedback +"/UsePattern",	keep_array=False)
	
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
				FeedbackTransp.Cut.Flux.Profiles[k].Name = h5_read(if_mesh, "/mesher/FeedbackTransp/Prof_Name{:d}".format(k+1), keep_array=False)
							
			for k in range(nFeedbackProfiles):
				FeedbackTransp.Cut.Flux.Profiles[k].xValues	= h5_read(if_mesh, cut+"/Flux_"+ FeedbackTransp.Cut.Flux.Profiles[k].Name+"_xValues", order='F')
				FeedbackTransp.Cut.Flux.Profiles[k].Values	= h5_read(if_mesh, cut+"/Flux_"+ FeedbackTransp.Cut.Flux.Profiles[k].Name+"_Values", order='F')
		
			Mesher.feedback_transp_OK = True
			if(DEBUG > 1): print("\tload_soledge_mesh_file: feedback = OK")
		except:
			FeedbackTransp		   = 0
			Mesher.feedback_transp_OK = False
			if(DEBUG > 1): print("\tload_soledge_mesh_file: feedback = NO")
	else:
		TranspCuts					= []
		FeedbackTransp				= 0
		Mesher.transp_prof_OK		= False
		Mesher.transp_values_OK 	= False
		Mesher.feedback_transp_OK	= False
		
	if_mesh.close()
	
	Mesher.X_points			= X_points
	Mesher.Frontiers		= Frontiers
	Mesher.Zones			= Zones
	Mesher.Megazones		= Megazones
	Mesher.PMegazones		= PMegazones
	Mesher.TranspCuts		= TranspCuts
	Mesher.FeedbackTransp	= FeedbackTransp

	if(Mesher.Mesh_OK): recompute_megazones(Mesher)

	return Mesher
