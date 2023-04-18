
# Function definition is here
import types
import scipy.io
import numpy as np
from matplotlib.path			import Path
from routines.reset_config		import reset_config
from routines.globals			import DEBUG

#============================================================
# This routine reads fields from MatLab files
#============================================================

def load_matlab_mesh_file(mat_file):

	if(DEBUG > 0):	print("load_matlab_mesh_file: reading file ",mat_file)

	try:
		mat = scipy.io.loadmat(mat_file)
	except OSError:
		print("\t\tNot found mat file: "+mat_file)
		return

#	Clear data

	Mesher			= reset_config()
	X_points		= []
	Zones			= []
	PMegazones		= []
	Megazones		= []
	
#	start reading

	if(True):
#		try:
		Mesher.r2			= mat['r2D']
		Mesher.z2			= mat['z2D']
		Mesher.Br2			= mat['Br2D']
		Mesher.Bz2			= mat['Bz2D']
		Mesher.Bphi2		= mat['Bphi2D']
		Mesher.flux2		= mat['flux2D']
		Mesher.OutLength2	= -1.
		Mesher.Smooth2		= 0.
		Mesher.Smooth2D		= 0.
		Mesher.nr			= Mesher.r2.shape[0]
		Mesher.nz			= Mesher.r2.shape[1]
		Mesher.in_equ_OK	= True
		
		Mesher.r2D		= Mesher.r2
		Mesher.z2D		= Mesher.z2
		Mesher.flux2D	= Mesher.flux2
		Mesher.Br2D		= Mesher.Br2
		Mesher.Bz2D		= Mesher.Bz2
		Mesher.Bphi2D	= Mesher.Bphi2
		Mesher.equ_OK		= True
		Mesher.new_equ_OK	= True
		
		if(DEBUG > 1): print("\tread_matalab_save: in_equ = OK")
		if(DEBUG > 1): print("\tread_matalab_save: equ = OK")
	else:
#	except:
		Mesher.in_equ_OK	= False
		if(DEBUG > 1): print("\tread_matalab_save: in_equ = NO")

	if(Mesher.in_equ_OK):
		if(True):
#		try:
			Mesher.extrapol_val	= mat['extrapol_val'][0,0]
			Mesher.raise_value	= 1.e-4

			nX_points = mat['X_points'][0,0]['num'][0,0]
			X_points = []
			for k in range(nX_points):
				X_points.append(types.SimpleNamespace())	
				X_points[-1].R		= mat['X_points'][0,0]['R'][0][k]	
				X_points[-1].Z		= mat['X_points'][0,0]['Z'][0][k]	
				X_points[-1].psi 	= mat['X_points'][0,0]['psi'][0][k]	
				X_points[-1].index	= mat['X_points'][0,0]['index'][0][k] - 1			 								#python index
				X_points[-1].sel		= False
	
				X_points[-1].cut = []
				for n in range(4):
					X_points[-1].cut.append(types.SimpleNamespace())
					X_points[-1].cut[-1].Rs		= mat['X_points'][0,0]['cut'][0,k]['R'][0,n]
					X_points[-1].cut[-1].Zs		= mat['X_points'][0,0]['cut'][0,k]['Z'][0,n]
					X_points[-1].cut[-1].psi		= mat['X_points'][0,0]['cut'][0,k]['psi'][0,n]
					X_points[-1].cut[-1].psilim	= mat['X_points'][0,0]['cut'][0,k]['psilim'][0,n]
					
					X_points[-1].cut[-1].R		= mat['X_points'][0,0]['cut'][0,k]['arc'][0,n]['R'][0]
					X_points[-1].cut[-1].Z		= mat['X_points'][0,0]['cut'][0,k]['arc'][0,n]['Z'][0]
					X_points[-1].cut[-1].type	= mat['X_points'][0,0]['cut'][0,k]['arc'][0,n]['type'][0,0]
	
				X_points[-1].branch = []
				for n in range(4):
					X_points[-1].branch.append(types.SimpleNamespace())
					X_points[-1].branch[-1].R		= mat['X_points'][0,0]['branch'][0,k]['R'][0,n]
					X_points[-1].branch[-1].Z		= mat['X_points'][0,0]['branch'][0,k]['Z'][0,n]
					X_points[-1].branch[-1].theta	= mat['X_points'][0,0]['branch'][0,k]['theta'][0,n]
					X_points[-1].branch[-1].arc		= mat['X_points'][0,0]['branch'][0,k]['arc'][0,n] - 1 				#python index
			
			Mesher.extrapol_OK	= True
			if(DEBUG > 1): print("\tread_matalab_save: extrapol = OK")
		else:
#		except:
			X_points		= []
			Mesher.extrapol_OK	= False
			if(DEBUG > 1): print("\tread_matalab_save: extrapol = NO")

		
	if(Mesher.extrapol_OK):
		if(True):
#		try:
			Mesher.psicore		= mat['psicore'][0,0]
			Mesher.psiout		= mat['psiout'][0,0]
			
			for k in range(nX_points):
				for n in range(4):
					X_points[k].cut[n].psiR	= mat['X_points'][0,0]['cut'][0,k]['arc'][0,n]['psiR'][0]
					X_points[k].cut[n].psiZ	= mat['X_points'][0,0]['cut'][0,k]['arc'][0,n]['psiZ'][0]	
			Mesher.xPoints_OK = True
			if(DEBUG > 1): print("\tread_matalab_save: xPoints = OK")
		else:
#		except:
			Mesher.xPoints_OK	= False
			if(DEBUG > 1): print("\tread_matalab_save: xPoints = NO")

		
	if(Mesher.xPoints_OK):
		if(True):
#		try:
			Mesher.Frontiers = []
			nFrontiers = mat['frontiers'][0,0]['num'][0,0]
			for k in range(nFrontiers):
				Mesher.Frontiers.append(types.SimpleNamespace())
				Mesher.Frontiers[-1].R		 = mat['frontiers'][0,0]['lim'][0,k]['R'][:,0]
				Mesher.Frontiers[-1].Z		 = mat['frontiers'][0,0]['lim'][0,k]['Z'][:,0]
				Mesher.Frontiers[-1].P1	 	 = mat['frontiers'][0,0]['lim'][0,k]['P1']
				Mesher.Frontiers[-1].psiA	 = mat['frontiers'][0,0]['lim'][0,k]['psiA'][0,0]
				Mesher.Frontiers[-1].psiB	 = mat['frontiers'][0,0]['lim'][0,k]['psiB'][0,0]
				Mesher.Frontiers[-1].sel	 = False
				Mesher.Frontiers[-1].psimin	 = np.min(Mesher.Frontiers[-1].P1[:,2])
				Mesher.Frontiers[-1].psimax  = np.max(Mesher.Frontiers[-1].P1[:,2])
				
			Mesher.Frontiers_OK	= True
			if(DEBUG > 1): print("\tread_matalab_save: Frontiers = OK")
		else:
#		except:
			Mesher.Frontiers	= []
			Mesher.Frontiers_OK	= False
			if(DEBUG > 1): print("\tread_matalab_save: Frontiers = NO")


	if(Mesher.Frontiers_OK):
#		try:
		if(True):
			shift_pcoord =  np.array([0,1,1,0])
			nZones = mat['zones'][0,0]['num'][0,0]
			Zones = []
			for k in range(nZones):
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
				
				Zones[k].east.R		= mat['zones'][0,0]['zone'][0,k]['east'][0,0]['R'][0]
				Zones[k].east.Z		= mat['zones'][0,0]['zone'][0,k]['east'][0,0]['Z'][0]
				Zones[k].west.R		= mat['zones'][0,0]['zone'][0,k]['west'][0,0]['R'][0]
				Zones[k].west.Z		= mat['zones'][0,0]['zone'][0,k]['west'][0,0]['Z'][0]
				Zones[k].north.R	= mat['zones'][0,0]['zone'][0,k]['north'][0,0]['R'][0]
				Zones[k].north.Z	= mat['zones'][0,0]['zone'][0,k]['north'][0,0]['Z'][0]
				Zones[k].south.R	= mat['zones'][0,0]['zone'][0,k]['south'][0,0]['R'][0]
				Zones[k].south.Z	= mat['zones'][0,0]['zone'][0,k]['south'][0,0]['Z'][0]
				Zones[k].coord		= mat['zones'][0,0]['zone'][0,k]['coord'][0] - 1 											#python index
				Zones[k].pA.coord	= mat['zones'][0,0]['zone'][0,k]['pA'][0,0]['coord'][0] - shift_pcoord 						#python index
				Zones[k].pB.coord	= mat['zones'][0,0]['zone'][0,k]['pB'][0,0]['coord'][0] - shift_pcoord 						#python index
				Zones[k].pC.coord	= mat['zones'][0,0]['zone'][0,k]['pC'][0,0]['coord'][0] - shift_pcoord 						#python index
				Zones[k].pD.coord	= mat['zones'][0,0]['zone'][0,k]['pD'][0,0]['coord'][0] - shift_pcoord 						#python index
				
				Zones[k].mz					= mat['zones'][0,0]['zone'][0,k]['mz'][0,0] - 1 									#python index
				Zones[k].pmz				= mat['zones'][0,0]['zone'][0,k]['pmz'][0,0] - 1 									#python index
				Zones[k].Xtype_east			= mat['zones'][0,0]['zone'][0,k]['Xtype_east'][0,0]
				Zones[k].Xtype_west			= mat['zones'][0,0]['zone'][0,k]['Xtype_west'][0,0]	
				
				Zones[k].Neighbour.east		= mat['zones'][0,0]['zone'][0,k]['Neighbour'][0,0]['east'][0,0] - 1 				#python index
				Zones[k].Neighbour.west		= mat['zones'][0,0]['zone'][0,k]['Neighbour'][0,0]['west'][0,0] - 1 				#python index
				Zones[k].Neighbour.north	= mat['zones'][0,0]['zone'][0,k]['Neighbour'][0,0]['north'][0,0] - 1 				#python index
				Zones[k].Neighbour.south	= mat['zones'][0,0]['zone'][0,k]['Neighbour'][0,0]['south'][0,0] - 1 				#python ind		
				
			nPMegazones = mat['Pmegazone'][0,0]['num'][0,0]
			PMegazones	= []
			for k in range(nPMegazones):
				PMegazones.append(types.SimpleNamespace())	
				PMegazones[k].list		= mat['Pmegazone'][0,0]['mz'][0,k]['list'][0] - 1 				#python index
				
				
			nMegazones = mat['megazone'][0,0]['num'][0,0]
			Megazones  = []
			for k in range(nMegazones):
				Megazones.append(types.SimpleNamespace())	
				Megazones[k].list		= mat['megazone'][0,0]['mz'][0,k]['list'][0] - 1 				#python index
				Megazones[k].isperiodic	= mat['megazone'][0,0]['mz'][0,k]['isperiodic'][0,0]
	
			Mesher.MagZones_OK	= True
			if(DEBUG): print("\tread_matalab_save: Limits = OK")
		else:
#		except:
			Zones				= []
			PMegazones			= []
			Megazones			= []
			Mesher.MagZones_OK	= False
			if(DEBUG > 1): print("\tread_matalab_save: Limits = NO")
		
#	Try if segments have been defined

	if(Mesher.MagZones_OK):
		if(True):
#		try:
			nPMegazones =len(PMegazones)
			for k in range(nPMegazones):
				PMegazones[k].isaligned = mat['Pmegazone'][0,0]['mz'][0,k]['isaligned'][0,0]
	
				PMegazones[k].refpoints = types.SimpleNamespace()		
				if(PMegazones[k].isaligned):
					PMegazones[k].subrefpoints		= [types.SimpleNamespace(), types.SimpleNamespace()]
					PMegazones[k].subrefpoints[0].R = mat['Pmegazone'][0,0]['mz'][0,k]['subrefpoints'][0,0]['R'][0]
					PMegazones[k].subrefpoints[0].Z = mat['Pmegazone'][0,0]['mz'][0,k]['subrefpoints'][0,0]['Z'][0]
					PMegazones[k].subrefpoints[1].R = mat['Pmegazone'][0,0]['mz'][0,k]['subrefpoints'][0,1]['R'][0]
					PMegazones[k].subrefpoints[1].Z = mat['Pmegazone'][0,0]['mz'][0,k]['subrefpoints'][0,1]['Z'][0]
					PMegazones[k].align_psimin		= mat['Pmegazone'][0,0]['mz'][0,k]['align_psimin'][0,0]
					PMegazones[k].align_psimax		= mat['Pmegazone'][0,0]['mz'][0,k]['align_psimax'][0,0]
		#				
				PMegazones[k].refpoints.R	= mat['Pmegazone'][0,0]['mz'][0,k]['refpoints'][0,0]['R'][0]
				PMegazones[k].refpoints.Z	= mat['Pmegazone'][0,0]['mz'][0,k]['refpoints'][0,0]['Z'][0]
				PMegazones[k].refpoints.nz	= mat['Pmegazone'][0,0]['mz'][0,k]['refpoints'][0,0]['nz'][0,0] - 1 				#python index
				PMegazones[k].refpoints.nzB = mat['Pmegazone'][0,0]['mz'][0,k]['refpoints'][0,0]['nzB'][0,0] - 1 				#python index
				
				PMegazones[k].ismeshed	  = True
	
				if(PMegazones[k].isaligned):																					#Mesher parameters in python grid-gen
					for n in range(2):
						PMegazones[k].subrefpoints[n].nPoints		= 0
						PMegazones[k].subrefpoints[n].RefineType	= 0
						PMegazones[k].subrefpoints[n].RefineSide	= 0
						PMegazones[k].subrefpoints[n].AdjustMode	= 0
						PMegazones[k].subrefpoints[n].ParamL		= 0
						PMegazones[k].subrefpoints[n].ParamR		= 0
				else:
					PMegazones[k].refpoints.nPoints		= 0
					PMegazones[k].refpoints.RefineType	= 0
					PMegazones[k].refpoints.RefineSide	= 0
					PMegazones[k].refpoints.AdjustMode	= 0
					PMegazones[k].refpoints.ParamL		= 0
					PMegazones[k].refpoints.ParamR		= 0
	
				PMegazones[k].meshchanged = True
				
			nMegazones =len(Megazones)
			for k in range(nMegazones):
				Megazones[k].refpoints		= types.SimpleNamespace()		
	
				Megazones[k].refpoints.R	= mat['megazone'][0,0]['mz'][0,k]['refpoints'][0,0]['R'][0]
				Megazones[k].refpoints.Z	= mat['megazone'][0,0]['mz'][0,k]['refpoints'][0,0]['Z'][0]
				Megazones[k].refpoints.psi	= mat['megazone'][0,0]['mz'][0,k]['refpoints'][0,0]['psi'][0]
				Megazones[k].ismeshed	  	= True
	
				Megazones[k].refpoints.nPoints		= 0												#Mesher parameters in python grid-gen
				Megazones[k].refpoints.RefineType	= 0
				Megazones[k].refpoints.RefineSide	= 0
				Megazones[k].refpoints.AdjustMode	= 0
				Megazones[k].refpoints.ParamL		= 0
				Megazones[k].refpoints.ParamR		= 0
	
			for k in range(len(Zones)):
				Zones[k].northaligned		= mat['zones'][0,0]['zone'][0,k]['northaligned'][0,0]
				Zones[k].southaligned		= mat['zones'][0,0]['zone'][0,k]['southaligned'][0,0]
	
				if(Zones[k].northaligned):
					Zones[k].subNorth = [types.SimpleNamespace(), types.SimpleNamespace()]
					for n in range(2):
						Zones[k].subNorth[n].R = mat['zones'][0,0]['zone'][0,k]['subNorth'][0,n]['R'][0]
						Zones[k].subNorth[n].Z = mat['zones'][0,0]['zone'][0,k]['subNorth'][0,n]['Z'][0]
						Zones[k].subNorth[n].ismeshed = True
						
				if(Zones[k].southaligned):
					Zones[k].subSouth = [types.SimpleNamespace(), types.SimpleNamespace()]
					for n in range(2):
						Zones[k].subSouth[n].R = mat['zones'][0,0]['zone'][0,k]['subSouth'][0,n]['R'][0]
						Zones[k].subSouth[n].Z = mat['zones'][0,0]['zone'][0,k]['subSouth'][0,n]['Z'][0]
						Zones[k].subSouth[n].ismeshed = True
				
				Zones[k].east.ismeshed  	= True
				Zones[k].west.ismeshed  	= True
				Zones[k].north.ismeshed		= True
				Zones[k].south.ismeshed		= True
				Zones[k].orthomeshchanged	= True
	
			Mesher.Segments_OK	= True
			if(DEBUG > 1): print("\tread_matalab_save: Segments = OK")
		else:
#		except:
			if(Mesher.MagZones_OK):
				for k in range(len(Zones)):
					Zones[k].northaligned	= False
					Zones[k].southaligned	= False
					Zones[k].east.ismeshed  = False
					Zones[k].west.ismeshed  = False
					Zones[k].north.ismeshed = False
					Zones[k].south.ismeshed = False
					
				for k in range(len(Megazones)):
					Megazones[k].ismeshed 			= False
				
				for k in range(len(PMegazones)):
					PMegazones[k].ismeshed 			= False
					PMegazones[k].isaligned 		= False
					PMegazones[k].subrefpoints		= [types.SimpleNamespace(), types.SimpleNamespace()]
					for n in range(2):
						PMegazones[k].subrefpoints[n].R = np.array([])
						PMegazones[k].subrefpoints[n].Z = np.array([])
	
			Mesher.Segments_OK	= False
			if(DEBUG > 1): print("\tread_matalab_save: Segments = NO")

	if(Mesher.Segments_OK):		
		if(True):
#		try:
			for k in range(len(Zones)):
				Zones[k].meshortho = mat['zones'][0,0]['zone'][0,k]['meshortho'][0,0]

			for k in range(len(Zones)):
				Zones[k].gridR	= mat['zones'][0,0]['zone'][0,k]['gridR']
				Zones[k].gridZ	= mat['zones'][0,0]['zone'][0,k]['gridZ']
				
			Mesher.MagGrid_OK	= True
			if(DEBUG): print("\tread_matalab_save: Grid = OK")
		else:
#		except:
			Mesher.MagGrid_OK	= False
			if(DEBUG > 1): print("\tread_matalab_save: Grid = NO")
		
#	Ballooning data
	
	TranspCuts			  = []
	Mesher.transp_prof_OK = False
	if(DEBUG > 1): print("\tread_matalab_save: transp_prof = NO")

	Mesher.TranspCuts = TranspCuts

#	Reading base (what is needed to soledge) mesh data	
	
	try:
		if(mat['Rwall'].shape[0] == 1):
			Mesher.Rwall 	= mat['Rwall'][0,:]
			Mesher.Zwall 	= mat['Zwall'][0,:]
		else:
			Mesher.Rwall 	= mat['Rwall'][:,0]
			Mesher.Zwall 	= mat['Zwall'][:,0]

		Mesher.wall_OK	= True
	except:
		Mesher.wall_OK	= False
	
#	try:
	if(True):

		for k in range(len(Zones)):
			if(len(Zones) <= k): Zones.append(types.SimpleNamespace())
	
			Zones[k].Nx = mat['zones'][0,0]['zone'][0,k]['Nx'][0,0]
			Zones[k].Nz = mat['zones'][0,0]['zone'][0,k]['Nz'][0,0]

			Zones[k].x	= mat['zones'][0,0]['zone'][0,k]['x']
			Zones[k].z	= mat['zones'][0,0]['zone'][0,k]['z']
		
			Zones[k].xb = mat['zones'][0,0]['zone'][0,k]['xb'] 
			Zones[k].zb	= mat['zones'][0,0]['zone'][0,k]['zb'] 
			

			
			"""	
			Zones[k].Neighbour		 = types.SimpleNamespace()
			Zones[k].Neighbour.east		= mat['zones'][0,0]['zone'][0,k]['Neighbour'][0,0]['east'][0,0] - 1 				#python index
			Zones[k].Neighbour.west		= mat['zones'][0,0]['zone'][0,k]['Neighbour'][0,0]['west'][0,0] - 1 				#python index
			Zones[k].Neighbour.north	= mat['zones'][0,0]['zone'][0,k]['Neighbour'][0,0]['north'][0,0] - 1 				#python index
			Zones[k].Neighbour.south	= mat['zones'][0,0]['zone'][0,k]['Neighbour'][0,0]['south'][0,0] - 1 				#python ind		
			"""
						
			Zones[k].MagNeighbour		= types.SimpleNamespace()
			Zones[k].MagNeighbour.north = mat['zones'][0,0]['zone'][0,k]['MagNeighbour'][0,0]['north'][0,0]
			Zones[k].MagNeighbour.south = mat['zones'][0,0]['zone'][0,k]['MagNeighbour'][0,0]['south'][0,0]
			Zones[k].MagNeighbour.east  = mat['zones'][0,0]['zone'][0,k]['MagNeighbour'][0,0]['east'][0,0]
			Zones[k].MagNeighbour.west  = mat['zones'][0,0]['zone'][0,k]['MagNeighbour'][0,0]['west'][0,0]
	
			Zones[k].SouthP			= types.SimpleNamespace()
			Zones[k].NorthP			= types.SimpleNamespace()
			Zones[k].WestP			= types.SimpleNamespace()
			Zones[k].EastP			= types.SimpleNamespace()
		
			Zones[k].Br				= mat['zones'][0,0]['zone'][0,k]['Br']
			Zones[k].SouthP.Br		= mat['zones'][0,0]['zone'][0,k]['SouthP'][0,0]['Br'][0] 
			Zones[k].NorthP.Br		= mat['zones'][0,0]['zone'][0,k]['NorthP'][0,0]['Br'][0] 
			Zones[k].WestP.Br		= mat['zones'][0,0]['zone'][0,k]['WestP'][0,0]['Br'][0] 
			Zones[k].EastP.Br		= mat['zones'][0,0]['zone'][0,k]['EastP'][0,0]['Br'][0] 
	
			Zones[k].Bz				= mat['zones'][0,0]['zone'][0,k]['Bz'] 
			Zones[k].SouthP.Bz		= mat['zones'][0,0]['zone'][0,k]['SouthP'][0,0]['Bz'][0] 
			Zones[k].NorthP.Bz		= mat['zones'][0,0]['zone'][0,k]['NorthP'][0,0]['Bz'][0] 
			Zones[k].WestP.Bz		= mat['zones'][0,0]['zone'][0,k]['WestP'][0,0]['Bz'][0] 
			Zones[k].EastP.Bz		= mat['zones'][0,0]['zone'][0,k]['EastP'][0,0]['Bz'][0] 
	
			Zones[k].Bphi			= mat['zones'][0,0]['zone'][0,k]['Bphi'] 
			Zones[k].SouthP.Bphi	= mat['zones'][0,0]['zone'][0,k]['SouthP'][0,0]['Bphi'][0] 
			Zones[k].NorthP.Bphi	= mat['zones'][0,0]['zone'][0,k]['NorthP'][0,0]['Bphi'][0] 
			Zones[k].WestP.Bphi		= mat['zones'][0,0]['zone'][0,k]['WestP'][0,0]['Bphi'][0] 
			Zones[k].EastP.Bphi		= mat['zones'][0,0]['zone'][0,k]['EastP'][0,0]['Bphi'][0] 
	
			Zones[k].gridRc			= mat['zones'][0,0]['zone'][0,k]['gridRc']
			Zones[k].SouthP.R		= mat['zones'][0,0]['zone'][0,k]['SouthP'][0,0]['R'][0] 
			Zones[k].NorthP.R		= mat['zones'][0,0]['zone'][0,k]['NorthP'][0,0]['R'][0] 
			Zones[k].WestP.R		= mat['zones'][0,0]['zone'][0,k]['WestP'][0,0]['R'][0] 
			Zones[k].EastP.R		= mat['zones'][0,0]['zone'][0,k]['EastP'][0,0]['R'][0] 
	
			Zones[k].gridZc			= mat['zones'][0,0]['zone'][0,k]['gridZc']
			Zones[k].SouthP.Z		= mat['zones'][0,0]['zone'][0,k]['SouthP'][0,0]['Z'][0]
			Zones[k].NorthP.Z		= mat['zones'][0,0]['zone'][0,k]['NorthP'][0,0]['Z'][0]
			Zones[k].WestP.Z		= mat['zones'][0,0]['zone'][0,k]['WestP'][0,0]['Z'][0] 
			Zones[k].EastP.Z		= mat['zones'][0,0]['zone'][0,k]['EastP'][0,0]['Z'][0]
	
		for k in range(nMegazones):	
			for iZone in Megazones[k].list: Zones[iZone].mz = k

		Mesher.use_penalization = False
		for k in range(len(Zones)):
			Zones[k].Chi = mat['zones'][0,0]['zone'][0,k]['chi']	
			if(Zones[k].Chi.max() > 0.): Mesher.use_penalization = True

		"""
		if(Mesher.use_penalization):
			nWall = Mesher.Rwall.shape[0]
			if((Mesher.Rwall[0] == Mesher.Rwall[-1]) and (Mesher.Zwall[0] == Mesher.Zwall[-1])):
				nWallPath = nWall
			else:
				nWallPath = nWall + 1
		
			WallVerts = np.empty((nWallPath,2),dtype='f8')
			WallCodes = np.empty(nWallPath,dtype='i')
		
			if((Mesher.Rwall[0] == Mesher.Rwall[-1]) and (Mesher.Zwall[0] == Mesher.Zwall[-1])):
				WallVerts[:,0] = Mesher.Rwall
				WallVerts[:,1] = Mesher.Zwall
			else:
				WallVerts[:-1,0] = Mesher.Rwall
				WallVerts[:-1,1] = Mesher.Zwall
				WallVerts[-1,:]  = WallVerts[0,:]
		
			WallCodes[1:-1] = Path.LINETO
			WallCodes[0]	= Path.MOVETO
			WallCodes[-1]	= Path.CLOSEPOLY
			Mesher.WallPath = Path(WallVerts, WallCodes)
			for k in range(len(Zones)):
	
				nTot	 	 = Zones[k].gridRc.size
				RZ1		 	 = np.empty((nTot,2), dtype='f8')
				RZ1[:,0]	 = Zones[k].gridRc.reshape(nTot)
				RZ1[:,1]	 = Zones[k].gridZc.reshape(nTot)
				
				IsIn		 = Mesher.WallPath.contains_points(RZ1)								#find point iside wall contour
				Zones[k].Chi = np.where(IsIn, np.zeros(nTot, dtype='f8'), np.ones(nTot, dtype='f8'))
				Zones[k].Chi = Zones[k].Chi.reshape(Zones[k].gridRc.shape)
			"""
					
		Mesher.Mesh_OK			= True
		if(DEBUG > 1): print("\tread_matalab_save: Mesh = OK")
#	except:
#		Mesher.Mesh_OK			= False
#		if(DEBUG): print("\tread_matalab_save: Mesh = NO")

	Mesher.transp_values_OK = False
	if(DEBUG > 1): print("\tread_matalab_save: transp_values = NO")

	"""
	Wall	   	  = types.SimpleNamespace()
	Wall.Material = []
	Wall.Pump	  = []
	Wall.Puff	  = []
	try:

		Eirene.R	= mat['eirene'][0,0]['R'][0]
		Eirene.Z	= mat['eirene'][0,0]['Z'][0]
		Direct		= mat['eirene'][0,0]['direct'][0,0]

		Rknots		= mat['eirene'][0,0]['Rknots'][0]
		Zknots		= mat['eirene'][0,0]['Zknots'][0]
		
		nTriangles	= mat['eirene'][0,0]['ntriangles'][0,0]
		Triangles = np.empty((nTriangles), dtype=[('k','i4'), ('i','i4'), ('j','i4'),('p1','i4'), ('p2','i4'), ('p3','i4'), 
													('neigh1','i4'), ('typeneigh1','i4'), ('BC1','i4'),
												('neigh2','i4'), ('typeneigh2','i4'), ('BC2','i4'),
												('neigh3','i4'), ('typeneigh3','i4'), ('BC3','i4')])

		Triangles.k		= mat['eirene'][0,0]['triangles'][0,0]['k'][0] - 1
		Triangles.i		= mat['eirene'][0,0]['triangles'][0,0]['i'][0] - 1
		Triangles.j		= mat['eirene'][0,0]['triangles'][0,0]['j'][0] - 1 
		Triangles.p1	= mat['eirene'][0,0]['triangles'][0,0]['p1'][0] - 1
		Triangles.p2	= mat['eirene'][0,0]['triangles'][0,0]['p2'][0] - 1
		Triangles.p3	= mat['eirene'][0,0]['triangles'][0,0]['p3'][0] - 1
		
		Eirene.nknots_new	= mat['eirene'][0,0]['nknots_new'][0]
		Eirene.tri_knots	= mat['eirene'][0,0]['tri_knots']

		nPump	= mat['eirene'][0,0]['numpump'][0,0]
		nPuffs	= mat['eirene'][0,0]['npuff'][0,0]

			wall: [1x1 struct]
			triR: [19870x3 double]
			triZ: [19870x3 double]
			knots_interp: [1x10244 struct]
			triangle: [1x460 struct]
	 	nummaterial: 1
		 numpump: 1
		material: [1x1 struct]
			pump: [1x1 struct]
			puff: [1x1 struct]	
			"""
	Mesher.X_points		= X_points
	Mesher.Zones		= Zones
	Mesher.Megazones	= Megazones
	Mesher.PMegazones	= PMegazones
	
	Eirene				= types.SimpleNamespace()
	
	return Mesher, Eirene
		
	