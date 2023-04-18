from tkinter 					import messagebox
from tkinter.messagebox 		import Message

import types
import numpy					as np
from routines.intersect_contour	import intersect_contour
from routines.part_contour		import part_contour
from routines.globals			import EXTERNAL_PLASMA_WALL, INTERNAL_PLASMA_WALL

def subdivide_selected_segment(Config, MagZones, zone_sel, side_sel): 

	if((side_sel != 0) and (side_sel == 1)):
		print("\tERROR: subdivide_selected_segment called with invalid side_sel=",side_sel)
		print("\t\t\tCHECK CODE!!!!")
		exit()

	success = 0
#	Arc of selected segment  north/south side (along flux surface, out/in)

	c1		= types.SimpleNamespace()
	c1.arc	= [types.SimpleNamespace()]
	c2		= types.SimpleNamespace()
	c2.arc	= [types.SimpleNamespace()]
	cutWall	= types.SimpleNamespace()

	if(side_sel == 0):								# north side 
		R		= MagZones[zone_sel].north.R		# Arc of selected segment
		Z		= MagZones[zone_sel].north.Z
	else:											# south side 
		R		= MagZones[zone_sel].south.R		# Arc of selected segment
		Z		= MagZones[zone_sel].south.Z

	c1.arc[0].x = R
	c1.arc[0].y = Z

#	Loop on walls to find intersection points
#	We also check whether end points of the selected segment are in the plasma or not
	cutWall.x		= np.empty(0, dtype='f8')
	cutWall.y 		= np.empty(0, dtype='f8')
	cutWall.iwall	= np.empty(0, dtype='i4')
	endsInPlasma	= [1,1]
	for iPWall in range(len(Config.iPwalls)):
		iWall = Config.iPwalls[iPWall]
		Wall  = Config.Walls[iWall]

		c2.arc[0].x = Wall.Rwall
		c2.arc[0].y = Wall.Zwall
		X			= intersect_contour(c1,c2)
		if(len(X) >= 1):
			cutWall.x 		= np.append(cutWall.x, np.array([X[k].x for k in range(len(X))]))
			cutWall.y		= np.append(cutWall.y, np.array([X[k].y for k in range(len(X))]))
			cutWall.iwall	= np.append(cutWall.iwall, iWall*np.ones(len(X), dtype='i4'))

#		First point
		InWall = Wall.WallPath.contains_point([R[0],Z[0]])
		if(((Wall.Type == EXTERNAL_PLASMA_WALL) and (not InWall)) or 
		   ((Wall.Type == INTERNAL_PLASMA_WALL) and      InWall)):	endsInPlasma[0] = 0

#		Last point
		InWall = Wall.WallPath.contains_point([R[-1],Z[-1]])
		if(((Wall.Type == EXTERNAL_PLASMA_WALL) and (not InWall)) or 
		   ((Wall.Type == INTERNAL_PLASMA_WALL) and      InWall)):	endsInPlasma[1] = 0


#	Define a reference extremety for the segment
	FindRef = False
	if(len(cutWall.x) > 0):											# There is something to do only if this segment instersects a wall
		if  ((endsInPlasma[0] == 1) and (endsInPlasma[1] == 0)):	# The first point is used as the reference for the distance
			xref = R[0]
			yref = Z[0]
			FindRef = True
		elif((endsInPlasma[0] == 0) and (endsInPlasma[1] == 1)):	# The last point is used as the reference for the distance
			xref = R[-1]
			yref = Z[-1]
			FindRef = True
		elif((endsInPlasma[0] == 1) and (endsInPlasma[1] == 1)):	# Both ends are in the plasma
#			SHITTY CASE WHERE WE DO NOT KNOW WHAT TO DO FOR THE MOMENT (cut in 3?)
			FindRef = False
		else:														# Both ends are in outside the plasma => private SOL => no need to align
			FindRef = False
		
#	Select the relevant intersection (ie the one the closest to the reference extremety)

	if(FindRef):
		d = np.sqrt((cutWall.x - xref)**2 + (cutWall.y-yref)**2)	# distance of intersection points to the reference extremety
		a = np.argmin(d)

		X	= [types.SimpleNamespace()]
		X[0].x = cutWall.x[a]
		X[0].y = cutWall.y[a]
		iAlignWall = cutWall.iwall[a]

#	Cuts the segment at the selected place

		cin			= types.SimpleNamespace()
		cin.x		= R
		cin.y		= Z
		p1			= types.SimpleNamespace()
		p1.x		= X[0].x
		p1.y		= X[0].y
		cout		= part_contour(cin,p1)
		subR1		= cout.arc[0].x
		subZ1		= cout.arc[0].y
		subR2		= cout.arc[1].x
		subZ2		= cout.arc[1].y

		if(side_sel == 0):				# dealing with a north side
			if(not hasattr(MagZones[zone_sel],'subNorth')): MagZones[zone_sel].subNorth = [types.SimpleNamespace(), types.SimpleNamespace()]

			MagZones[zone_sel].iAlignWall			= iAlignWall
			MagZones[zone_sel].subNorth[0].R 		= subR1
			MagZones[zone_sel].subNorth[0].Z		= subZ1
			MagZones[zone_sel].subNorth[1].R		= subR2
			MagZones[zone_sel].subNorth[1].Z		= subZ2
			MagZones[zone_sel].subNorth[0].ismeshed	= False
			MagZones[zone_sel].subNorth[1].ismeshed	= False
			if(MagZones[zone_sel].Neighbour.north > -1):
				nz = MagZones[zone_sel].Neighbour.north
				if(not hasattr(MagZones[nz],'subSouth')):		MagZones[nz].subSouth			= [types.SimpleNamespace(), types.SimpleNamespace()]

				MagZones[nz].iAlignWall				= iAlignWall
				MagZones[nz].subSouth[0].R			= subR1
				MagZones[nz].subSouth[0].Z			= subZ1
				MagZones[nz].subSouth[1].R			= subR2
				MagZones[nz].subSouth[1].Z			= subZ2
				MagZones[nz].subSouth[0].ismeshed	= False
				MagZones[nz].subSouth[1].ismeshed	= False

		if(side_sel == 1):				# dealing with a south side
			if(not hasattr(MagZones[zone_sel],'subSouth')):	MagZones[zone_sel].subSouth				= [types.SimpleNamespace(), types.SimpleNamespace()]
			MagZones[zone_sel].subSouth[0].ismeshed	= False
			MagZones[zone_sel].subSouth[0].R			= subR1
			MagZones[zone_sel].subSouth[0].Z			= subZ1
			MagZones[zone_sel].subSouth[1].R			= subR2
			MagZones[zone_sel].subSouth[1].Z			= subZ2
			MagZones[zone_sel].iAlignWall				= iAlignWall
			MagZones[zone_sel].subSouth[1].ismeshed	= False
			if(MagZones[zone_sel].Neighbour.south > -1):
				nz = MagZones[zone_sel].Neighbour.south
				if(not hasattr(MagZones[nz],'subNorth')):		MagZones[nz].subNorth				= [types.SimpleNamespace(), types.SimpleNamespace()]
				MagZones[nz].iAlignWall				= iAlignWall
				MagZones[nz].subNorth[0].R			= subR1
				MagZones[nz].subNorth[0].Z			= subZ1
				MagZones[nz].subNorth[1].R			= subR2
				MagZones[nz].subNorth[1].Z			= subZ2
				MagZones[nz].subNorth[0].ismeshed	= False
				MagZones[nz].subNorth[1].ismeshed	= False
		success = 1
	else:
		success = 0
		messagebox.showwarning("Subdivide Segmeent", "Selected segment cannot be aligned with wall")

	return success