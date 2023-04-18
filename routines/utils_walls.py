import types
import numpy as np
from matplotlib.path				import Path
from math 							import sqrt

from routines.globals				import *
from routines.get_in_out			import get_in_out
from routines.intersect_contour		import intersect_contour
from routines.find_closest_segment	import find_closest_segment


def walls_define_path(Config):

#	define boundary path

	Config.iPwalls	 = []
	Config.iEwalls	 = []
	Config.iExtWalls  = [0]
	Config.iIntEWalls = []
	for iWall in range(len(Config.Walls)):
		Wall = Config.Walls[iWall]
		CloseWall = True
		if(Wall.Type == EXTERNAL_PLASMA_WALL): 
			Config.iExtWalls[0] = iWall
			Config.iPwalls.append(iWall)
		elif(Wall.Type == INTERNAL_PLASMA_WALL): Config.iPwalls.append(iWall)
		elif(Wall.Type == EXTERNAL_EIRENE_WALL): Config.iExtWalls.append(iWall)
		elif(Wall.Type == INTERNAL_EIRENE_WALL): Config.iIntEWalls.append(iWall)
		else: 									 CloseWall = False

		if(CloseWall):
			Config.Walls[iWall].Closed = True
			if((Wall.Rwall[0] != Wall.Rwall[-1]) or (Wall.Zwall[0] != Wall.Zwall[-1])):					#If needed close wall
				Wall.Rwall = np.append(Wall.Rwall, Wall.Rwall[0])
				Wall.Zwall = np.append(Wall.Zwall, Wall.Zwall[0])

			Area = np.sum((Wall.Zwall[1:]+Wall.Zwall[:-1])*(Wall.Rwall[1:]-Wall.Rwall[0:-1]))
			if(Area > 0.):																				#if needed order counter clockwise
				Wall.Rwall = Wall.Rwall[::-1]
				Wall.Zwall = Wall.Zwall[::-1]

			Wall.Clockwise = -1
			Wall.WallPath = Path(np.array([Wall.Rwall ,Wall.Zwall]).T, closed=True)						#define path
		else:
			if((Wall.Rwall[0] == Wall.Rwall[-1]) and (Wall.Zwall[0] == Wall.Zwall[-1])): Wall.Closed = True
			else:																		 Wall.Closed = False


def wall_define_eirene_walls(Config):
	EWalls = []
	for iWall in range(len(Config.Walls)):
		CWall = Config.Walls[iWall]

		if(CWall.Type == EXTERNAL_EIRENE_WALL):
			EWalls.append(types.SimpleNamespace())
			EWalls[-1].iWall = iWall
			EWalls[-1].Wall  = CWall
			EWalls[-1].ESubZones 	 		 = [types.SimpleNamespace()]
			EWalls[-1].ESubZones[0].dr_knots = (CWall.Rwall.max()-CWall.Rwall.min())/20
			EWalls[-1].ESubZones[0].dz_knots = (CWall.Zwall.max()-CWall.Zwall.min())/20
			EWalls[-1].ESubZones[0].Rbord	 = np.copy(CWall.Rwall) 
			EWalls[-1].ESubZones[0].Zbord	 = np.copy(CWall.Zwall)
	return EWalls

def wall_fix_eirene_walls(Config, EWalls):
	if(len(EWalls) > 0):
		for i in range(len(EWalls)):
			iWall = EWalls[i].iWall
			EWalls[i].Wall = Config.Walls[iWall]
			if(len(EWalls[i].ESubZones) == 1):
				EWalls[i].ESubZones[0].Rbord = np.copy(Config.Walls[iWall].Rwall)
				EWalls[i].ESubZones[0].Zbord = np.copy(Config.Walls[iWall].Zwall)
	elif(len(Config.Walls) > 1):
		EWalls =  wall_define_eirene_walls(Config)

	return EWalls

def walls_copy(Walls):

	NewWalls = []
	for k in range(len(Walls)):
		NewWalls.append(wall_copy(Walls[k]))

	return NewWalls

def wall_copy(Wall):
	NewWall				= types.SimpleNamespace()
	NewWall.Rwall 		= np.copy(Wall.Rwall)													#Wall contour
	NewWall.Zwall 		= np.copy(Wall.Zwall)
	NewWall.Type		= Wall.Type
	NewWall.LineType	= np.copy(Wall.LineType)
	NewWall.Changed		= Wall.Changed
	NewWall.Closed		= Wall.Closed

	return NewWall


def walls_esubzones_copy(ESubZones):
	NewESubZones = []
	for k in range(len(ESubZones)):
		NewESubZones.append(wall_esubzone_copy(ESubZones[k]))

	return NewESubZones

def wall_esubzone_copy(ESubZone):
	NewESubZone			 = types.SimpleNamespace()
	NewESubZone.Rbord	 = np.copy(ESubZone.Rbord) 
	NewESubZone.Zbord	 = np.copy(ESubZone.Zbord)
	NewESubZone.dr_knots = ESubZone.dr_knots
	NewESubZone.dz_knots = ESubZone.dz_knots

	return NewESubZone


def plot2d_walls(Ax, Walls, plasma_wall=True, eirene_wall=True, extra_wall=True, marker=None, end_marker=False):

	for iWall in range(len(Walls)):
		Wall = Walls[iWall]

		LineColor = LINE_COLORS[Wall.LineType[0] % len(LINE_COLORS)]
		LineStyle = LINE_STYLES[Wall.LineType[1] % len(LINE_STYLES)]
		LineWidth = eval(LINE_WIDTH[Wall.LineType[2] % len(LINE_WIDTH)])
		if  ((Wall.Type == INTERNAL_PLASMA_WALL) or (Wall.Type == EXTERNAL_PLASMA_WALL)): Plot_Wall = plasma_wall
		elif((Wall.Type == INTERNAL_EIRENE_WALL) or (Wall.Type == EXTERNAL_EIRENE_WALL)): Plot_Wall = eirene_wall
		else:																			  Plot_Wall = extra_wall

		if(Plot_Wall):
			if(hasattr(Wall,'ESubZones')):
				for iSubEZone in range(len(Wall.ESubZones)):
					SubZone = Wall.ESubZones[iSubEZone]
					if(marker != None):	Ax.plot(SubZone.Rbord, SubZone.Zbord, color="gray", linestyle=LineStyle, linewidth=LineWidth, marker=marker)
					else:				Ax.plot(SubZone.Rbord, SubZone.Zbord, color="gray", linestyle=LineStyle, linewidth=LineWidth)

					if(Wall.Type == EXTERNAL_EIRENE_WALL): Ax.text(np.mean(SubZone.Rbord), np.mean(SubZone.Zbord), "({:d},{:d})".format(iWall+1,iSubEZone+1), color='m', ha="center",fontsize=10)
						
			if(marker != None):	Ax.plot(Wall.Rwall, Wall.Zwall, color=LineColor, linestyle=LineStyle, linewidth=LineWidth, marker=marker)
			else:				Ax.plot(Wall.Rwall, Wall.Zwall, color=LineColor, linestyle=LineStyle, linewidth=LineWidth)
			if(end_marker):
				Ax.plot(Wall.Rwall[0],  Wall.Zwall[0], "bo")
				Ax.plot(Wall.Rwall[-1], Wall.Zwall[-1], "bo")
	return  


#	This routine provides array for plot_and_ask interface
#	==============================

def set_XYlines_walls(Config, wall_line="k-"):

	Xarrs = []
	Yarrs = []
	lines = []
	for iWall in Config.iPwalls:
		Xarrs.append(Config.Walls[iWall].Rwall)
		Yarrs.append(Config.Walls[iWall].Zwall)
		lines.append(wall_line)

	return Xarrs, Yarrs, lines


#	This routine provides min & max of walls
#	========================

def get_min_max_walls(Config, plasma_wall=True, eirene_wall=False, extra_wall=False):

	Rwalls = []
	Zwalls = []
	for iWall in range(len(Config.Walls)):
		Wall = Config.Walls[iWall]

		if  ((((Wall.Type == INTERNAL_PLASMA_WALL) or (Wall.Type == EXTERNAL_PLASMA_WALL)) and plasma_wall) or \
		     (((Wall.Type == INTERNAL_EIRENE_WALL) or (Wall.Type == EXTERNAL_EIRENE_WALL)) and eirene_wall) or \
			  ((Wall.Type == DRAWING_WALL) and extra_wall)):
			Rwalls.append(Wall.Rwall)
			Zwalls.append(Wall.Zwall)

	r_min = np.min(np.array([np.min(Rwall) for Rwall in Rwalls]))
	r_max = np.max(np.array([np.max(Rwall) for Rwall in Rwalls]))
	z_min = np.min(np.array([np.min(Zwall) for Zwall in Zwalls]))
	z_max = np.max(np.array([np.max(Zwall) for Zwall in Zwalls]))

	return r_min, r_max, z_min, z_max

#	This routine returns the maximum distance of selected walls from one point
#	===========================================

def get_dmax_points_walls(Config, Rpt, Zpt, plasma_wall=True, eirene_wall=False, extra_wall=False):
	
	d_max = 0.
	for iWall in range(len(Config.Walls)):
		Wall = Config.Walls[iWall]
		if  ((((Wall.Type == INTERNAL_PLASMA_WALL) or (Wall.Type == EXTERNAL_PLASMA_WALL)) and plasma_wall) or \
		     (((Wall.Type == INTERNAL_EIRENE_WALL) or (Wall.Type == EXTERNAL_EIRENE_WALL)) and eirene_wall) or \
			  ((Wall.Type == DRAWING_WALL) and extra_wall)):
			d_max_new = np.max(np.sqrt((Wall.Rwall - Rpt)**2 + (Wall.Zwall - Zpt)**2))
			if(d_max_new > d_max): d_max = d_max_new

	return d_max

#	This routine returns the intersection between a countour and walls
#	===========================================

def intersect_walls(Config, cLines):

	cWalls 	   		= types.SimpleNamespace()
	cWalls.arc 		= []
	for iWall in range(len(Config.Walls)):
		Wall = Config.Walls[iWall]

		if  ((((Wall.Type == INTERNAL_PLASMA_WALL) or (Wall.Type == EXTERNAL_PLASMA_WALL)) and plasma_wall) or \
		     (((Wall.Type == INTERNAL_EIRENE_WALL) or (Wall.Type == EXTERNAL_EIRENE_WALL)) and eirene_wall) or \
			  ((Wall.Type == DRAWING_WALL) and extra_wall)):
			cWalls.arc.append(types.SimpleNamespace())
			cWalls.arc[-1].x	= Wall.Rwall
			cWalls.arc[-1].y	= Wall.Zwall
			

	X = intersect_contour(cWalls, cLines)

	return X

#	This routine returns points inside or outside all selected walls
#	===========================================

def get_in_out_walls(Config, PtsR, PtsZ, plasma_wall=True, eirene_wall=False, extra_wall=False):

	pts_inout = np.ones(len(PtsR), dtype = 'f4')
	for iWall in range(len(Config.Walls)):
		Wall = Config.Walls[iWall]
		InOut = 0
		if  (((Wall.Type == INTERNAL_PLASMA_WALL) and plasma_wall) or \
		     ((Wall.Type == INTERNAL_EIRENE_WALL) and eirene_wall)):
			InOut = -1
		elif(((Wall.Type == EXTERNAL_PLASMA_WALL) and plasma_wall) or \
		     ((Wall.Type == EXTERNAL_EIRENE_WALL))and eirene_wall):
			InOut =  1

		if(InOut != 0):
			wall_pts_in, wall_pts_out = get_in_out(Wall.Rwall, Wall.Zwall, PtsR, PtsZ)
			if(InOut == 1): pts_inout[wall_pts_out] = 0
			else:			pts_inout[wall_pts_in]  = 0

	pts_in	  = np.where(pts_inout == 1)[0]
	pts_out   = np.where(pts_inout == 0)[0]

	return pts_in, pts_out


#	===========================================


def get_spline_min_walls(Config, spline, plasma_wall=True, eirene_wall=False, extra_wall=False):
	
	psi_min= 1.e50
	for iWall in range(len(Config.Walls)):
		Wall = Config.Walls[iWall]

		if  ((((Wall.Type == INTERNAL_PLASMA_WALL) or (Wall.Type == EXTERNAL_PLASMA_WALL)) and plasma_wall) or \
		     (((Wall.Type == INTERNAL_EIRENE_WALL) or (Wall.Type == EXTERNAL_EIRENE_WALL)) and eirene_wall) or \
			  ((Wall.Type == DRAWING_WALL) and extra_wall)):
			psi_min_new = np.min(spline.ev(Wall.Rwall, Wall.Zwall))

			if(psi_min_new < psi_min): psi_min = psi_min_new

	return psi_min

def walls_find_closest_segment(Config, Rpt, Zpt, plasma_wall=True, eirene_wall=True, extra_wall=True, skip_wall=-1):

	dMin=1.e50
	for iWall in range(len(Config.Walls)):
		Wall = Config.Walls[iWall]
		if  ((((Wall.Type == INTERNAL_PLASMA_WALL) or (Wall.Type == EXTERNAL_PLASMA_WALL)) and plasma_wall) or \
		     (((Wall.Type == INTERNAL_EIRENE_WALL) or (Wall.Type == EXTERNAL_EIRENE_WALL)) and eirene_wall) or \
			  ((Wall.Type == DRAWING_WALL) and extra_wall)):
			d, i = find_closest_segment(Rpt, Zpt, Wall.Rwall, Wall.Zwall)
			if(d < dMin):
				kWall  = iWall
				iPtSeg = i
				dMin   = d

	return kWall, iPtSeg


def walls_find_closest_point(Config, Rpt, Zpt, plasma_wall=True, eirene_wall=True, extra_wall=True):

	dMin=1.e50
	for iWall in range(len(Config.Walls)):
		Wall = Config.Walls[iWall]
		if  ((((Wall.Type == INTERNAL_PLASMA_WALL) or (Wall.Type == EXTERNAL_PLASMA_WALL)) and plasma_wall) or \
		     (((Wall.Type == INTERNAL_EIRENE_WALL) or (Wall.Type == EXTERNAL_EIRENE_WALL)) and eirene_wall) or \
			  ((Wall.Type == DRAWING_WALL) and extra_wall)):
			d = np.sqrt((Wall.Rwall-Rpt)**2 + (Wall.Zwall-Zpt)**2)
			imin = np.argmin(d)
			if(d[imin] < dMin):
				kWall	= iWall
				iPt		= imin
				dMin	= d[imin]

	return kWall, iPt
