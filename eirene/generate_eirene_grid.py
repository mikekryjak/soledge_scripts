
import os
import types
import numpy as np

from routines.globals									import *

from eirene.find_center_position						import find_center_position
from eirene.generate_knots								import generate_knots
from eirene.find_inplasma								import find_inplasma
from eirene.direct_or_indirect							import zones_direct_or_indirect
from eirene.detect_cut									import detect_cut
from eirene.detect_aligned_mesh							import detect_aligned_mesh
from eirene.remove_knots_in_solid						import remove_knots_in_solid
from eirene.generate_knots_on_wall						import generate_knots_on_wall

from eirene.finalization_remove							import finalization_remove
from eirene.generate_triangles							import generate_triangles_in_the_plasma, generate_triangles_on_wall, generate_triangles_center
from eirene.triangles_routines							import KnotsInterp_reshape


#from soledge.clibs.eirene_clib							import find_neighbors
from eirene.find_neighbors								import find_neighbors

from eirene.find_connexity								import find_connexity
from eirene.check_triangles_orientation_size			import check_triangles_orientation_size
from eirene.remove_unused_triangles						import remove_unused_triangles
from eirene.generate_data_for_interpolation_pass1_2 	import generate_data_for_interpolation_pass1_2
from eirene.generate_data_for_interpolation_pass3_4 	import generate_data_for_interpolation_pass3_4
from eirene.check_interpolation_triangles				import check_interpolation_triangles

from eirene.generate_wall_triangles						import generate_wall_triangles
from eirene.generate_vacuum_triangles					import generate_vacuum_triangles
from eirene.generate_vacuum_wall_triangles				import generate_vacuum_wall_triangles

from eirene.generate_interpolation_flux_on_wall			import generate_interpolation_flux_on_wall
from eirene.assign_triangles_to_trans_mesh				import assign_triangles_to_trans_mesh
from eirene.compute_weights								import compute_weights
from eirene.remapping_para								import remapping_para
from eirene.assign_triangles_to_trans_mesh_perp			import assign_triangles_to_trans_mesh_perp
from eirene.compute_weights_perp						import compute_weights_perp
from eirene.remapping_perp								import remapping_perp
from eirene.check_remapping_weights						import check_remapping_weights
from eirene.wall_segments								import wall_segments
from eirene.find_sequence								import find_sequence
from eirene.get_ext_plasma_triseq						import get_ext_plasma_triseq
from eirene.assign_non_penalized_cell					import assign_non_penalized_cell

from interfaces.progressbar 							import ProgressBar


# Function definition is here
#=======================================================================================================================================

#=========================================================
# This routine generate eirene grid
#=========================================================
#

def generate_eirene_grid(Root, Config, EWalls, ToTheCenter, plot_check=0):

	print("generate_eirene_grid")

	Zones	 = Config.Zones
	Walls	 = Config.Walls
	X_points = Config.X_points

	if(ToTheCenter): Rcenter, Zcenter = find_center_position(Zones)

#	Note:	In matlab version there is move_wall_inward routine which remove duplicated points and move:
#			 external wall (wall.type=1) inward
#			internal wall (wall.type=-1) outward
#			Inward/outward movement is done to avoid problem at the strike points with alligne mesh!?!
#			We manage differently knots on the wall and so we do not need move wall in/out but we clean duplicated points from the wall when we load it!
	
	Nsteps = 12
	prog_bar = ProgressBar(Root, Title="Eirene grid generation",Label="Generating eirene knots", Value=0.) 

	nKnots, RKnots, ZKnots = generate_knots(Zones, X_points)
	InPWallsKnots = find_inplasma(Config, Zones, RKnots, ZKnots)

	Direct = zones_direct_or_indirect(Zones, RKnots, ZKnots)						#Direct (0=indirect,1=direct)
	
	print("\tknots generated in the plasma")

	prog_bar.Update(1/Nsteps,Label="Detecting wall intersections")

	detect_cut(Root, Config, Zones, RKnots, ZKnots)							#% refinement of Zones[k].iscrossed[i,j,1:4,iWall] (one value for each 4 segments of each quadrangle, zero if no crossing with iWall)

	CheckAlligned = False													#Before setting to true detect_aligned_mesh must be fixed
	detect_aligned_mesh(Zones, Walls, CheckAlligned)						#This routine is not updated to the last IsCrossed shape

	prog_bar.Update(2/Nsteps,Label="Removing knots in the wall")

#	All knots that are inside a material are removed:  output in ToKnotsPlasma, RKnotsPlasma,ZKnotsPlasma 
	ToKnotsPlasma, RKnotsPlasma, ZKnotsPlasma	= remove_knots_in_solid(Zones, RKnots, ZKnots, InPWallsKnots)

	prog_bar.Update(3/Nsteps,Label="Generating knots on the wall")

#	All knots in the wall: output in ToKnotsWall, RKnotsWall, ZKnotsWall 
	ToKnotsWall, RKnotsWall, ZKnotsWall			= generate_knots_on_wall(Root, Config, Zones, RKnots, ZKnots)

#	Remove all unnecessary knots
	nKnots, RKnots, ZKnots						= finalization_remove(Zones, ToKnotsPlasma, RKnotsPlasma, ZKnotsPlasma, ToKnotsWall, RKnotsWall, ZKnotsWall)

	print("\tknots generated on the wall")

	ToKnotsPlasma	= 0
	RKnotsPlasma	= 0
	ZKnotsPlasma	= 0
	ToKnotsWall		= 0
	RKnotsWall		= 0
	ZKnotsWall		= 0

#	adding center
	if(ToTheCenter):
		nKnotsCenter	= nKnots.max()+1;
		nKnots			= np.append(nKnots, nKnotsCenter)
		RKnots			= np.append(RKnots, Rcenter)
		ZKnots			= np.append(ZKnots, Zcenter)
		if(DEBUG > 1): print('\tknots generated at the center')

		
	prog_bar.Update(4/Nsteps,Label="Generating triangles in the plasma")
	
	Triangles, TriKnots = generate_triangles_in_the_plasma(Zones, Direct)
	
	print("\ttriangles generated in the plasma")

	prog_bar.Update(5/Nsteps,Label="Generating triangles on the wall")
	Triangles, TriKnots = generate_triangles_on_wall(Zones,Triangles, TriKnots, Direct)
	print("\ttriangles generated on the wall")
	
	if(ToTheCenter): Triangles, TriKnots = generate_triangles_center(Zones, nKnotsCenter, Triangles, TriKnots, Direct)

	prog_bar.Update(6/Nsteps,Label="Finding triangle neighbors")
	Triangles = find_neighbors(Triangles, TriKnots, RKnots, ZKnots) 
	
	prog_bar.Update(7/Nsteps,Label="check triangles")

	Triangles = check_triangles_orientation_size(Root, RKnots, ZKnots, Triangles)
	nKnots, RKnots, ZKnots, Triangles, TriKnots = remove_unused_triangles(Zones, RKnots, ZKnots, Triangles, TriKnots)

	find_connexity(Root, Triangles, RKnots, ZKnots)
	nKnots, RKnots, ZKnots, Triangles, TriKnots = remove_unused_triangles(Zones, RKnots, ZKnots, Triangles, TriKnots)

	prog_bar.Update(8/Nsteps,Label="Generating soledge-eirene interpolations")

	KnotsInterp = generate_data_for_interpolation_pass1_2(Zones, nKnots,  RKnots, ZKnots, Triangles, TriKnots)

	generate_data_for_interpolation_pass3_4(Zones, RKnots, ZKnots, KnotsInterp)

	if(ToTheCenter):	
		KnotsInterp.assp[-1] = 5
		KnotsInterp.nsol[-1] = 0
		KnotsInterp.neir[-1] = 0

	check_interpolation_triangles(Zones, RKnots, ZKnots, KnotsInterp)

	prog_bar.Update(9/Nsteps,Label="Determining parallel wall flux interpolations")

	WallTriangles = generate_wall_triangles(Triangles)
	generate_interpolation_flux_on_wall(Zones)

	prog_bar.Update(10/Nsteps,Label="Find correspondence between quadrangle and triangles")

	status = assign_triangles_to_trans_mesh(Root, WallTriangles, RKnots, ZKnots)
	if(not status): return															#exit due to error in assign_triangles_to_trans_mesh

	compute_weights(Zones, WallTriangles, RKnots, ZKnots)
	remapping_para(Zones, WallTriangles, RKnots, ZKnots)

	prog_bar.Update(11/Nsteps,Label="Determining perpendicular wall flux interpolations")

	status = assign_triangles_to_trans_mesh_perp(Root, WallTriangles, RKnots, ZKnots)
	if(not status): return															#exit due to error in assign_triangles_to_trans_mesh_perp

	compute_weights_perp(Zones, WallTriangles, RKnots, ZKnots)
	remapping_perp(Zones, WallTriangles, RKnots, ZKnots)

	check_remapping_weights(Root, Zones, WallTriangles)

	EireneWall = types.SimpleNamespace()
	EireneWall.TriSequences, EireneWall.TriFaces, EireneWall.R12, EireneWall.Z12 = wall_segments(WallTriangles, RKnots, ZKnots)
	EireneWall.EWalls = EWalls

#	Now define vacuum triangles
#	=================

	if(len(EWalls) > 0):
		Triangles, RKnots, ZKnots = generate_vacuum_triangles(Root, Config, EireneWall, Triangles, WallTriangles, RKnots, ZKnots)

		KnotsInterp = KnotsInterp_reshape(KnotsInterp, nKnotsInterp=len(ZKnots))									#Extend KnotsInterp

		TriKnots	= np.array([Triangles.p1,Triangles.p2,Triangles.p3]).T
		Triangles	= find_neighbors(Triangles, TriKnots, RKnots, ZKnots)											#Fix neighbors

		WallTriangles	= generate_vacuum_wall_triangles(Triangles, WallTriangles)									#Fix wall triangles

		EireneWall.TriSequences, EireneWall.TriFaces, EireneWall.R12, EireneWall.Z12 = wall_segments(WallTriangles, RKnots, ZKnots)		#Update wall

	get_ext_plasma_triseq(Config, EireneWall)																		#get/set external triangle sequence

	if(Root.FileType == FILE_TYPE_SOLEDGE3X): assign_non_penalized_cell(Zones, Triangles)							#Non penalized triangles for S3X


	Eirene		  			= types.SimpleNamespace()
	Eirene.ToTheCenter		= ToTheCenter
#	Eirene.Direct			= Direct							#No more used

#	Eirene.Knots			= nKnots									#in matlab is nknots_
	Eirene.RKnots 			= RKnots							#saved in save_eirene_triangles 	(in matlab is R_ (in m.))
	Eirene.ZKnots			= ZKnots							#saved in save_eirene_triangles		(in matlab is Z_ (in m.))
	Eirene.KnotsInterp  	= KnotsInterp						#saved in save_eirene_triangles		(in matlab is knots_interp
	Eirene.WallTriangles	= WallTriangles						#saved in save_eirene_triangles		(in matlab is triangle)
	Eirene.Triangles 		= Triangles							#saved in save_eirene_triangles (p1/2/3,BC/1/2/3,k,i,j) in save_neighbors (neigh1/2/3, typeneigh1/2/3)		in matlab is triangles_
	Eirene.Wall				= EireneWall

#	Wall definition

	Eirene.Wall.TypeMat = []
	Eirene.Wall.IsPump  = []
	for iPwall in range(len(Eirene.Wall.TriSequences)):
		Eirene.Wall.TypeMat.append(np.zeros((len(Eirene.Wall.TriSequences[iPwall])), dtype='i4'))
		Eirene.Wall.IsPump.append(np.zeros((len(Eirene.Wall.TriSequences[iPwall])), dtype='bool'))

	Eirene.Wall.Material			= []
	Eirene.Wall.Material.append(types.SimpleNamespace())
	Eirene.Wall.Material[-1].Name 	= "W"

	Eirene.Wall.Pump		 		= []
	Eirene.Wall.Puff				= []

	Eirene.Surfaces					= []

	if(DEBUG > 1):
		print("\tgenerated:")
		print("\t\t", len(Eirene.RKnots)," nodes ")
		print("\t\t", len(Eirene.Triangles)," triangles ")
		print("\t\t", len(Eirene.WallTriangles)," wall triangles ")
		
	if(DEBUG > 1):print("\tEirene grid generated")
	
	prog_bar = 0

	return Eirene
