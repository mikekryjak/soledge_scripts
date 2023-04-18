import types
import numpy							as np

from mesh.generate_soledge_grid			import generate_soledge_grid
from mesh.generate_soledge_xz			import generate_soledge_xz
from mesh.detect_core					import detect_core
from mesh.neigh_props					import neigh_props
from interfaces.progressbar				import ProgressBar


def generate_soledge_mesh(Root, Config):

	prog_bar = ProgressBar(Root, Label="Generating soledge mesh....", Value=0.) 
	
	MagZones		= Config.MagZones
	MagMegazones	= Config.MagMegazones
	MagPMegazones	= Config.MagPMegazones
	
#	Copy from unsplitted data to splitted ones

	Zones  = []
	for k in range(len(MagZones)):
		Zones.append(types.SimpleNamespace())
		
		Zones[k].gridR	= np.copy(MagZones[k].gridR)
		Zones[k].gridZ	= np.copy(MagZones[k].gridZ)
		
		Zones[k].Neighbour	= types.SimpleNamespace()
		Zones[k].Neighbour.north	= MagZones[k].Neighbour.north
		Zones[k].Neighbour.south	= MagZones[k].Neighbour.south
		Zones[k].Neighbour.east		= MagZones[k].Neighbour.east
		Zones[k].Neighbour.west		= MagZones[k].Neighbour.west
		Zones[k].mz					= MagZones[k].mz
		Zones[k].pmz				= MagZones[k].pmz
		Zones[k].magz				= np.array([k,0,0], dtype='i4')
		MagZones[k].list			= np.array([k], dtype='i4')
		
	Megazones = []
	for k in range(len(MagMegazones)):
		Megazones.append(types.SimpleNamespace())
		Megazones[k].list  		= np.copy(MagMegazones[k].list)
		Megazones[k].isperiodic	= MagMegazones[k].isperiodic	

	PMegazones = []
	for k in range(len(MagPMegazones)):
		PMegazones.append(types.SimpleNamespace())
		PMegazones[k].list  		= np.copy(MagPMegazones[k].list)
		
	Config.Zones 	 	= Zones
	Config.Megazones	= Megazones
	Config.PMegazones	= PMegazones

#	Compute mesh
#======================================================	
	
	Step   = 0
	nSteps = 4
	generate_soledge_grid(Config)
	Step   += 1
	prog_bar.Update(Step/nSteps)
	
	generate_soledge_xz(Config)
	Step   += 1
	prog_bar.Update(Step/nSteps)

	detect_core(Config)
	Step   += 1
	prog_bar.Update(Step/nSteps)
	
	neigh_props(Config)
	prog_bar = 0

#	Copy mesh to unsplitted zones
#=======================================================

	for k in range(len(MagZones)):												#duplicate some array for convenience

		MagZones[k].gridRc	=  np.copy(Zones[k].gridRc)
		MagZones[k].gridZc	=  np.copy(Zones[k].gridZc)
		MagZones[k].Br	 	=  np.copy(Zones[k].Br)
		MagZones[k].Bz	 	=  np.copy(Zones[k].Bz)
		MagZones[k].Bphi 	=  np.copy(Zones[k].Bphi)
		MagZones[k].x 	  	= np.copy(Zones[k].x)
		MagZones[k].z	  	= np.copy(Zones[k].z)

		MagZones[k].Nx		= MagZones[k].gridRc.shape[0]
		MagZones[k].Nz		= MagZones[k].gridRc.shape[1]

		MagZones[k].Neighbour.north	= Zones[k].Neighbour.north					#Update Neighbour
		MagZones[k].Neighbour.south	= Zones[k].Neighbour.south
		MagZones[k].Neighbour.east	= Zones[k].Neighbour.east
		MagZones[k].Neighbour.west	= Zones[k].Neighbour.west

	
	Root.set_Mesh_OK(Config, True)

	
	return
