import types
import numpy							as np
import scipy.interpolate				as interp
from interfaces.plot_and_ask			import plot_and_ask
from routines.intersect_contour			import intersect_2contours
from interfaces.progressbar				import ProgressBar
from routines.utils_walls				import set_XYlines_walls
from routines.globals					import DEBUG
from mesh.compute_mesh_intersections	import extend_RZ_east_west
	
def generate_feedback_profiles(Root, Config):

	if(DEBUG > 0): print("generate_feedback_profiles")
	
	X_points		= Config.X_points
	MagZones		= Config.MagZones
	MagMegazones	= Config.MagMegazones	
	Zones			= Config.Zones
	FeedbackTransp	= Config.FeedbackTransp
	
	pro_bar = ProgressBar(Root, Label="Wait Generation")
	
#	Root.update()

	NxTot	= 0
	for k in range(len(FeedbackTransp.Cut.Flux.nz)):						#Loop on Flux cutted zones
		iMagZone = FeedbackTransp.Cut.Flux.nz[k]
		NxTot	+= MagZones[iMagZone].gridRc.shape[0]

	if(len(X_points) > 0):
		iSep = 0
		i	 = 0
		while (MagMegazones[MagZones[FeedbackTransp.Cut.Flux.nz[i]].mz].isperiodic):
			iSep += MagZones[FeedbackTransp.Cut.Flux.nz[i]].gridRc.shape[0]
			i	 += 1
	else:
		iSep = 100000

	FeedbackTransp.Data.IndexSep = iSep -	1

	FeedbackTransp.Data.NxTot		= NxTot
		
	FeedbackTransp.Data.Profiles = np.zeros((len(FeedbackTransp.Cut.Flux.Profiles)+1, NxTot), dtype='f8')				#One more profile for coordinate
	FeedbackTransp.Data.RefMin   = np.zeros((len(FeedbackTransp.Cut.Flux.Profiles), 2), dtype='f8')
	 
#	Find poloidal position from intersection of first cutted zone

	iMagZone		= FeedbackTransp.Cut.Flux.nz[0]

	gridR, gridZ = extend_RZ_east_west(MagZones, iMagZone, 1)
	Ri, Zi, is1, is2 = intersect_2contours(FeedbackTransp.Cut.Flux.r12[0,:], FeedbackTransp.Cut.Flux.z12[0,:], gridR, gridZ)						#Intersection with second flux contour
	if(len(Ri) != 1):
		print("generate_feedback_profiles.1: Error in number of intersections=", len(Ri))
		
		set_XYlines_walls(Config,"k.-")
		Xarrs.extend([FeedbackTransp.Cut.Flux.r12[0,:], gridR])
		Yarrs.extend([FeedbackTransp.Cut.Flux.z12[0,:], gridZ])
		lines.extend(["r-","g.-"])

		LinesData = [Xarrs, Yarrs, lines]
		choiche = plot_and_ask(Root, LinesData=LinesData, title="Error in number of intersections")
		exit()	
			
	jMagZ1 = np.argmin((MagZones[iMagZone].gridRc[1,:] - Ri[0])**2 + (MagZones[iMagZone].gridZc[1,:] - Zi[0])**2)
	
	gridR, gridZ = extend_RZ_east_west(MagZones, iMagZone, -2)
	Ri, Zi, is1, is2 = intersect_2contours(FeedbackTransp.Cut.Flux.r12[0,:], FeedbackTransp.Cut.Flux.z12[0,:], gridR, gridZ) 						#Intersection with second-last flux contour
	if(len(Ri) != 1):
		Xarrs, Yarrs, lines= set_XYlines_walls(Config, wall_line="k.-")
		Xarrs.extend([FeedbackTransp.Cut.Flux.r12[0,:], gridR])
		Yarrs.extend([FeedbackTransp.Cut.Flux.z12[0,:], gridZ])
		lines.extend(["r-","g.-"])

		LinesData = [Xarrs, Yarrs, lines]
		choiche = plot_and_ask(Root, LinesData=LinesData, title="Error in number of intersections")
		print("generate_feedback_profiles.-2: Error in number of intersections=", len(Ri))
		exit()

	jMagZ2		= np.argmin((MagZones[iMagZone].gridRc[-2,:]- Ri[0])**2 + (MagZones[iMagZone].gridZc[-2,:] - Zi[0])**2)

#	 Find poloidal intersection and intersected Zones from MagZone and Mag intersection

	jMagZ	= int(0.5*(jMagZ1+jMagZ2))

	north   = -1
	for sz_k in MagZones[FeedbackTransp.Cut.Flux.nz[0]].list:
		sz_sj	= Zones[sz_k].magz[2]
		if((sz_sj <= jMagZ) and (sz_sj + Zones[sz_k].gridRc.shape[1] > jMagZ)):
			FeedbackTransp.Data.jZ = jMagZ - sz_sj
			FeedbackTransp.Data.nZ = np.array([sz_k])
			north = Zones[sz_k].Neighbour.north
			break

	if(north == -1):
		print("generate_feedback_profiles.: Not found Zone for intersected MagZone")
		exit()
		
	while(north >=0):
		FeedbackTransp.Data.nZ = np.append(FeedbackTransp.Data.nZ, north)
		north = Zones[north].Neighbour.north	
	
	FeedbackTransp.Data.input_r  = np.zeros(NxTot, dtype='f8')
	FeedbackTransp.Data.x		  = np.zeros(NxTot, dtype='f8')
	FeedbackTransp.Data.index	  = np.arange(NxTot, dtype='f8')

#	Prepare x and input_r arrays

	NxTot	= 0
	for k in range(len(FeedbackTransp.Cut.Flux.nz)):
		iMagZone	 = FeedbackTransp.Cut.Flux.nz[k]
		Nx		 = MagZones[iMagZone].gridRc.shape[0]
		FeedbackTransp.Data.x[NxTot:NxTot+Nx] 	  = MagZones[iMagZone].x[:,0]
		FeedbackTransp.Data.input_r[NxTot:NxTot+Nx] = MagZones[iMagZone].gridRc[:,jMagZ]
		NxTot	+= Nx

	iProg	= 0
	NxTot	= 0
	for k in range(len(FeedbackTransp.Cut.Flux.nz)):					#Loop on Flux cutted zones
		iMagZone		= FeedbackTransp.Cut.Flux.nz[k]

		nGrid = MagZones[iMagZone].gridR.shape[0]
		rz = np.empty((nGrid, 2), dtype='f8')
		rz[0, 0]  = FeedbackTransp.Cut.Flux.r12[k,0]					
		rz[0, 1]  = FeedbackTransp.Cut.Flux.z12[k,0]					
		rz[-1, 0] = FeedbackTransp.Cut.Flux.r12[k,1]					
		rz[-1, 1] = FeedbackTransp.Cut.Flux.z12[k,1]	
						
		for iGrid in range(1,nGrid-1):							#grid contours
			gridR, gridZ = extend_RZ_east_west(MagZones, iMagZone, iGrid)
			Ri, Zi, is1, is2  = intersect_2contours(FeedbackTransp.Cut.Flux.r12[k,:], FeedbackTransp.Cut.Flux.z12[k,:], gridR, gridZ)
			if(len(Ri) != 1):
				Xarrs, Yarrs, lines= set_XYlines_walls(Config, wall_line="k.-")
				Xarrs.extend([FeedbackTransp.Cut.Flux.r12[k,:], gridR])
				Yarrs.extend([FeedbackTransp.Cut.Flux.z12[k,:], gridZ])
				lines.extend(["r-","g.-"])

				LinesData = [Xarrs, Yarrs, lines]
				choiche = plot_and_ask(Root, LinesData=LinesData, title="Error in number of intersections")
				print("generate_feedback_profiles.k: Error in number of intersections=", len(Ri))
				exit()

			rz[iGrid, 0] = Ri[0]					
			rz[iGrid, 1] = Zi[0]
								
		dGrid = np.sqrt((rz[:,0] - FeedbackTransp.Cut.Flux.r12[k,0])**2 + (rz[:,1] - FeedbackTransp.Cut.Flux.z12[k,0])**2) + FeedbackTransp.Cut.Flux.d12[k,0]
		dGrid = 0.5*(dGrid[:-1] + dGrid[1:])
		
		for iProf in range(len(FeedbackTransp.Cut.Flux.Profiles)):
			FeedbackTransp.Data.Profiles[iProf, NxTot : NxTot + len(dGrid)] = np.interp(dGrid, FeedbackTransp.Cut.Flux.Profiles[iProf].xValues,  FeedbackTransp.Cut.Flux.Profiles[iProf].Values)

		FeedbackTransp.Data.Profiles[-1, NxTot : NxTot + len(dGrid)] = dGrid

		NxTot	+= len(dGrid)
	
		iProg	+= 1
		pro_bar.Update(iProg/len(FeedbackTransp.Cut.Flux.nz))			
	
	pro_bar = 0

	if(DEBUG > 0): print("generate_feedback_profiles: completed")

	return
