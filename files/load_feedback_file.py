import os
import types
import numpy						as np
import scipy.interpolate			as interp
from routines.globals				import DEBUG
from routines.intersect_contour		import intersect_contour

#====================================
# This routine save profiles as text
#====================================

def load_feedback_file(FeedFile, FeedbackFileType, Config, Cut):

	if(DEBUG > 0): print("load_feedback_file: Reading from ",FeedFile)

	FileData = np.loadtxt(FeedFile)	

	if(FeedbackFileType == 0):
		Headers  = ["xv", "Ni", "Ti", "Te", "In_Ni", "In_Te", "In_Ti"]		#feedProf xv, Ni, Ti, Te
	else:
		Headers  = ["xv", "D", "Chi", "Chie"]									#diffusion xv, D, Chi, Chie
#		Chi 	 = 0.5*(FileData[:,2] + FileData[:,3])
#		FileData = np.append(FileData, Chi.reshape(Chi.shape[0],1), axis=1)
		
	MagZones	 = Config.MagZones
	X_points = Config.X_points
	Profiles = Cut.Flux.Profiles

	dCut	 = np.empty(0, dtype='f8')
	xCut	 = np.empty(0, dtype='f8')

	c1		= types.SimpleNamespace()
	c1.arc	= [types.SimpleNamespace()]
	c2		= types.SimpleNamespace()
	c2.arc	= [types.SimpleNamespace()]
	for k in range(len(Cut.Flux.nz)):						#Loop on Flux cutted zones
		c1.arc[0].x = Cut.Flux.r12[k,:]
		c1.arc[0].y = Cut.Flux.z12[k,:]
		iZone		= Cut.Flux.nz[k]

		nGrid = MagZones[iZone].gridRc.shape[0]
		rz = np.empty((nGrid, 2), dtype='f8')
		rz[0, 0] = Cut.Flux.r12[k,0]					
		rz[0, 1] = Cut.Flux.z12[k,0]					
		rz[-1, 0] = Cut.Flux.r12[k,1]					
		rz[-1, 1] = Cut.Flux.z12[k,1]	
						
		for iGrid in range(1,nGrid-1):
			c2.arc[0].x = MagZones[iZone].gridRc[iGrid,:]			#grid contour
			c2.arc[0].y = MagZones[iZone].gridZc[iGrid,:]						
			X  = intersect_contour(c1,c2)
			if(len(X) != 1):
				print("\tError in load_SOLEDGE_feedbak_file: number of intersections=", len(X))
				exit()	
			rz[iGrid, 0] = X[0].x					
			rz[iGrid, 1] = X[0].y
								
		dGrid = np.sqrt((rz[:,0] - Cut.Flux.r12[k,0])**2 + (rz[:,1] - Cut.Flux.z12[k,0])**2) + Cut.Flux.d12[k,0]				


		dCut = np.append(dCut,dGrid)
		xCut = np.append(xCut, MagZones[iZone].x[:,0])
			

	for i in range(len(Profiles)):

		ix = indexes_upper(Headers, "xv")
		iv = indexes_upper(Headers, Profiles[i].Name)
		if((len(ix) == 1) and (len(iv) == 1)):
			ix = ix[0]
			iv = iv[0]

			dFileData		 = np.interp(FileData[:,ix] , xCut, dCut)

			InRange = np.where((dFileData  > Profiles[i].xValues[0]) &  (dFileData  < Profiles[i].xValues[-1]))[0]
			iFirst			 = np.argmin(np.abs(dFileData - Profiles[i].xValues[0]))
			iLast			 = np.argmin(np.abs(dFileData - Profiles[i].xValues[-1]))
			xValuesNew		 = np.empty((len(InRange)+2), dtype='f8')
			ValuesNew		 = np.empty((len(InRange)+2), dtype='f8')
										
			xValuesNew[1:-1] = dFileData[InRange]
			ValuesNew[1:-1]  = FileData[InRange,iv]
			
			xValuesNew[0]	 = Profiles[i].xValues[0]
			xValuesNew[-1]	 = Profiles[i].xValues[-1]
			
			ValuesNew[0]	 = FileData[iFirst,iv]
			ValuesNew[-1]	 = FileData[iLast,iv]
			
			Profiles[i].xValues	 = xValuesNew
			Profiles[i].Values	 = ValuesNew
			if(DEBUG > 1): print("Loaded: ",Profiles[i].Name)

	if(DEBUG > 0): print("load_feedback_file: Completed")

	return

def indexes_upper(list_strs, str):
	return  [i for i, j in enumerate(list_strs) if j.upper() == str.upper()]
