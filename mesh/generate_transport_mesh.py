import types
import numpy as np
from math							import sqrt
from matplotlib.path 				import Path
from routines.intersect_contour		import intersect_contour
from routines.globals				import *
from mesh.get_rz_core_sep			import get_rz_core_sep
from interfaces.progressbar			import ProgressBar
	
def generate_transport_mesh(MainWindow, Config, TranspCuts):
	c1		= types.SimpleNamespace()
	c1.arc	= [types.SimpleNamespace()]
	c2		= types.SimpleNamespace()
	c2.arc	= [types.SimpleNamespace()]

	MagZones = Config.MagZones
	Zones	 = Config.Zones

	IsBallooningMode = False
	for i in range(len(TranspCuts)):
		if(TranspCuts[i].BallooningMode != 0):
			IsBallooningMode = True
			break
	if(IsBallooningMode):
		Rcore, Zcore, CoreMegazone = get_rz_core_sep(Config, core_and_sep=False)
		CorePath = Path(np.array([Rcore, Zcore]).T, closed=True)
		InCore   = np.where(CorePath.contains_points(np.array([Config.r2D.reshape(-1),Config.z2D.reshape(-1)]).T))[0]
		iAxis	 = InCore[np.argmin(Config.flux2D.reshape(-1)[InCore])]

		BrAxis	 = Config.Br2D.reshape(-1)[iAxis]
		BzAxis	 = Config.Bz2D.reshape(-1)[iAxis]
		BphiAxis = Config.Bphi2D.reshape(-1)[iAxis]
		InCore	 = 0
		CorePath = 0
		
		

	Config.flux2D
	pro_bar = ProgressBar(MainWindow, Label="Wait Generation")
#	MainWindow.update()
	nTot	= len(MagZones)
	iProg	= 0
	for i in range(len(TranspCuts)):									#Loop on Cuts
		for k in range(len(TranspCuts[i].Flux.nz)):						#Loop on Flux cutted zones
			c1.arc[0].x = TranspCuts[i].Flux.r12[k,:]
			c1.arc[0].y = TranspCuts[i].Flux.z12[k,:]
			iMagZone	= TranspCuts[i].Flux.nz[k]
			if(not MagZones[iMagZone].IsTranspValues):
				nGrid = MagZones[iMagZone].gridR.shape[0]
				dGrid = np.empty(nGrid, dtype='f8')
		#		Xj    = np.empty(nGrid, dtype='i4')
				for iGrid in range(1,nGrid-1):
					c2.arc[0].x = MagZones[iMagZone].gridR[iGrid,:]			#grid contour
					c2.arc[0].y = MagZones[iMagZone].gridZ[iGrid,:]						
					X  = intersect_contour(c1,c2)
					if(len(X) != 1):
						print("\tError in generate_transport_mesh: number of intersections=", len(X))
						exit()					
					dGrid[iGrid] = np.sqrt((X[0].x - TranspCuts[i].Flux.r12[k,0])**2 + (X[0].y - TranspCuts[i].Flux.z12[k,0])**2)					
		#			Xj[iGrid]    = np.argmin((c2.arc[0].x - X[0].x)**2 + (c2.arc[0].y - X[0].y)**2)	

				if(dGrid[-2] > dGrid[1]):
		#			Xj[0]     = np.argmin((MagZones[iMagZone].gridR[0,:]  - c1.arc[0].x[0] )**2 + (MagZones[iMagZone].gridZ[0,:]  - c1.arc[0].y[0] )**2)	
		#			Xj[-1]    = np.argmin((MagZones[iMagZone].gridR[-1,:] - c1.arc[0].x[-1])**2 + (MagZones[iMagZone].gridZ[-1,:] - c1.arc[-1].y[0])**2)	
					dGrid[0]  = TranspCuts[i].Flux.d12[k,0]
					dGrid[-1] = TranspCuts[i].Flux.d12[k,1]
				else:
		#			Xj[0]     = np.argmin((MagZones[iMagZone].gridR[-1,:] - c1.arc[0].x[0] )**2 + (MagZones[iMagZone].gridZ[-1,:] - c1.arc[0].y[0] )**2)	
		#			Xj[-1]    = np.argmin((MagZones[iMagZone].gridR[0,:]  - c1.arc[0].x[-1])**2 + (MagZones[iMagZone].gridZ[0,:]  - c1.arc[-1].y[0])**2)	
					dGrid[0]  = TranspCuts[i].Flux.d12[k,1]
					dGrid[-1] = TranspCuts[i].Flux.d12[k,0]
					
				dGrid[1:-1]  += TranspCuts[i].Flux.d12[k,0]				
				dFluxes		= 0.5*(dGrid[1:] + dGrid[:-1])

#				define gridZc
#				gridZc = 0.25*(MagZones[iMagZone].gridZ[:-1,:-1] + MagZones[iMagZone].gridZ[:-1,1:] + MagZones[iMagZone].gridZ[1:,1:] + MagZones[iMagZone].gridZ[1:,:-1])

				"""
				if(TranspCuts[i].BallooningMode != 0):
		
					Xj		= (0.5*(Xj[1:] + Xj[:-1])).astype(np.int)
					Xi		= np.arange(MagZones[iMagZone].gridZ.shape[0]-1, dtype='i4')

					if(TranspCuts[i].BallooningMode	== 1):				#Ballooning 1/Bt
						BAtCut = MagZones[iMagZone].Bphi[Xi,Xj[Xi]]
					elif(TranspCuts[i].BallooningMode == 2):			#Ballooning 1/B
						BAtCut = np.sqrt(MagZones[iMagZone].Br[Xi,Xj[Xi]]**2 + MagZones[iMagZone].Bz[Xi,Xj[Xi]]**2 + MagZones[iMagZone].Bphi[Xi,Xj[Xi]]**2)
					else:
						print("\tError in generate_transport_mesh: invalid BallooningMode=", TranspCuts[i].BallooningMode)
						exit()					
				"""

				for kMagZone in TranspCuts[i].Flux.MagZones[k]:				#Loop on Megazone
					iThetaZone = np.where(TranspCuts[i].Theta.nz == kMagZone); iThetaZone = iThetaZone[0]
					if(len(iThetaZone) != 1):
						print("\tError in generate_transport_mesh: not found zone=", len(iThetaZone))
						exit()
					else:
						iThetaZone = iThetaZone[0]

					FactdTheta = (TranspCuts[i].Theta.d12[iThetaZone,1]-TranspCuts[i].Theta.d12[iThetaZone,0])/MagZones[kMagZone].gridZ.shape[1]	
					dThetas	   = FactdTheta*(np.arange(MagZones[kMagZone].gridZ.shape[1]-1, dtype='f8') + 0.5)+TranspCuts[i].Theta.d12[iThetaZone,0]

					if(TranspCuts[i].BallooningMode != 0):
						if(TranspCuts[i].BallooningMode	== 1):				#Ballooning 1/Bt
							BallFact   = BphiAxis/MagZones[kMagZone].Bphi

						elif(TranspCuts[i].BallooningMode == 2):			#Ballooning 1/B
							BallFact =  sqrt(BrAxis**2 + BzAxis**2 + BphiAxis**2)/ \
										np.sqrt(MagZones[kMagZone].Br**2 + MagZones[kMagZone].Bz**2 + MagZones[kMagZone].Bphi**2)

						elif(TranspCuts[i].BallooningMode	== 3):			#Ballooning 1/Bt^2
							BallFact   = BphiAxis/MagZones[kMagZone].Bphi**2

						elif(TranspCuts[i].BallooningMode == 4):			#Ballooning 1/B^2
							BallFact =  (BrAxis**2 + BzAxis**2 + BphiAxis**2)/ \
										(MagZones[kMagZone].Br**2 + MagZones[kMagZone].Bz**2 + MagZones[kMagZone].Bphi**2)

						"""
						if(MagZones[iMagZone].mz == MagZones[kMagZone].mz):
							Bfact = BAtCut
						else:					
							Bfact = np.interp(MagZones[kMagZone].x[:,0], MagZones[iMagZone].x[:,0], BAtCut)


						BallFact = np.outer(Bfact, np.ones(MagZones[kMagZone].gridZ.shape[1]-1, dtype='f8'))/BallFact
						"""

					MagZones[kMagZone].Ballooning = []
					for iTransp in range(len(TranspCuts[i].Flux.Profiles)):
						vFluxes = np.interp(dFluxes, TranspCuts[i].Flux.Profiles[iTransp].xValues,  TranspCuts[i].Flux.Profiles[iTransp].Values)
						
						if(MagZones[iMagZone].mz != MagZones[kMagZone].mz):
							vFluxes = np.interp(MagZones[kMagZone].x[:,0], MagZones[iMagZone].x[:,0], vFluxes)
						vThetas = np.interp(dThetas, TranspCuts[i].Theta.Profiles[iTransp].xValues, TranspCuts[i].Theta.Profiles[iTransp].Values)
							
						MagZones[kMagZone].Ballooning.append(np.outer(vFluxes, vThetas))			#Multiply flux dependence to theta dependence
						if(TranspCuts[i].BallooningMode != 0): 
							MagZones[kMagZone].Ballooning[iTransp] *= BallFact							#Multiply by balloonig factor
					
					iProg	+= 1
					pro_bar.Update(iProg/nTot)
#					MainWindow.update()
					
					MagZones[kMagZone].IsTranspValues = True

							
#	Transfer Ballooning values to splitted zones

	if(not (Zones is MagZones)):																	#Transfer for true MagZones
		for k in range(len(MagZones)):
			for sz_k in MagZones[k].list:
				if(not Zones[sz_k].IsTranspValues):
					sz_si	= Zones[sz_k].magz[1]
					sz_sj	= Zones[sz_k].magz[2]
					sz_ei	= sz_si + Zones[sz_k].Br.shape[0]
					sz_ej	= sz_sj + Zones[sz_k].Br.shape[1]
					Zones[sz_k].Ballooning = []
					for iTransp in range(len(MagZones[k].Ballooning)):
						Zones[sz_k].Ballooning.append(MagZones[k].Ballooning[iTransp][sz_si:sz_ei, sz_sj:sz_ej])

					Zones[sz_k].IsTranspValues	= True

			MagZones[k].Ballooning 		= []
					
	pro_bar = 0
	return