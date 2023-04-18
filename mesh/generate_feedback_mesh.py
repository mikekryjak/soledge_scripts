import types
import numpy as np
import numpy.matlib as mat
from interfaces.progressbar			import ProgressBar
from routines.globals				import DEBUG, BALLOONING_NAMES

	
def generate_feedback_mesh(MainWindow, Config):

	if(DEBUG > 0): print("generate_feedback_mesh")
	
	pro_bar = ProgressBar(MainWindow, Label="Wait Generation")
	iProg = 0

	MagZones		= Config.MagZones
	MagMegazones	= Config.MagMegazones
	Zones			= Config.Zones
	FeedbackTransp	= Config.FeedbackTransp

	FeedbackTransp.Data.RefMin[0:2,0] = FeedbackTransp.Data.NRef			# Ne, Ni
	FeedbackTransp.Data.RefMin[0:2,1] = FeedbackTransp.Data.Nmin*0.1
	
	FeedbackTransp.Data.RefMin[2:4,0] = FeedbackTransp.Data.TRef			#Te & Ti
	FeedbackTransp.Data.RefMin[2:4,1] = FeedbackTransp.Data.Tmin

	FeedbackTransp.Data.RefMin[4:,0] = 1.									#D & Chi
	FeedbackTransp.Data.RefMin[4:,1] = FeedbackTransp.Data.Dmin
	
	for nm in range(len(MagMegazones)):
		isperiodic = MagMegazones[nm].isperiodic
		for k in MagMegazones[nm].list:
			Nx = MagZones[k].Nx
			Nz = MagZones[k].Nz
			
			x  = MagZones[k].x[:,0]
			
			MagZones[k].FeedbackValues    = []
			for i in range(len(FeedbackTransp.Cut.Flux.Profiles)):
				MagZones[k].FeedbackValues.append(interpolate_feedback(x, i, isperiodic, Nz, FeedbackTransp))

			MagZones[k].IsFeedbackValues = True
			
		iProg	+= 1
		pro_bar.Update(iProg/len(MagMegazones))
#		MainWindow.update()

	if(FeedbackTransp.Data.UsePattern != "No"):						#Ballooning definition based on 1/B
		if(not Config.transp_prof_OK):
			TranspCuts					= [types.SimpleNamespace()]
			TranspCuts[0].Flux			= types.SimpleNamespace()
			TranspCuts[0].Flux.Profiles	= []
			for i in range(len(BALLOONING_NAMES)):
				TranspCuts[0].Flux.Profiles.append(types.SimpleNamespace())
				TranspCuts[0].Flux.Profiles[-1].Name = BALLOONING_NAMES[i]

			Config.TranspCuts		= TranspCuts
			Config.transp_values_OK = True

		Bsum = 0.
		nTot = 0
		for k in range(len(MagZones)):								#Compute averge B
			Bsum += np.sum(np.sqrt(MagZones[k].Br**2 + MagZones[k].Bz**2 + MagZones[k].Bphi**2))
			nTot += MagZones[k].Br.size
		Bave = Bsum/nTot

		for k in range(len(MagZones)):								#Ballooning as Bave/B
			MagZones[k].Ballooning = []
			MagZones[k].Ballooning.append(Bave/np.sqrt(MagZones[k].Br**2 + MagZones[k].Bz**2 + MagZones[k].Bphi**2))

			for i in range(1,len(BALLOONING_NAMES)):			#the same for D, Mu, Chi
				MagZones[k].Ballooning.append(np.copy(MagZones[k].Ballooning[0]))

			for i in range(len(FeedbackTransp.Cut.Flux.Profiles)-2,len(FeedbackTransp.Cut.Flux.Profiles)):
				MagZones[k].FeedbackValues[i] *= MagZones[k].Ballooning[0]

							
#	Transfer data to splitted zones
	for k in range(len(MagZones)):								#Ballooning as Bave/B
		for sz in range(len(MagZones[k].list)):
			sz_k	= MagZones[k].list[sz]
			sz_si	= Zones[sz_k].magz[1]
			sz_sj	= Zones[sz_k].magz[2]
			sz_ei	= sz_si + Zones[sz_k].Nx
			sz_ej	= sz_sj + Zones[sz_k].Nz

			Zones[sz_k].FeedbackValues = []
			for iFeed in range(len(MagZones[k].FeedbackValues)):
				Zones[sz_k].FeedbackValues.append(MagZones[k].FeedbackValues[iFeed][sz_si:sz_ei, sz_sj:sz_ej])
			Zones[sz_k].IsFeedbackValues = True

			if(FeedbackTransp.Data.UsePattern != "No"):						#Ballooning definition based on 1/B
				Zones[sz_k].Ballooning = []
				for iBall in range(len(MagZones[k].Ballooning)):
					Zones[sz_k].Ballooning.append(MagZones[k].Ballooning[iBall][sz_si:sz_ei, sz_sj:sz_ej])

		MagZones[k].FeedbackValues = []
		if(FeedbackTransp.Data.UsePattern != "No"): MagZones[k].Ballooning = []

	pro_bar = 0

	if(DEBUG > 0): print("generate_feedback_mesh: completed")

	return

def interpolate_feedback(x, iType, isperiodic, Nz, FeedbackTransp):

	if(isperiodic or (iType > 3)):
		Values =  np.interp(x, FeedbackTransp.Data.x,  FeedbackTransp.Data.Profiles[iType, :])
	
	else:
		index		= np.interp(x, FeedbackTransp.Data.x, FeedbackTransp.Data.index)
		index		= index.astype(np.int)
		
		iExp		= np.where(index <  FeedbackTransp.Data.IndexSep); iExp 	= iExp[0]
		iInterp		= np.where(index >= FeedbackTransp.Data.IndexSep); iInterp = iInterp[0]
		
		Values		= np.empty(len(x), dtype='f8')
		if(len(iExp) > 0   ): Values[iExp]	 = FeedbackTransp.Data.Profiles[iType, FeedbackTransp.Data.IndexSep]* \
												np.exp((FeedbackTransp.Data.input_r[index[iExp]] - FeedbackTransp.Data.input_r[FeedbackTransp.Data.IndexSep])/((5./2.)*FeedbackTransp.Data.Lambda))
		if(len(iInterp) > 0): Values[iInterp] = np.interp(x[iInterp], FeedbackTransp.Data.x,  FeedbackTransp.Data.Profiles[iType, :])

	Values  = np.where(Values > FeedbackTransp.Data.RefMin[iType, 1], Values, FeedbackTransp.Data.RefMin[iType, 1])
	Values /=  FeedbackTransp.Data.RefMin[iType, 0]
	
	return mat.repmat(Values,Nz,1).T