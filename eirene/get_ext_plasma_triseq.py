import numpy as np

#	This routine find the Eirene Wall Triangles sequence corresponding to the external plasma wall 

def get_ext_plasma_triseq(Config, EireneWall):

	CWalls		 = Config.Walls
	EWalls		 = EireneWall.EWalls
	TriSequences = EireneWall.TriSequences
	R12			 = EireneWall.R12
	Z12			 = EireneWall.Z12

#	Find EIRENE wall corresponding to external plasma wall
#	================================

	if(len(TriSequences) > 1):									#more than one plasma wall, find the oen correspondind to the external one
		CWall = CWalls[Config.iExtWalls[0]]						#external plasma wall
		Rwall = CWall.Rwall
		Zwall = CWall.Zwall	
		iRZMaxMin = np.array([np.argmax(Rwall),np.argmin(Rwall),np.argmax(Zwall), np.argmin(Zwall)])
		DistTot = np.zeros(len(TriSequences), dtype='f8')
		for iEWall in range(len(TriSequences)):
			for i in iRZMaxMin:
				DistTot[iEWall] += np.min((R12[iEWall][:,0]-Rwall[i])**2 + (Z12[iEWall][:,0]-Zwall[i])**2)

		EireneWall.iTriSeqExtPlasma = np.argmin(DistTot)												#Eirene plasma walll
	else:
		EireneWall.iTriSeqExtPlasma = 0

	return EireneWall.iTriSeqExtPlasma