import numpy as np
from interfaces.plot_and_ask	import plot_and_ask
from eirene.polygon_routines 	import polygon_count_in_out
from routines.globals			import DEBUG
# Function definition is here
#=======================================================================================================================================

# CAVEAT!!!!!
# This function assumes that no cell will be crossed by several walls at
# the same time

# First lap to count number of cuts between each cell and a wall =>
# zonesE.zone(k).iscrossed(i,j,:,:) stored the wall number crossed by each
# segment of each quadrangle (0 if no crossing)

#=========================================================
# This routine detect cut
#=========================================================
#

def detect_cut(Root, Config, Zones, RKnots, ZKnots):

	if(DEBUG > 0):	print("detect_cut")
	
	nres = 10																# raffinement
	nZones = len(Zones)
	for k in range(nZones):
		for iPwall in range(len(Config.iPwalls)):
			Wall = Config.Walls[Config.iPwalls[iPwall]]
			WallPath = Wall.WallPath

			IsCrossed = np.sum(Zones[k].IsCrossed[:,:,:,iPwall], axis=2)

			ii, jj = np.where(IsCrossed >= 3)
			if((DEBUG > 1) and (len(ii) > 0)):	print("\traffinement on ",len(ii)," points on zone ",k+1)
				
			for l in range(len(ii)):
				i = ii[l]
				j = jj[l]
		
				iKnodeA = Zones[k].KnotA[i,j]
				iKnodeB = Zones[k].KnotB[i,j]
				iKnodeC = Zones[k].KnotC[i,j]
				iKnodeD = Zones[k].KnotD[i,j]
		
		#		segment1
				if(polygon_count_in_out( RKnots[iKnodeA], ZKnots[iKnodeA], RKnots[iKnodeB], ZKnots[iKnodeB], nres, WallPath) == 1):
					Zones[k].IsCrossed[i,j,0,iPwall] = 1
				else:
					Zones[k].IsCrossed[i,j,0,iPwall] = 0
					
		#		segment2
				if(polygon_count_in_out( RKnots[iKnodeB], ZKnots[iKnodeB], RKnots[iKnodeC], ZKnots[iKnodeC], nres, WallPath) == 1):
					Zones[k].IsCrossed[i,j,1,iPwall] = 1
				else:
					Zones[k].IsCrossed[i,j,1,iPwall] = 0

		
		#		segment3
				if(polygon_count_in_out( RKnots[iKnodeC], ZKnots[iKnodeC], RKnots[iKnodeD], ZKnots[iKnodeD], nres, WallPath) == 1):
					Zones[k].IsCrossed[i,j,2,iPwall] = 1
				else:
					Zones[k].IsCrossed[i,j,2,iPwall] = 0

		
		#		segment4
				if(polygon_count_in_out( RKnots[iKnodeD], ZKnots[iKnodeD], RKnots[iKnodeA], ZKnots[iKnodeA], nres, WallPath) == 1):
					Zones[k].IsCrossed[i,j,2,iPwall] = 1
				else:
					Zones[k].IsCrossed[i,j,2,iPwall] = 0
				

				if(np.sum(Zones[k].IsCrossed[i,j,:,iPwall]) != 2):
	#				manual
					Xarrs  = [Wall.Rwall, np.array([RKnots[iKnodeA],RKnots[iKnodeB],RKnots[iKnodeC],RKnots[iKnodeD],RKnots[iKnodeA]])]
					Yarrs  = [Wall.Zwall, np.array([ZKnots[iKnodeA],ZKnots[iKnodeB],ZKnots[iKnodeC],ZKnots[iKnodeD],ZKnots[iKnodeA]])]
					lines  = ["k-","b-"]
					LinesData = [Xarrs, Yarrs, lines]
					
					xPos   = np.array([RKnots[iKnodeA],RKnots[iKnodeB],RKnots[iKnodeC],RKnots[iKnodeD]])
					yPos   = np.array([ZKnots[iKnodeA],ZKnots[iKnodeB],ZKnots[iKnodeC],ZKnots[iKnodeD]])
					Texts  = ["A","B","C","D"]
					colors = ["g","g","g","g"]
					sizes  = np.ones(len(Texts), dtype="i4")*18

					TextsData = [xPos, yPos, Texts, colors, sizes]
					AskTexts =  ["Segment 1: [AB=1], [BC=2], [CD=3], [AD=4]?", "Segment 2: [AB=1], [BC=2], [CD=3], [AD=4]?"]
					choice = plot_and_ask(Root, LinesData=LinesData, TextsData=TextsData, AskTexts=AskTexts, title="Select segments to cut")

					if(choice[0] > 0):
						choice1 = max(1,min(4,eval(choice[1])))
						choice2 = max(1,min(4,eval(choice[2])))

						Zones[k].IsCrossed[i,j,:,iPwall] = 0
						Zones[k].IsCrossed[i,j,choice1-1,iPwall] = 1
						Zones[k].IsCrossed[i,j,choice2-1,iPwall] = 1


	if(DEBUG > 0):	print("detect_cut: Completed")
