import types
import numpy			as np
from math				import floor, sqrt
from tkinter			import messagebox
from routines.globals	import DEBUG

#	pA/pB/pC/pD.coord		= [psi, theta, branch, type]	see find_cotes.py

def find_zones(self, Config):

	if(DEBUG > 0): print("\nfind_zones\n")
	
	X_points	= Config.X_points
	Cotes    	= Config.Cotes
	Segments	= Config.Segments
	nMagZones		= 0
	MagZones 		= []

	nCotes = len(Cotes)
	Nseg   = len(Segments)
	thresh = 1.e-10
	for k in range(nCotes):
		p1 = Cotes[k].p1.coord
		p2 = Cotes[k].p2.coord
		comp1 = np.array([p1[1],p1[2],p1[3],p2[1],p2[2],p2[3]])
		list1 = np.array([], dtype='i4')
		list2 = np.array([], dtype='i4')
		for k2 in range(nCotes):
			if(k2 != k):
				p1_ 	= Cotes[k2].p1.coord
				p2_		= Cotes[k2].p2.coord
				comp2	= np.array([p1_[1],p1_[2],p1_[3],p2_[1],p2_[2],p2_[3]])
				somme	= 0
				for k3 in range(6):
					if(comp1[k3] == comp2[k3]):											somme += 1
					else:
						if((k3 == 1) or (k3 == 4)):
							if(((comp1[k3] % 1) == 0.5) or ((comp2[k3] % 1) == 0.5)): 	somme += 1				#joker

				if((somme == 6) and (len(np.where(list2==k2)[0]) == 0)): list1 = np.append(list1, k2)

				comp2 = np.array([p2_[1],p2_[2],p2_[3],p1_[1],p1_[2],p1_[3]])
				somme = 0
				for k3 in range(6):
					if(comp1[k3] == comp2[k3]):											somme += 1
					else:
						if((k3 == 1) or (k3 == 4)):
							if(((comp1[k3] % 1) == 0.5) or ((comp2[k3] % 1) == 0.5)): 	somme += 1					#joker

				if((somme == 6) and (len(np.where(list1==k2)[0]) == 0)): list2 = np.append(list2, k2)

		sens = np.append(np.ones(len(list1), dtype='i4'),2*np.ones(len(list2), dtype='i4'))
		list1 = np.append(list1, list2)
		d = np.zeros(len(list1))
		for m in range(len(list1)): d[m] = Cotes[k].p1.coord[0] - Cotes[list1[m]].p1.coord[0]		

		a = np.where(d < 0)[0]
		
		if(len(a) > 0):
			
			list1	= list1[a]
			sens	= sens[a]
			d		= np.abs(d[a])
			b		= np.where(d == d.min())[0]
			list1	= list1[b]
			sens	= sens[b]

#			We have to select between the remaining cotes
			if (len(list1) > 1):
				dist = np.zeros(len(list1), dtype = 'f8')
				for iclist in range(len(list1)):
					for ipt in range(len(Cotes[k].R)):
						dist[iclist] += np.min((Cotes[k].R[ipt]-Cotes[list1[iclist]].R)**2+(Cotes[k].Z[ipt]-Cotes[list1[iclist]].Z)**2)

				choice	= np.where(dist == dist.min())[0]
				list1	= list1[choice]
				sens	= sens[choice]

#			Check sens for zones with one/two common point

			dist = sqrt((Cotes[k].R[-1]-Cotes[k].R[0])**2 + (Cotes[k].Z[-1]-Cotes[k].Z[0])**2)
			sum  = np.sum((Cotes[k].R[-2]-Cotes[k].R[1])*(Cotes[k].Z - Cotes[k].Z[1]) - (Cotes[k].Z[-2]-Cotes[k].Z[1])*(Cotes[k].R - Cotes[k].R[1]))
			for m in range(len(list1)):
				k2 	  = list1[m]
				dist1 = sqrt((Cotes[k2].R[-1]-Cotes[k2].R[0])**2 + (Cotes[k2].Z[-1]-Cotes[k2].Z[0])**2)
				if((dist < 1e-6) or (dist1 < 1e-6)):
					if(sum != 0):
						sum1 = np.sum((Cotes[k2].R[-2]-Cotes[k2].R[1])*(Cotes[k2].Z - Cotes[k2].Z[1]) - (Cotes[k2].Z[-2]-Cotes[k2].Z[1])*(Cotes[k2].R - Cotes[k2].R[1]))
						if(sum*sum1 > 0):	sens[m] = 1
						else:				sens[m] = 2
					else:
						direction =  (Cotes[k].R[-2]-Cotes[k].R[1])*(Cotes[k2].R[-2]-Cotes[k2].R[1]) + (Cotes[k].Z[-2]-Cotes[k].Z[1])*(Cotes[k2].Z[-2]-Cotes[k2].Z[1])
						if(direction > 0):	sens[m] = 1
						else:			 	sens[m] = 2

			for m in range(len(list1)):
	
				MagZones.append(types.SimpleNamespace())
				MagZones[-1].pA		= types.SimpleNamespace()
				MagZones[-1].pB		= types.SimpleNamespace()
				MagZones[-1].pC		= types.SimpleNamespace()
				MagZones[-1].pD		= types.SimpleNamespace()
				MagZones[-1].south		= types.SimpleNamespace()
				MagZones[-1].north		= types.SimpleNamespace()
				MagZones[-1].east		= types.SimpleNamespace()
				MagZones[-1].west		= types.SimpleNamespace()
					
				MagZones[-1].coord		= np.array([list1[m],k,-1,-1])		# [north, south, east, west]
				MagZones[-1].pA.coord	= Cotes[k].p1.coord
				MagZones[-1].pD.coord	= Cotes[k].p2.coord
				MagZones[-1].south.R	= Cotes[k].R
				MagZones[-1].south.Z	= Cotes[k].Z
				if(sens[m] == 1):
					MagZones[-1].pB.coord = Cotes[list1[m]].p1.coord
					MagZones[-1].pC.coord = Cotes[list1[m]].p2.coord
					MagZones[-1].north.R  = Cotes[list1[m]].R
					MagZones[-1].north.Z  = Cotes[list1[m]].Z
				else:
					MagZones[-1].pC.coord = Cotes[list1[m]].p1.coord
					MagZones[-1].pB.coord = Cotes[list1[m]].p2.coord
					MagZones[-1].north.R  = Cotes[list1[m]].R[::-1]
					MagZones[-1].north.Z  = Cotes[list1[m]].Z[::-1]
				if(DEBUG > 1):	
					print("nZone =           ", len(MagZones))
					print("MagZones[{:d}].pA.coord".format(len(MagZones)), np.array(MagZones[-1].pA.coord) + [0,1,1,0])
					print("MagZones[{:d}].pB.coord".format(len(MagZones)), np.array(MagZones[-1].pB.coord) + [0,1,1,0])
					print("MagZones[{:d}].pC.coord".format(len(MagZones)), np.array(MagZones[-1].pC.coord) + [0,1,1,0])
					print("MagZones[{:d}].pD.coord".format(len(MagZones)), np.array(MagZones[-1].pD.coord) + [0,1,1,0])

#					print("MagZones[{:d}].south.R".format(len(MagZones)),MagZones[-1].south.R)
#					print("MagZones[{:d}].south.Z".format(len(MagZones)),MagZones[-1].south.Z)
#					print("MagZones[{:d}].north.R".format(len(MagZones)),MagZones[-1].north.R)
#					print("MagZones[{:d}].north.Z".format(len(MagZones)),MagZones[-1].north.Z)


				R=(MagZones[-1].south.R[floor(MagZones[-1].south.R.shape[0]/2)]+MagZones[-1].north.R[floor(MagZones[-1].north.R.shape[0]/2)])*0.5
				Z=(MagZones[-1].south.Z[floor(MagZones[-1].south.Z.shape[0]/2)]+MagZones[-1].north.Z[floor(MagZones[-1].north.Z.shape[0]/2)])*0.5
				self.Ax.text(R,Z,"Z_{:d}".format(len(MagZones)), color='green', size = "x-small")
			
#				figure[0]
#				hold on
#				plot(MagZones[-1].south.R,MagZones[-1].south.Z,'r-','LineWidth',2)
#				plot(MagZones[-1].north.R,MagZones[-1].north.Z,'b--','LineWidth',2)
#				R=(MagZones[-1].south.R(floor(end/2))+MagZones[-1].north.R(floor(end/2)))*0.5
#				Z=(MagZones[-1].south.Z(floor(end/2))+MagZones[-1].north.Z(floor(end/2)))*0.5
#				text(R,Z,num2str(lend(MagZones)))
	nMagZones = len(MagZones)
	for k in range(nMagZones):
		FoundWest = False
		comp1 = np.append(MagZones[k].pA.coord,MagZones[k].pB.coord)									#complete west
		for n in range(Nseg):
			comp2 = np.array([Segments[n].psimin, Segments[n].nX, Segments[n].nB, Segments[n].type, 
							  Segments[n].psimax, Segments[n].nX, Segments[n].nB, Segments[n].type])
			somme = 0
			for k3 in range(8):
				if(abs(comp1[k3] - comp2[k3]) < thresh):							 somme += 1
				else:
					if((k3 == 2) or (k3 == 6)):
						if(((comp1[k3] % 1) == 0.5) or ((comp2[k3] % 1) == 0.5)): somme += 1			#joker


			if(somme == 8):
				if(FoundWest): print("\tError in find_zones: alrady found West for zone =",k+1)

				MagZones[k].coord[3] = n
				MagZones[k].west.R	  = Segments[n].R
				MagZones[k].west.Z	  = Segments[n].Z
				FoundWest = True

		if(not  FoundWest): print("\tError in find_zones:not found West for zone =",k+1)


		FoundEast = False
		comp1 = np.append(MagZones[k].pD.coord,MagZones[k].pC.coord)									#complete est
		for n in range(Nseg):
			comp2 =  np.array([Segments[n].psimin, Segments[n].nX, Segments[n].nB, Segments[n].type,
							   Segments[n].psimax, Segments[n].nX, Segments[n].nB, Segments[n].type])
			somme = 0
			for k3 in range(8):
				if(abs(comp1[k3] - comp2[k3]) < thresh):								somme += 1
				else:
					if((k3 == 2) or (k3 == 6)):
						if(((comp1[k3] % 1) == 0.5) or ((comp2[k3] % 1) == 0.5)): somme += 1		#joker

			if(somme == 8):
				if(FoundEast): print("\tError in find_zones: alrady found East for zone =",k+1)

				MagZones[k].coord[2] = n
				MagZones[k].east.R	  = Segments[n].R
				MagZones[k].east.Z   = Segments[n].Z
				FoundEast = True

		if(not  FoundEast): print("\tError in find_zones:not found East for zone =",k+1)
		if(DEBUG > 1):	
			print("\tnZone =   ", k+1)
			print("\tcoord     ", np.array(MagZones[k].coord) + [1,1,1,1])
			print("\tCote&Seg  ", np.array([Cotes[MagZones[k].coord[0]].iCote,Cotes[MagZones[k].coord[1]].iCote,Segments[MagZones[k].coord[2]].iSeg,Segments[MagZones[k].coord[3]].iSeg]) + [1,1,1,1])

#		%     plot(MagZones[k].west.R,MagZones[k].west.Z,'g-','LineWidth',2)
#		%     plot(MagZones[k].east.R,MagZones[k].east.Z,'k--','LineWidth',2)

#	complete cote with X-points if necessary

	"""
	for k in range(nMagZones):
		if(MagZones[k].pA.coord[3] == 1):															#Xpoint
			nX = int(MagZones[k].pA.coord[1])
			if(MagZones[k].pA.coord[0]== X_points[nX].psi):
				if((MagZones[k].south.R[0] != X_points[nX].R) or (MagZones[k].south.Z[0] != X_points[nX].Z)):
					MagZones[k].south.R = np.append(X_points[nX].R,MagZones[k].south.R)
					MagZones[k].south.Z = np.append(X_points[nX].Z,MagZones[k].south.Z)

		if(MagZones[k].pB.coord[3] == 1):															#Xpoint
			nX = int(MagZones[k].pB.coord[1])
			if(MagZones[k].pB.coord[0] == X_points[nX].psi):
				if((MagZones[k].north.R[0] != X_points[nX].R) or (MagZones[k].north.Z[0] != X_points[nX].Z)):
					MagZones[k].north.R = np.append(X_points[nX].R,MagZones[k].north.R)
					MagZones[k].north.Z = np.append(X_points[nX].Z,MagZones[k].north.Z)

		if(MagZones[k].pC.coord[3] == 1):															#Xpoint
			nX = int(MagZones[k].pC.coord[1])
			if(MagZones[k].pC.coord[0] == X_points[nX].psi):
				if((MagZones[k].north.R[-1] != X_points[nX].R) or (MagZones[k].north.Z[-1] != X_points[nX].Z)):
					MagZones[k].north.R = np.append(MagZones[k].north.R,X_points[nX].R)
					MagZones[k].north.Z = np.append(MagZones[k].north.Z,X_points[nX].Z)

		if(MagZones[k].pD.coord[3] == 1):															#Xpoint
			nX = int(MagZones[k].pD.coord[1])
			if(MagZones[k].pD.coord[0] == X_points[nX].psi):
				if((MagZones[k].south.R[-1] != X_points[nX].R) or (MagZones[k].south.Z[-1] != X_points[nX].Z)):
					MagZones[k].south.R = np.append(MagZones[k].south.R,X_points[nX].R)
					MagZones[k].south.Z = np.append(MagZones[k].south.Z,X_points[nX].Z)
	"""

	for k in range(nMagZones):
		MagZones[k].Xtype_east	   = Segments[MagZones[k].coord[2]].Xtype
		MagZones[k].Xtype_west	   = Segments[MagZones[k].coord[3]].Xtype
		MagZones[k].northaligned   = False
		MagZones[k].southaligned   = False
		MagZones[k].east.ismeshed  = False
		MagZones[k].west.ismeshed  = False
		MagZones[k].north.ismeshed = False
		MagZones[k].south.ismeshed = False
		MagZones[k].meshortho	   = False

	Config.MagZones = MagZones

	self.Fig.canvas.draw()
	
	if(DEBUG > 0): print("\nfind_zones: completed\n")

	return True

#	messagebox.showinfo("Grd-Gen", "Domain decomposed in {:d} zones".format(len(MagZones)))
