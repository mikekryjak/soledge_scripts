import types
import numpy as np
from tkinter					import messagebox
from tkinter.messagebox 		import Message
from mesh.reverse_megazone		import reverse_megazone
from routines.globals			import DEBUG

def analyze_topology(self, Config):

	if(DEBUG > 0): print("analyze_topology")

	MagZones  = Config.MagZones
	remaining = np.arange(len(MagZones), dtype='i4')

	print("\tlen(MagZones)=",len(MagZones))

	iMagMegazone  = 0
	MagMegazones  = []
	thresh		  = 1.e-10;
	while(len(remaining)  >  0):

		MagMegazones.append(types.SimpleNamespace())
		iMagMegazone		+= 1
		ref					= remaining[0]
		MagMegazones[-1].list	= np.array([ref], dtype='i4')

#		look east
		cur  = ref
		num  = -1
		type = MagZones[cur].pC.coord[3]
		while((type != 2) and (num != ref)):
			east = np.append(MagZones[cur].pD.coord, MagZones[cur].pC.coord)
			for k in range(len(remaining)):
				num   = remaining[k]
				west  = np.append(MagZones[num].pA.coord, MagZones[num].pB.coord)
				somme = 0
				for i in range(8):
					if(abs(west[i] - east[i]) < thresh):					  somme += 1
					elif((i == 2) or (i == 6)):
						if(((west[i] % 1) == 0.5) or ((east[i] % 1) == 0.5)): somme += 1			# joker

				if(somme == 8):
					MagMegazones[-1].list = np.append(MagMegazones[-1].list, num)
					cur  = num
					type = MagZones[num].pC.coord[3]
					if(DEBUG > 1):	
						print("\tLook east: direct zone {:d} add to MagMegazones[{:d}].list ".format(num+1,len(MagMegazones)),MagMegazones[-1].list+1)					
					break

#				Reverse search

				if(num != cur):
					west  = np.append(MagZones[num].pD.coord, MagZones[num].pC.coord)
					somme = 0
					for i in range(8):
						if(abs(west[i] - east[i]) < thresh):						  somme += 1
						elif((i == 2) or (i == 6)):
								if(((west[i] % 1) == 0.5) or ((east[i] % 1) == 0.5)): somme += 1	# joker

					if(somme == 8):
						MagMegazones[-1].list	= np.append(MagMegazones[-1].list, num)
						cur						= num
						type					= MagZones[num].pB.coord[3]
#						%return zone
						temp_R					= MagZones[num].east.R
						temp_Z					= MagZones[num].east.Z
						MagZones[num].east.R	= MagZones[num].west.R
						MagZones[num].east.Z	= MagZones[num].west.Z
						MagZones[num].west.R	= temp_R
						MagZones[num].west.Z	= temp_Z
						MagZones[num].north.R	= MagZones[num].north.R[::-1]
						MagZones[num].north.Z	= MagZones[num].north.Z[::-1]
						MagZones[num].south.R	= MagZones[num].south.R[::-1]
						MagZones[num].south.Z	= MagZones[num].south.Z[::-1]

						pA						= MagZones[num].pA
						pB 						= MagZones[num].pB
						MagZones[num].pA		= MagZones[num].pD
						MagZones[num].pB		= MagZones[num].pC
						MagZones[num].pC 		= pB
						MagZones[num].pD		= pA
						if(DEBUG > 1):	
							print("\tLook east: reverted zone {:d} add to MagMegazones[{:d}].list ".format(num+1,len(MagMegazones)),MagMegazones[-1].list+1)					
						break
			if(somme != 8):
				print("\tError in analyze_topology: somme=", somme, " it should be 8!!! ")
				print("Megazone        =", iMagMegazone+1)
				print("Zone            =", cur+1)
				print("Remaining Zones =", remaining+1)
				messagebox.showwarning("analyze_topology", "Error somme={:d} it should be 8!!!\nMegazone={:d}\nZone={:d}\n".format(somme,iMagMegazone+1,cur+1))
				return False
#				exit()
#		%     if(([MagMegazones[-1.list[-1] == MagMegazones[-1].list[0]) and (
#		%             (~MagZones[MagMegazones[-1].list[0]].pA.coord[3] == 2))
#		%         MagMegazones[-1].isperiodic = true
#		%         MagMegazones[-1].list = MagMegazones[-1].list[1:end-1)
#		%
#		%     else
#		%         MagMegazones[-1].isperiodic = false

#		%look west
		cur  = ref
		num  = -1
		type = MagZones[cur].pB.coord[3]
		while((type != 2) and (num != ref)):
			west = np.append(MagZones[cur].pA.coord, MagZones[cur].pB.coord)
			for k in range(len(remaining)):
				num   = remaining[k]
				east  = np.append(MagZones[num].pD.coord, MagZones[num].pC.coord)
				somme = 0
				for i in range(8):
					if(abs(west[i] - east[i]) < thresh):					  somme += 1
					elif((i == 2) or (i == 6)):
						if(((west[i] % 1) == 0.5) or ((east[i] % 1) == 0.5)): somme += 1	# joker

				if(somme == 8):
					MagMegazones[-1].list = np.append(num, MagMegazones[-1].list)
					cur  = num
					type = MagZones[num].pB.coord[3]
					if(DEBUG > 1):	
						print("\tLook west: direct zone {:d} add to MagMegazones[{:d}].list ".format(num+1,len(MagMegazones)),MagMegazones[-1].list+1)
					break

				if(num != cur):
					east  = np.append(MagZones[num].pA.coord, MagZones[num].pB.coord)
					somme = 0
					for i in range(8):
						if(abs(west[i] - east[i]) < thresh):					  somme += 1
						elif((i == 2) or (i == 6)):
							if(((west[i] % 1) == 0.5) or ((east[i] % 1) == 0.5)): somme += 1	# joker

					if(somme == 8):
						MagMegazones[-1].list = np.append(num, MagMegazones[-1].list)
						cur 					= num
						type					= MagZones[num].pC.coord[3]							#python indexes
#						%return zone
						temp_R					= MagZones[num].east.R
						temp_Z					= MagZones[num].east.Z
						MagZones[num].east.R	= MagZones[num].west.R
						MagZones[num].east.Z	= MagZones[num].west.Z
						MagZones[num].west.R	= temp_R
						MagZones[num].west.Z	= temp_Z
						MagZones[num].north.R	= MagZones[num].north.R[::-1]
						MagZones[num].north.Z	= MagZones[num].north.Z[::-1]
						MagZones[num].south.R	= MagZones[num].south.R[::-1]
						MagZones[num].south.Z	= MagZones[num].south.Z[::-1]
						pA						= MagZones[num].pA
						pB						= MagZones[num].pB
						MagZones[num].pA		= MagZones[num].pD
						MagZones[num].pB		= MagZones[num].pC
						MagZones[num].pC		= pB
						MagZones[num].pD		= pA
						if(DEBUG > 1):	
							print("Look west: reverting zone {:d} add to MagMegazones[{:d}].list ".format(num+1,len(MagMegazones)),MagMegazones[-1].list+1)
						break
			if(somme != 8):
				print("Error in analyze_topology: somme=", somme, " it should be 8!!! ")
				print("Megazone        =", iMagMegazone+1)
				print("Zone            =", cur+1)
				print("Remaining Zones =", remaining+1)
				messagebox.showwarning("analyze_topology", "Error somme={:d} it should be 8!!!\nMegazone={:d}\nZone={:d}\n".format(somme,iMagMegazone+1,cur+1))
				return False

#		remove from remaining

		remaining  = np.setdiff1d(remaining, MagMegazones[-1].list)

	for k in range(len(MagMegazones)):
		MagMegazones[k].isperiodic = False
		ref  = MagMegazones[k].list[0]
		list = np.array([ref],dtype='i4')
		for k2 in range(1,len(MagMegazones[k].list)):
			if(MagMegazones[k].list[k2] != ref):
				list = np.append(list, MagMegazones[k].list[k2])
			else:
				MagMegazones[k].isperiodic = True
				break
		MagMegazones[k].list = list
		for nz in range(len(MagMegazones[k].list)):
			MagZones[MagMegazones[k].list[nz]].mz = k
		
		if(DEBUG > 1):	
			print("\tMagMegazones[{:d}].isperiodic  ".format(k+1), MagMegazones[k].isperiodic)
			print("\tMagMegazones[{:d}].list        ".format(k+1), MagMegazones[k].list + 1)

#	%megazone1 used as a reference for east west direction
	done = np.array([], dtype='i4')
	todo = np.array([2], dtype='i4')
	while(len(done) < len(MagMegazones)):

		ref = todo[0]
		if(DEBUG > 1):
			print("\tDef Magzones: MagMegazones[{:d}].list=".format(ref+1),MagMegazones[ref].list+1)

		for k in range(len(MagMegazones[ref].list)):
			MagZones[MagMegazones[ref].list[k]].Neighbour	    = types.SimpleNamespace()
			MagZones[MagMegazones[ref].list[k]].Neighbour.north = -1
			MagZones[MagMegazones[ref].list[k]].Neighbour.south = -1

#			%north scan
			north= np.append(MagZones[MagMegazones[ref].list[k]].pB.coord, MagZones[MagMegazones[ref].list[k]].pC.coord)
			for k2 in range(len(MagZones)):
				south = np.append(MagZones[k2].pA.coord, MagZones[k2].pD.coord)

				somme = 0
				for i in range(8):
					if(north[i] == south[i]): somme += 1
#					elif((i == 2) or (i == 6)):								#python indexes
#						if((west[i] == -2) or (east[i] == -2)): somme += 1		#By P. Innocente removed because -2 replaced by fractionary values
#									%joker

				if(somme == 8):
#					%check sens
					sum = 0
					vx = MagZones[MagMegazones[ref].list[k]].north.R[-1] - MagZones[MagMegazones[ref].list[k]].north.R[0]
					vy = MagZones[MagMegazones[ref].list[k]].north.Z[-1] - MagZones[MagMegazones[ref].list[k]].north.Z[0]
					for k5 in range(len(MagZones[MagMegazones[ref].list[k]].north.R)):
						vx1 = MagZones[MagMegazones[ref].list[k]].north.R[k5] - MagZones[MagMegazones[ref].list[k]].north.R[0]
						vy1 = MagZones[MagMegazones[ref].list[k]].north.Z[k5] - MagZones[MagMegazones[ref].list[k]].north.Z[0]
						sum = sum+vx*vy1-vy*vx1

					sum1 = 0
					vx = MagZones[k2].south.R[-1] - MagZones[k2].south.R[0]
					vy = MagZones[k2].south.Z[-1] - MagZones[k2].south.Z[0]
					for k5 in range(len(MagZones[k2].south.R)):
						vx1  = MagZones[k2].south.R[k5]- MagZones[k2].south.R[0]
						vy1  = MagZones[k2].south.Z[k5]- MagZones[k2].south.Z[0]
						sum1 = sum1+vx*vy1-vy*vx1

					if(sum1*sum >= 0):
						MagZones[MagMegazones[ref].list[k]].Neighbour.north = k2
						ii = np.where(done == MagZones[k2].mz)[0]
						if(len(ii) == 0): todo = np.append(todo, MagZones[k2].mz)
						break

#				%reverse
				south = np.append(MagZones[k2].pD.coord, MagZones[k2].pA.coord)
				somme = 0
				for i in range(8):
					if(north[i] == south[i]): somme += 1
#					elif((i == 2) or (i == 6)):											#python indexes
#						if((west[i] == -2) or (east[i] == -2)): somme += 1					#By P. Innocente removed because meaningless
#									%joker

				if(somme == 8):
#					check sens
					sum = 0
					vx  = MagZones[MagMegazones[ref].list[k]].north.R[-1] - MagZones[MagMegazones[ref].list[k]].north.R[0]
					vy  = MagZones[MagMegazones[ref].list[k]].north.Z[-1] - MagZones[MagMegazones[ref].list[k]].north.Z[0]
					for k5 in range(len(MagZones[MagMegazones[ref].list[k]].north.R)):
						vx1 = MagZones[MagMegazones[ref].list[k]].north.R[k5] - MagZones[MagMegazones[ref].list[k]].north.R[0]
						vy1 = MagZones[MagMegazones[ref].list[k]].north.Z[k5] - MagZones[MagMegazones[ref].list[k]].north.Z[0]
						sum = sum+vx*vy1-vy*vx1

					sum1 = 0
					vx = MagZones[k2].south.R[0] - MagZones[k2].south.R[-1]
					vy = MagZones[k2].south.Z[0] - MagZones[k2].south.Z[-1]
					for k5 in range(len(MagZones[k2].south.R)):
						vx1  = MagZones[k2].south.R[k5] - MagZones[k2].south.R[-1]
						vy1  = MagZones[k2].south.Z[k5] - MagZones[k2].south.Z[-1]
						sum1 = sum1+vx*vy1-vy*vx1
					if(sum1*sum >= 0):
						MagZones[MagMegazones[ref].list[k]].Neighbour.north = k2
						reverse_megazone(MagZones, MagMegazones, MagZones[k2].mz)
						ii = np.where(done == MagZones[k2].mz)[0]
						if(len(ii) == 0): 	todo = np.append(todo, MagZones[k2].mz)
						break


#			%south scan
			south = np.append(MagZones[MagMegazones[ref].list[k]].pA.coord, MagZones[MagMegazones[ref].list[k]].pD.coord)
			for k2 in range(len(MagZones)):
				north = np.append(MagZones[k2].pB.coord, MagZones[k2].pC.coord)
				somme = 0
				for i in range(8):
					if(north[i] == south[i]): somme += 1
#					elif((i == 2) or (i == 6)):											#python indexes
#							if((north[i] == -2) or (north[i] == -2)): somme += 1				#By P. Innocente removed because meaningless
#									%joker

				if(somme == 8):
#					%check sens
					sum = 0
					vx  = MagZones[MagMegazones[ref].list[k]].south.R[-1] - MagZones[MagMegazones[ref].list[k]].south.R[0]
					vy  = MagZones[MagMegazones[ref].list[k]].south.Z[-1] - MagZones[MagMegazones[ref].list[k]].south.Z[0]
					for k5 in range(len(MagZones[MagMegazones[ref].list[k]].south.R)):
						vx1 = MagZones[MagMegazones[ref].list[k]].south.R[k5] - MagZones[MagMegazones[ref].list[k]].south.R[0]
						vy1 = MagZones[MagMegazones[ref].list[k]].south.Z[k5] - MagZones[MagMegazones[ref].list[k]].south.Z[0]
						sum = sum+vx*vy1-vy*vx1

					sum1	= 0
					vx		= MagZones[k2].north.R[-1] - MagZones[k2].north.R[0]
					vy		= MagZones[k2].north.Z[-1] - MagZones[k2].north.Z[0]
					for k5 in range(len(MagZones[k2].north.R)):
						vx1  = MagZones[k2].north.R[k5] - MagZones[k2].north.R[0]
						vy1  = MagZones[k2].north.Z[k5] - MagZones[k2].north.Z[0]
						sum1 = sum1+vx*vy1-vy*vx1

					if(sum1*sum >= 0):
						MagZones[MagMegazones[ref].list[k]].Neighbour.south = k2
						ii = np.where(done == MagZones[k2].mz)[0]
						if(len(ii) == 0): 	todo = np.append(todo, MagZones[k2].mz)
						break

#				reverse
				north = np.append(MagZones[k2].pC.coord, MagZones[k2].pB.coord)
				somme = 0
				for i in range(8):
					if(north[i] == south[i]): somme += 1
#					elif((i == 2) or (i == 6)):												#python indexes
#						if((north[i] == -2) or (south[i] == -2)):	somme += 1					#By P. Innocente removed because meaningless
#									joker
				if(somme == 8):
#					check sens
					sum = 0
					vx  = MagZones[MagMegazones[ref].list[k]].south.R[-1] - MagZones[MagMegazones[ref].list[k]].south.R[0]
					vy  = MagZones[MagMegazones[ref].list[k]].south.Z[-1] - MagZones[MagMegazones[ref].list[k]].south.Z[0]
					for k5 in range(len(MagZones[MagMegazones[ref].list[k]].south.R)):
						vx1 = MagZones[MagMegazones[ref].list[k]].south.R[k5] - MagZones[MagMegazones[ref].list[k]].south.R[0]
						vy1 = MagZones[MagMegazones[ref].list[k]].south.Z[k5] - MagZones[MagMegazones[ref].list[k]].south.Z[0]
						sum = sum+vx*vy1-vy*vx1
					sum1 = 0
					vx   = MagZones[k2].north.R[0] - MagZones[k2].north.R[-1]
					vy   = MagZones[k2].north.Z[0] - MagZones[k2].north.Z[-1]
					for k5 in range(len(MagZones[k2].north.R)):
						vx1  = MagZones[k2].north.R[k5] - MagZones[k2].north.R[-1]
						vy1  = MagZones[k2].north.Z[k5] - MagZones[k2].north.Z[-1]
						sum1 = sum1+vx*vy1-vy*vx1
					if(sum1*sum >= 0):
						MagZones[MagMegazones[ref].list[k]].Neighbour.south = k2
						reverse_megazone(MagZones, MagMegazones, MagZones[k2].mz)
						ii = np.where(done == MagZones[k2].mz)[0]
						if(len(ii) == 0): 	todo = np.append(todo, MagZones[k2].mz)
						break
		done  = np.append(done, ref)

#		remove from todo

		todo  = np.setdiff1d(todo, done)

	for k in range(len(MagMegazones)):
		for nz in range(len(MagMegazones[k].list)-1):
			MagZones[MagMegazones[k].list[nz  ]].Neighbour.east = MagMegazones[k].list[nz+1]
			MagZones[MagMegazones[k].list[nz+1]].Neighbour.west = MagMegazones[k].list[nz]
		if(MagMegazones[k].isperiodic):
			MagZones[MagMegazones[k].list[0 ]].Neighbour.west = MagMegazones[k].list[-1]
			MagZones[MagMegazones[k].list[-1]].Neighbour.east = MagMegazones[k].list[0]
		else:
			MagZones[MagMegazones[k].list[0 ]].Neighbour.west = -1
			MagZones[MagMegazones[k].list[-1]].Neighbour.east = -1

#	weak check
	for k in range(len(MagZones)):

		if(MagZones[k].Neighbour.north  >  -1):
			N = MagZones[k].Neighbour.north
			if (MagZones[N].Neighbour.south != k):
				print("\nanalyze_topology: Problem with neighbors: Zone ",k)
#				exit()

		if(MagZones[k].Neighbour.south  >  -1):
			S = MagZones[k].Neighbour.south
			if (MagZones[S].Neighbour.north != k):
				print("\nanalyze_topology:Problem with neighbors: Zone ",k)
#				exit()

		if(MagZones[k].Neighbour.east  >  -1):
			E = MagZones[k].Neighbour.east
			if (MagZones[E].Neighbour.west != k):
				print("\nanalyze_topology:Problem with neighbors: Zone ",k)
#				exit()

		if(MagZones[k].Neighbour.west  >  -1):
			W = MagZones[k].Neighbour.west
			if (MagZones[W].Neighbour.east != k):
				print("\nanalyze_topology:Problem with neighbors: Zone ",k)
#				exit()

	for k in range(len(MagMegazones)): MagMegazones[k].ismeshed = False

	iPMegazone = 0
	MagPMegazones  = []
	for k in range(len(MagZones)):
		"""
		print("MagZone = ",k+1)
		print("south.R = ",MagZones[k].south.R)
		print("south.Z = ",MagZones[k].south.Z)
		print("west.R  = ",MagZones[k].west.R)
		print("west.Z  = ",MagZones[k].west.Z)
		print("north.R = ",MagZones[k].north.R)
		print("north.Z = ",MagZones[k].north.Z)
		print("east.R  = ",MagZones[k].east.R)
		print("east.Z  = ",MagZones[k].east.Z)
		"""

		if(MagZones[k].Neighbour.south == -1):
			MagPMegazones.append(types.SimpleNamespace())
			MagPMegazones[-1].list = np.array([k], dtype='i4')
			MagZones[k].pmz		= iPMegazone
			num					= k
			while(MagZones[num].Neighbour.north != -1):
				num					= MagZones[num].Neighbour.north
				MagZones[num].pmz		= iPMegazone
				MagPMegazones[-1].list = np.append(MagPMegazones[-1].list, num)
				
			iPMegazone 			+= 1
			if(DEBUG > 1):
				print("\tPMegazone              ", len(MagPMegazones))
				print("\tMagPMegazones[-1].list ", MagPMegazones[-1].list+1)

	for k in range(len(MagPMegazones)):
		MagPMegazones[k].ismeshed 			= False
		MagPMegazones[k].isaligned 			= False
		MagPMegazones[k].subrefpoints		= [types.SimpleNamespace(), types.SimpleNamespace()]
		MagPMegazones[k].subrefpoints[0].R = np.array([])
		MagPMegazones[k].subrefpoints[0].Z = np.array([])
		MagPMegazones[k].subrefpoints[1].R = np.array([])
		MagPMegazones[k].subrefpoints[1].Z = np.array([])

	Config.MagMegazones  = MagMegazones
	Config.MagPMegazones = MagPMegazones

	print("\tlen(MagMegazones)=",len(MagMegazones))
	print("\tlen(MagPMegazones)=",len(MagPMegazones))
	
	if(DEBUG > 0): print("\nanalyze_topology: completed\n")

	return True
