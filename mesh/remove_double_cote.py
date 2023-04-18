import types
import numpy as np
from routines.globals				import DEBUG

def remove_double_cote(Mesher, Tol = 1.e-10, Thresh = 1.e-10):

	if(DEBUG > 0): print("\nremove_double_cote")
	
	OldCotes	= Mesher.Cotes
	NewCotes	= []
	m			= 0
	for nc in range(len(OldCotes)):
		inCote	= False
		for k in range(m):
			difftest1 = np.sum(np.abs(list[k,:] - np.append(OldCotes[nc].p1.coord[[0,1,3]],OldCotes[nc].p2.coord[[0,1,3]])))
			difftest2 = np.sum(np.abs(list[k,:] - np.append(OldCotes[nc].p2.coord[[0,1,3]],OldCotes[nc].p1.coord[[0,1,3]])))		
			if((difftest1 < 6*Thresh) or (difftest2 < 6*Thresh)):
				d = np.empty(len(OldCotes[nc].R))
				for k2 in range(len(OldCotes[nc].R)):
					d[k2] =	np.min(np.sqrt((NewCotes[k].R - OldCotes[nc].R[k2])**2 + (NewCotes[k].Z-OldCotes[nc].Z[k2])**2))

				dm = np.sum(np.where(d < Tol, 1, 0))
				if(len(d) == 2):
					if(dm == 2): 
						inCote = True
						break
				elif(dm > 2):
					inCote = True 							#more than 2 points overlap
					break

		d = np.sum(np.sqrt((OldCotes[nc].R[1:]-OldCotes[nc].R[:-1])**2+(OldCotes[nc].Z[1:]-OldCotes[nc].Z[:-1])**2)) 
		if(d == 0): zerol = True							#remove zero lenght cotes
		else:		zerol = False

		coord12 = np.array([np.append(OldCotes[nc].p1.coord[[0,1,3]],OldCotes[nc].p2.coord[[0,1,3]])])
		if((not inCote) and (not zerol)):
			if(m == 0): list = coord12
			else:		list = np.append(list, coord12, axis=0)
			m	= list.shape[0]
			p	= list.shape[1]
			
			NewCotes.append(types.SimpleNamespace())
			NewCotes[-1].p1		= OldCotes[nc].p1
			NewCotes[-1].p2		= OldCotes[nc].p2
			NewCotes[-1].R		= OldCotes[nc].R
			NewCotes[-1].Z		= OldCotes[nc].Z
			NewCotes[-1].iCote	= OldCotes[nc].iCote
			if(DEBUG > 1):	
				print("\tncotes=             ",len(NewCotes))
				print("\tnc=                 ",nc+1)
				print("\tCotes[-1].p1.coord= ",np.array(NewCotes[-1].p1.coord) + [0,1,1,0])
				print("\tCotes[-1].p2.coord= ",np.array(NewCotes[-1].p2.coord) + [0,1,1,0])

	list = np.array([], dtype='i4')
	for nc in range(len(NewCotes)):
		if(((NewCotes[nc].p1.coord[3]==1) and (NewCotes[nc].p1.coord[2]==-1)) and \
		   ((NewCotes[nc].p2.coord[3]==1) and (NewCotes[nc].p2.coord[2]==-1))): list = np.append(list,nc)			#boucle

#	dist = np.zeros(len(list), dtype='f8')
	if(len(list) > 2):
		print('should look at remove_double_cote.py - might be a problem')


	Mesher.Cotes = NewCotes
	
	if(DEBUG > 0): print("remove_double_cote: completed\n")

	return True

