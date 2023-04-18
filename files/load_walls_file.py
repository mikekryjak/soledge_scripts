
import numpy as np
import os
import types
from routines.globals			import DEBUG, DRAWING_WALL

#=================================================
# 
#=================================================

def load_walls_file(WallFile):

	if(DEBUG > 0):	print("load_walls_file")

	FileExt = os.path.splitext(os.path.split(WallFile)[1])[1]

	if(FileExt == ".csv"):	return load_walls_file_csv(WallFile)
	else:					return load_walls_file_txt(WallFile)

def load_walls_file_txt(WallFile):
	fid = open(WallFile,'r')

	Walls = []
	for full_line in fid:
		line  = full_line[:-1]
		if(line.upper().find("WALL =") > -1): continue
		iType = line.upper().find("TYPE =")
		if(iType > -1):
			if(len(Walls) > 0):
				RZ_mat = np.array(RZ_mat, dtype='f8')
				Walls[-1].Rwall = np.copy(RZ_mat[:,0])
				Walls[-1].Zwall = np.copy(RZ_mat[:,1])
				if((Walls[-1].Rwall[0] == Walls[-1].Rwall[-1]) and (Walls[-1].Zwall[0] == Walls[-1].Zwall[-1])): Walls[-1].Closed = True
				else:																		 					 Walls[-1].Closed = False

			Walls.append(types.SimpleNamespace())
			Walls[-1].Changed = True
			RZ_mat = []

			iLineType = line.upper().find("LINETYPE")
			if(iLineType > -1):																					#Found LineType
				iEq = line[iType+4:iLineType].find("=")
				Walls[-1].Type = int(eval(line[iType+4+iEq+1:iLineType]))
				iEq = line[iLineType+8:].find("=")
				Walls[-1].LineType = np.array(eval(line[iLineType+8+iEq+1:]), dtype='i4')
			else:																								#Set default Linetype
				iEq = line[iType+4:].find("=")
				Walls[-1].Type = int(eval(line[iType+4+iEq+1:]))
				Walls[-1].LineType = np.array([0,0,0], dtype='i4')

		elif(len(Walls) == 0):
			Walls.append(types.SimpleNamespace())
			Walls[-1].Changed = True
			RZ_mat = []

			Walls[-1].Type = 1
			Walls[-1].LineType = np.array([0,0,0], dtype='i4')
		else:
			s_line = line.split(",")
			if(len(s_line) == 2): 
				RZ = eval(line)
			else:
				s_line = line.split()
				if(len(s_line) == 2):
					RZ = eval(s_line[0] + "," + s_line[1])
			RZ_mat.append(RZ)

	fid.close()

#	fix last data

	RZ_mat = np.array(RZ_mat, dtype='f8')
	Walls[-1].Rwall = np.copy(RZ_mat[:,0])
	Walls[-1].Zwall = np.copy(RZ_mat[:,1])
	if((Walls[-1].Rwall[0] == Walls[-1].Rwall[-1]) and (Walls[-1].Zwall[0] == Walls[-1].Zwall[-1])): Walls[-1].Closed = True
	else:																		 					 Walls[-1].Closed = False

	if(DEBUG > 0):	print("load_walls_file: Completed")
	
	return Walls


def load_walls_file_csv(WallFile):

	Walls = []
	try:
		file_data = np.loadtxt(WallFile, delimiter=",", skiprows = 0)
	except:
		file_data = np.loadtxt(WallFile, delimiter=",", skiprows = 1)

	if(file_data.shape[1] == 2):									#Single wall R,Z file
		Walls.append(types.SimpleNamespace())
		Walls[-1].Type = DRAWING_WALL
		Walls[-1].LineType = np.array([0,0,0], dtype='i4')
		Walls[-1].Changed = True
		Walls[-1].Rwall = np.copy(file_data[:,0])
		Walls[-1].Zwall = np.copy(file_data[:,1])
		if((Walls[-1].Rwall[0] == Walls[-1].Rwall[-1]) and (Walls[-1].Zwall[0] == Walls[-1].Zwall[-1])): Walls[-1].Closed = True
		else:																		 					 Walls[-1].Closed = False
		Walls[-1].Clockwise = 0
		Walls[-1].Changed	 = True

	elif(file_data.shape[1] == 4):									#Segments Rs,Zs,Re,Ze file
		MaxDist2 = 1e-8
		FreeSegs = np.arange(file_data.shape[0])
		while(len(FreeSegs) > 0):
			Walls.append(types.SimpleNamespace())
			Walls[-1].Type = DRAWING_WALL
			Walls[-1].LineType = np.array([0,0,0], dtype='i4')
			Walls[-1].Changed = True
			Rwall = np.array([file_data[FreeSegs[0],0], file_data[FreeSegs[0],2]])
			Zwall = np.array([file_data[FreeSegs[0],1], file_data[FreeSegs[0],3]])
			FreeSegs = FreeSegs[1:]
			for dir in range(2):																								#Two directions: Forward and reverser
				Found = True
				while (Found  and (len(FreeSegs) > 0)):
					Found = False
					d2		= (file_data[FreeSegs,0]-Rwall[-1])**2+(file_data[FreeSegs,1]-Zwall[-1])**2							#Check first point of segments
					id2Min	= np.argmin(d2)
					if(d2[id2Min] < MaxDist2):
						Rwall = np.append(Rwall, file_data[FreeSegs[id2Min],2])
						Zwall = np.append(Zwall, file_data[FreeSegs[id2Min],3])
						Found = True
					else:
						d2		= (file_data[FreeSegs,2]-Rwall[-1])**2+(file_data[FreeSegs,3]-Zwall[-1])**2						#Check second point of segments
						id2Min	= np.argmin(d2)
						if(d2[id2Min] < MaxDist2):
							Rwall = np.append(Rwall, file_data[FreeSegs[id2Min],0])
							Zwall = np.append(Zwall, file_data[FreeSegs[id2Min],1])
							Found = True
					if(Found):
						if(id2Min < len(FreeSegs)-1):	FreeSegs = np.append(FreeSegs[:id2Min], FreeSegs[id2Min+1:])
						else:							FreeSegs = FreeSegs[:id2Min]

				Rwall = Rwall[::-1]																								#Reverse curve direction
				Zwall = Zwall[::-1]

			Walls[-1].Rwall = np.copy(Rwall)
			Walls[-1].Zwall = np.copy(Zwall)
			if((Walls[-1].Rwall[0] == Walls[-1].Rwall[-1]) and (Walls[-1].Zwall[0] == Walls[-1].Zwall[-1])): Walls[-1].Closed = True
			else:																		 					 Walls[-1].Closed = False
			Walls[-1].Clockwise = 0
			Walls[-1].Changed	 = True
				

	if(DEBUG > 0):	print("load_walls_file: Completed")
	
	return Walls