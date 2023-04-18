
import numpy as np
import os
import types

from eirene.set_wall_material_pump			import set_wall_material_pump
from eirene.set_wall_puff					import set_wall_puff
from routines.globals						import DEBUG

#=================================================
# This routine write soledge2D.elemente file
#=================================================

def load_wall_material_pump_puff(WallFile, Eirene):

	if(DEBUG > 0):	print("load_wall_material_pump_puff")

	fid = open(WallFile,'r')
	fid.readline()
	fid.readline()
	fid.readline()
	fid.readline()
	fid.readline()
	
	R = np.empty(3, dtype = 'f8')
	Z = np.empty(3, dtype = 'f8')

#	Material section
	
	fid.readline()
	Line = fid.readline()
	nMat = int(eval(Line.split()[1]))
	fid.readline()

	Eirene.Wall.Material = []
	for n in range(nMat):
		Eirene.Wall.Material.append(types.SimpleNamespace())
		Values = fid.readline().split("\t")
		Eirene.Wall.Material[-1].Name = Values[1]
		nSegs = int(eval(Values[2]))
		
		for k in range(nSegs):
			for i in range(3):
				Values = fid.readline().split("\t")
				R[i] = eval(Values[2])
				Z[i] = eval(Values[3])		
			set_wall_material_pump(Eirene, False, n, R, Z)

#	Pump section
			
	fid.readline()
	fid.readline()
	fid.readline()
	Line = fid.readline()
	nPump = int(eval(Line.split()[1]))
	fid.readline()
	
	Eirene.Wall.Pump = []
	for n in range(nPump):
		Values = fid.readline().split("\t")
		
		Eirene.Wall.Pump.append(types.SimpleNamespace())
		Eirene.Wall.Pump[-1].Name = Values[1]
		nSegs = int(eval(Values[2]))
		
		for k in range(nSegs):
			for i in range(3):
				Values = fid.readline().split("\t")
				R[i] = eval(Values[2])
				Z[i] = eval(Values[3])
			set_wall_material_pump(Eirene, True, n, R, Z)

	fid.readline()
	fid.readline()
	fid.readline()
	Line = fid.readline()
	nPuff = int(eval(Line.split()[1]))
	fid.readline()

# Puff section

	Eirene.Wall.Puff = []
	for n in range(nPuff):
		Values = fid.readline().split("\t")
		set_wall_puff(Eirene, -1, eval(Values[2]), eval(Values[3]))
		Eirene.Wall.Puff[-1].Name = Values[1]
	fid.close()

	if(DEBUG > 0):	print("load_wall_material_pump_puff: Completed")
	
		