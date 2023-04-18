
import numpy as np
import os
from routines.globals			import DEBUG

#=================================================
# This routine write soledge2D.elemente file
#=================================================

def save_wall_material_pump_puff(WallFile, Eirene):

	if(DEBUG > 0):	print("save_wall_material_pump_puff")

	fid = open(WallFile,'w')
	fid.write("##################################################################\n")
	fid.write("####   File Wall material, pump and puff positions        ########\n")
	fid.write("##################################################################\n")
	fid.write("\n\n")
	
	fid.write("Number of Materials:\n")
	fid.write("Nmat= {:d}\n".format(len(Eirene.Wall.Material)))
	fid.write("\n")
	
	for n in range(len(Eirene.Wall.Material)):
		R = []
		Z = []
		for iWall in range(len(Eirene.Wall.IsPump)):
			iMats = np.where(~Eirene.Wall.IsPump[iWall])[0]
			if(len(iMats) > 0):
				Rw, Zw = compute_rz_123(Eirene, iWall, iMats, n)
				R.extend(Rw)
				Z.extend(Zw)

		fid.write('{:d}\t{:s}\t{:d}\n'.format(n+1, Eirene.Wall.Material[n].Name, len(R)))
		for k in range(len(R)):
			fid.write('\tstart:\t{:f}\t{:f}\n'.format(R[k][0], Z[k][0]))
			fid.write('\tmid:  \t{:f}\t{:f}\n'.format(R[k][1], Z[k][1]))
			fid.write('\tend:  \t{:f}\t{:f}\n'.format(R[k][2], Z[k][2]))

	fid.write("\n\n")
	fid.write("Number of Pumps:\n")
	fid.write("Npump= {:d}\n".format(len(Eirene.Wall.Pump)))
	fid.write("\n")
	
	for n in range(len(Eirene.Wall.Pump)):
		R = []
		Z = []
		for iWall in range(len(Eirene.Wall.IsPump)):
			iMats = np.where(Eirene.Wall.IsPump[iWall])[0]
			if(len(iMats) > 0):
				Rw, Zw = compute_rz_123(Eirene, iWall, iMats, n)
				R.extend(Rw)
				Z.extend(Zw)

		if(len(R) > 0):
			fid.write('{:d}\t{:s}\t{:d}\n'.format(n+1, Eirene.Wall.Pump[n].Name, len(R)))
			for k in range(len(R)):
				fid.write('\tstart:\t{:f}\t{:f}\n'.format(R[k][0], Z[k][0]))
				fid.write('\tmid:  \t{:f}\t{:f}\n'.format(R[k][1], Z[k][1]))
				fid.write('\tend:  \t{:f}\t{:f}\n'.format(R[k][2], Z[k][2]))

	fid.write("\n\n")
	fid.write("Number of puff:\n")
	fid.write("Npuff= {:d}\n".format(len(Eirene.Wall.Puff)))
	fid.write("\n")

	for n in range(len(Eirene.Wall.Puff)):
		fid.write("{:d}\t{:s}\t{:f}\t{:f}\n".format(n+1, Eirene.Wall.Puff[n].Name, Eirene.Wall.Puff[n].R, Eirene.Wall.Puff[n].Z))

	fid.close()

	if(DEBUG > 0):	print("save_wall_material_pump_puff: Completed")
	
	
def compute_rz_123(Eirene, iWall, iMats, n):
	
	TypeMat = Eirene.Wall.TypeMat[iWall]
	R12		= Eirene.Wall.R12[iWall]
	Z12		= Eirene.Wall.Z12[iWall]
	IsPump	= Eirene.Wall.IsPump[iWall]
	iMat = np.where(TypeMat[iMats] == n)[0]
	if(len(iMat) > 0):
		iMat  	= iMats[iMat]
		diMat	= iMat[1:] - iMat[:-1]
		iEnd	= np.where(diMat > 1)[0]
		if(len(iEnd) > 0):
			iEnd	= np.append(iEnd, len(diMat))
			iStart  = np.append(0,iEnd[:-1]+1)
		else:
			iStart = np.array([0])
			iEnd   = np.array([len(iMat)-1])
				
		if((len(iStart) > 1) and (iMat[iStart[0]] == 0) and (iMat[iEnd[-1]] == IsPump.shape[0]-1)):
			iStart[0] = iStart[-1]
			iStart    = iStart[:-1]
			iEnd      = iEnd[:-1]
		
		iStart = iMat[iStart]
		iEnd   = iMat[iEnd]
			
		Rw = []
		Zw = []
		R = np.empty((3), dtype = 'f8')
		Z = np.empty((3), dtype = 'f8')
		FactOffset = 0.1
		for iSeg in range(len(iStart)):
			R[0] = R12[iStart[iSeg],0] + FactOffset*(R12[iStart[iSeg],1] - R12[iStart[iSeg],0]) 
			Z[0] = Z12[iStart[iSeg],0] + FactOffset*(Z12[iStart[iSeg],1] - Z12[iStart[iSeg],0]) 
			R[2] = R12[iEnd[iSeg],1]   - FactOffset*(R12[iEnd[iSeg],1]   - R12[iEnd[iSeg],0])
			Z[2] = Z12[iEnd[iSeg],1]   - FactOffset*(Z12[iEnd[iSeg],1]   - Z12[iEnd[iSeg],0]) 
				
			if(iStart[iSeg] < iEnd[iSeg]):
				kMid = int(0.5*(iStart[iSeg] + iEnd[iSeg])+0.5) 
			else:
				kMid = int(0.5*(iStart[iSeg] + iEnd[iSeg] + IsPump.shape[0]) + 0.5) 
				if(kMid >= IsPump.shape[0]): kMid -= IsPump.shape[0]

			R[1] = 0.5*(R12[kMid,0] + R12[kMid,1]) 
			Z[1] = 0.5*(Z12[kMid,0] + Z12[kMid,1])
			Rw.append(np.copy(R))
			Zw.append(np.copy(Z))

		return Rw, Zw
	else:
		return [],[]
		