import numpy as np
import os
from routines.globals					import BC_CORE, DEBUG, FILE_TYPE_SOLEDGE2D, FILE_TYPE_SOLEDGE3X
#=================================================
# This routine write .....
#=================================================

def save_puffs(Dir, Eirene, FileType=FILE_TYPE_SOLEDGE2D):

	if(DEBUG > 0):	print("save_puffs")

	try:
		os.mkdir(Dir)
	except OSError:
		pass

	if(DEBUG > 1):	print("\tsaving to: "+ Dir + "puffs")

	fid = open(Dir + 'puffs','w')
		
	for n in range(len(Eirene.Wall.Puff)):
		fid.write('{:d}\t'.format(n+1))									#python index
		fid.write('{:s}\t'.format(Eirene.Wall.Puff[n].Name))
		fid.write('{:d}\t'.format(Eirene.Wall.Puff[n].ntri+1))				#python index
		fid.write('{:d}\t'.format(Eirene.Wall.Puff[n].side+1))				#python index
		fid.write('{:f}\t'.format(Eirene.Wall.Puff[n].R))
		fid.write('{:f}\n'.format(Eirene.Wall.Puff[n].Z))

	fid.close()

	if(DEBUG > 0):	print("save_puffs: Completed")
