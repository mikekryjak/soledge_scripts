# Function definition is here

import os
from routines.globals			import DEBUG, FILE_TYPE_SOLEDGE2D, FILE_TYPE_SOLEDGE3X

#=================================================
# This routine save XXXX.surfaces file
#=================================================

def save_surfaces(Dir, Eirene, FileType=FILE_TYPE_SOLEDGE2D):

	if(len(Eirene.Surfaces) == 0): return

	if(DEBUG > 0):	print("save_surfaces")

	try:
		os.mkdir(Dir)
	except OSError:
		pass

	if(FileType == FILE_TYPE_SOLEDGE2D):	prefix = "soledge2D"
	elif(FileType == FILE_TYPE_SOLEDGE3X):	prefix = "raptorX"
	else:
		print("\tERROR: I am sorry invalid FileType=",FileType," Check code!")
		exit()

	if(DEBUG > 1):	print("\tsaving to: "+Dir + prefix + ".surfaces")

	fid = open(Dir + prefix + ".surfaces",'w')

#	Count number of elementary segments in Surfaces

	Surfaces = Eirene.Surfaces
	nSurf = 0
	for n in range(len(Surfaces)):
		nSurf += len(Surfaces[n].R) - 1

	fid.write('{:d}\n'.format(nSurf))
	for n in range(len(Surfaces)):
		Surface = Surfaces[n]
		for k in range(len(Surface.R)-1):
			fid.write('{:12.7e}\t'.format(Surface.R[k]))
			fid.write('{:12.7e}\t'.format(Surface.Z[k]))
			fid.write('{:12.7e}\t'.format(Surface.R[k+1]))
			fid.write('{:12.7e}\t'.format(Surface.Z[k+1]))
			fid.write('{:12.7e}\n'.format(Surface.Rcoef))

	fid.close()

	if(DEBUG > 0):	print("save_surfaces: Completed")
