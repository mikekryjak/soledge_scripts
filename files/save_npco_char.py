# Function definition is here

import os
from routines.globals		import DEBUG, FILE_TYPE_SOLEDGE2D, FILE_TYPE_SOLEDGE3X

#=================================================
# This routine soledge2D.npco_char file
#=================================================

def save_npco_char(Dir, Eirene, FileType=FILE_TYPE_SOLEDGE2D):

	if(DEBUG > 0):	print("save_npco_char")

	try:
		os.mkdir(Dir)
	except OSError:
		pass

	if(FileType == FILE_TYPE_SOLEDGE2D):	prefix = "soledge2D"
	elif(FileType == FILE_TYPE_SOLEDGE3X):	prefix = "raptorX"
	else:
		print("\tERROR: I am sorry invalid FileType=",FileType," Check code!")
		exit()

	if(DEBUG > 1):	print("\tsaving to: "+Dir + prefix + ".npco_char")

	fid = open(Dir + prefix + ".npco_char",'w')

	RKnots = Eirene.RKnots
	ZKnots = Eirene.ZKnots
	nKnots = len(RKnots)
	fid.write('{:d}\n'.format(nKnots))
	for n in range(nKnots):
		fid.write('{:d}\t'.format(n+1))
		fid.write('{:12.7e}\t'.format(100.*RKnots[n]))
		fid.write('{:12.7e}\n'.format(100.*ZKnots[n]))

	fid.close()

	if(DEBUG > 0):	print("save_npco_char: Completed")
