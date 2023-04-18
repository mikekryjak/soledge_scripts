import os
import types
import numpy				as np
from routines.globals		import DEBUG

#=========================================================
# This routine read a file with header
#=========================================================

def load_chords_file(ChordsdFile):

	if(DEBUG > 0): print("load_chords_file: Reading from: ",ChordsdFile)

	FileExt = os.path.splitext(os.path.split(ChordsdFile)[1])[1]

	if(FileExt == ".chords"):
		FileName = ChordsdFile

		if(not os.path.isfile(FileName)):
			if(len(os.path.dirname(FileName)) == 0):				#If file does not exist and no path try to find on path
				iEnd = FileName.find("_")
				if(iEnd > -1): base_path = FileName[:iEnd]
				else:		   base_path = os.path.splitext(os.path.split(ExpFile)[1])[0]
				try:
					base_path = os.environ[base_path+"_path"]
					FileName  = base_path+"/"+FileName
					if(not os.path.isfile(FileName)): FileName  = ""
				except:
					FileName  = ""
			else:
				FileName  = ""

		fid  = open(FileName,'r')
		head = fid.readline()
		fid.close

		Headers  = head.split("'")
		
		Chords	 			 = types.SimpleNamespace()
		Chords.SubDiags		 = [types.SimpleNamespace()]
		Chords.label				 		= Headers[1]
		Chords.SubDiags[0].los_name 		= Headers[1]
		Chords.SubDiags[0].n_samps_lines	= eval(Headers[2])

		FileData = np.loadtxt(FileName, skiprows = 1)

		Chords.SubDiags[0].data		   = False
		Chords.SubDiags[0].lines_start = FileData[:,0:3].T
		Chords.SubDiags[0].lines_end   = FileData[:,3:6].T
		Chords.SubDiags[0].x_values    = FileData[:,6]
		Chords.SubDiags[0].n_chans     = np.arange(1,len(Chords.SubDiags[0].x_values)+1)

		print("Chords.SubDiags[0].lines_start.shape=",Chords.SubDiags[0].lines_start.shape)

	if(DEBUG > 0): print("load_chords_file: completed")

	return Chords

