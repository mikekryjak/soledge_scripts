import os
import types
import numpy					as np
from files.load_chords_file		import load_chords_file
from routines.globals			import DEBUG

#=========================================================
# This routine read a file with header
#=========================================================

def load_diag_file(DiagFile):

	if(DEBUG > 0): print("load_diag_file: Reading from: ",DiagFile)

	FileExt = os.path.splitext(os.path.split(DiagFile)[1])[1]
	
	if(len(FileExt) == 0): FileName = DiagFile + ".diag"
	else:				   FileName = DiagFile

	if(not os.path.isfile(FileName)):
		if(len(os.path.dirname(FileName)) == 0):				#If file does not exist and no path try to find on path
			iEnd = FileName.find("_")
			if(iEnd > -1): base_path = FileName[:iEnd]
			else:		   base_path = os.path.splitext(os.path.split(FileName)[1])[0]
			try:
				base_path = os.environ[base_path+"_path"]
				FileName  = base_path+"/"+FileName
				if(not os.path.isfile(FileName)): FileName  = ""
			except:
				FileName  = ""
		else:
			FileName  = ""


	print("\tReading diag file from = ", FileName)
	fid   = open(FileName,'r')
	head  = fid.readline()
	datal = fid.readline()													#read next line to check if there are data lines 
	fid.close

	if(FileExt == ".csv"):	delimiter = ","
	else:					delimiter = None

	Headers = head.split(delimiter)

	DiagDir  = os.path.dirname(DiagFile)
	if((Headers[1][0] != "/") and (Headers[1][0] != "$") and (len(DiagDir) > 0)):
		ChordsFile = DiagDir + "/" + Headers[1]
	elif(Headers[1][0] == "$"):
		sFile  = Headers[1].find("/")
		ChordsFile = os.environ[Headers[1][1:sFile]] + Headers[1][sFile]
	else:
		ChordsFile = Headers[1]
	
	Diag = load_chords_file(ChordsFile)
	Diag.par_name  = Headers[0]
	Diag.marker	   = "o"
	Diag.color	   = "g"

	if(len(datal) > 0): 
		try:
			FileData = np.loadtxt(DiagFile, delimiter=delimiter, skiprows = 1)
		except:
			FileData = None
	else:
		FileData = None

	if((FileData is not None) and (FileData.shape[0] > 1)):
		Diag.SubDiags[0].data   = True
		Diag.SubDiags[0].x_data = FileData[:,0]
		Diag.SubDiags[0].y_data	= FileData[:,1]
		if(FileData.shape[0] > 2): Diag.SubDiags[0].y_err = FileData[:,2]
		else:					   Diag.SubDiags[0].y_err = np.empty(0,dtype='f8')
	else:
		Diag.SubDiags[0].data   = False

	return Diag
	
