import os
from tkinter				import messagebox
import numpy				as np
from routines.globals		import DEBUG

#=========================================================
# This routine read a file with header
#=========================================================

def load_text_data(LoadFile):

	if(DEBUG > 0): print("load_text_data: Reading from: ",LoadFile)

	FileExt = os.path.splitext(os.path.split(LoadFile)[1])[1]

	ToLoad	   = True
	WithHeader = False
	CommaDelim = False

	if(FileExt == ".csv"):
		try:									#try no header comma delimiter
			FileData = np.loadtxt(LoadFile, delimiter=",")
			ToLoad	   = False
			Delimiter  = ","
		except:
			pass
		
		if(ToLoad):
			try:									#try with header comma delimiter
				Delimiter  = ","
				FileData = np.loadtxt(LoadFile, delimiter=Delimiter, skiprows = 1)
				ToLoad	   = False
				WithHeader = True
			except:
				pass

	if(ToLoad):
		try:									#try no header any space delimiter
			FileData = np.loadtxt(LoadFile)
			ToLoad = False
			Delimiter  = None
		except:
			pass
		
		if(ToLoad):
			try:								#try with header any space delimiter
				FileData = np.loadtxt(LoadFile, skiprows = 1)
				ToLoad	   = False
				WithHeader = True
				Delimiter  = None
			except:
				pass
					
	if(ToLoad):
		messagebox.showerror("Load from file", "Failed load data from file:{:s}\n".format(LoadFile))
		return [], 0
			
	if(WithHeader):								#read header
		fid  = open(LoadFile,'r')
		head = fid.readline()
		fid.close

		Headers = head.split(Delimiter)														#try splitting with space or tab
		
		if(len(Headers) != FileData.shape[1]):
			messagebox.showerror("Load from file", "Headers does not match columns in file:{:s}\n".format(LoadFile))
			return [], 0

		for i in range(len(Headers)): Headers[i] = Headers[i].strip()

	else:
		Headers = []
		for k in range(FileData.shape[1]):	Headers.append("{:d}".format(k+1))

	return Headers, FileData

