#===================================================
# Routine to read an integer value checchig it is ok
#===================================================
import numpy as np
import sys

def input_int(Prompt=''):
	
	is_ok = 0
	while is_ok != 1:
		aa = input(Prompt)
		if len(aa) > 0:
			try:	
				value = int(eval(aa))
				is_ok = 1
			except:
				try:	
					value = np.array(eval(aa), dtype='i')
					is_ok = 1
				except:
					pass
					
	
	return value

#===================================================
# Routine to read a float value checchig it is ok
#===================================================

def input_float(Prompt=''):
	
	is_ok = 0
	while is_ok != 1:
		aa = input(Prompt)

		if len(aa) > 0:
			try:	
				value = float(eval(aa))
				is_ok = 1
			except:
				try:	
					value = np.array(eval(aa), dtype='f8')
					is_ok = 1
				except:
					pass
					
				
	
	return value

#===================================================
# Routine scan a string to get data
#===================================================

def scan_array(String, Dtype, Nvalues=1, Separator=","):
	if(Separator != " "):
		NewStrings = String.replace(Separator," ")
		NewStrings = NewStrings.split()
	else:
		NewStrings = String.split()

	if Nvalues == 1:
		Data	= np.array([NewStrings[0]], dtype=Dtype)
	else:
		Data	= np.array([NewStrings[0:Nvalues]], dtype=Dtype)

#	Data = Data.astype(Dtype)

	return Data


#===================================================
# Routine to count file lines
#===================================================

def file_get_nlines(filename):

	fid = open(filename,'r')

	n=0
	while fid.readline() : n +=1

	fid.close()

	return n
