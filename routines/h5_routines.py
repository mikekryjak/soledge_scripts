#import h5py
import types
import numpy 			as np
from routines.globals	import DEBUG

#====================================================================================================
# write and read data in a h5df file
#====================================================================================================

def h5_write(file_id, Name, Data, order='F', delete = False, return_dset = False):

	if(delete): h5_delete(file_id, Name)

	if(isinstance(Data, np.ndarray)):
		if (order != 'C') and (len(Data.shape) > 1):
			dset  = file_id.create_dataset(Name, Data.T.shape, Data.dtype)
			dset[()]  = Data.T
		else:
			dset  = file_id.create_dataset(Name, Data.shape, Data.dtype)
			dset[()]  = Data
	else:	
		if(isinstance(Data, list)):
			Data2 = np.array(list)
			dset  = file_id.create_dataset(Name, Data2.shape, Data2.dtype)
			dset[()]  = Data2
		elif(isinstance(Data, bool) or isinstance(Data, np.bool_)):
			dset  = file_id.create_dataset(Name, (1,), dtype='bool')								#write scalars as array fro compatibility with Matlab version
			dset[()]  = Data
		elif(isinstance(Data, int) or isinstance(Data, np.int64)):
			dset  = file_id.create_dataset(Name, (1,), dtype='i8')
			dset[()]  = Data
		elif(isinstance(Data, float) or isinstance(Data, np.float64)):
			dset  = file_id.create_dataset(Name, (1,), dtype='f8')
			dset[()]  = Data
		elif(isinstance(Data, np.int32)):
			dset  = file_id.create_dataset(Name, (1,), dtype='i4')
			dset[()]  = Data
		elif(isinstance(Data, np.float32)):
			dset  = file_id.create_dataset(Name, (1,), dtype='f4')
			dset[()]  = Data
		elif(isinstance(Data, str)):
			dset  = file_id.create_dataset(Name, (1,), dtype='S{:d}'.format(len(Data)))
			dset[()]  = np.string_(Data)
		else:
			print("\th5_write: unknow data type from field:", Name)
	
	if(return_dset): return dset

#====================================================================================================
# read a data froma a h5df dataset
#====================================================================================================

def h5_read(file_id, Name, order='F', keep_array=True, messages=1, raise_err = 1):
	try:
		Data = file_id[Name][...]
	except:
		if((DEBUG > 0) and (messages > 0)): print("\th5_read: error reading: ",Name)
		if(raise_err != 0): raise
		else:				return

	if (order != 'C') and (len(Data.shape) > 1):
		return Data.T
	else:
		if(str(Data.dtype)[1]=="S"):
			Data = Data.astype(np.str)

		if(len(Data.shape) == 0): 
			return Data.tolist()														#Manage scalar data
		elif(keep_array or (len(Data.shape) > 1) or (Data.shape[0] > 1)):		
			return Data
		else:
			return Data[0]

#====================================================================================================
# delete a node if exists
#====================================================================================================

def h5_delete(file_id, Node):
	e = Node in file_id
	if(e): del file_id[Node]

#====================================================================================================
# check if a node exists
#====================================================================================================

def h5_exist(file_id, Node):
	e = Node in file_id
	return e

#====================================================================================================
# write signal description and units
#====================================================================================================

def	h5_write_signal_descr_units(file_id, Node, Descr="", Units=""):
	if(len(Descr) > 0): h5_write(file_id, Node+"/descr", Descr, delete=True)
	if(len(Units) > 0): h5_write(file_id, Node+"/units", Units, delete=True)

#====================================================================================================
# write signal_0d
#====================================================================================================

def h5_write_signal_0d(file_id, Node, Values, Descr="", Units=""):
	if(not  h5_exist(file_id, Node)): file_id.create_group(Node)

	h5_write(file_id, Node+"/dtype",    0, delete=True)
	h5_write(file_id, Node+"/values", Values, delete=True)

	h5_write_signal_descr_units(file_id, Node, Descr, Units)

#====================================================================================================
# write signal_1d
#====================================================================================================

def h5_write_signal_1d(file_id, Node, Values, Times, Descr="", Units=""):
	if(not  h5_exist(file_id, Node)): file_id.create_group(Node)

	h5_write(file_id, Node+"/dtype",    1, delete=True)
	h5_write(file_id, Node+"/values", Values, delete=True)
	h5_write(file_id, Node+"/times",  Times, delete=True)

	h5_write_signal_descr_units(file_id, Node, Descr, Units)

#====================================================================================================
# write signal_2d
#====================================================================================================

def h5_write_signal_2d(file_id, Node, Values, Times, xValues, Descr="", Units=""):

	if(not  h5_exist(file_id, Node)): file_id.create_group(Node)

	dset = h5_write(file_id, Node+"/values",   Values, delete=True, return_dset = True)
	vsz_range = np.array([Times[0], xValues[0], Times[-1], xValues[-1]])
	dset.attrs.create("vsz_range", vsz_range)

	h5_write(file_id, Node+"/dtype",    2, delete=True)
	h5_write(file_id, Node+"/times",    Times, delete=True)
	h5_write(file_id, Node+"/xvalues",  xValues, delete=True)

	h5_write_signal_descr_units(file_id, Node, Descr, Units)


#====================================================================================================
# write signal
#====================================================================================================

def h5_write_signal(file_id, Node, Sig):

	if(Sig.type == 0):
		h5_write_signal_0d(file_id, Node, Sig.values, Descr=Sig.descr, Units=Sig.units)
	elif(Sig.type == 1):
		h5_write_signal_1d(file_id, Node, Sig.values, Sig.times, Descr=Sig.descr, Units=Sig.units)
	elif(Sig.type == 2):
		h5_write_signal_2d(file_id, Node, Sig.values, Sig.times, Sig.xvalues, Descr=Sig.descr, Units=Sig.units)
	else:
		if(DEBUG > 0): print("\th5_write_signal: not allowed type: ",Sig.type)


#====================================================================================================
# read signal description and units
#====================================================================================================

def h5_read_signal_descr_units(file_id, Node):

	if(h5_exist(file_id, Node+"/descr")): Descr = h5_read(file_id, Node+"/descr")
	else:									Descr = ""

	if(h5_exist(file_id, Node+"/units")): Units = h5_read(file_id, Node+"/units")
	else:									Units = ""

	return Descr, Units

#====================================================================================================
# read signal
#====================================================================================================

def h5_read_signal(file_id, Node):

	sig = types.SimpleNamespace()

	if(h5_exist(file_id, Node+"/dtype")): dtype = h5_read(file_id, Node+"/dtype", keep_array=False)
	else: 
		if(h5_exist(file_id, Node)):	
			dtype = -1
		else:
			sig.type = -1
			if(DEBUG > 0): print("\th5_read_signal: error reading: ",Node)
			return sig


	if(dtype > -1):
		sig.type   = dtype
		sig.values = h5_read(file_id, Node+"/values")
		sig.descr, sig.units = h5_read_signal_descr_units(file_id, Node)
	elif(dtype > -2):
		sig.type   = 0
		sig.values = h5_read(file_id, Node, keep_array=True)
		sig.descr  = ""
		sig.units  = ""
	else:
		sig.type   = -1
		
	if(dtype > 0): sig.times   = h5_read(file_id, Node+"/times")
	if(dtype > 1): sig.xvalues = h5_read(file_id, Node+"/xvalues")
	
	return sig
