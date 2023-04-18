import os
import h5py
from h5_routines 		import *

# 	##########################
#	writes diagnostic data to a hdf5 pulse-file
# 	##########################

def save_exp_data_to_h5(exp_name, shot, diag_name, values, use_signals = False):

	full_name = os.getenv(exp_name+"_path",'.')+"/"+exp_name
	if_data = h5py.File(full_name+"_{:d}.h5".format(shot), "a")

	if(not h5_exist(if_data, diag_name)): if_data.create_group(diag_name)

	if(use_signals):
		for i in range(len(values.node_names)):
			if(values.__dict__[values.data_names[i]].type == 0):
				h5_write_signal_0d(if_data, diag_name + "/" +values.node_names[i], values.__dict__[values.data_names[i]].values,
											values.__dict__[values.data_names[i]].descr, values.__dict__[values.data_names[i]].units)
			elif(values.__dict__[values.data_names[i]].type == 1):
				h5_write_signal_1d(if_data, diag_name + "/" +values.node_names[i], values.__dict__[values.data_names[i]].values, values.__dict__[values.data_names[i]].times, 
											values.__dict__[values.data_names[i]].descr, values.__dict__[values.data_names[i]].units)
			elif(values.__dict__[values.data_names[i]].type == 2):
				h5_write_signal_2d(if_data, diag_name + "/" +values.node_names[i], values.__dict__[values.data_names[i]].values, values.__dict__[values.data_names[i]].times, values.__dict__[values.data_names[i]].xvalues,
											values.__dict__[values.data_names[i]].descr, values.__dict__[values.data_names[i]].units)
	else:
		for i in range(len(values.node_names)):
			if(values.data_type[i] == 0):
				h5_write_signal_0d(if_data, diag_name + "/" +values.node_names[i], values.__dict__[values.data_names[i]], values.data_descr[i], values.data_units[i])
			elif(values.data_type[i] == 1):
				h5_write_signal_1d(if_data, diag_name + "/" +values.node_names[i], values.__dict__[values.data_names[i]], values.__dict__[values.time_name], values.data_descr[i], values.data_units[i])
			elif(values.data_type[i] == 2):
				h5_write_signal_2d(if_data, diag_name + "/" +values.node_names[i], values.__dict__[values.data_names[i]], values.__dict__[values.time_name], values.__dict__[values.pos_name], values.data_descr[i], values.data_units[i])


	if_data.close()

	return

# 	##########################
#	reads diagnostic data from a hdf5 pulse-file
# 	##########################

def read_exp_data_from_h5(exp_name, shot, diag_name, values, use_signals = False):

	full_name = os.getenv(exp_name+"_path",'.')+"/"+exp_name
	if_data = h5py.File(full_name+"_{:d}.h5".format(shot), "r")

	if(not h5_exist(if_data, diag_name)):
		print("ERROR: ",diag_name," does not exist in:\n\tpulse file :",exp_name," for shot: {:d}".format(shot)) 
		return

	if(not hasattr(values,'data_descr')): values.data_descr = []
	if(not hasattr(values,'data_units')): values.data_units = []
	if(len(values.data_descr) < len(values.node_names)): for i in range(len(values.data_descr), len(values.node_names)) : values.data_descr.append("")
	if(len(values.data_units) < len(values.node_names)): for i in range(len(values.data_units), len(values.node_names)) : values.data_units.append("")

	for i in range(len(values.node_names)):
		if(values.data_type[i] == 0):
			values.__dict__[values.data_names[i]], values.data_descr[i], values.data_units[i] = h5_read_signal_0d(if_data, diag_name + "/" +values.node_names[i])
		elif(values.data_type[i] == 1):
			values.__dict__[values.data_names[i]], values.__dict__[values.time_name], values.data_descr[i], values.data_units[i] = h5_read_signal_1d(if_data, diag_name + "/" +values.node_names[i])
		elif(values.data_type[i] == 2):
			values.__dict__[values.data_names[i]], values.__dict__[values.time_name], values.__dict__[values.pos_name], values.data_descr[i], values.data_units[i] = h5_read_signal_2d(if_data, diag_name + "/" +values.node_names[i])

	if_data.close()

	return
