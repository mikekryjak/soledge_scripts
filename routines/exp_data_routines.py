import os
import h5py
from numpy								import *
import numpy.matlib						as mat

from math										import sqrt
from files.load_chords_file				import load_chords_file
from routines.globals 					import DEBUG, LENGTH_TYPES
from routines.h5_routines 				import *
from routines.interp2d					import interp2d
from routines.set_profile_from_filedata	import interpolate_on_length

# 	##########################
#	writes diagnostic data to a hdf5 pulse-file
# 	##########################

def save_exp_data_to_h5(exp_name, shot, diag_name, values, use_signals = False):

	full_name = os.getenv(exp_name+"_path",'.')+"/"+exp_name
	if_data = h5py.File(full_name+"_{:d}.h5".format(shot), "a")

	if(not h5_exist(if_data, diag_name)): if_data.create_group(diag_name)

	if(use_signals):
		for i in range(len(values.node_names)):
			h5_write_signal(if_data, diag_name + "/" +values.node_names[i], values.__dict__[values.data_names[i]])
	else:
		for i in range(len(values.node_names)):
			if(values.data_type[i] == 0):
				h5_write_signal_0d(if_data, diag_name + "/" +values.node_names[i], values.__dict__[values.data_names[i]], Descr = values.data_descr[i], Units = values.data_units[i])
			elif(values.data_type[i] == 1):
				h5_write_signal_1d(if_data, diag_name + "/" +values.node_names[i], values.__dict__[values.data_names[i]], values.__dict__[values.time_name], Descr = values.data_descr[i], Units = values.data_units[i])
			elif(values.data_type[i] == 2):
				h5_write_signal_2d(if_data, diag_name + "/" +values.node_names[i], values.__dict__[values.data_names[i]], values.__dict__[values.time_name], values.__dict__[values.pos_name], Descr = values.data_descr[i], Units = values.data_units[i])


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
	if(len(values.data_descr) < len(values.node_names)):
		for i in range(len(values.data_descr), len(values.node_names)) : values.data_descr.append("")
	if(len(values.data_units) < len(values.node_names)):
		for i in range(len(values.data_units), len(values.node_names)) : values.data_units.append("")

	for i in range(len(values.node_names)):
		if(values.data_type[i] == 0):
			values.__dict__[values.data_names[i]], values.data_descr[i], values.data_units[i] = h5_read_signal_0d(if_data, diag_name + "/" +values.node_names[i])
		elif(values.data_type[i] == 1):
			values.__dict__[values.data_names[i]], values.__dict__[values.time_name], values.data_descr[i], values.data_units[i] = h5_read_signal_1d(if_data, diag_name + "/" +values.node_names[i])
		elif(values.data_type[i] == 2):
			values.__dict__[values.data_names[i]], values.__dict__[values.time_name], values.__dict__[values.pos_name], values.data_descr[i], values.data_units[i] = h5_read_signal_2d(if_data, diag_name + "/" +values.node_names[i])

	if_data.close()

	return


# 	#################
#	reads points diagnostic data
# 	#################

def get_exp_data_point_help():
	print("\nGlobals signal parameters for points on cut type diagnostic")
	print("\tsig_times             <signal name> with times")
	print("\tsig_rho_pol           <signal name> x in rho_pol coordinate")
	print("\tsig_psi_pol           <signal name> x in psi_pol coordinate")
	print("\tsig_rho_v             <signal name> x in rho from volume coordinate")
	print("\tsig_length            <signal name> x in cut lenght coordinate")
	print("\tsig_dist              <signal name> x in distance from separatrix")
	print("\tsig_xm                <signal name> x in xm mesh coordinate")
	print("\nSignal name definition and signal specific parameters")
	print("\tsig_name               <signal name> where name data name (Ex. sig_dense, sig_tempi)")		
	print("\tsig_name_times         <signal name> with times")		
	print("\tsig_name_err           <signal name> error (name is data name)")		
	print("\tsig_name_skip_chans    <array> of channels to skip (Ex. [1,4,7])")
	print("\tsig_name_nchans        <signal name> with channel numbers")
	print("\tsig_name_x_fun         <function> of to modify x (Ex. 'x-5')")
	print("\tsig_name_data_fun      <function> to modify data (Ex. '5*y')")
	print("\tsig_name_data_min      <value> minimum acceptedd data values")
	print("\tsig_name_data_max      <value> maximum accepted data values")
	print("\tsig_name_err_max       <value> maximum accepted error on data")		

def get_exp_data_point(xName, yName, exp_descr, Lengths, OutxOffset=0., InxOffset=0., NameExt="", DiagType = "POINT"):

	if(DEBUG > 0): print("get_exp_data_point")
	if(DEBUG > 1):
		print("\txName=", xName)
		print("\tyName=", yName)

	diag_data = []

	full_name = os.getenv(exp_descr.pulse_name+"_path",'.')+"/"+exp_descr.pulse_name+"_{:d}.h5".format(int(exp_descr.shot))
	try:
		if_data = h5py.File(full_name, "r")
	except:
		print("ERROR: does not exist pulse file :"+full_name) 
		return diag_data

	OutxOff	 =  0.
	for k in range(len(LENGTH_TYPES)):
		if(xName.upper() == LENGTH_TYPES[k]):
			OutXType = k-1
			break
	if(OutXType == -1):
		OutXType = 0
		OutxOff	 = OutxOffset

	for id in range(len(exp_descr.diags)):
		if(exp_descr.diags[id].type == DiagType):
			par_names	= [exp_descr.diags[id].signals[k].par_name for k in range(len(exp_descr.diags[id].signals))]
			found_names = [False for k in range(len(par_names))]

			InXType = -2
			InxOff  = 0. 
			for k in range(len(LENGTH_TYPES)):
				ix = indexes_upper(par_names, yName+"_" + LENGTH_TYPES[k], found=found_names)
				if((len(NameExt) > 0) and (len(ix) == 0)): ix = indexes_upper(par_names, LENGTH_TYPES[k]+"_" + NameExt, found=found_names)
				if(len(ix) == 0): ix = indexes_upper(par_names, LENGTH_TYPES[k], found=found_names)
				if(len(ix) != 0):
					InXType = k -1
					break

				if(InXType == -1): 
					InXType = 0
					if(InxOffset != 0.): InxOff	= InxOffset
					else:  				 InxOff	= OutxOffset

			if(InXType == -1): 
				InXType = 0
				if(InxOffset != 0.): InxOff	= InxOffset
				else:  				 InxOff	= OutxOffset

			
			iv_err = []
			iv	   = []
			if(len(NameExt) > 0): 
				iv	   = indexes_upper(par_names, yName+"_" + NameExt, found=found_names)
				if(not exp_descr.diags[id].average): iv_err = indexes_upper(par_names, yName+"_" + NameExt+"_err", found=found_names)
			if(len(iv) == 0):
				iv	   = indexes_upper(par_names, yName, found=found_names)
				iv_err = indexes_upper(par_names, yName+"_err", found=found_names)

			if(len(NameExt) > 0): it = indexes_upper(par_names, yName+"_"+ NameExt+"_times" , found=found_names)
			else:				  it = []
			if(len(it) == 0): 	  it = indexes_upper(par_names, yName+"_times" , found=found_names)	
			if(len(it) == 0): 	  it = indexes_upper(par_names, "times", found=found_names)	

			if((len(ix) == 1) and (len(iv) == 1)):
				if(DEBUG > 1):
					print("\tFound:             ", yName)
					print("\tFor x reading:     ", exp_descr.diags[id].signals[ix[0]].sig_name)

				sData = h5_read_signal(if_data, exp_descr.diags[id].signals[ix[0]].sig_name)
				if(sData.type == -1): continue
				xData = sData.values

				if(DEBUG > 1): print("\tFor y reading:     ", exp_descr.diags[id].signals[iv[0]].sig_name)
				sData = h5_read_signal(if_data, exp_descr.diags[id].signals[iv[0]].sig_name)
				if(sData.type == -1): continue
				yData = sData.values

				if(len(it) == 1):
					if(DEBUG > 1): print("\tFor t reading:     ", exp_descr.diags[id].signals[it[0]].sig_name)
					sData = h5_read_signal(if_data, exp_descr.diags[id].signals[it[0]].sig_name)
					if(sData.type == -1): continue
					tData = sData.values
				elif(sData.type > 0):
					tData = sData.times
				else:
					continue

#				check for y_err

				if(len(iv_err) == 1): 
					if(DEBUG > 1):print("\tFor y_err reading: ", exp_descr.diags[id].signals[iv_err[0]].sig_name)
					sData = h5_read_signal(if_data, exp_descr.diags[id].signals[iv_err[0]].sig_name)
					if(sData.type >= 0): yData_err = sData.values
					else:				 iv_err = []

#				Fix for single time	

				if((len(tData) == 1) and (len(yData.shape) == 1) and (len(xData.shape) == 1)):
					yData = yData.reshape(1,yData.shape[0])
					xData = xData.reshape(1,xData.shape[0])
					if((len(iv_err) == 1) and (len(yData_err.shape) == 1)): yData_err = yData_err.reshape(1,yData_err.shape[0])

#				check for nchans

				iv_nchans = indexes_upper(par_names, yName+"_" + NameExt+"_nchans", found=found_names)
				if(len(iv_nchans) == 1):
					if(DEBUG > 1):print("\tReading: ", exp_descr.diags[id].signals[iv_nchans[0]].sig_name," as nchans for: ",name)
					sData = h5_read_signal(if_data, exp_descr.diags[id].signals[iv_nchans[0]].sig_name)
					if(sData.type >= 0): NumChans  = sData.values
					else:				 iv_nchans = []

					if(NumChans.shape[0] != yData.shape[1]):
						if(DEBUG > 0): print("\tMismatch between nChans signal and data signal for: ",name)
						continue

				if(len(iv_nchans) != 1): NumChans = np.linspace(1,yData.shape[1],num=yData.shape[1])
				
				if(len(NameExt) == 0):	yFullName = yName
				else:					yFullName = yName + "_" + NameExt

#				check for skip_chans

				iv_skip_chans = indexes_upper(par_names, yFullName+"_skip_chans", found=found_names)
				if(len(iv_skip_chans) == 1):
					skip_chans = eval(exp_descr.diags[id].signals[iv_skip_chans[0]].sig_name)
					iKeepChans = np.linspace(0,NumChans.shape[0]-1,num=NumChans.shape[0]-1).np.int
					KeepChans  = NumChans[iKeepChans]
					for iChan in skip_chans:
						iKeepChans = iKeepChans[np.where(iChan != KeepChans)[0]]
						KeepChans  = NumChans[iKeepChans]
						
					if(len(KeepChans) < len(NumChans)):
						yData = yData[:,iKeepChans]
						if(len(iv_err) == 1): yData_err = yData_err[:,iKeepChans]

#				check for data_fun and apply it

				iv_data_fun = indexes_upper(par_names, yFullName+"_data_fun", found=found_names)
				if(len(iv_data_fun) == 1):
					data_fun = exp_descr.diags[id].signals[iv_data_fun[0]].sig_name
					x	  = xData
					y	  = yData
					yData = eval(data_fun)
					x	  = 0
					y	  = 0

#				check for x_fun

				iv_x_fun = indexes_upper(par_names, yFullName+"_x_fun", found=found_names)
				if(len(iv_x_fun) == 1):
					x_fun = exp_descr.diags[id].signals[iv_x_fun[0]].sig_name
					x	  = xData
					xData = eval(x_fun)
					x	  = 0
				
#				check for err_max

				if(len(iv_err) == 1): 
					iv_err_max = indexes_upper(par_names, yFullName+"_err_max", found=found_names)
					if(len(iv_err_max) == 1): err_max = exp_descr.diags[id].signals[iv_err_max[0]].sig_name
				else: iv_err_max =[]

#				check for x_min

				iv_x_min = indexes_upper(par_names, yFullName+"_x_min", found=found_names)
				if(len(iv_x_min) == 1): x_min = exp_descr.diags[id].signals[iv_x_min[0]].sig_name

#				check for x_max

				iv_x_max = indexes_upper(par_names, yFullName+"_x_max", found=found_names)
				if(len(iv_x_max) == 1): x_max = exp_descr.diags[id].signals[iv_x_max[0]].sig_name

#				check for data_min

				iv_data_min = indexes_upper(par_names, yFullName+"_data_min", found=found_names)
				if(len(iv_data_min) == 1): data_min = exp_descr.diags[id].signals[iv_data_min[0]].sig_name

#				check for data_max

				iv_data_max = indexes_upper(par_names, yFullName+"_data_max", found=found_names)
				if(len(iv_data_max) == 1): data_max = exp_descr.diags[id].signals[iv_data_max[0]].sig_name
				
				NotFound = 0
				for k in range(len(found_names)):
					if(not found_names[k]): NotFound += 1

				if(NotFound > 0):
					print("\tWARNING:")
					print("\tIn experimental descrition following parameters are not used or unknow")
					for k in range(len(found_names)):
						if(not found_names[k]): print("\t\t",par_names[k])
					
				sData = 0

				if((tData.max() < exp_descr.tstart) or (tData.min() > exp_descr.tend)):
					if(DEBUG > 1):
						print("\tNo data points in requested time intervall:",exp_descr.tstart," - ", exp_descr.tend)
						print("\t\tdata available in range:",tData.min()," - ",tData.max())

					continue				#No data points in requested intervall

				if(exp_descr.tstart < exp_descr.tend):
					iTimes	= np.where((tData >= exp_descr.tstart) & (tData <= exp_descr.tend))[0]	
					nx		 = yData.shape[1]
					xData	 = xData[iTimes, :].reshape(-1)
					yData	 = yData[iTimes, :].reshape(-1)
					if((exp_descr.diags[id].average) and (len(iTimes) > 1)):
						tData	 = tData[iTimes]
						tData    = mat.repmat(tData,nx,1).T.reshape(-1)
						nTimes	 = len(tData)

					if(len(iv_err) == 1): yData_err = yData_err[iTimes, :].reshape(-1)

					if(len(iv_err_max) == 1):
						iOk   	  = np.where(yData_err < err_max)[0]
						xData	  = xData[iOk]
						yData	  = yData[iOk]
						if(len(iv_err) == 1): yData_err = yData_err[iOk]
						if((exp_descr.diags[id].average) and (len(iTimes) > 1)): tData	  = tData[iOk]

					if(len(iv_x_min) == 1):
						iOk   	  = np.where(xData > x_min)[0]
						xData	  = xData[iOk]
						yData	  = yData[iOk]
						if((exp_descr.diags[id].average) and (len(iTimes) > 1)): tData	  = tData[iOk]
						if(len(iv_err) == 1): yData_err = yData_err[iOk]

					if(len(iv_x_max) == 1):
						iOk   	  = np.where(xData < x_max)[0]
						xData	  = xData[iOk]
						yData	  = yData[iOk]
						if((exp_descr.diags[id].average) and (len(iTimes) > 1)): tData	  = tData[iOk]
						if(len(iv_err) == 1): yData_err = yData_err[iOk]

					if(len(iv_data_min) == 1):
						iOk   	  = np.where(yData > data_min)[0]
						xData	  = xData[iOk]
						yData	  = yData[iOk]
						if((exp_descr.diags[id].average) and (len(iTimes) > 1)): tData	  = tData[iOk]
						if(len(iv_err) == 1): yData_err = yData_err[iOk]

					if(len(iv_data_max) == 1):
						iOk   	  = np.where(yData < data_max)[0]
						xData	  = xData[iOk]
						yData	  = yData[iOk]
						if((exp_descr.diags[id].average) and (len(iTimes) > 1)): tData	  = tData[iOk]
						if(len(iv_err) == 1): yData_err = yData_err[iOk]

					if(exp_descr.diags[id].remove_zero):
						iNonZero = np.where(yData > 0.)[0]
						xData	 = xData[iNonZero]
						yData	 = yData[iNonZero]
						if((exp_descr.diags[id].average) and (len(iTimes) > 1)): tData	  = tData[iNonZero]
						if(len(iv_err) == 1): yData_err = yData_err[iNonZero]

					if((exp_descr.diags[id].average) and (len(iTimes) > 1)):
						nyPts	 = 40
						len_arr  = nyPts
						xDataTmp = np.copy(xData)
						tmpTimes, xData, yData = interp2d(tData, xDataTmp, yData, nTimes, nyPts, 'none')

						yData_err = np.std(yData,axis=0)[1:-1]
						iv_err	  = [0]
						xData	  = np.mean(xData,axis=0)[1:-1]
						yData	  = np.mean(yData,axis=0)[1:-1]
						tmpTimes = 0.
					
					iTimes = 0
				else:
					time_ave = 0.5*(exp_descr.tstart + exp_descr.tend)
					iTimes = np.argmin(np.abs(tData - time_ave))
					xData	 = xData[iTimes, :]
					yData	 = yData[iTimes, :]
					if(len(iv_err) == 1): yData_err = yData_err[iTimes, :]

					if(len(iv_err_max) == 1):
						iOk   	  = np.where(yData_err < err_max)[0]
						xData	  = xData[iOk]
						yData	  = yData[iOk]
						yData_err = yData_err[iOk]

					if(len(iv_data_min) == 1):
						iOk   	  = np.where(yData > data_min)[0]
						xData	  = xData[iOk]
						yData	  = yData[iOk]
						if(len(iv_err) == 1): yData_err = yData_err[iOk]

					if(len(iv_data_max) == 1):
						iOk   	  = np.where(yData < data_max)[0]
						yData	  = yData[iOk]
						if(len(iv_err) == 1): yData_err = yData_err[iOk]

					if(exp_descr.diags[id].remove_zero):
						iNonZero = np.where(yData > 0.)[0]
						xData	  = xData[iOk]
						yData	 = yData[iNonZero]
						if(len(iv_err) == 1): yData_err = yData_err[iNonZero]
					
#				Define data


				xValues, yValues = interpolate_on_length(xData, yData, Lengths[:,OutXType], Lengths, OutXType, InXType, OutxOff, extrapolate = False, yName=yName)

#				for ii in range(len(xData)):   print("xData,yData    =",xData[ii],yData[ii])
#				for ii in range(len(xValues)): print("xValues,yValues=",xValues[ii],yValues[ii])

				if((len(xValues) > 0) and (len(yValues) > 0)):

					diag_data.append(types.SimpleNamespace())

					diag_data[-1].type	  = exp_descr.diags[id].type
					diag_data[-1].xValues = xValues
					diag_data[-1].yValues = yValues

					if(len(iv_err) == 1): 
						diag_data[-1].errors = True
						diag_data[-1].xValues, diag_data[-1].yErrors = interpolate_on_length(xData, yData_err, Lengths[:,OutXType], Lengths, OutXType, InXType, OutxOff, extrapolate = False, yName=yName+"_err")
					else:
						diag_data[-1].errors = False
#				

					diag_data[-1].xName		= par_names[ix[0]]
					diag_data[-1].yName		= par_names[iv[0]]
					diag_data[-1].label		= exp_descr.diags[id].label
					diag_data[-1].marker	= exp_descr.diags[id].marker
					diag_data[-1].color		= exp_descr.diags[id].color
					diag_data[-1].NumChans	= NumChans

	if_data.close()

	return diag_data


# 	######################
#	reads lines integrated diagnostic data
# 	######################


def get_exp_data_rz_help():
	print("\nSpecific global diagnostic for rx_point type diagnostic")
	print("\tr_min          <value> min r-axis 2D plot(D=auto)")
	print("\tr_max          <value> max r-axis 2D plot (D=auto)")
	print("\tz_min          <value> min z-axis 2D plot (D=auto)")
	print("\tz_max          <value> max z-axis 2D plot (D=auto)")
	print("\tv_min          <value> min v-axis 2D plot (D=auto)")
	print("\tv_max          <value> max v-axis 2D plot (D=auto)")
	print("\tv_label        <value> v-axis label 2D plot(D=None)")
	print("\tv_log          <true/false> to set log v-axis 2D plot(D=False)")

	print("\nGlobal signal parameters for rx_point type diagnostic")
	print("\tsig_times             <signal name> with times")
	print("\tsig_r                 <signal name> R coordinate")
	print("\tsig_z                 <signal name> Z coordinate")
	print("\tsig_r_as_x            <true/false> to use r numer as x in plot")
	print("\tsig_z_as_x            <true/false> to use z numer as x in plot")
	print("\tsig_chan_as_x         <true/false> to use channel numer as x in plot")
	print("\tsig_skip_chans        <array> of channels to skip (Ex. [1,4,7])")
	print("\tsig_nchans            <signal name> with channel numbers")
	print("\tsig_r_fun             <function> of to modify r (Ex. 'r-5')")
	print("\tsig_z_fun             <function> of to modify r (Ex. 'z-5')")
	print("\tsig_data_fun          <function> to modify data (Ex. '5*y')")
	print("\tsig_data_min          <value> minimum accepted data values")
	print("\tsig_data_max          <value> maximum accepted data values")
	print("\tsig_err_max           <value> maximum accepted error on data")

	print("\nSignal name definition and signal specific parameter")
	print("\tsig_name_data          <signal name> with data, name will be the name of parameter (Ex. sig_hrts_data)")		
	print("\tsig_name_times         <signal name> with times")		
	print("\tsig_name_err           <signal name> error, name is name in sig_data_name (Ex. sig_hrts_err)")		
	print("\tsig_name_global        where name is name sig_name_data name and global is one of previous globals options")

def get_exp_data_rz(exp_descr, DiagType = "RZ_POINTS"):

	if(DEBUG > 0): print("get_exp_data_rz")

	diag_data = []
	if_data	  = None
	for exp_descr_diag in exp_descr.diags:
		if(exp_descr_diag.type == "RZ_POINTS"):
			if(exp_descr_diag.points_file != None):
				print("\nRZ_POINTS FROM TEXT FILE NOT IMPLEMENTED\n")
				if(if_data == None): return
#				Points = load_points_file(exp_descr_diag.points_file)					
#				diag_data.append(Points)
#				FoundDiag = True
			else:
				if_data, diag_data, FoundDiag = get_exp_data_rz_from_pulse_file(if_data, diag_data, exp_descr, exp_descr_diag)
				if(if_data == None): return

			if(FoundDiag):
				diag_data[-1].type		= exp_descr_diag.type
				diag_data[-1].par_name  = exp_descr_diag.par_name
				diag_data[-1].label     = exp_descr_diag.label
				diag_data[-1].title		= exp_descr_diag.title
				diag_data[-1].marker    = exp_descr_diag.marker
				diag_data[-1].color     = exp_descr_diag.color
				diag_data[-1].x_min		= exp_descr_diag.x_min
				diag_data[-1].x_max		= exp_descr_diag.x_max
				diag_data[-1].y_min		= exp_descr_diag.y_min
				diag_data[-1].y_max		= exp_descr_diag.y_max
				diag_data[-1].x_label	= exp_descr_diag.x_label
				diag_data[-1].y_label	= exp_descr_diag.y_label
				diag_data[-1].y_factor	= exp_descr_diag.y_factor
				diag_data[-1].x_log		= exp_descr_diag.x_log
				diag_data[-1].y_log		= exp_descr_diag.y_log

				diag_data[-1].r_min		= exp_descr_diag.r_min
				diag_data[-1].r_max		= exp_descr_diag.r_max
				diag_data[-1].z_min		= exp_descr_diag.z_min
				diag_data[-1].z_max		= exp_descr_diag.z_max
				diag_data[-1].v_min		= exp_descr_diag.v_min
				diag_data[-1].v_max		= exp_descr_diag.v_max
				diag_data[-1].v_label	= exp_descr_diag.v_label
				diag_data[-1].v_log		= exp_descr_diag.v_log

	if(if_data != None): if_data.close()

	return diag_data

#	####################
#	Read RZ data from pulse-file
#	####################

def get_exp_data_rz_from_pulse_file(if_data, diag_data, exp_descr, exp_descr_diag):

	full_name = os.getenv(exp_descr.pulse_name+"_path",'.')+"/"+exp_descr.pulse_name+"_{:d}.h5".format(exp_descr.shot)
	try:
		if_data = h5py.File(full_name.format(exp_descr.shot), "r")
	except:
		print("ERROR: does not exist in pulse file :"+full_name) 
		return if_data, diag_data, True

	print("\tFor RZ_POINTS found exp_descr_diag.type=",exp_descr_diag.type)

	if(exp_descr_diag.type == "RZ_POINTS"):
		par_names = [exp_descr_diag.signals[k].par_name for k in range(len(exp_descr_diag.signals))]
		found_names = [False for k in range(len(par_names))]
			
		itimes		 	= indexes_upper(par_names, "times", found=found_names)
		ir_data		  	= indexes_upper(par_names, "r", found=found_names)
		iz_data		  	= indexes_upper(par_names, "z", found=found_names)
		ir_as_x			= indexes_upper(par_names, "r_as_x", found=found_names)
		iz_as_x			= indexes_upper(par_names, "z_as_x", found=found_names)
		ichan_as_x		= indexes_upper(par_names, "chan_as_x", found=found_names)
		iskip_chans		= indexes_upper(par_names, "skip_chans", found=found_names)
		inchans			= indexes_upper(par_names, "nchans", found=found_names)
		ixchans		 	= indexes_upper(par_names, "xchans", found=found_names)
		ir_fun			= indexes_upper(par_names, "r_fun", found=found_names)
		iz_fun			= indexes_upper(par_names, "z_fun", found=found_names)
		idata_fun		= indexes_upper(par_names, "data_fun", found=found_names)
		idata_min	 	= indexes_upper(par_names, "data_min", found=found_names)
		idata_max	 	= indexes_upper(par_names, "data_max", found=found_names)
		ierr_max		= indexes_upper(par_names, "err_max", found=found_names)

		sig_names = [] 
		for name in par_names:
			try:
				dindex = name.upper().index("_DATA")
				if(len(name) - dindex == len("_DATA")): sig_names.append(name[:dindex])
			except:
				pass

		IsFirstName = True
		for name in sig_names:

#			Read data and parameters
#			================

			iv     = indexes_upper(par_names, name+"_data", found=found_names)

			if(len(iv) == 1):
				iv_times  		 = indexes_upper(par_names, name+"_times", found=found_names)
				iv_r_data		 = indexes_upper(par_names, name+"_r", found=found_names)
				iv_z_data		 = indexes_upper(par_names, name+"_z", found=found_names)
				iv_r_as_x		 = indexes_upper(par_names, name+"r_as_x", found=found_names)
				iv_z_as_x		 = indexes_upper(par_names, name+"z_as_x", found=found_names)
				iv_chan_as_x	 = indexes_upper(par_names, name+"chan_as_x", found=found_names)
				iv_skip_chans	 = indexes_upper(par_names, name+"_skip_chans", found=found_names)
				iv_nchans	 	 = indexes_upper(par_names, name+"_nchans", found=found_names)
				iv_xchans 		 = indexes_upper(par_names, name+"_xchans", found=found_names)
				iv_r_fun 		 = indexes_upper(par_names, name+"_r_fun", found=found_names)
				iv_z_fun 		 = indexes_upper(par_names, name+"_z_fun", found=found_names)
				iv_data_fun		 = indexes_upper(par_names, name+"_data_fun", found=found_names)
				iv_data_min	 	 = indexes_upper(par_names, name+"_data_min")
				iv_data_max	 	 = indexes_upper(par_names, name+"_data_max", found=found_names)
				iv_err			 = indexes_upper(par_names, name+"_err", found=found_names)
				iv_err_max		 = indexes_upper(par_names, name+"_err_max", found=found_names)

#				reading signal data

				if(DEBUG > 1): print("\tReading: ", exp_descr_diag.signals[iv[0]].sig_name," as data for: ",name)
				sData = h5_read_signal(if_data, exp_descr_diag.signals[iv[0]].sig_name)
				if(sData.type == -1):
					if(DEBUG > 0): print("\tInvalid data signal for: ",name)
					continue
				yData = sData.values

#				reading time signal

				if(len(iv_times) != 1): iv_times = itimes
				if(len(iv_times) == 1):
					if(DEBUG > 1): print("\tReading: ",exp_descr_diag.signals[iv_times[0]].sig_name," as times for: ",name)
					tData = h5_read_signal(if_data, exp_descr_diag.signals[iv_times[0]].sig_name)
					if(tData.type == -1):
						if(DEBUG > 0): print("\tInvalid times signal for: ",name)
						continue
					tData = tData.values
				elif(sData.type > 0):
					tData = sData.times
				else:
					if(DEBUG > 0): print("\tTimes data not available for: ",name)
					continue

				if((tData.max() < exp_descr.tstart) or (tData.min() > exp_descr.tend)):
					if(DEBUG > 1): print("\tNo data points in requested time intervall for: ",name)
					continue


#				check for y_err and reading

				if(len(iv_err) == 1): 
					if(DEBUG > 1):print("\tReading: ", exp_descr_diag.signals[iv_err[0]].sig_name," as error for: ",name)
					sData = h5_read_signal(if_data, exp_descr_diag.signals[iv_err[0]].sig_name)
					if(sData.type >= 0): yData_err = sData.values
					else:				 iv_err = []
					if((yData_err.shape[0] != yData.shape[0]) or (yData_err.shape[1] != yData.shape[1])):
						if(DEBUG > 0): print("\tMismatch between error signal and data signal for: ",name)
						continue

#				Get parameter to use for x in plot

				if(len(iv_r_as_x) != 1):					iv_r_as_x = ir_as_x
				if(len(iv_r_as_x) == 1): 					r_as_x = eval(exp_descr_diag.signals[iv_r_as_x[0]].sig_name)
				elif(exp_descr_diag.r_as_x != None):	r_as_x = exp_descr_diag.r_as_x
				else:					 					r_as_x = False

				if(len(iv_z_as_x) != 1):					iv_z_as_x = iz_as_x
				if(len(iv_z_as_x) == 1):					z_as_x = eval(exp_descr_diag.signals[iv_z_as_x[0]].sig_name)
				elif(exp_descr_diag.z_as_x != None):	z_as_x = exp_descr_diag.z_as_x
				else:										z_as_x = False

				if(len(iv_chan_as_x) != 1):						iv_chan_as_x = ichan_as_x
				if(len(iv_chan_as_x) == 1):						chan_as_x = eval(exp_descr_diag.signals[iv_chan_as_x[0]].sig_name)
				elif(exp_descr_diag.chan_as_x != None):	chan_as_x = exp_descr_diag.chan_as_x
				else:											chan_as_x = False

#				if none of r, z and nchans requested as x try to read or get x chans signal

				if((not r_as_x) and (not z_as_x) and (not chan_as_x)):
					if(len(iv_xchans) != 1): iv_xchans = ixchans
					if(len(iv_xchans) == 1):
						if(DEBUG > 1): print("\tReading: ",exp_descr_diag.signals[iv_xchans[0]].sig_name," as x chans or: ",name)
						xData = h5_read_signal(if_data, exp_descr_diag.signals[iv_xchans[0]].sig_name)
						if(xData.type == -1):
							if(DEBUG > 0): print("\tInvalid x chans signal for: ",name)
							continue
						xData = xData.values
						if(xData.shape[0] != yData.shape[1]):
							if(DEBUG > 0): print("\tMismatch xchans and data signal for: ",name)
							continue

					elif(sData.type > 0):
						xData = sData.xvalues
					else:
						if(DEBUG > 0): print("\tx chans data not available for: ",name)
						continue
				else:
					xData = np.linspace(1,yData.shape[1],num=yData.shape[1])

#				reading r signal

				if(len(iv_r_data) != 1): iv_r_data = ir_data
				if(len(iv_r_data) == 1):
					if(DEBUG > 1): print("\tReading: ",exp_descr_diag.signals[iv_r_data[0]].sig_name," as r chans or: ",name)
					rData = h5_read_signal(if_data, exp_descr_diag.signals[iv_r_data[0]].sig_name)
					if(rData.type == -1):
						if(DEBUG > 0): print("\tInvalid r chans signal for: ",name)
						continue
					rData = rData.values
					if(rData.shape[0] != yData.shape[1]):
						if(DEBUG > 0): print("\tMismatch r and data signal for: ",name)
						continue
				else:
					if(DEBUG > 0): print("\tr chans data not available for: ",name)
					continue
								
#				reading z signal

				if(len(iv_z_data) != 1): iv_z_data = iz_data
				if(len(iv_z_data) == 1):
					if(DEBUG > 1): print("\tReading: ",exp_descr_diag.signals[iv_z_data[0]].sig_name," as r chans or: ",name)
					zData = h5_read_signal(if_data, exp_descr_diag.signals[iv_z_data[0]].sig_name)
					if(zData.type == -1):
						if(DEBUG > 0): print("\tInvalid z chans signal for: ",name)
						continue
					zData = zData.values
					if(zData.shape[0] != yData.shape[1]):
						if(DEBUG > 0): print("\tMismatch z and data signal for: ",name)
						continue
				else:
					if(DEBUG > 0): print("\tz chans data not available for: ",name)
					continue
								
#				check for nchans and read

				if(len(iv_nchans) != 1): iv_nchans = inchans
				if(len(iv_nchans) == 1):
					if(DEBUG > 1):print("\tReading: ", exp_descr_diag.signals[iv_nchans[0]].sig_name," as nchans for: ",name)
					sData = h5_read_signal(if_data, exp_descr_diag.signals[iv_nchans[0]].sig_name)
					if(sData.type >= 0): NumChans  = sData.values
					else:				 iv_nchans = []

					if(NumChans.shape[0] != yData.shape[1]):
						if(DEBUG > 0): print("\tMismatch between nChans signal and data signal for: ",name)
						continue

				if(len(iv_nchans) != 1): NumChans = np.linspace(1,yData.shape[1],num=yData.shape[1]).astype(int)
					
#				check for skip_chans

				if(len(iv_skip_chans) != 1): iv_skload_chords_fileip_chans = iskip_chans
				if(len(iv_skip_chans) == 1):
					skip_chans = eval(exp_descr_diag.signals[iv_skip_chans[0]].sig_name)
					iKeepChans = np.linspace(0,NumChans.shape[0]-1,num=NumChans.shape[0]-1).astype(int)
					KeepChans  = NumChans[iKeepChans]
					for iChan in skip_chans:
						iKeepChans = iKeepChans[np.where(iChan != KeepChans)[0]]
						KeepChans  = NumChans[iKeepChans]

					if(len(KeepChans) < len(NumChans)):
						NumChans	= KeepChans
						if(len(rData.shape) > 1):
							rData = rData[:,iKeepChans]
							zData = zData[:,iKeepChans]
						else:
							rData = rData[iKeepChans]
							zData = zData[iKeepChans]
						yData = yData[:,iKeepChans]
								
						if(len(iv_err) == 1): yData_err = yData_err[:,iKeepChans]
					
#				check for err_max

				if(len(iv_err_max) != 1): iv_err_max = ierr_max
				if(len(iv_err_max) == 1): err_max = exp_descr_diag.signals[iv_err_max[0]].sig_name

#				check for data_fun and apply it

				if(len(iv_data_fun) != 1): iv_data_fun= idata_fun
				if(len(iv_data_fun) == 1):
					data_fun = exp_descr_diag.signals[iv_data_fun[0]].sig_name
					r	  = rData
					z	  = zData
					y	  = yData
					yData = eval(data_fun)
					r = z = y= 0

#				check for r_fun

				if(len(iv_r_fun) != 1): iv_r_fun= ir_fun
				if(len(iv_r_fun) == 1):
					r_fun = exp_descr_diag.signals[iv_r_fun[0]].sig_name
					r	  = zData
					rData = eval(r_fun)
					r	  = 0

#				check for z_fun

				if(len(iv_z_fun) != 1): iv_z_fun= iz_fun
				if(len(iv_z_fun) == 1):
					r_fun = exp_descr_diag.signals[iv_z_fun[0]].sig_name
					z	  = zData
					zData = eval(z_fun)
					z	  = 0

#				check for data_min

				if(len(iv_data_min) != 1): iv_data_min = idata_min
				if(len(iv_data_min) == 1): data_min = exp_descr_diag.signals[iv_data_min[0]].sig_name

#				check for data_max

				if(len(iv_data_max) != 1): iv_data_max = idata_max
				if(len(iv_data_max) == 1): data_max = exp_descr_diag.signals[iv_data_max[0]].sig_name			
				
				NotFound = 0
				for k in range(len(found_names)):
					if(not found_names[k]): NotFound += 1

				if(NotFound > 0):
					print("\tWARNING:")
					print("\tIn experimental description following parameters are not used or unknow")
					for k in range(len(found_names)):
						if(not found_names[k]): print("\t\t",par_names[k])

				sData = 0

#				Process data
#				========

				if(chan_as_x):	xData = np.copy(NumChans)
				if(r_as_x):			xData = np.copy(rData)
				if(z_as_x):			xData = np.copy(zData)

#				Average or select data
				if(exp_descr.tstart < exp_descr.tend):
					iTimes	= np.where((tData >= exp_descr.tstart) & (tData <= exp_descr.tend))[0]	
					yData	= yData[iTimes, :]
					if(len(iv_err) == 1): yData_err    = yData_err[iTimes, :]
					if(len(xData.shape) > 1): xData = xData[iTimes,:]
					if(len(rData.shape) > 1):
						rData = rData[iTimes,:]
						zData = zData[iTimes,:]

					n_chans  = np.empty(0, dtype='i4')
					x_values  = np.empty(0, dtype='f8')
					r_values  = np.empty(0, dtype='f8')
					z_values  = np.empty(0, dtype='f8')

					x_data  = np.empty(0, dtype='f8')						#rz coordinate used by experimental data
					r_data  = np.empty(0, dtype='f8')
					z_data  = np.empty(0, dtype='f8')
					y_data  = np.empty(0, dtype='f8')
					y_data_err = np.empty(0, dtype='f8')

					for iChan in range(len(NumChans)):
						yvalues = yData[:, iChan]
						iOk = np.arange(len(yvalues))

						if((len(iv_err) == 1) and (len(iv_err_max) == 1)):  iOk = iOk[np.where(yData_err[iOk, iChan] < err_max)[0]]
						if(len(iv_data_min) == 1):							iOk = iOk[np.where(yvalues[iOk] > data_min)[0]]
						if(len(iv_data_max) == 1):							iOk = iOk[np.where(yvalues[iOk] < data_max)[0]]
						if(exp_descr_diag.remove_zero):				iOk = iOk[np.where(yvalues[iOk] > 0.)[0]]

						yvalues = yvalues[iOk]
						if(len(iv_err) == 1): yerr = yData_err[iOk, iChan]

						if((exp_descr_diag.average) and (len(rData.shape) == 1)):
							SetAve = True
							if(len(yvalues) > 0):
								if(len(iv_err) == 1): 
									weights	= 1./yerr**2
									avwh	= np.average(yvalues, weights=weights, returned=True)
									yvalue= avwh[0]
									yerr   = sqrt(1./(avwh[1]/len(yerr)))
								else:
									yvalue = np.mean(yvalues)
									if(abs(yvalue) > 1.):	yerr  = np.std(yvalues/yvalue)*yvalue
									else:							yerr   = np.std(yvalues)
									if((len(iv_err_max) == 1) and (yerr > err_max)): SetAve = False

								if(SetAve):
									n_chans = np.append(n_chans,iChan)
									x_values = np.append(x_values,xData[iChan])
									r_values = np.append(r_values,rData[iChan])
									z_values = np.append(z_values,zData[iChan])

									x_data = np.append(x_data,xData[iChan])
									r_data = np.append(r_data,rData[iChan])
									z_data= np.append(z_data,zData[iChan])
									y_data = np.append(y_data,yvalue)
									y_data_err = np.append(y_data_err,yerr)
						else:
							if(len(xData.shape) == 1):
								x_values = np.append(x_values, xData[iChan])
								x_data  = np.append(x_data, np.ones(len(yvalues),dtype='f8')*xData[iChan])
							else:
								x_data  = np.append(x_data, xData[:, iChan])

							if(len(rData.shape) == 1):
								n_chans = np.append(n_chans,iChan)
								r_values = np.append(r_values, rData[iChan])
								z_values = np.append(z_values, zData[iChan])

								r_data  = np.append(r_data, np.ones(len(yvalues),dtype='f8')*rData[iChan])
								z_data  = np.append(z_data, np.ones(len(yvalues),dtype='f8')*zData[iChan])
							else:
								r_data  = np.append(r_data, np.ones(len(yvalues),dtype='i4')*iChan)
								r_data  = np.append(r_data, rData[:, iChan])
								z_data  = np.append(z_data, zData[:, iChan])
							y_data   = np.append(y_data, yvalues)
							if(len(iv_err) == 1): y_data_err = np.append(y_data_err, yerr)

					if(len(xData.shape) != 1):	x_values = x_data
					if(len(rData.shape) != 1):	r_values = r_data, z_values = z_data
				else:
					time_ave = exp_descr.tstart
					iTimes = np.argmin(np.abs(tData - time_ave))
					y_data  = yData[iTimes, :]

					iOk = np.arange(len(y_data))
					if(len(iv_err_max) == 1):						iOk = iOk[np.where(yData_err[iTimes,iOk] < err_max)[0]]
					if(len(iv_data_min) == 1):						iOk = iOk[np.where(y_data[iOk] > data_min)[0]]
					if(len(iv_data_max) == 1):					iOk = iOk[np.where(y_data[iOk] < data_max)[0]]
					if(exp_descr_diag.remove_zero):	iOk = iOk[np.where(y_data[iOk] > 0.)[0]]

					if(len(xData.shape) > 1):	x_data = xData[iTimes, iOk]
					else:						x_data = xData[iOk]
					if(len(rData.shape) > 1):
						r_data = rData[iTimes, iOk]
						z_data = zData[iTimes, iOk]
					else:
						r_data = rData[iOk]
						z_data = zData[iOk]

					y_data = y_data[iOk]
					if(len(iv_err) == 1):	y_data_err   = yData_err[iTimes,iOk]

					if(len(iv_vfactor) == 1):
						y_data	  *= values_facto
						if(len(iv_err) == 1): y_data_err *= values_factor

					n_chans = NumChans[iOk]
					x_values = x_data
					r_values = r_data
					z_values = z_data

				if(IsFirstName):
					IsFirstName = False
					
#					Define data
					diag_data.append(types.SimpleNamespace())
					diag_data[-1].SubDiags = []

				diag_data[-1].SubDiags.append(types.SimpleNamespace())
				SubDiag = diag_data[-1].SubDiags[-1]

				SubDiag.type		= exp_descr_diag.type
				SubDiag.data		= True

				SubDiag.n_chans		= n_chans
				SubDiag.x_values	= x_values
				SubDiag.r_values	= r_values
				SubDiag.z_values	= z_values
				SubDiag.x_data		= x_data
				SubDiag.r_data		= r_data
				SubDiag.z_data		= z_data
				SubDiag.y_data		= y_data
				if((len(iv_err) == 1) or exp_descr_diag.average): 
					SubDiag.err_data	= True
					SubDiag.y_err		= y_data_err
				else:
					SubDiag.err_data	= False
					SubDiag.y_err		= np.empty(0,dtype='f8')


	return if_data, diag_data, True


# 	######################
#	reads lines integrated diagnostic data
# 	######################


def get_exp_data_los_help():
	print("\nSpecific global diagnostic for line type diagnostic")
	print("\tr_min          <value> min r-axis 2D plot(D=auto)")
	print("\tr_max          <value> max r-axis 2D plot (D=auto)")
	print("\tz_min          <value> min z-axis load_chords_file2D plot (D=auto)")
	print("\tz_max          <value> max z-axis 2D plot (D=auto)")
	print("\tv_min          <value> min v-axis 2D plot (D=auto)")
	print("\tv_max          <value> max v-axis 2D plot (D=auto)")
	print("\tv_label        <value> v-axis label 2D plot(D=None)")
	print("\tv_log          <true/false> to set log v-axis 2D plot(D=False)")

	print("\nGlobal signal parameters for line type diagnostic")
	print("\tsig_times             <signal name> with times")
	print("\tsig_lines_start_end   <signal name> with start_end points lines of sigth [6,nlines] signal")
	print("\tsig_lines_start       <signal name> start lines points of sigth [3,nlines] signal")
	print("\tsig_lines_end         <signal name> start lines points of sigth [3,nlines] signal")
	print("\tsig_chan_as_x         <true/talse> to use channel numer as x in plot")
	print("\tsig_skip_chans        <array> of channels to skip (Ex. [1,4,7])")
	print("\tsig_nchans            <signal name> with channel numbers")
	print("\tsig_xchans            <signal name> with x of los")
	print("\tsig_x_fun             <function> of to modify x (Ex. 'x-5')")
	print("\tsig_data_fun          <function> to modify data (Ex. '5*y')")
	print("\tsig_data_min          <value> minimum accepted data values")
	print("\tsig_data_max          <value> maximum accepted data values")
	print("\tsig_err_max           <value> maximum accepted error on data")

	print("\nSignal                name definition and signal specific parameter")
	print("\tsig_name_data         <signal name> with data, name will be the name of parameter(Ex. sig_bolov_data)")		
	print("\tsig_name_err          <signal name> erorr, name is name in sig_data_name (Ex. sig_bolov_err)")		
	print("\tsig_name_global       where name is name sig_name_data name and global is one of globals options")

def get_exp_data_los(exp_descr, DiagType = "LINES"):

	if(DEBUG > 0): print("get_exp_data_los")

	diag_data = []
	if_data	  = None
	for exp_descr_diag in exp_descr.diags:
		if(exp_descr_diag.type == "LINES"):

			if(exp_descr_diag.chords_file != None):
				Chords = load_chords_file(exp_descr_diag.chords_file)
				diag_data.append(Chords)
				FoundDiag = True
			else:
				if_data, diag_data, FoundDiag = get_exp_data_los_from_pulse_file(if_data, diag_data, exp_descr, exp_descr_diag)
				if(if_data == None): return

			if(FoundDiag):
				diag_data[-1].type		= exp_descr_diag.type
				diag_data[-1].par_name  = exp_descr_diag.par_name
				diag_data[-1].label     = exp_descr_diag.label
				diag_data[-1].title		= exp_descr_diag.title
				diag_data[-1].marker    = exp_descr_diag.marker
				diag_data[-1].color     = exp_descr_diag.color
				diag_data[-1].x_min		= exp_descr_diag.x_min
				diag_data[-1].x_max		= exp_descr_diag.x_max
				diag_data[-1].y_min		= exp_descr_diag.y_min
				diag_data[-1].y_max		= exp_descr_diag.y_max
				diag_data[-1].x_label	= exp_descr_diag.x_label
				diag_data[-1].y_label	= exp_descr_diag.y_label
				diag_data[-1].y_factor	= exp_descr_diag.y_factor
				diag_data[-1].x_log		= exp_descr_diag.x_log
				diag_data[-1].y_log		= exp_descr_diag.y_log

				diag_data[-1].r_min		= exp_descr_diag.r_min
				diag_data[-1].r_max		= exp_descr_diag.r_max
				diag_data[-1].z_min		= exp_descr_diag.z_min
				diag_data[-1].z_max		= exp_descr_diag.z_max
				diag_data[-1].v_min		= exp_descr_diag.v_min
				diag_data[-1].v_max		= exp_descr_diag.v_max
				diag_data[-1].v_label	= exp_descr_diag.v_label
				diag_data[-1].v_log		= exp_descr_diag.v_log

	if(if_data != None): if_data.close()

	return diag_data

#	======================
#	Read line of sight data from pulse-file
#	======================

def get_exp_data_los_from_pulse_file(if_data, diag_data, exp_descr, exp_descr_diag):

	if(if_data == None):																					#Moved hera to manage non data diagnostic
		full_name = os.getenv(exp_descr.pulse_name+"_path",'.')+"/"+exp_descr.pulse_name+"_{:d}.h5".format(exp_descr.shot)
		try:
			if_data = h5py.File(full_name.format(exp_descr.shot), "r")
		except:
			print("ERROR: does not exist in pulse file :"+full_name) 
			return None, diag_data, False

	par_names = [exp_descr_diag.signals[k].par_name for k in range(len(exp_descr_diag.signals))]
	found_names = [False for k in range(len(par_names))]
		
	itimes		 	= indexes_upper(par_names, "times", found=found_names)
	ilos_st_end  	= indexes_upper(par_names, "lines_start_end", found=found_names)
	ilos_st	     	= indexes_upper(par_names, "lines_start", found=found_names)
	ilos_end     	= indexes_upper(par_names, "lines_end", found=found_names)
	ichan_as_x		= indexes_upper(par_names, "chan_as_x", found=found_names)
	iskip_chans		= indexes_upper(par_names, "skip_chans", found=found_names)
	inchans			= indexes_upper(par_names, "nchans", found=found_names)
	ixchans		 	= indexes_upper(par_names, "xchans", found=found_names)
	ix_fun			= indexes_upper(par_names, "x_fun", found=found_names)
	idata_fun		= indexes_upper(par_names, "data_fun", found=found_names)
	idata_min	 	= indexes_upper(par_names, "data_min", found=found_names)
	idata_max	 	= indexes_upper(par_names, "data_max", found=found_names)
	ierr_max		= indexes_upper(par_names, "err_max", found=found_names)

	sig_names = [] 
	for name in par_names:
		try:
			dindex = name.upper().index("_DATA")
			if(len(name) - dindex == len("_DATA")): sig_names.append(name[:dindex])
		except:
			pass

	IsFirstName = True
	for name in sig_names:

#		Read data and parameters
#		================

		iv     = indexes_upper(par_names, name+"_data", found=found_names)

		if(len(iv) == 1):
			iv_times  		 = indexes_upper(par_names, name+"_times", found=found_names)
			iv_los_st_end	 = indexes_upper(par_names, name+"_lines_start_end", found=found_names)
			iv_los_st		 = indexes_upper(par_names, name+"_lines_start", found=found_names)
			iv_los_end		 = indexes_upper(par_names, name+"_lines_end", found=found_names)
			iv_skip_chans	 = indexes_upper(par_names, name+"_skip_chans", found=found_names)
			iv_nchans	 	 = indexes_upper(par_names, name+"_nchans", found=found_names)
			iv_xchans 		 = indexes_upper(par_names, name+"_xchans", found=found_names)
			iv_x_fun 		 = indexes_upper(par_names, name+"_x_fun", found=found_names)
			iv_data_fun		 = indexes_upper(par_names, name+"_data_fun", found=found_names)
			iv_data_min	 	 = indexes_upper(par_names, name+"_data_min")
			iv_data_max	 	 = indexes_upper(par_names, name+"_data_max", found=found_names)
			iv_err			 = indexes_upper(par_names, name+"_err", found=found_names)
			iv_err_max		 = indexes_upper(par_names, name+"_err_max", found=found_names)

#			reading signal data

			if(DEBUG > 1): print("\tReading: ", exp_descr_diag.signals[iv[0]].sig_name," as data for: ",name)
			sData = h5_read_signal(if_data, exp_descr_diag.signals[iv[0]].sig_name)
			if(sData.type == -1):
				if(DEBUG > 0): print("\tInvalid data signal for: ",name)
				continue
			yData = sData.values

#			reading time signal

			if(len(iv_times) != 1): iv_times = itimes
			if(len(iv_times) == 1):
				if(DEBUG > 1): print("\tReading: ",exp_descr_diag.signals[iv_times[0]].sig_name," as times for: ",name)
				tData = h5_read_signal(if_data, exp_descr_diag.signals[iv_times[0]].sig_name)
				if(tData.type == -1):
					if(DEBUG > 0): print("\tInvalid times signal for: ",name)
					continue
				tData = tData.values
			elif(sData.type > 0):
				tData = sData.times
			else:
				if(DEBUG > 0): print("\tTimes data not available for: ",name)
				continue

			if((tData.max() < exp_descr.tstart) or (tData.min() > exp_descr.tend)):
				if(DEBUG > 1): print("\tNo data points in requested time intervall for: ",name)
				continue


#			check for y_err and read signal

			if(len(iv_err) == 1): 
				if(DEBUG > 1):print("\tReading: ", exload_chords_filep_descr_diag.signals[iv_err[0]].sig_name," as error for: ",name)
				sData = h5_read_signal(if_data, exp_descr_diag.signals[iv_err[0]].sig_name)
				if(sData.type >= 0): yData_err = sData.values
				else:				 iv_err = []

				if((yData_err.shape[0] != yData.shape[0]) or (yData_err.shape[1] != yData.shape[1])):
					if(DEBUG > 0): print("\tMismatch between error signal and data signal for: ",name)
					continue

#			reading x chans signal

			if((len(ichan_as_x) == 0) or (not exp_descr_diag.signals[ichan_as_x[0]])):
				if(len(iv_xchans) != 1): iv_xchans = ixchans
				if(len(iv_xchans) == 1):
					if(DEBUG > 1): print("\tReading: ",exp_descr_diag.signals[iv_xchans[0]].sig_name," as x chans or: ",name)
					xData = h5_read_signal(if_data, exp_descr_diag.signals[iv_xchans[0]].sig_name)
					if(xData.type == -1):
						if(DEBUG > 0): print("\tInvalid x chans signal for: ",name)
						continue
					xData = xData.values
					if(xData.shape[0] != yData.shape[1]):
						if(DEBUG > 0): print("\tMismatch xchans and data signal for: ",name)
						continue

				elif(sData.type > 0):
					xData = sData.xvalues
				else:
					if(DEBUG > 0): print("\tx chans data not available for: ",name)
					continue
			else:
				xData = np.linspace(1,yData.shape[1],num=yData.shape[1])
								
#			check for nchans and read

			if(len(iv_nchans) != 1): iv_nchans = inchans
			if(len(iv_nchans) == 1):
				if(DEBUG > 1):print("\tReading: ", exp_descr_diag.signals[iv_nchans[0]].sig_name," as nchans for: ",name)
				sData = h5_read_signal(if_data, exp_descr_diag.signals[iv_nchans[0]].sig_name)
				if(sData.type >= 0): NumChans  = sData.values
				else:				 iv_nchans = []

				if(NumChans.shape[0] != yData.shape[1]):
					if(DEBUG > 0): print("\tMismatch between nChans signal and data signal for: ",name)
					continue

			if(len(iv_nchans) != 1): NumChans = np.linspace(1,yData.shape[1],num=yData.shape[1]).astype(int)

#			lines of sight coordinates data

			if(len(iv_los_st_end) != 1): iv_los_st_end	= ilos_st_end
			if(len(iv_los_st) != 1):	 iv_los_st		= ilos_st
			if(len(iv_los_end) != 1):	 iv_los_end		= ilos_end
			if(len(iv_los_st_end) == 1): 
				if(DEBUG > 1): print("\tReading: ", exp_descr_diag.signals[iv_los_st_end[0]].sig_name," as los for: ",name)
				sData = h5_read_signal(if_data, exp_descr_diag.signals[iv_los_st_end[0]].sig_name)
				if(sData.type == -1):
					if(DEBUG > 0): print("\tInvalid lines of sight coordinates signal for: ",name)
					continue
				lines_start = sData.values[:3,:]
				lines_end   = sData.values[3:,:]

			elif((len(iv_los_st) == 1) and (len(iv_los_end) == 1)):
				if(DEBUG > 1): print("\tReading: ", exp_descr_diag.signals[iv_los_st[0]].sig_name," as los start for: ",name)
				sData = h5_read_signal(if_data, exp_descr_diag.signals[iv_los_st[0]].sig_name)
				if(sData.type == -1):
					if(DEBUG > 0): print("\tInvalid lines of sight start coordinates signal for: ",name)
					continue
				lines_start = sData.values

				if(DEBUG > 1): print("\tReading: ", exp_descr_diag.signals[iv_los_end[0]].sig_name," as los end for: ",name)
				sData = h5_read_signal(if_data, exp_descr_diag.signals[iv_los_end[0]].sig_name)
				if(sData.type == -1):
					if(DEBUG > 0): print("\tInvalid lines of sight end coordinates signal for: ",name)
					continue
				lines_end = sData.values
			else:
				if(DEBUG > 0): print("\tNot available lines of sight coordinates for: ",name)
				continue
						
			if((len(lines_start.shape) != 2)  or (len(lines_end.shape) != 2)):
				if(DEBUG > 0): print("\tInvalid lines of sight coordinates signal shape for: ",name)
				continue

			if((lines_start.shape[1] != yData.shape[1]) or (lines_end.shape[1] != yData.shape[1])):
				if(DEBUG > 0): print("\tMismatch between lines of sight coordinates and data signal for: ",name)
				continue
					
#			check for skip_chans

			if(len(iv_skip_chans) != 1): iv_skip_chans = iskip_chans
			if(len(iv_skip_chans) == 1):
				skip_chans = eval(exp_descr_diag.signals[iv_skip_chans[0]].sig_name)
				iKeepChans = np.linspace(0,NumChans.shape[0]-1,num=NumChans.shape[0]-1).astype(int)
				KeepChans  = NumChans[iKeepChans]
				for iChan in skip_chans:
					iKeepChans = iKeepChans[np.where(iChan != KeepChans)[0]]
					KeepChans  = NumChans[iKeepChans]

				if(len(KeepChans) < len(NumChans)):
					NumChans	= KeepChans
					lines_start = lines_start[:,iKeepChans]
					lines_end	= lines_end[:,iKeepChans]
					xData = xData[iKeepChans]
					yData = yData[:,iKeepChans]
					if(len(iv_err) == 1): yData_err = yData_err[:,iKeepChans]
				
#			check for err_max

			if(len(iv_err_max) != 1): iv_err_max = ierr_max
			if(len(iv_err_max) == 1): err_max = exp_descr_diag.signals[iv_err_max[0]].sig_name

#			check for data_fun and apply it

			if(len(iv_data_fun) != 1): iv_data_fun= idata_fun
			if(len(iv_data_fun) == 1):
				data_fun = exp_descr_diag.signals[iv_data_fun[0]].sig_name
				x	  = xData
				y	  = yData
				yData = eval(data_fun)
				y	  = 0
				x	  = 0

#			check for x_fun and apply it

			if(len(iv_x_fun) != 1): iv_x_fun= ix_fun
			if(len(iv_x_fun) == 1):
				x_fun = exp_descr_diag.signals[iv_x_fun[0]].sig_name
				x	  = xData
				xData = eval(x_fun)
				x	  = 0

#			check for data_min

			if(len(iv_data_min) != 1): iv_data_min = idata_min
			if(len(iv_data_min) == 1): data_min = exp_descr_diag.signals[iv_data_min[0]].sig_name

#			check for data_max

			if(len(iv_data_max) != 1): iv_data_max = idata_max
			if(len(iv_data_max) == 1): data_max = exp_descr_diag.signals[iv_data_max[0]].sig_name
					
				
			NotFound = 0
			for k in range(len(found_names)):
				if(not found_names[k]): NotFound += 1

			if(NotFound > 0):
				print("\tWARNING:")
				print("\tIn experimental descrition following parameters are not used or unknow")
				for k in range(len(found_names)):
					if(not found_names[k]): print("\t\t",par_names[k])

			sData = 0

#			Process data
#			========

			if(exp_descr_diag.chan_as_x): xData = NumChans

#			Average or select data

			x_lines = xData
			if(exp_descr.tstart < exp_descr.tend):
				iTimes	= np.where((tData >= exp_descr.tstart) & (tData <= exp_descr.tend))[0]	
				nChans  = yData.shape[1]
				yData	= yData[iTimes, :]

				y_data		= np.empty(0, dtype='f8')
				x_data		= np.empty(0, dtype='f8')
				y_data_err	= np.empty(0, dtype='f8')
				for iChan in range(nChans):
					yvalues = yData[:, iChan]

					iOk = np.arange(len(yvalues))
					if((len(iv_err) == 1) and (len(iv_err_max) == 1)):	iOk = iOk[np.where(yData_err[iOk, iChan] < err_max)[0]]
					if(len(iv_data_min) == 1):							iOk = iOk[np.where(yvalues > data_min)[0]]
					if(len(iv_data_max) == 1):							iOk = iOk[np.where(yvalues[iOk] < data_max)[0]]
					if(exp_descr_diag.remove_zero):						iOk = iOk[np.where(yvalues[iOk] > 0.)[0]]

					if(len(iOk) > 0):
						yvalues = yvalues[iOk]
						if(len(iv_err) == 1): yerrs = yData_err[iOk, iChan]
						if(exp_descr_diag.average):
							SetAve = True
							if(len(iv_err) == 1): 
								weights	= 1./yerrs**2
								avwh	= np.average(yvalues, weights=weights, returned=True)
								yvalue = avwh[0]
								yerr   = sqrt(1./(avwh[1]/len(yerrs)))
							else:
								yvalue = np.mean(yvalues)
								if(abs(yvalue) > 1.):	yerr   = np.std(yvalues/yvalue)*yvalue
								else:					yerr   = np.std(yvalues)
								if((len(iv_err_max) == 1) and (yerr > err_max)): SetAve = False
							if(SetAve):
								x_data		= np.append(x_data,xData[iChan])
								y_data		= np.append(y_data, yvalue)
								y_data_err	= np.append(y_data_err, yerr)
						else:
							x_data  = np.append(x_data, np.ones(len(yvalues),dtype='f8')*xData[iChan])
							y_data   = np.append(y_data, yvalues)
							if(len(iv_err) == 1): y_data_err = np.append(y_data_err, yerr)

			else:
				time_ave = exp_descr.tstart
				iTimes = np.argmin(np.abs(tData - time_ave))
				y_data  = yData[iTimes, :]

				iOk = np.arange(len(y_data))
				if(len(iv_err_max) == 1):				iOk = iOk[np.where( yData_err[iTimes, iOk] < err_max)[0]]
				if(len(iv_data_min) == 1):				iOk = iOk[np.where(yvalues > data_min)[0]]
				if(len(iv_data_max) == 1):				iOk = iOk[np.where(yvalues < data_max)[0]]
				if(exp_descr_diag.remove_zero):	iOk = iOk[np.where(y_data > 0.)[0]]

				x_data = xData[iOk]
				y_data = y_data[iOk]
				if(len(iv_err) == 1): yerr = yData_err[iTimes, iOk]
				if(len(iv_vfactor) == 1):
					yvalues	  *= values_factor
					if(len(iv_err) == 1): yerr *= values_factor

			if(IsFirstName):
				IsFirstName = False
					
#				Define data

				diag_data.append(types.SimpleNamespace())
				diag_data[-1].SubDiags = []

			diag_data[-1].SubDiags.append(types.SimpleNamespace())
			SubDiag = diag_data[-1].SubDiags[-1]

			SubDiag.data		= True
			SubDiag.los_name	= name

			SubDiag.n_chans   	= NumChans
			SubDiag.x_values	= x_lines
			SubDiag.x_data		= x_data
			SubDiag.y_data		= y_data

			if((len(iv_err) == 1) or exp_descr_diag.average): 
				SubDiag.err_data	= True
				SubDiag.y_err		= y_data_err
			else:
				SubDiag.err_data	= False
				SubDiag.y_err		= np.empty(0,dtype='f8')

			SubDiag.lines_start	= lines_start
			SubDiag.lines_end	= lines_end

	return if_data, diag_data, True

#	========================================================000

# 	#################################
#	reads spectra and compare with lines integrated data
# 	#################################


def get_exp_data_spectra_help():
	print("\nSpecific global diagnostic for spectra type diagnostic")
	print("\tadas_files     <adas_folder2/adas_file2>,<adas_folder2/adas_file2> adas files")
	print("\tuse_wavelength <value> flag to plot in function of wavelength or 1/wavelength (D=true)")
	print("\tlabel_lines    <value> option value to plot lines label (0=none,1=ion name, 2=ion name+wavelenght)")
	print("\tthreshold      <value> threshold value to plot lines")
	print("\tfwhm           <value> FWHM of instrumentation in Amstrong")
	print("\tw_min          <value> min wavelengh(A)(D=auto)")
	print("\tw_max          <value> max wavelengh(A)(D=auto)")
	print("\tr_min          <value> min r-axis 2D plot(D=auto)")
	print("\tr_max          <value> max r-axis 2D plot (D=auto)")
	print("\tz_min          <value> min z-axis 2D plot (D=auto)")
	print("\tz_max          <value> max z-axis 2D plot (D=auto)")

	print("\nGlobal signal parameters for spectra type diagnostic")
	print("\tsig_times             <signal name> with times")
	print("\tsig_lines_start_end   <signal name> with start_end points lines of sigth [6,nlines] signal")
	print("\tsig_lines_start       <signal name> start lines points of sigth [3,nlines] signal")
	print("\tsig_lines_end         <signal name> start lines points of sigth [3,nlines] signal")
	print("\tsig_chan_as_x         <true/talse> to use channel numer as x in plot")
	print("\tsig_skip_chans        <array> of channels to skip (Ex. [1,4,7])")
	print("\tsig_nchans            <signal name> with channel numbers")
	print("\tsig_x_fun             <function> of to modify x (Ex. 'x-5')")
	print("\tsig_data_fun          <function> to modify data (Ex. '5*y')")
	print("\tsig_data_min          <value> minimum accepted data values")
	print("\tsig_data_max          <value> maximum accepted data values")
	print("\tsig_err_max           <value> maximum accepted error on data")

	print("\nSpecific global diagnostic for spectra type diagnostic")
	print("\tr_min                 <value> min r-axis 2D plot(D=auto)")
	print("\tr_max                 <value> max r-axis 2D plot (D=auto)")
	print("\tz_min                 <value> min z-axis 2D plot (D=auto)")
	print("\tz_max                 <value> max z-axis 2D plot (D=auto)")
	print("\tv_min                 <value> min v-axis 2D plot (D=auto)")
	print("\tv_max                 <value> max v-axis 2D plot (D=auto)")
	print("\tv_label               <value> v-axis label 2D plot(D=None)")
	print("\tv_log                 <true/false> to set log v-axis 2D plot(D=False)")

	print("\nGlobal signal parameters for spectra type diagnostic")
	print("\tsig_times             <signal name> with times")
	print("\tsig_x                 <signal name> x coordinate")
	print("\tsig_skip_chans        <array> of channels to skip (Ex. [1,4,7])")
	print("\tsig_nchans            <signal name> with channel numbers")
	print("\tsig_x_fun             <function> of to modify r (Ex. 'r-5')")
	print("\tsig_data_fun          <function> to modify data (Ex. '5*y')")
	print("\tsig_data_min          <value> minimum accepted data values")
	print("\tsig_data_max          <value> maximum accepted data values")
	print("\tsig_err_max           <value> maximum accepted error on data")

	print("\nSignal                name definition and signal specific parameter")
	print("\tsig_name_data         <signal name> with data, name will be the name of parameter(Ex. sig_bolov_data)")		
	print("\tsig_name_times        <signal name> with times")		
	print("\tsig_name_x            <signal name> with x coordinate, name will be the name of parameter(Ex. sig_bolov_x)")		
	print("\tsig_name_err          <signal name> erorr, name is name in sig_data_name (Ex. sig_bolov_err)")		
	print("\tsig_name_global       where name is name sig_name_data name and global is one of globals options")

def get_exp_data_spectra(exp_descr, DiagType = "SPECTRA"):

	if(DEBUG > 0): print("get_exp_data_spectra")

	diag_data = []
	if_data	  = None
	for exp_descr_diag in exp_descr.diags:
		if(exp_descr_diag.type == "SPECTRA"):
			if(exp_descr_diag.chords_file != None):
				Chords = load_chords_file(exp_descr_diag.chords_file)
				diag_data.append(Chords)
				FoundDiag = True
			else:
				print("SPECTRA FROM PULSE-FILE NOT IMPLEMENTED")
#				if_data, diag_data, FoundDiag = get_exp_data_spectra_from_pulse_file(if_data, diag_data, exp_descr, exp_descr_diag)
				if(if_data == None): return

			if(FoundDiag):
				diag_data[-1].type				= exp_descr_diag.type
				diag_data[-1].par_name 			= exp_descr_diag.par_name
				diag_data[-1].adas_files		= exp_descr_diag.adas_files
				diag_data[-1].use_wavelength	= exp_descr_diag.use_wavelength
				diag_data[-1].label_lines		= exp_descr_diag.label_lines	
				diag_data[-1].label_threshold	= exp_descr_diag.label_threshold
				diag_data[-1].threshold_min		= exp_descr_diag.threshold_min
				diag_data[-1].threshold_max		= exp_descr_diag.threshold_max
				diag_data[-1].fwhm				= exp_descr_diag.fwhm
				diag_data[-1].w_min				= exp_descr_diag.w_min
				diag_data[-1].w_max				= exp_descr_diag.w_max
				diag_data[-1].label				= exp_descr_diag.label
				diag_data[-1].title				= exp_descr_diag.title
				diag_data[-1].marker			= exp_descr_diag.marker
				diag_data[-1].color				= exp_descr_diag.color
				diag_data[-1].x_min				= exp_descr_diag.x_min
				diag_data[-1].x_max				= exp_descr_diag.x_max
				diag_data[-1].y_min				= exp_descr_diag.y_min
				diag_data[-1].y_max				= exp_descr_diag.y_max
				diag_data[-1].x_label			= exp_descr_diag.x_label
				diag_data[-1].y_label			= exp_descr_diag.y_label
				diag_data[-1].y_factor			= exp_descr_diag.y_factor
				diag_data[-1].x_log				= exp_descr_diag.x_log
				diag_data[-1].y_log				= exp_descr_diag.y_log

				diag_data[-1].r_min			= exp_descr_diag.r_min
				diag_data[-1].r_max			= exp_descr_diag.r_max
				diag_data[-1].z_min			= exp_descr_diag.z_min
				diag_data[-1].z_max			= exp_descr_diag.z_max
				diag_data[-1].v_min			= exp_descr_diag.v_min
				diag_data[-1].v_max			= exp_descr_diag.v_max
				diag_data[-1].v_label		= exp_descr_diag.v_label
				diag_data[-1].v_log			= exp_descr_diag.v_log

	if(if_data != None): if_data.close()

	return diag_data


def indexes_upper(list_strs, str, found=None):
	iFound = []
	for i, list_st in enumerate(list_strs):
		if(list_st.upper() == str.upper()):
			iFound.append(i)
			if(found is not None): found[i] = True
	return iFound
