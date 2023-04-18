import os
import types
import numpy					as np
from files.load_f90nml			import load_f90nml
from routines.globals			import DEBUG

#=========================================================
# This routine read a experimental description
#=========================================================

def load_exp_descr_help():
	print("\nGlobal parameters")
	print("\tpulse_type     <name> pulse file type ('H5')")
	print("\tpulse_name     <name> pulse name (Ex. 'jet')")
	print("\tshot           <value> pulse numeber")
	print("\ttstart         <value> start time data interval")
	print("\ttend           <value> end time data interval")
	print("\taverage        <true/false> to average data")
	print("\tremove_zero    <true/flase> to remove zero values")
	print("\tmarker         <name> marker symbol (Ex. 'o')")
	print("\tcolor          <name> color name (Ex. 'r')")
	print("\nGlobal diagnostic parameters")
	print("\taverage        <true/false> to average data")
	print("\tremove_zero    <true/false> to remove zero values")
	print("\tmarker         <name> marker symbol (Ex. 'o')")
	print("\tcolor          <name> color name (Ex. 'r')")
	print("\tlabel          <name> label string")
	print("\tx_min          <value> min x-axis (D=auto)")
	print("\tx_max          <value> max x-axis (D=auto)")
	print("\ty_min          <value> min y-axis (D=auto)")
	print("\ty_max          <value> max y-axis (D=auto)")
	print("\tx_label        <value> x-axis label (D=None)")
	print("\ty_label        <value> y-axis label (D=None)")
	print("\ty_factor       <value> factor on values (D=1.)")
	print("\tx_log          <true/false> to set log x-axis (D=False)")
	print("\ty_log          <true/false> to set log y-axis (D=False)")

def load_exp_descr(ExpFile):

	FileExt = os.path.splitext(os.path.split(ExpFile)[1])[1]

	if(len(FileExt) == 0): FileName = ExpFile + ".exp"
	else:				   FileName = ExpFile
	
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

	if(len(FileName) == 0):
		print("load_exp_descr: ERROR not available file: ",ExpFile)
		exit()
		
	exp_file_descr = load_f90nml(FileName)
	exp_file_keys  = exp_file_descr.keys_to_strings()

	exp_descr = types.SimpleNamespace()

	exp_descr.pulse_type	= get_from_dict(exp_file_descr, 'globals', 'pulse_type', "H5")
	exp_descr.pulse_name	= get_from_dict(exp_file_descr, 'globals', 'pulse_name', "SOLEDGE")
	exp_descr.shot			= get_from_dict(exp_file_descr, 'globals', 'shot', 0.)
	exp_descr.tstart		= get_from_dict(exp_file_descr, 'globals', 'tstart', 0.)
	exp_descr.tend			= get_from_dict(exp_file_descr, 'globals', 'tend', 0.)
	exp_descr.average		= get_from_dict(exp_file_descr, 'globals', 'average', False)
	exp_descr.remove_zero	= get_from_dict(exp_file_descr, 'globals', 'remove_zero', False)
	exp_descr.marker		= get_from_dict(exp_file_descr, 'globals', 'marker', "o")
	exp_descr.color  		= get_from_dict(exp_file_descr, 'globals', 'color', "r")
	exp_descr.label  		= get_from_dict(exp_file_descr, 'globals', 'label', "")
	exp_descr.x_min  		= get_from_dict(exp_file_descr, 'globals', 'x_min', None)
	exp_descr.x_max  		= get_from_dict(exp_file_descr, 'globals', 'x_max', None)
	exp_descr.y_min  		= get_from_dict(exp_file_descr, 'globals', 'y_min', None)
	exp_descr.y_max  		= get_from_dict(exp_file_descr, 'globals', 'y_max', None)
	exp_descr.x_label  		= get_from_dict(exp_file_descr, 'globals', 'x_label', "")
	exp_descr.y_label  		= get_from_dict(exp_file_descr, 'globals', 'y_label', "")
	exp_descr.y_factor 		= get_from_dict(exp_file_descr, 'globals', 'y_factor', None)
	exp_descr.x_log 		= get_from_dict(exp_file_descr, 'globals', 'x_log', None)
	exp_descr.y_log 		= get_from_dict(exp_file_descr, 'globals', 'y_log', None)
	
	exp_descr.diags = []
	for exp_key in exp_file_keys:

		if(exp_key != "globals"):
			is_on = get_from_dict(exp_file_descr, exp_key, 'is_on', True)
			if(not is_on): continue

			exp_descr.diags.append(types.SimpleNamespace())
			exp_descr.diags[-1].type		= get_from_dict(exp_file_descr, exp_key, 'type', "POINT").upper()
			exp_descr.diags[-1].marker		= get_from_dict(exp_file_descr, exp_key, 'marker', exp_descr.marker)
			exp_descr.diags[-1].color  		= get_from_dict(exp_file_descr, exp_key, 'color', exp_descr.color)
			exp_descr.diags[-1].label		= get_from_dict(exp_file_descr, exp_key, 'label', exp_descr.label)
			exp_descr.diags[-1].title	  	= get_from_dict(exp_file_descr, exp_key, 'title', "")
			exp_descr.diags[-1].remove_zero = get_from_dict(exp_file_descr, exp_key, 'remove_zero', exp_descr.remove_zero)
			exp_descr.diags[-1].average		= get_from_dict(exp_file_descr, exp_key, 'average', exp_descr.average)
			exp_descr.diags[-1].x_min  		= get_from_dict(exp_file_descr, exp_key, 'x_min', exp_descr.x_min)
			exp_descr.diags[-1].x_max  		= get_from_dict(exp_file_descr, exp_key, 'x_max', exp_descr.x_max)
			exp_descr.diags[-1].y_min  		= get_from_dict(exp_file_descr, exp_key, 'y_min', exp_descr.y_min)
			exp_descr.diags[-1].y_max  		= get_from_dict(exp_file_descr, exp_key, 'y_max', exp_descr.y_max)
			exp_descr.diags[-1].x_label  	= get_from_dict(exp_file_descr, exp_key, 'x_label', exp_descr.x_label)
			exp_descr.diags[-1].y_label  	= get_from_dict(exp_file_descr, exp_key, 'y_label', exp_descr.y_label )
			exp_descr.diags[-1].y_factor 	= get_from_dict(exp_file_descr, exp_key, 'y_factor', exp_descr.y_factor)
			exp_descr.diags[-1].x_log  		= get_from_dict(exp_file_descr, exp_key, 'x_log', exp_descr.x_log)
			exp_descr.diags[-1].y_log  		= get_from_dict(exp_file_descr, exp_key, 'y_log', exp_descr.y_log)

			sub_exp_keys					= exp_file_descr[exp_key].keys_to_strings()

#			load diags specific paramters

			if((exp_descr.diags[-1].type == "LINES") or (exp_descr.diags[-1].type == "SPECTRA") or (exp_descr.diags[-1].type == "RZ_POINTS")):

				ii = indexes_upper(sub_exp_keys, "par_name")
				if(len(ii) == 1):
					exp_descr.diags[-1].par_name = exp_file_descr[exp_key]['par_name']
					del sub_exp_keys[ii[0]]
				else:
					print("ERROR: par_name is needed for ",exp_descr.diags[-1].type," diagnostic: ",exp_key)
					exit()

				exp_descr.diags[-1].chan_as_x	= get_from_dict(exp_file_descr, exp_key, 'chan_as_x', False)
				exp_descr.diags[-1].r_min		= get_from_dict(exp_file_descr, exp_key, 'r_min', None)
				exp_descr.diags[-1].r_max		= get_from_dict(exp_file_descr, exp_key, 'r_max', None)
				exp_descr.diags[-1].z_min		= get_from_dict(exp_file_descr, exp_key, 'z_min', None)
				exp_descr.diags[-1].z_max		= get_from_dict(exp_file_descr, exp_key, 'z_max', None)
				exp_descr.diags[-1].v_min		= get_from_dict(exp_file_descr, exp_key, 'v_min', None)
				exp_descr.diags[-1].v_max		= get_from_dict(exp_file_descr, exp_key, 'v_max', None)
				exp_descr.diags[-1].v_label 	= get_from_dict(exp_file_descr, exp_key, 'v_label', "")
				exp_descr.diags[-1].v_log		= get_from_dict(exp_file_descr, exp_key, 'v_log', False)

			if(exp_descr.diags[-1].type == "RZ_POINTS"):
				exp_descr.diags[-1].points_file	= get_from_dict(exp_file_descr, exp_key, 'points_file', None)
				exp_descr.diags[-1].r_as_x	= get_from_dict(exp_file_descr, exp_key, 'r_as_x', None)
				exp_descr.diags[-1].z_as_x	= get_from_dict(exp_file_descr, exp_key, 'z_as_x', None)

			if((exp_descr.diags[-1].type == "LINES") or (exp_descr.diags[-1].type == "SPECTRA")):
				exp_descr.diags[-1].chords_file	= get_from_dict(exp_file_descr, exp_key, 'chords_file', None)

			if(exp_descr.diags[-1].type == "SPECTRA"):
				exp_descr.diags[-1].adas_files			= get_from_dict(exp_file_descr, exp_key, 'adas_files', None)
				exp_descr.diags[-1].fwhm				= get_from_dict(exp_file_descr, exp_key, 'fwhm', None)
				exp_descr.diags[-1].w_min				= get_from_dict(exp_file_descr, exp_key, 'w_min', None)
				exp_descr.diags[-1].w_max				= get_from_dict(exp_file_descr, exp_key, 'w_max', None)
				exp_descr.diags[-1].use_wavelength		= get_from_dict(exp_file_descr, exp_key, 'use_wavelength', True)
				exp_descr.diags[-1].label_lines			= get_from_dict(exp_file_descr, exp_key, 'label_lines', 0)
				exp_descr.diags[-1].label_threshold		= get_from_dict(exp_file_descr, exp_key, 'label_threshold', 0.)
				exp_descr.diags[-1].threshold_min		= get_from_dict(exp_file_descr, exp_key, 'threshold_min', 0.)
				exp_descr.diags[-1].threshold_max		= get_from_dict(exp_file_descr, exp_key, 'threshold_max', 0.)
				if(len(exp_descr.diags[-1].adas_files[0]) == 1): exp_descr.diags[-1].adas_files = [exp_descr.diags[-1].adas_files]


#			Load signals
				
			exp_descr.diags[-1].signals = []
			for sub_exp_key in sub_exp_keys:
				if(sub_exp_key[:4] == "sig_"):
					exp_descr.diags[-1].signals.append(types.SimpleNamespace())
					exp_descr.diags[-1].signals[-1].par_name = sub_exp_key[4:]
					exp_descr.diags[-1].signals[-1].sig_name = exp_file_descr[exp_key][sub_exp_key]
						
	return exp_descr


def get_from_dict(dict_data, section, name, def_data):
	try:
		return dict_data[section][name]
	except:
		return def_data

def indexes_upper(list_strs, str):
	return  [i for i, j in enumerate(list_strs) if j.upper() == str.upper()]
