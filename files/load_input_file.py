import types
import numpy as np
from routines.globals			import DEBUG

#=========================================================
# This load  soledge input file
#=========================================================

def load_input_file(Path="."):

	if((len(Path) > 0) and (Path[-1] != "/")): Path = Path + "/"


	if(DEBUG > 1): print("load_soledge_input_file: reading file ",Path+"input.txt")
	try:
		fid = open(Path+"input.txt","r")
	except OSError:
		print("\tNot found file: ",Path+"input.txt")
		return 0

	input = types.SimpleNamespace()
	input.flags					= types.SimpleNamespace()
	input.global_parameters 	= types.SimpleNamespace()
	input.boundary_conditions	= types.SimpleNamespace()
	input.transport_parameters	= types.SimpleNamespace()
	input.ballooning_parameters = types.SimpleNamespace()
	input.flags					= types.SimpleNamespace()

	skip_lines(fid,1)													#Read file version
	try:
		input.file_ver = read_val(fid,"=", 1) 
	except:
		input.file_ver = 1
	skip_lines(fid,1)

	read_simulation_parameters_block(fid, input)
	read_chemistry_block(fid, input)
	read_boundary_conditions_block(fid, input)
	read_transport_parameters_block(fid, input)
	if(input.file_ver > 4): read_threshold_parameters_block(fid, input)
	read_reference_parameters_block(fid, input)

	fid.close()

	return input

#===========================================================

def read_simulation_parameters_block(fid, input):

	flags			  = input.flags
	global_parameters = input.global_parameters

	skip_lines(fid,7)
	flags.is_SLAB = read_val(fid,"=", True) 
	skip_lines(fid,2)
	flags.is_to_the_centre = read_val(fid,"=", True) 
	skip_lines(fid,2)
	flags.is_Pen = read_val(fid,"=", True) 
	skip_lines(fid,2)
	flags.solve_temperature = read_val(fid,"=", True)
	skip_lines(fid,2)

	if(input.file_ver > 1):
		flags.equalTi = read_val(fid,"=", True)
		skip_lines(fid,2)
		flags.forceNoZhdanov = read_val(fid,"=", True)
		skip_lines(fid,2)

	flags.neutral_model = read_val(fid,"=",1)
	skip_lines(fid,2)
	flags.turbulence_model = read_val(fid,"=",1)
	skip_lines(fid,2)
	flags.restart = read_val(fid,"=",True)
	skip_lines(fid,2)
	if(input.file_ver < 5):
		global_parameters.N_iterations = read_val(fid,"=",1) 
		skip_lines(fid,2)
		global_parameters.N_save = read_val(fid,"=",1) 
		skip_lines(fid,2)
		global_parameters.CFL = read_val(fid,"=",1.) 
	else:
		global_parameters.CFL = read_val(fid,"=",1.) 
		skip_lines(fid,2)
		global_parameters.N_iterations = read_val(fid,"=",1) 
		skip_lines(fid,2)
		if(input.file_ver > 5):
			global_parameters.N_save_progress = read_val(fid,"=",1) 
			skip_lines(fid,2)
		global_parameters.N_save_evolutions = read_val(fid,"=",1) 
		skip_lines(fid,2)
		global_parameters.N_cycles_points = read_val(fid,"=",1) 
		skip_lines(fid,2)
		global_parameters.N_cycles_profiles= read_val(fid,"=",1) 
	if(input.file_ver > 3):
		skip_lines(fid,2)  
		global_parameters.decomposition = read_val(fid,"=",True) 
	skip_lines(fid,2)  

#===========================================================

def read_chemistry_block(fid, input):

	global_parameters = input.global_parameters

	skip_lines(fid,5)
	global_parameters.N_species =  read_val(fid,"=",1) 

	skip_lines(fid,2)
	global_parameters.Species = fid.readline()
#	global_parameters.Species = fid.readline().split(",")
#	check_element_list()
#	compute_ion_number()
#	affect_element_and_charge_to_species()
	skip_lines(fid,2)

#===========================================================

def read_boundary_conditions_block(fid, input):

	boundary_conditions	= input.boundary_conditions

	skip_lines(fid,6)
	boundary_conditions.BCn_model = read_values(fid,"=",1)

	skip_lines(fid,4)
	boundary_conditions.BCn = read_values(fid,"=",1.)

	skip_lines(fid,3)
	boundary_conditions.BCT_model = read_values(fid,"=",1)

	skip_lines(fid,4)
	boundary_conditions.BCTe = read_val(fid,"=",2) 

	boundary_conditions.BCTi = read_values(fid,"=",1.)
	skip_lines(fid,2)

	input.boundary_conditions = boundary_conditions

#===========================================================

def read_transport_parameters_block(fid, input):

	transport_parameters  = input.transport_parameters
	ballooning_parameters = input.ballooning_parameters
	flags				  = input.flags
	skip_lines(fid,5)
	transport_parameters.Dn_p = read_values(fid,"=",1.)

	skip_lines(fid,2)
	transport_parameters.Dn_t = read_values(fid,"=",1.)

	skip_lines(fid,2)
	transport_parameters.nu_p = read_values(fid,"=",1.)

	skip_lines(fid,2)
	transport_parameters.nu_t = read_values(fid,"=",1.)

	skip_lines(fid,2)
	transport_parameters.chie_p = read_val(fid,"=",1.) 

	skip_lines(fid,2)
	transport_parameters.chie_t = read_val(fid,"=",1.)
	skip_lines(fid,2)
	transport_parameters.chii_p = read_values(fid,"=",1.)

	skip_lines(fid,2)
	transport_parameters.chii_t = read_values(fid,"=",1.)

	skip_lines(fid,2)
	transport_parameters.pinch_model = read_values(fid,"=",1)

	skip_lines(fid,2)
	transport_parameters.v_pinch = read_values(fid,"=",1.)

	skip_lines(fid,2)
	ballooning_parameters.ballooning_model = read_val(fid,"=",1)

	if(ballooning_parameters.ballooning_model == 3): flags.radialFeedback = True
	else:											 flags.radialFeedback = False

	skip_lines(fid,2)
	ballooning_parameters.zbal = read_val(fid,"=",1.) 

	skip_lines(fid,1)
	ballooning_parameters.minmaxbal = read_val(fid,"=",1.) 

	skip_lines(fid,1)
	ballooning_parameters.sigmabal = read_val(fid,"=",1) 

	skip_lines(fid,8)
	transport_parameters.Coulomb_log = read_val(fid,"=",1.) 

	skip_lines(fid,2)
	transport_parameters.delta_e = read_val(fid,"=",1.) 

	if(input.file_ver > 1):
		skip_lines(fid,3)
		transport_parameters.pseudoConv_coef= read_val(fid,"=",1.) 

	skip_lines(fid,2)
	transport_parameters.gamma_i = read_val(fid,"=",1.) 

	skip_lines(fid,3)
	transport_parameters.Flux_limiter0 = read_val(fid,"=",1.)
	transport_parameters.Flux_limiter  = read_values(fid,"=",1.)
	transport_parameters.Flux_limiter  = np.append(transport_parameters.Flux_limiter0, transport_parameters.Flux_limiter)

	if(input.file_ver > 1):
		transport_parameters.Flux_limiter_Zh  = read_values(fid,"=",1.)
		transport_parameters.Flux_limiter_Zh  = np.append(transport_parameters.Flux_limiter_Zh[0], transport_parameters.Flux_limiter_Zh)

	transport_parameters.Flux_limiter_nu = read_values(fid,"=",1) 

	skip_lines(fid,2)
	flags.non_zero_parallel_viscosity = read_val(fid,"=", True)

	if(input.file_ver == 3):
		skip_lines(fid,2)
		transport_parameters.R_Coupling_Fact  = read_values(fid,"=",1.)
		transport_parameters.R_Coupling_Fact  = np.append(transport_parameters.R_Coupling_Fact[0], transport_parameters.R_Coupling_Fact)

	skip_lines(fid,2)

#===========================================================

def read_threshold_parameters_block(fid, input):
	input.threshold_parameters			= types.SimpleNamespace()

	skip_lines(fid,5)
	input.threshold_parameters.coupling_min_density	= read_values(fid,"=",1.) 
	skip_lines(fid,2)
	input.threshold_parameters.min_density	= read_values(fid,"=",1.) 
	skip_lines(fid,2)
	input.threshold_parameters.min_temperature	= read_val(fid,"=",1.) 
	skip_lines(fid,2)

#===========================================================

def read_reference_parameters_block(fid, input):
	input.reference_parameters			= types.SimpleNamespace()
	input.reference_parameters.fields	= types.SimpleNamespace()

	skip_lines(fid,5)
	input.reference_parameters.fields.n0	= read_val(fid,"=",1.) 

	skip_lines(fid,2)
	input.reference_parameters.fields.T0eV	= read_val(fid,"=",1.) 


#===========================================================
#===========================================================

def skip_lines(fid, nLines):
	for k in range(nLines): 
		tmp = fid.readline()

def read_val(fid, Sep, DType):
	tmp   = fid.readline()
	value = tmp.split(Sep)[1]
	if(isinstance(DType,bool)):
		if(value == "T"):	return True
		else:				return False

	elif (isinstance(DType,int)):
		return int(eval(value))

	elif (isinstance(DType,float)):
		return float(eval(value))

	elif (isinstance(DType,str)):
		return value

	else:
		print("\tunknown reuested data type")


def read_values(fid, Sep, DType):
	tmp = fid.readline()
	values = tmp.split(Sep)[1]
	values = values.split(",")
	try:
		if(isinstance(DType,bool)):
			if(value == "T"):	return True
			else:				return False

		elif (isinstance(DType,int)):
			return np.array([int(eval(k)) for k in values], dtype='i4')

		elif (isinstance(DType,float)):
			return np.array([float(eval(k)) for k in values], dtype='f8')

		elif (isinstance(DType,str)):
			return values

		else:
			print("\tunknown requested data type")
	except:
			print("\terror parsing: ",tmp)
