from routines.f90nml_parser		import Parser

#	quick routine for reading a fortran namelist file
 
def load_f90nml(nml_fname):

	return Parser().read(nml_fname)