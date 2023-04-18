import os
import numpy								as np
from routines.set_profile_from_filedata		import set_cut_profile_from_filedata
from routines.globals						import DEBUG,LENGTH_TYPES
from files.load_text_data					import load_text_data


#====================================
# This routine save profiles as text
#====================================

def load_profiles_as_text(ProfFile, Cut, xName, OutxOffset=0., InxOffset=0.):

	if(DEBUG > 0):
		print("load_profiles_as_text: Reading from ",ProfFile)
		print("\tValid names of x coordinates: ", LENGTH_TYPES)
			

	Headers, FileData = load_text_data(ProfFile)

	if(len(Headers) == 0): return

#	print("Headers=",Headers)

	Profiles = Cut.Profiles

	if(hasattr(Cut, "Lengths")):
		Lengths  = Cut.Lengths
	else:
		Lengths  = np.empty((0,0), dtype='f8')

	for i in range(len(Profiles)):
		if(DEBUG > 0): print("\tSearch for following data:    ", Profiles[i].Name)
		set_cut_profile_from_filedata(Headers, FileData, Profiles[i], xName, Lengths, OutxOffset=OutxOffset, InxOffset=InxOffset)

	if(DEBUG > 0): print("load_profiles_as_text: Completed")

	return
