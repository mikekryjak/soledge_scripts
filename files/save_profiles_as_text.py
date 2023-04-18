import os
import numpy				as np
from routines.globals		import DEBUG

#====================================
# This routine save profiles as text
#====================================

def save_profiles_as_text(ProfFile, Cut, xName):

	if(DEBUG > 0): print("save_profiles_as_text: Saving to ",ProfFile)

	Profiles = Cut.Profiles

	FileExt = os.path.splitext(os.path.split(ProfFile)[1])[1]

	if(FileExt == ".csv"):
		Delimiter = ", "
	elif(FileExt == ".tsv"):
		Delimiter = "\t"
	else:
		Delimiter = " "

	nMax = 0
	Header = "" 
	for i in range(len(Profiles)):
		nMax = max(nMax, len(Profiles[i].xValues))
		Header = Header + xName + "_" + Profiles[i].Name + Delimiter +  Profiles[i].Name + Delimiter

	Header = Header[:-len(Delimiter)]

	ProfilesData = np.zeros((nMax, 2*len(Profiles)), dtype='f8')
	for i in range(len(Profiles)):
		ProfilesData[:len(Profiles[i].xValues), 2*i]   = Profiles[i].xValues
		ProfilesData[:len(Profiles[i].xValues), 2*i+1] = Profiles[i].Values
		if(len(Profiles[i].xValues) < nMax):
			ProfilesData[len(Profiles[i].xValues):, 2*i]   = np.nan
			ProfilesData[len(Profiles[i].xValues):, 2*i+1] = np.nan
	

	np.savetxt(ProfFile, ProfilesData, header=Header, delimiter=Delimiter, fmt="%15.7e", comments="")
