import numpy			as np
from routines.globals	import DEBUG,LENGTH_TYPES

def set_cut_profile_from_filedata(Headers, FileData, Profile, xName, Lengths, OutxOffset=0., InxOffset=0., extrapolate = True):

#	xOffset shift of xValues not transfered to Lengths

	Profile.xValues, ValuesNew = set_profile_from_filedata(Headers, FileData, xName, Profile.Name, Profile.xValues, Lengths, OutxOffset=OutxOffset, InxOffset=InxOffset, extrapolate = extrapolate)
	if(ValuesNew.size == Profile.xValues.size): Profile.Values =  ValuesNew
	return

def set_profile_from_filedata(Headers, FileData, xName, yName, xInValues, Lengths, OutxOffset=0., InxOffset=0., extrapolate = False, NameExt=""):

	if(DEBUG > 0): print("set_profile_from_filedata")

	Empty = np.empty(0, dtype='i4')												#Empty array to return when unable to find a right array

	OutxOff	 =  0.
	for k in range(len(LENGTH_TYPES)):
		if(xName.upper() == LENGTH_TYPES[k]):
			OutXType = k-1
			break
	if(OutXType == -1):
		OutXType = 0
		OutxOff	 = OutxOffset


	InXType = -2
	InxOff  = 0. 
	for k in range(len(LENGTH_TYPES)):
		ix = indexes_upper(Headers, LENGTH_TYPES[k]+"_" + yName)
		if((len(NameExt) > 0) and (len(ix) == 0)): ix = indexes_upper(Headers, LENGTH_TYPES[k]+"_" + NameExt)
		if(len(ix) == 0): ix = indexes_upper(Headers, LENGTH_TYPES[k])
		if(len(ix) != 0):
			InXType = k -1
			break

	if(InXType == -1): 
		InXType = 0
		if(InxOffset != 0.): InxOff	= InxOffset
		else:  				 InxOff	= OutxOffset

	if(xInValues.size == 0): xDefValues = Lengths[:,OutXType]
	else:					 xDefValues = xInValues - OutxOff												#Transform Dist in Length

	if(len(NameExt) > 0): iv = indexes_upper(Headers, yName+"_" + NameExt)
	else:				  iv = []
	if(len(iv) == 0):	  iv = indexes_upper(Headers, yName)

	if((len(ix) == 1) and (len(iv) == 1)):
		ix = ix[0]
		iv = iv[0]

		okRange = np.where(np.isfinite(FileData[:,ix]) & np.isfinite(FileData[:,iv]))[0]
		xData = np.copy(FileData[okRange,ix]) - InxOff 													#Copy and transform Dist in Length

		vData = np.copy(FileData[okRange,iv])

		return interpolate_on_length(xData, vData, xDefValues, Lengths, OutXType, InXType, OutxOff, extrapolate = extrapolate, yName=yName)
	else:
		print("\tWARNING!!: in Headers no fields for")
		if(len(ix) == 0): print("\t\tx: ", LENGTH_TYPES," or with ###_"+yName)
		if(len(iv) == 0): print("\t\ty: ",yName)
		print("\tHeaders: ",Headers)

		return  xDefValues + OutxOff, Empty

	if(DEBUG > 0): print("set_profile_from_filedata: Completed")

#	#########################################

def interpolate_on_length(xData, vData, xDefValues, Lengths, OutXType, InXType, OutxOff, extrapolate = True, yName=""):

	Empty = np.empty(0, dtype='i4')												#Empty array to return when unable to find a right array

	if(OutXType != InXType):
		if(Lengths.shape[0] > 1):
			ixData = np.where(np.isfinite(Lengths[:,InXType]) & np.isfinite(Lengths[:,OutXType]))[0]
			if(len(ixData) > 1):
				InConvRange	 = np.where((xData  >= Lengths[ixData,InXType].min()) &  (xData <= Lengths[ixData,InXType].max()))[0]
				if(DEBUG > 2): 
					print("length_types[InXType]        =",LENGTH_TYPES[InXType+1])
					print("min(Lengths[ixData,InXType]) =",Lengths[ixData,InXType].min())
					print("max(Lengths[ixData,InXType]) =",Lengths[ixData,InXType].max())
					print("length_types[OutXType]       =",LENGTH_TYPES[OutXType+1])
					print("min(Lengths[ixData,OutXType])=",Lengths[ixData,OutXType].min())
					print("max(Lengths[ixData,OutXType])=",Lengths[ixData,OutXType].max())
				if(len(InConvRange) > 0):
					xData = np.copy(xData[InConvRange])
					vData = np.copy(vData[InConvRange])

					xData = np.interp(xData, Lengths[ixData,InXType], Lengths[ixData,OutXType])						#If needed transform x
					if(DEBUG > 1): print("\tconverted xdata for ", yName)
				else:
					if(DEBUG > 0):
						print("\tWARNING!!: No data in conversion range for " + yName)
						print("\t\tInput type  "+LENGTH_TYPES[InXType+1])
						print("\t\tOutput type "+LENGTH_TYPES[OutXType+1])
						print("\t\tMin(Lengths[ixData,InXType]) =",Lengths[ixData,InXType].min())
						print("\t\tMax(Lengths[ixData,InXType]) =",Lengths[ixData,InXType].max())
						print("\t\tMin(xData)                   =",xData.min())
						print("\t\tMax(xData)                   =",xData.max())
						print("\t\t================================")
						print("\t\t\txData                       =",xData)
						print("\t\t================================")
						print("\t\tLengths[ixData,InXType]      =",Lengths[ixData,InXType])
							
					return xDefValues + OutxOff, Empty

			else:
				if(DEBUG > 0): print("\tERROR!!!: No data from mesh for requested parameter. Check code")
				return xDefValues + OutxOff, Empty
		else:
			if(DEBUG > 0): print("\tERROR!!!: No data from mesh for requested parameter. Check code")
			return xDefValues + OutxOff, Empty

	InRange	= np.where((xData  >= xDefValues[0]) &  (xData <= xDefValues[-1]))[0]
	if(len(InRange) > 0):
		if(extrapolate):
			xValuesNew		 = np.empty((len(InRange)+2), dtype='f8')
			ValuesNew		 = np.empty((len(InRange)+2), dtype='f8')
											
			xValuesNew[1:-1] = xData[InRange]
			ValuesNew[1:-1]  = vData[InRange]
				
			xValuesNew[0]	 = xDefValues[0]
			xValuesNew[-1]	 = xDefValues[-1]
				
			if(InRange[0] > 0):
				ValuesNew[0] = vData[InRange[0]-1] + (vData[InRange[0]] - vData[InRange[0]-1])*(xDefValues[0] - xData[InRange[0]-1])/(xData[InRange[0]] - xData[InRange[0]-1])
			else:
				ValuesNew[0] = vData[0]

			if(InRange[-1] < len(xData)-1):
				ValuesNew[-1] = vData[InRange[-1]] + (vData[InRange[-1]+1] - vData[InRange[-1]])*(xDefValues[-1] - xData[InRange[-1]])/(xData[InRange[-1]+1] - xData[InRange[-1]])
			else:
				ValuesNew[-1]	 = vData[-1]

			if(abs(xValuesNew[1] - xValuesNew[0]) < 1.e-6*abs(xValuesNew[-1] - xValuesNew[0])):
				xValuesNew = np.append(xValuesNew[0],xValuesNew[2:])
				ValuesNew  = np.append(ValuesNew[0],  ValuesNew[2:])

			if(abs(xValuesNew[-1] - xValuesNew[-2]) < 1.e-6*abs(xValuesNew[-1] - xValuesNew[0])):
				xValuesNew = np.append(xValuesNew[:-2],xValuesNew[-1])
				ValuesNew  = np.append(ValuesNew[:-2],  ValuesNew[-1])

		else:
			xValuesNew = np.copy(xData[InRange])
			ValuesNew  = np.copy(vData[InRange])
		return xValuesNew + OutxOff, ValuesNew
	else:
		if(DEBUG > 1):
			print("\tINFO: No data in range for " + yName)
			print("\txDefValues[0]  =",xDefValues[0])
			print("\txDefValues[-1] =",xDefValues[-1])
			print("\tMin(xData)     =",xData.min())
			print("\tMax(xData)     =",xData.max())
		return  xDefValues + OutxOff, Empty


def indexes_upper(list_strs, str):
	return  [i for i, j in enumerate(list_strs) if j.upper() == str.upper()]
