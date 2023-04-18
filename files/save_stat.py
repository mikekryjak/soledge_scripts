import os

def save_stat(filename, header, data_array, format_array):

	if(not os.path.isfile(filename)):
		fid = open(filename,"w")
		fid.write(header+"\n")
	else:
		fid = open(filename,"a")
		
	for i  in range(len(data_array)): fid.write(format_array[i].format(data_array[i]))
	fid.write("\n")
	
	fid.close
