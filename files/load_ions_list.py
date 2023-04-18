
def load_ions_list(Path="."):

	if((len(Path) > 0) and (Path[-1] != "/")): Path = Path + "/"

	fid = open(Path+'ions_list','r')
	ions =[]
	while (True):
		ion = fid.readline()
		if(len(ion) == 0): break
		ions.append(ion[1:-1])

	fid.close()

	return ions

	
