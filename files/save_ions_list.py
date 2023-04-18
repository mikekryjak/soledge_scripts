
def save_ions_list(path=".", ions=["e-","D+"]):

	fid = open(path+'/ions_list','w')

	for ion in ions: fid.write(" "+ion+"\n")

	fid.close()
