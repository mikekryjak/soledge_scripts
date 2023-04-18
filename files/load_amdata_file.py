import os
import types
import numpy				as np
from routines.globals		import DEBUG

#=========================================================
# This routine read amdata
#=========================================================


def load_amdata_file(specie, data_type="ionization"):

	if(data_type == "ionization"):
		FileName="scd_"+specie+".dat"
	elif(data_type == "recombination"):
		FileName="acd_"+specie+".dat"
	elif(data_type == "line_excitation"):
		FileName="plt_"+specie+".dat"
	elif(data_type == "line_recombination"):
		FileName="prb_"+specie+".dat"
	else:
		print("\tERROR: Unrecognised data_tye=",data_type)
		print("\t\tallowed values are: ionization, recombination,line_excitation,line_recombination")

	fid  = open(FileName,'r')

	amdata = types.SimpleNamespace()

	skip_lines(fid,7)
	amdata.Element = read_val(fid,":","") 
	amdata.Z = read_val(fid,"=",1) 

	line = fid.readline()
	amdata.izmin  = eval(line[14:16])
	amdata.izmax  = eval(line[24:27])
#	print("Element,Z,izmin,izmax=",amdata.Element,amdata.Z,amdata.izmin,amdata.izmax)

	skip_lines(fid,1)
	amdata.degree  = read_val(fid,"=",1) 
	skip_lines(fid,1)

	amdata.nNe = read_val(fid,"=",1) 
	amdata.Ne_min = read_val(fid,"=",1.) 
	amdata.Ne_max = read_val(fid,"=",1.) 

#	print("NNe,Ne_min,Ne_max=",amdata.NNe,amdata.Ne_min,amdata.Ne_max)

	ncol=5
	nl=int(amdata.nNe/ncol)
	if(nl*ncol < amdata.nNe): nl += 1
	skip_lines(fid,1+nl)

#	now read temperature data
	amdata.iz_data=[]
	for iz in range(amdata.izmin,amdata.izmax+1):
		skip_lines(fid,3)
		amdata.iz_data.append(types.SimpleNamespace())
		amdata.iz_data[-1].nTe = read_val(fid,"=",1) 
		amdata.iz_data[-1].Te_min = read_val(fid,"=",1.) 
		amdata.iz_data[-1].Te_max = read_val(fid,"=",1.) 

#		print("iz,nTe,Te_min,Te_max=",iz,amdata.iz_data[-1].nTe,amdata.iz_data[-1].Te_min,amdata.iz_data[-1].Te_max)

		nl=int(amdata.iz_data[-1].nTe/ncol)
		if(nl*ncol < amdata.iz_data[-1].nTe): nl += 1
		skip_lines(fid,1+nl+1)

#		number of fitting coefficients

		nf=int((amdata.degree+1)*amdata.degree/2)
		nl=int(nf/ncol)
		if(nl*ncol < nf): nl += 1
		len_values = 19
		amdata.iz_data[-1].coefficients=np.empty(nf,dtype='f8')
		for il in range(nl):
			line = fid.readline()
			n_line = min(5,nf-il*ncol)
			values = []
			for k in range(n_line): values.append(eval(line[len_values*k:len_values*(k+1)]))
			amdata.iz_data[-1].coefficients[il*ncol:il*ncol+n_line] =values
#		print("coefficients",amdata.iz_data[-1].coefficients)

		skip_lines(fid,4)

	fid.close

	return amdata

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
