#!/usr/bin/python

# Function definition is here

import numpy as np

#=========================================================
# This parse comands to get value
#=========================================================

def cli_get_value(qualifier, args, default=0):

	for i in range(len(args)):
		if(args[i] == qualifier): 
			if(i < len(args) -1):
				if(args[i+1] != "-"): return cli_cvt_value(args[i+1], default)

			print("cli_get_value: option without value: ",qualifier)
			return default

	return default

#=========================================================
# This check if a qualifier is present
#=========================================================

def cli_present(qualifier, args):

	for arg in args:
		if(arg == qualifier): return 1

	return 0



#===============================================================
# This routine try to convert values to the most appropiate one
# -		List are converted in numpy array
# -		string containing number are converted in number
# -		float equal to int are converted in int
#
#==============================================================

def cli_cvt_value(value, default=0):

	if(isinstance(default, str)): return value
	elif(isinstance(default, list)):
		if(len(default) == 0):
			try:
				value =  eval(value)
			except:
				pass
			if(isinstance(value, list)): return value
			else:						 return [value]
		else:
			if(isinstance(default[0], str)):
				if((value[0] == "'") or (value[0] == '"') or (value[0] == '[')):
					return value[1:-1].split(",")
				else:
					return value.split(",")
			elif(isinstance(default[0], float)):
				value =  eval(value)
				try:
					fvalue = []
					for f in value: fvalue.append(float(f))
					return fvalue
				except TypeError:
					return value
			elif(isinstance(default[0], int)):
				value =  eval(value)
				try:
					ivalue = []
					for i in value: ivalue.append(int(i))
					return ivalue
				except TypeError:
					return value
			else:
				print("cli_cvt_value: not manage value class=",value.__class__)
				return value
				
	elif(isinstance(default, float)):
		try:
			return float(value)
		except TypeError:
			return value
	elif (isinstance(default, int)):
		try:
			return int(value)
		except TypeError:
			return value
	else:
		print("cli_cvt_value: not manage value class=",value.__class__)
		return value


	
