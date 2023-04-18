#!/usr/bin/python

# Function definition is here

import csv
import numpy as np

#====================================================================================================
# Read a csv file
#====================================================================================================

def cvs_read(file, DataType='f8'):
	f=open(file)
	fdata=csv.reader(f,quotechar='"',quoting=csv.QUOTE_NONNUMERIC)
	for row in fdata:
		next(fdata)

	n_lines = fdata.line_num-1
	f.closed

	f=open(file)
	fdata=csv.reader(f,quotechar='"',quoting=csv.QUOTE_NONNUMERIC)
	titles = next(fdata)

	data = np.zeros((n_lines, len(titles)), dtype=DataType)
	for i in range(n_lines):
		data[i,:] = next(fdata)

	f.closed

	return titles, data
