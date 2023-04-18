#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	from neigh_to_h5			import neigh_to_h5

	neigh_to_h5()
	exit()

import numpy as np
import h5py

def neigh_to_h5():

	f2 = open('soledge2D.neighbors','r')
	ntri=int(f2.readline())
	neigh_tri=np.zeros((ntri,3))
	neigh_face=np.zeros((ntri,3))
	cnt=0
	for line in f2:
		line = line.strip()
		columns = line.split()
		neigh_tri[cnt,0]=int(columns[1])
		neigh_face[cnt,0]=int(columns[2])
		neigh_tri[cnt,1]=int(columns[4])
		neigh_face[cnt,1]=int(columns[5])
		neigh_tri[cnt,2]=int(columns[7])
		neigh_face[cnt,2]=int(columns[8])
		cnt=cnt+1
	f2.close()

	f=h5py.File('triangles.h5','a')

	print("writing triangles.h5/triangles/neigh_tri")
	f.create_dataset('/triangles/neigh_tri',data=np.transpose(neigh_tri))
	print("writing triangles.h5/triangles/neigh_face")
	f.create_dataset('/triangles/neigh_face',data=np.transpose(neigh_face))
	f.close()
