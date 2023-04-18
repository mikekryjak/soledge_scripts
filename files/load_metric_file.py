import numpy as np
import types
import os
import h5py
from routines.h5_routines	 import h5_read
from routines.globals		 import DEBUG

#=========================================================
# This routine read metric file
#=========================================================

def load_metric_file(Path="", nZones=0): 

	if((len(Path) > 0) and (Path[-1] != "/")): Path = Path + "/"

	if(nZones == 0):
		if_mesh = h5py.File(Path+"../mesh.h5", "r")
		nZones = h5_read(if_mesh,"NZones", keep_array=False)
		if_mesh.close()

#	Metric file
	
	if(DEBUG > 0): print("load_metric_file: Loading: ",Path+"metric")
	try:
		metric = []
		if_metric = h5py.File(Path+"metric", "r")
	except:
		if(DEBUG > 0): print("Metric file not available")
		return metric

	for k in range(nZones):

		metric.append(types.SimpleNamespace())
		zone = "zone{:d}".format(k+1)
			
		metric[k].cpp	  		= h5_read(if_metric, zone+ '/cpp', order = 'F')					#[Nx+2,Nz+2]
		metric[k].cpt	  		= h5_read(if_metric, zone+ '/cpt', order = 'F')					#[Nx+2,Nz+2]
		metric[k].ctt	  		= h5_read(if_metric, zone+ '/ctt', order = 'F')					#[Nx+2,Nz+2]

		metric[k].c_pp	  		= h5_read(if_metric, zone+ '/c_pp', order = 'F')					#[Nx+2,Nz+2]
		metric[k].c_pt	  		= h5_read(if_metric, zone+ '/c_pp', order = 'F')					#[Nx+2,Nz+2]
		metric[k].c_tt	  		= h5_read(if_metric, zone+ '/c_tt', order = 'F')					#[Nx+2,Nz+2]

		metric[k].Jac	  		= h5_read(if_metric, zone+ '/Jac', order = 'F')					#[Nx+2,Nz+2]
		metric[k].G	  			= h5_read(if_metric, zone+ '/G', order = 'F')						#[Nx+2,Nz+2]

		metric[k].dPdR	  		= h5_read(if_metric, zone+ '/dPdR', order = 'F')					#[Nx+2,Nz+2]
		metric[k].dPdZ	  		= h5_read(if_metric, zone+ '/dPdZ', order = 'F')					#[Nx+2,Nz+2]
		metric[k].dTdR	  		= h5_read(if_metric, zone+ '/dTdR', order = 'F')					#[Nx+2,Nz+2]
		metric[k].dTdZ	  		= h5_read(if_metric, zone+ '/dTdZ', order = 'F')					#[Nx+2,Nz+2]
		metric[k].dRdP	  		= h5_read(if_metric, zone+ '/dRdP', order = 'F')					#[Nx+2,Nz+2]
		metric[k].dZdP	  		= h5_read(if_metric, zone+ '/dZdP', order = 'F')					#[Nx+2,Nz+2]
		metric[k].dRdT	  		= h5_read(if_metric, zone+ '/dRdT', order = 'F')					#[Nx+2,Nz+2]
		metric[k].dZdT	  		= h5_read(if_metric, zone+ '/dZdT', order = 'F')					#[Nx+2,Nz+2]

		metric[k].dS_north	  	= h5_read(if_metric, zone+ '/dS_north', order = 'F')				#[Nx,Nz]
		metric[k].dS_south	  	= h5_read(if_metric, zone+ '/dS_south', order = 'F')				#[Nx,Nz]
		metric[k].dS_east	  	= h5_read(if_metric, zone+ '/dS_east', order = 'F')				#[Nx,Nz]
		metric[k].dS_west	  	= h5_read(if_metric, zone+ '/dS_west', order = 'F')				#[Nx,Nz]
		metric[k].dvol		  	= h5_read(if_metric, zone+ '/dvol', order = 'F')					#[Nx,Nz]

		metric[k].sinepitch_east  	= h5_read(if_metric, zone+ '/sinepitch_east', order = 'F')		#[Nx,Nz]
		metric[k].sinepitch_west  	= h5_read(if_metric, zone+ '/sinepitch_west', order = 'F')		#[Nx,Nz]

	if_metric.close()

	return metric
