#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot2d_dt import plot2d_dt
	from routines.cli_routines import cli_present,cli_get_value

	if(cli_present("-h",sys.argv)):
		print("\nThis script plot soledge dt on mesh\n")
		print("plot2d_dt options")
		print("\t-path        Directory with simulation [d='./']")
		print("\t-evolution   If > 0 plot evolution file [d=0]")
		print("\t-d_only      plot only e- and D+ [d=false]")
		print("\t-dt_range    Range dt plot scale (0.=auto) [D=[[0.,0.]]")
		print("\t-marker_val  Marker factor of dt_min [D=1.1]")
		print("\t-log_scale   Use log scale for z colors [d=false]")
		print("\t-extra_walls Flag to show extra walls [d=false]")
		print("\t-no_shade    no shade 2d plot [d=false]")
		print("\t-one_plot    One plot on each figure [d=false]")
		print("\t-save        Save figures on files, values=none/png/ps/eps/pdf [D='none']")
		print()
		exit()

	path		= cli_get_value("\t-path",			sys.argv,"./")
	evolution	= cli_get_value("\t-evolution",		sys.argv,	0)
	d_only		= cli_present("-d_only",			sys.argv)
	dt_range	= cli_get_value("\t-dt_range",		sys.argv,	[0.,0.])
	marker_val	= cli_get_value("\t-marker_val",	sys.argv,	1.1)
	log_scale	= cli_present("\t-log_scale",		sys.argv)
	one_plot	= cli_present("\t-one_plot",		sys.argv)
	extra_walls	= cli_present("-extra_walls",	sys.argv)
	no_shade	= cli_present("-no_shade",		sys.argv)
	save		= cli_get_value("\t-save",			sys.argv, "none")
	plot2d_dt(path=path, evolution=evolution, d_only=d_only, dt_range=dt_range, log_scale=log_scale, one_plot=one_plot, save=save)
	exit()

#=======================================

# Function definition is here

import h5py
import os
import numpy							as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot 				as pyp
from matplotlib.colors 					import LogNorm
from matplotlib.backends.backend_pdf	import PdfPages
from routines.h5_routines 				import h5_read
from routines.compute_dt 				import compute_dt
from routines.utils_walls				import plot2d_walls

from files.load_refpar_file				import load_refpar_file
from files.load_plasma_files			import load_plasma_files
from files.load_soledge_mesh_file		import load_soledge_mesh_file
from files.load_metric_file		import load_metric_file
from files.load_input_file				import load_input_file

#=========================================================
# This routine to plot M, ne/ni and Te/Ti results
#=========================================================

def plot2d_dt(path="./", evolution=0, d_only=0, dt_range=[0.,0.], marker_val=1.1, log_scale=0, one_plot=0, extra_walls=0, no_shade=0, save="none"):

	print("plot2d_dt")

	shading = 'gouraud'
	if(no_shade == 1): shading = 'flat'

#	Preparing for data reading

	if((len(path) > 0) and (path[-1] != "/")): path = path + "/"
	input = load_input_file(path)

#	Read simulation time
	Config	= load_soledge_mesh_file(path+"mesh.h5")
	Zones	= Config.Zones
	nZones  = len(Zones)
	
	RefPar	= load_refpar_file(path+ "Results/")

	Plasmas = load_plasma_files(path, Evolution=evolution, DeNorm=False)

	if(d_only != 0): Plasmas = Plasmas[:2]

	Metric = load_metric_file(path+"Results/", nZones=nZones)

# Load or compute transport

	nPlasmas = len(Plasmas)
	Dn_p = []
	for iIon in range(nPlasmas-1):
		if(iIon == 0): iSpec = 0
		else:					iSpec=1
		Dn_p.append([])
		for iZone in range(len(Config.Zones)):
			if((input.ballooning_parameters.ballooning_model == 2) and (Config.transp_values_OK)):
				Dn_p[-1].append(Config.Zones[iZone].Ballooning[1]*input.transport_parameters.Dn_p[iSpec])
			else:
				Dn_p[-1].append(input.transport_parameters.Dn_p[iSpec])

	CFL		= input.global_parameters.CFL
	Dt		= compute_dt(Config, Plasmas, Metric, RefPar, Dn_p, CFL=CFL, Dn_p=Dn_p)
	
#	Preparing for plot

	if(save == "pdf"):	pdf = PdfPages("plot2d_dt_t={:.3f}.".format(RefPar.time)+save)   #pdf in one file only

	Labels = ["$dt_{adv}\ (ns)$", "$dt_{diff}\ (ns)$"]
	Fig = []
	Ax  = []
	Im  = []

	for iIon in range(len(Plasmas)-1):
		if(one_plot != 0):
			for i in range(2):
				Fig.append(pyp.figure())
				Fig[-1].patch.set_facecolor('white')
				Ax.append(Fig[-1].add_subplot(1,1,1))
#				Fig[-1].tight_layout()
		else:
			Fig.append(pyp.figure(figsize=(12, 8)))
			Fig[-1].patch.set_facecolor('white')
			for i in range(2):
				Ax.append(Fig[-1].add_subplot(1,2,i+1))
#				Fig[-1].tight_layout()

	for i in range(len(Ax)):
		Ax[i].set_aspect(1.)
		Ax[i].autoscale(enable=True, axis='both', tight=True)
		Ax[i].set_xlabel("$R\ (m)$")
		Ax[i].set_ylabel("$Z\ (m)$")
		Ax[i].tick_params(axis='both', which='major', labelsize=8)
		Ax[i].autoscale(axis='both', tight=True)
#		Ax[i].set_axis_bgcolor("#bdb76b")

	for iIon in range(nPlasmas-1):
		iPlasma = iIon + 1

		DtAdvMin  = np.empty(nZones, dtype='f8')
		DtAdvMax  = np.empty(nZones, dtype='f8')
		DtDiffMin = np.empty(nZones, dtype='f8')
		DtDiffMax = np.empty(nZones, dtype='f8')
		len_in    = np.zeros(nZones, dtype='i4')
		i = 0
		for k in range(nZones):
			Dt[iIon][k].Adv  *= 1e9
			Dt[iIon][k].Diff *= 1e9
			ii, jj  = np.where(Zones[k].Chi == 0)
			if(len(ii) > 0):
				DtAdvMin[i]  = np.min(Dt[iIon][k].Adv[ii,jj])
				DtAdvMax[i]  = np.max(Dt[iIon][k].Adv[ii,jj])
				DtDiffMin[i] = np.min(Dt[iIon][k].Diff[ii,jj])
				DtDiffMax[i] = np.max(Dt[iIon][k].Diff[ii,jj])
				len_in[k]	 = len(ii)
				i			+= 1

		DtAdvMin  = np.min(DtAdvMin[:i])
		DtAdvMax  = np.max(DtAdvMax[:i])
		DtDiffMin = np.min(DtDiffMin[:i])
		DtDiffMax = np.max(DtDiffMax[:i])
	
		if(dt_range[0] > 0.):
			DtAdvMin  =  dt_range[0]
			DtDiffMin =  dt_range[0]

		if(dt_range[1] > 0.):
			DtAdvMax  =  dt_range[1]
			DtDiffMax =  dt_range[1]
		elif(dt_range[1] < -1.):
			DtAdvMax  =  -dt_range[1]*DtAdvMin
			DtDiffMax =  -dt_range[1]*DtDiffMin
		else:
			DtAdvMax  =  5.*DtAdvMin
			DtDiffMax =  5.*DtDiffMin

		iAxAdv = iIon*2
		Ax[iAxAdv].set_title("$"+Plasmas[iPlasma][0].ion + "\ Advection:\ dt_{min}="+"{:.3f}\ (ns)\ @\ CFL={:.2f}$".format(DtAdvMin,CFL))
		iAxDiff = iIon*2+1
		Ax[iAxDiff].set_title("$"+Plasmas[iPlasma][0].ion + "\ Diffusion:\ dt_{min}="+"{:.3f}\ (ns)\ @\ CFL={:.2f}$".format(DtDiffMin,CFL))

		Im.extend([0,0])
		for k in range(nZones):
			if(len_in[k] > 0):
				DtAdv   = np.ma.masked_where(Zones[k].Chi == 1., Dt[iIon][k].Adv)
				DtDiff  = np.ma.masked_where(Zones[k].Chi == 1., Dt[iIon][k].Diff)
		
				DtAdv = np.append(DtAdv, DtAdv[:,-1].reshape(DtAdv.shape[0],1), axis=1)
				DtAdv = np.append(DtAdv, DtAdv[-1,:].reshape(1,DtAdv.shape[1]), axis=0)
		
				DtDiff = np.append(DtDiff, DtDiff[:,-1].reshape(DtDiff.shape[0],1), axis=1)
				DtDiff = np.append(DtDiff, DtDiff[-1,:].reshape(1,DtDiff.shape[1]), axis=0)
				if(log_scale == 0):
					Im[iAxAdv] = Ax[iAxAdv].pcolormesh(Zones[k].gridR, Zones[k].gridZ, DtAdv,  vmin = DtAdvMin,  vmax = DtAdvMax, cmap=pyp.get_cmap('hot'), shading=shading)
					Im[iAxDiff] = Ax[iAxDiff].pcolormesh(Zones[k].gridR, Zones[k].gridZ, DtDiff, vmin = DtDiffMin,  vmax = DtDiffMax, cmap=pyp.get_cmap('hot'), shading=shading)
				else:
					Im[iAxAdv] = Ax[iAxAdv].pcolormesh(Zones[k].gridR, Zones[k].gridZ, DtAdv,  norm=LogNorm(vmin = DtAdvMin,  vmax = DtAdvMax), cmap=pyp.get_cmap('hot'), shading=shading)
					Im[iAxDiff] = Ax[iAxDiff].pcolormesh(Zones[k].gridR, Zones[k].gridZ, DtDiff, norm=LogNorm(vmin = DtDiffMin,  vmax = DtDiffMax), cmap=pyp.get_cmap('hot'), shading=shading)

				if(marker_val > 1.):
					ii,jj = np.where((Dt[iIon][k].Adv < DtAdvMin*marker_val) & (Zones[k].Chi == 0.))
					if(len(ii) > 0):
						Ax[iAxAdv].plot(Zones[k].gridRc[ii,jj], Zones[k].gridZc[ii,jj], 'ro',markersize=10)

					ii,jj = np.where((Dt[iIon][k].Diff < DtDiffMin*marker_val) & (Zones[k].Chi == 0.))
					if(len(ii) > 0):
						Ax[iAxDiff].plot(Zones[k].gridRc[ii,jj], Zones[k].gridZc[ii,jj], 'ro',markersize=10)

		for i in range(2): plot2d_walls(Ax[iIon*2+i], Config.Walls, extra_wall=extra_walls)

		if(one_plot != 0):
			for i in range(len(Labels)):
				cb = Fig[iIon*2+i].colorbar(Im[iIon*2+i], ax=Ax[iIon*2+i])
				cb.set_label(Labels[i])
		else:
			for i in range(len(Labels)):
				cb = Fig[iIon].colorbar(Im[iIon*2+i], ax=Ax[iIon*2+i])
				cb.set_label(Labels[i])

	if(save != "none"):
		for i in range(len(Fig)):
			if(one_plot != 1): Fig[i].set_size_inches(20.,15.)
			if(save == "pdf"):
				pdf.savefig(Fig[i])
			else:
				Fig[i].savefig("plot2d_dt_{:d}.".format(i+1)+save)
		pyp.show(block=False)
		pyp.close()
	else:
		pyp.show()

	if(save == "pdf"):	pdf.close()

	print("plot2d_dt: Completed")

	return
