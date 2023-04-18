#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot_eqdsk_data				import plot_eqdsk_data
	from routines.cli_routines			import *

	if(cli_present("-h",sys.argv)):
		print("\nThis script plot eqdsk file equilibrium parameters\n")
		print("plot_eqdsk_data options")
		print("\t-file           eqdsk filename [mag.eqdsk]")
		print("\t-in_contours    plot in wall countours")
		print("\t-reverse        If set reverse x axes")
		print("\t-contours       If set plot contours")
		print("\t-one_plot       One plot on each figure [d=false]")
		print("\t-save           Save figures on files, values=none/png/ps/eps/pdf [D='none']")
		print()
		exit()

	file			= cli_get_value("-file",			sys.argv, "mag.eqdsk")
	in_contours		= cli_present("-in_contours",	sys.argv)
	reverse			= cli_present("-reverse",			sys.argv)
	one_plot		= cli_present("-one_plot",		sys.argv)
	save			= cli_get_value("-save",			sys.argv, "none")
	plot_eqdsk_data(file, in_contours=in_contours, reverse=reverse, save=save, one_plot=one_plot)
	exit()

#=======================================

# Function definition is here

import numpy 							as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot 				as pyp
from matplotlib.backends.backend_pdf	import PdfPages
from files.eqdsk_routines				import load_eqdsk_file
from matplotlib.path		import Path
try:
	from matplotlib	import _cntr as cntr
except:
	from legacycontour  import _cntr as cntr

#=========================================================
# This routine plot magnetic fields as provided
#=========================================================
#

def plot_eqdsk_data(file="mag.eqdsk", in_contours=0, reverse=0, save="none", one_plot=0):

	print("plot_eqdsk_data")
	print("\treading magnetic file")

	eqdsk = load_eqdsk_file(file)

#	rz_c1, rz_c2, rz_c3, lengths = eqdsk_get_separatrix_contours(eqdsk)
	r_arr = np.arange(eqdsk.nw, dtype='f8')*eqdsk.rdim/(eqdsk.nw-1) + eqdsk.rleft
	z_arr = np.arange(eqdsk.nh, dtype='f8')*eqdsk.zdim/(eqdsk.nh-1) + (eqdsk.zmid - 0.5*eqdsk.zdim)
	eqdsk.psirz += eqdsk.sibry
	r_mat, z_mat = np.meshgrid(r_arr, z_arr, indexing='ij')

	Psi = np.arange(eqdsk.nw, dtype='f8')/(eqdsk.nw-1)*(eqdsk.sibry - eqdsk.simag) + eqdsk.simag

	print("\tplotting magnetic field on original mesh")

	levels = [40]
	linestyles = 1
	if(len(levels) > 1):
		sorted_levels = np.sort(levels)
		ind = np.where(np.diff(sorted_levels) > 0.); ind = ind[0]
		if(len(ind) > 0):
			sorted_levels = sorted_levels[np.insert(ind,0,0)]
		else:
			sorted_levels = int(sorted_levels[0])
	else:
		n_levels		= int(levels[0])
		min_flux		= np.min(eqdsk.psirz)
		max_flux		= np.max(eqdsk.psirz)
		n_levels_out	= int(n_levels*(min_flux-eqdsk.sibry)/(max_flux-min_flux))
		sorted_levels = (np.arange(n_levels,dtype='f8') + n_levels_out)/(n_levels)*(max_flux-min_flux)+eqdsk.sibry
	lines_tyles=['solid' , 'dashed' , 'dashdot' , 'dotted' ]

#	Prepare for pdf

	if(save == "pdf"):	pdf = PdfPages(os.path.splitext(os.path.split(MeshFile)[1])[0]+"_data.pdf")   				#pdf in one file only

	titles2d = ["$Flux\ [Wb/rad]$"]
	ylabel1d = ["$F=RB_t\ [T*m]$",  "$P\ [nt/m]$",		 "$FF'\ [(mT^2)/(Wb*rad)]$", "$P'\ [(nt/m^2)/(Wb/rad)]$", "$q$"]
	xlabel1d = ["$\Psi\ [Wb/rad]$", "$\Psi\ [Wb/rad]$", "$\Psi\ [Wb/rad]$",			"$\Psi\ [Wb/rad]$",				"$\Psi\ [Wb/rad]$"]	
	fig = []
	ax1d = []
	ax2d = []
	pl1d = []
	im2d = []
	if(one_plot != 1):
		fig.append(pyp.figure())
		for i in range(1):
			ax2d.append(fig[0].add_subplot(1,1,1))

		fig.append(pyp.figure())
		for i in range(len(xlabel1d)):
			ax1d.append(fig[1].add_subplot(3,2,i+1))
	else:
		fig = []
		for i in range(1):
			fig.append(pyp.figure())
			ax2d.append(fig[i].add_subplot(111))
		for i in range(1,len(xlabel1d)+1):
			fig.append(pyp.figure())
			ax1d.append(fig[i].add_subplot(111))

	im = []
	for i in range(len(ax2d)):
		ax2d[i].set_aspect(1.)
		ax2d[i].autoscale(enable=True, axis='both', tight=True)
		ax2d[i].set_xlabel("$R\ (m)$")
		ax2d[i].set_ylabel("$Z\ (m)$")
		if((one_plot == 1) or (i == 0)): ax2d[i].set_title(file)
		ax2d[i].autoscale(enable=True, axis='both', tight=True)


		im2d.append(ax2d[i].pcolormesh(r_mat,   z_mat, eqdsk.psirz.T,   shading='gouraud'))
		if(in_contours == 0):
			ax2d[i].contour(r_mat, z_mat, eqdsk.psirz.T, sorted_levels, linestyles=lines_tyles[linestyles], colors='k')
		else:
			WallPath = Path(np.array([eqdsk.rlim,eqdsk.zlim]).T, closed=True)									#define wall path

			cs = cntr.Cntr(r_mat, z_mat, eqdsk.psirz.T)															#define flux surfaces
			for level in sorted_levels:
				res = cs.trace(level)																			#get all contours for level
				nsegs = len(res) // 2
				segments, codes = res[:nsegs], res[nsegs:]
				for k in range(nsegs-1,-1,-1):
					x_count = segments[k][:,0]																	#define contour
					y_count = segments[k][:,1]

					IsIn = np.where(WallPath.contains_points(np.array([x_count,y_count]).T))[0]					#find countour points inside wall
					if(len(IsIn) > 0):
						ax2d[i].plot(x_count[IsIn], y_count[IsIn], linestyle=lines_tyles[linestyles], color='b')

		ax2d[i].plot(eqdsk.rlim,eqdsk.zlim,'k-')									# plot wall
		ax2d[i].plot(eqdsk.rbbbs,eqdsk.zbbbs,'r-')									# plot boundary
#		ax2d[i].plot(rz_c1[0,:], rz_c1[1,:], 'g-')
#		ax2d[i].plot(rz_c2[0,:], rz_c2[1,:], 'b-')
#		ax2d[i].plot(rz_c3[0,:], rz_c3[1,:], 'y-')

	for i in range(len(ax1d)):
		ax1d[i].autoscale(enable=True, axis='both', tight=True)
		ax1d[i].set_xlabel(xlabel1d[i])
		ax1d[i].set_ylabel(ylabel1d[i])
		if((one_plot == 1) or (i == 0)): ax1d[i].set_title(file)
		ax1d[i].grid()

		pl1d.append(ax1d[0].plot(Psi, eqdsk.fpol))
		pl1d.append(ax1d[1].plot(Psi, eqdsk.pres))
		pl1d.append(ax1d[2].plot(Psi, eqdsk.ffprim))
		pl1d.append(ax1d[3].plot(Psi, eqdsk.pprime))
		pl1d.append(ax1d[4].plot(Psi, eqdsk.qpsi))

	if(reverse == 1):												#Reverse x axis
		for i in range(len(ax1d)):
			xlim = ax1d[i].get_xlim()
			ax1d[i].set_xlim(xlim[1], xlim[0])

	cb = []
	if(one_plot != 1):
		for i in range(len(im2d)):
			cb.append(fig[0].colorbar(im2d[i],   ax=ax2d[i]))
	else:
		for i in range(len(im2d)):
			cb.append(fig[i].colorbar(im2d[i]))

	for i in range(len(cb)): cb[i].set_label(titles2d[i])

	if(save != "none"):
		for i in range(len(fig)):
			if(one_plot != 1): fig[i].set_size_inches(20.,15.)
			if(save == "pdf"):
				pdf.savefig(Fig[i])
			else:
				fig[i].savefig(os.path.splitext(os.path.split(MeshFile)[1])[0]+"_data_{:d}.".format(i+1)+save)

		pyp.show(block=False)
		pyp.close()
	else:
		pyp.show()

	print("plot_eqdsk_data: Completed")

	return

