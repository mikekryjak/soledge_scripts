#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot2d_eqdsk_fields	import plot2d_eqdsk_fields
	from routines.cli_routines	import *

	if(cli_present("-h",sys.argv)):
		print("\nThis script plot fileds and flux from an eqdsk file\n")
		print("plot2d_eqdsk_fields options")
		print("\t-file   <name>    eqdsk filename [mag.eqdsk]")
		print("\t-nr     <value>   if > than in eqdsk grid sinze in r [0]")
		print("\t-nz     <value>   if > than in eqdsk grid sinze in z [0]")
		print("\t-smooth <value>   smooth factor (>= 0.5) used if nr or nz > grid size [1]")
		print("\t-one_plot  One plot on each figure [d=false]")
		print("\t-save     Save figures on files, values=none/png/ps/eps/pdf [D='none']")
		print()
		exit()

	file		= cli_get_value("-file",		sys.argv, "mag.eqdsk")
	nr			= cli_get_value("-nr",			sys.argv, 0)
	nz			= cli_get_value("-nz",			sys.argv, 0)
	smooth		= cli_get_value("-smooth",		sys.argv, 0.)
	one_plot	= cli_present("-one_plot",	sys.argv)
	save		= cli_get_value("-save",		sys.argv, "none")
	plot2d_eqdsk_fields(file, nr=nr, nz=nz, smooth=smooth, save=save, one_plot=one_plot)
	exit()

#=======================================

# Function definition is here

import numpy 							as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot 				as pyp
from matplotlib.backends.backend_pdf	import PdfPages
from files.eqdsk_routines				import load_eqdsk_file, eqdsk_compute_fields

#=========================================================
# This routine plot magnetic fields as provided
#=========================================================
#

def plot2d_eqdsk_fields(file="mag.eqdsk", nr=0, nz=0, smooth=0., save="none", one_plot=0):

	print("eqdsk_plot_fields")
	print("\treading magnetic file")

	eqdsk = load_eqdsk_file(file)

	print("\tcomputing magnetic fields")
	rz_b_psi  = eqdsk_compute_fields(eqdsk, nr=nr, nz=nz, smooth=smooth)

	print("\tplotting magnetic field on original mesh")

#	Prepare for pdf

	if(save == "pdf"):	pdf = PdfPages(os.path.splitext(os.path.split(MeshFile)[1])[0]+"_fields.pdf")   				#pdf in one file only


	ax	   = []
	if(one_plot != 1):
		fig    = pyp.figure()
		for i in range(1,5): ax.append(fig.add_subplot(2,2,i))
	else:
		fig = []
		for i in range(1,5):
			fig.append(pyp.figure())
			ax.append(fig[i-1].add_subplot(111))

	im = []
	for i in range(len(ax)):
		ax[i].set_aspect(1.)
		ax[i].set_xlabel("$R\ (m)$")
		ax[i].set_ylabel("$Z\ (m)$")
		ax[i].set_title(file)
		ax[i].autoscale(enable=True, axis='both', tight=True)
#		im.append(ax[i].pcolormesh(rz_b_psi[:,:,0], rz_b_psi[:,:,1], rz_b_psi[:,:,2+i], shading='gouraud'))
		im.append(ax[i].pcolormesh(rz_b_psi[:,:,0], rz_b_psi[:,:,1], rz_b_psi[:,:,2+i], shading='flat'))
		ax[i].contour(rz_b_psi[:,:,0], rz_b_psi[:,:,1], rz_b_psi[:,:,2+i], linestyles='dashed', colors='k')
		ax[i].contour(rz_b_psi[:,:,0], rz_b_psi[:,:,1], rz_b_psi[:,:,2+i], linestyles='solid', colors='b', levels=[eqdsk.sibry])
		ax[i].plot(eqdsk.rlim,eqdsk.zlim,'k-')													# plot wall
		ax[i].plot(eqdsk.rbbbs,eqdsk.zbbbs,'r-')												# plot boundary

	cb = []
	if(one_plot != 1):
		for i in range(len(im)):
			cb.append(fig.colorbar(im[i],   ax=ax[i]))
	else:
		for i in range(len(im)):
			cb.append(fig[i].colorbar(im[i]))

	labels = ["$B_R\ (T)$", "$B_Z\ (T)$", "$B_Y\ (T)$", "$Flux\ (Wb/rad)$"]
	for i in range(len(cb)): cb[i].set_label(labels[i])

	if(save != "none"):
		for i in range(len(fig)):
			if(one_plot != 1): fig[i].set_size_inches(20.,15.)
			if(save == "pdf"):
				pdf.savefig(Fig[i])
			else:
				fig[i].savefig(os.path.splitext(os.path.split(MeshFile)[1])[0]+"_fields_{:d}.".format(i+1)+save)

		pyp.show(block=False)
		pyp.close()
	else:
		pyp.show()


	pyp.show()

	print("plot_mag_fields: Completed")

