#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from plot_rates import plot_rates
	from routines.cli_routines import *

	if(cli_present("-h",sys.argv)):
		print("\nThis script plot plasma parameters on soledge mesh\n")
		print("plot_rates options")
		print("\t-specie <value>   Impurity specie (Li, Be, C, N, Ne, Ar) [d='C']")
		print("\t-rate <value>     Rate type (0=ioniz., 1=recomb.,2=line excit,3= line recomb [d=0]")
		print("\t-ne <values>      electron density values [d=[]]")
		print("\t-te <values>      electron temperature values [d=[]]")
		print("\t-log_scale        Use log scale for z colors [d=false]")
		print("\t-save             Save figures on files, values=none/png/ps/eps/pdf [D='none']")
		print("\t-no_shade         no shade 2d plot [d=false]")
		print("\t-one_plot         All figure on one plot [d=false]")
		print()
		exit()

	specie	 	= cli_get_value("-specie",	sys.argv,  "C")
	rate	 	= cli_get_value("-rate",	sys.argv, 0)
	ne	= cli_get_value("-ne",	sys.argv,	[])
	te		= cli_get_value("-te",		sys.argv, [])
	log_scale	= cli_present("-log_scale",		sys.argv)
	no_shade	= cli_present("-no_shade",	sys.argv)
	save		= cli_get_value("-save",		sys.argv, "none")
	one_plot	= cli_present("-one_plot",	sys.argv)
	plot_rates(specie=specie, rate=rate, ne=ne, te=te, log_scale=log_scale, one_plot=one_plot, no_shade=no_shade, save=save)
	exit()

#=======================================

# Function definition is here

from math import log, log10, floor
import os
import numpy							as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot 				as pyp
from matplotlib.colors 					import LogNorm
from matplotlib.backends.backend_pdf	import PdfPages
from files.load_amdata_file		import load_amdata_file

#=========================================================
# 
#=========================================================

def plot_rates(specie="C", rate=0, ne=[], te=[], log_scale=0, one_plot=0, no_shade=0, save="none"):

	print("plot_rates")

	shading = 'gouraud'
	if(no_shade == 1): shading = 'flat'

#	Read references data
	data_types=["ionization","recombination","line_excitation","line_recombination"]

	amdata = load_amdata_file(specie, data_type=data_types[rate])

#	Preparing for plot

	if(save == "pdf"):	pdf = PdfPages("{:s}_{:s}_rates.{:s}".format(specie,data_types[rate],save))   #pdf in one file only
	else:				i_plot_file = 0

	nPLots = amdata.izmax-amdata.izmin+1

	nRows		= 2
	nCols		= 3
	PlotPerFig	= nRows*nCols
	nFigs		= int(nPLots/PlotPerFig)
	if(nFigs*PlotPerFig < nPLots): nFigs += 1

	Fig = []
	Ax  = []
	Im  = []
	iFig = []
	if(one_plot != 1):
		for iF in range(nFigs):	
			Fig.append(pyp.figure())
			for i in range(min(PlotPerFig,nPLots-iF*PlotPerFig)): 
				Ax.append(Fig[-1].add_subplot(nRows,nCols,i+1))
				Ax[-1].locator_params(axis='x',nbins=4)
				Im.append(0)
				iFig.append(Fig[-1])

			Fig[-1].tight_layout(pad=2., w_pad=3., h_pad=3.)
	else:
		for i in range(nPLots):
			Fig.append(pyp.figure())
			Ax.append(Fig[i].add_subplot(111))
			Im.append(0)

	for figure in Fig:  figure.patch.set_facecolor('white')


	for i in range(nPLots):
		Ax[i].autoscale(enable=True, axis='both', tight=True)
		Ax[i].set_xlabel("$n_e\ (m^{-3})$")
		Ax[i].set_ylabel("$T_e\ (eV)$")
		Ax[i].set_xscale("log")
		Ax[i].set_yscale("log")
		if(amdata.izmin+i > 1): Ax[i].set_title("{:s} {:s} {:d}+".format(data_types[rate], specie, amdata.izmin+i-1))
		else:								Ax[i].set_title("{:s} {:s}".format(data_types[rate], specie))
		
#	Plot parameters

		Ne2D, Te2D, rate_coefficients = compute_rate_coefficients(amdata,i)
		if(log_scale == 0):
#			Im[i] = Ax[i].pcolormesh(logNe2D, logTe2D, rate_coefficient, shading=shading,  vmin = Vmin,  vmax = Vmax)
			Im[i] = Ax[i].pcolormesh(Ne2D, Te2D, rate_coefficients, shading=shading)
		else:
			rate_coefficients = np.where(rate_coefficients < 1., 1., rate_coefficients)
			Im[i] = Ax[i].pcolormesh(Ne2D, Te2D, rate_coefficients, shading=shading,norm=LogNorm(vmin=rate_coefficients.min(), vmax=rate_coefficients.max()))

		cb = iFig[i].colorbar(Im[i], ax=Ax[i])
		cb.set_label("$rate (m^3/s)$")

	if(save != "none"):
		for iF in range(len(Labels)):
			if(one_plot != 1): FigPl[iF][0].set_size_inches(10.05,7.44)
			if(save == "pdf"):
				pdf.savefig(Fig[i])
			else:
				i_plot_file += 1
				Fig[i].savefig("{:s}_{:s}_rates.{:s}".format(specie,data_types[rate],save))
		pyp.show(block=False)
		pyp.close()
	else:
		pyp.show()


	if(save == "pdf"):	pdf.close()

	print("plot_rates: Completed")

	return


def compute_rate_coefficients(amdata, iz, nNe=0, nTe=0):
	if(nNe == 0): nNe = amdata.nNe
	if(nTe == 0): nTe = amdata.iz_data[iz].nTe

	logNe1D = np.linspace(log10(amdata.Ne_min), log10(amdata.Ne_max), nNe, dtype='f8')
	logTe1D = np.linspace(log10(amdata.iz_data[iz].Te_min), log10(amdata.iz_data[iz].Te_max), nTe, dtype='f8')

	logNe2D, logTe2D		= np.meshgrid(logNe1D, logTe1D)

	pow_logNe = np.empty((amdata.degree, nTe, nNe), dtype='f8')
	pow_logTe = np.empty((amdata.degree, nTe, nNe), dtype='f8')
	pow_logNe[0,:,:] = 1.
	pow_logTe[0,:,:] = 1.

	for k in range (1,amdata.degree):
		pow_logNe[k,:,:] =pow_logNe[k-1,:,:]*logNe2D
		pow_logTe[k,:,:] = pow_logTe[k-1,:,:]*logTe2D

	rate_coefficients = np.zeros((nTe, nNe), dtype='f8')
	ind=0
	for m in range(amdata.degree):
		for l in range (amdata.degree-m):
			rate_coefficients = rate_coefficients+amdata.iz_data[iz].coefficients[ind]*pow_logNe[m,:,:]*pow_logTe[l,:,:]
			ind += 1
	Ne2D = 10**logNe2D
	Te2D = 10**logTe2D
	rate_coefficients = np.exp(log(10.)*rate_coefficients)*Ne2D

	return Ne2D, Te2D, rate_coefficients

