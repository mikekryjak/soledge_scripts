#!/usr/bin/env python

# ----- Script section -----

if __name__ == "__main__":
	import sys
	from cvs_cvs_plot2d_on_tri		import cvs_cvs_plot2d_on_tri
	from routines.cli_routines	import *

	if(cli_present("-h",sys.argv)):
		print("\nThis script plot plasma parameters on eirene mesh\n")
		print("cvs_plot2d_on_tri options")
		print("\t-tri_knots  <value>  file with knots indexes [d='TriKnots.txt']")
		print("\t-knots_rz   <value>  file with knots (R,Z) coordinate [d='']")
		print("\t-knots_r    <value>  file with knots (R) coordinates [d='Knots_R.txt']")
		print("\t-knots_z    <value>  file with knots (Z) coordinates [d='Knots_Z.txt']")
		print("\t-v_file     <value>  file parameter on triangles [d='']")
		print("\t-plot_name  <value>  title on plot [d='']")
		print("\t-r_label    <value>  label on R axis [d='$R\ (m)$']")
		print("\t-z_label    <value>  label on Z axis [d='$Z\ (m)$']")
		print("\t-v_label    <value>  value name [d='']")
		print("\t-skip_first          skip first line of files [d=false]")
		print("\t-rmin      <value>   R min on plot [D=]")
		print("\t-rmax      <value>   R max on plot [D=]")
		print("\t-zmin      <value>   Z min on plot [D=]")
		print("\t-zmax      <value>   Z max on plot [D=]")
		print("\t-vmin      <value>   Value min on plot [D=]")
		print("\t-vmax      <value>   Value max on plot [D=]")
		print("\t-log_scale           Use log scale for z colors [d=false]")
		print("\t-save      <value>   Save figures on files, values=none/stat/png/ps/eps/pdf [D='none']")
		print()
		exit()

	tri_knots	= cli_get_value("-tri_knots",		sys.argv,"TriKnots.txt")
	knots_rz	= cli_get_value("-knots_rz",			sys.argv,"")
	knots_r	 	= cli_get_value("-knots_r",			sys.argv,"Knots_R.txt")
	knots_z	 	= cli_get_value("-knots_z",			sys.argv,"Knots_Z.txt")
	v_file		= cli_get_value("-v_file",			sys.argv,"")
	plot_name	= cli_get_value("-plot_name",		sys.argv,"")
	v_label		= cli_get_value("-v_label",			sys.argv,"")
	r_label		= cli_get_value("-r_label",			sys.argv,"$R\ (m)$'")
	z_label	 	= cli_get_value("-z_label",			sys.argv,"$Z\ (m)$'")
	skip_first	= cli_present("-skip_first",		sys.argv)

	rmin	 	= cli_get_value("-rmin",			sys.argv,-1e.100)
	rmax	 	= cli_get_value("-rmax",			sys.argv, 1.e100)
	zmin	 	= cli_get_value("-zmin",			sys.argv,-1.e100)
	zmax	 	= cli_get_value("-vmax",			sys.argv, 1.e100)
	vmin	 	= cli_get_value("-vmin",			sys.argv,-1.e100)
	zmax	 	= cli_get_value("-vmax",			sys.argv, 1.e100)
	save		= cli_get_value("-save",			sys.argv, "none")
	log_scale	= cli_present("-log_scale",		sys.argv)
	cvs_plot2d_on_tri(tri_knots, knots_rz, knots_r, knots_z, v_file, plot_name=plot_name, v_label=v_label, r_label=r_label, z_label=z_label, \
						skip_first=skip_first, log_scale=log_scale, save=save, vmin=vmin, vmax=vmax, rmin=rmin, rmax=rmax, zmin=zmin, zmax=zmax)
	exit()

#=======================================

# Function definition is here

import h5py
import os
import numpy as np
import matplotlib.pyplot				as pyp
import matplotlib.tri					as tri
from matplotlib.backends.backend_pdf	import PdfPages
from matplotlib.colors 					import LogNorm
from matplotlib.path 					import Path
from routines.utils_walls				import plot2d_walls
from routines.h5_routines				import h5_read
from mesh.get_rz_core_sep				import get_rz_core_sep
from files.load_refpar_file				import load_refpar_file
from files.load_plasma_files			import load_plasma_files
from files.load_soledge_mesh_file		import load_soledge_mesh_file
from files.load_ions_list				import load_ions_list
from files.save_stat					import save_stat

#==============================================================================
# This routine plots parameters on eirene mesh
#==============================================================================

	cvs_plot2d_on_tri(tri_knots, knots_rz, knots_r, knots_z, v_file, plot_name="", v_label="", r_label="", z_label="", \
						skip_first=0, log_scale=0, save="none", vmin=-1e100, vmax=1e100, rmin=-1e100, rmax=1e100, zmin=-1e100, zmax=1e100):

	print("cvs_plot2d_on_tri")

#	load mesh

	if(skip_first): skiprows=1
	else:			skiprows=0

	TriKnots = np.loadtxt(tri_knots, delimiter=",", skiprows=skiprows, delimiter=",")
	if(len(knots_rz) > 0):
		rz = np.loadtxt(knots_rz, delimiter=",", skiprows=skiprows, delimiter=",")
		r = rz[:,0]
		z = rz[:,1]
		rz=none
	else:
		r=np.loadtxt(knots_r, delimiter=",", skiprows=skiprows, delimiter=",")
		z=np.loadtxt(knots_z, delimiter=",", skiprows=skiprows, delimiter=",")

#	load values
	values = np.loadtxt(v_file, delimiter=",", skiprows=skiprows, delimiter=",")
	
	Fig = pyp.figure()
	Fig.patch.set_facecolor('white')
	Ax  = Fig.add_subplot(111)

	Ax.set_aspect(1.)
	Ax.autoscale(enable=True, axis='both', tight=True)
	Ax.set_xlabel("$R\ (m)$")
	AxPl[iF][i].set_ylabel("$Z\ (m)$")

	Im  = []
					Fig.append(pyp.figure())
					Fig[-1].patch.set_facecolor('white')
					for i in range(len(Labels[iF])): 
						Ax.append(Fig[-1].add_subplot(NrowNcol[iF][1],NrowNcol[iF][0],i+1))
						AxPl[-1].append(Ax[-1])
						FigNum[-1].append(len(Fig))
						PosNum[-1].append(i+1)
						FigPl[-1].append(Fig[-1])

					Fig[-1].tight_layout(pad=2., w_pad=3., h_pad=3.)
			else:
				for iF in range(len(Labels)):
					AxPl.append([])
					FigPl.append([])
					FigNum.append([])
					PosNum.append([])
					for i in range(len(Labels[iF])):
						Fig.append(pyp.figure())
						Fig[-1].patch.set_facecolor('white')
						Ax.append(Fig[-1].add_subplot(111))

						AxPl[-1].append(Ax[-1])
						FigNum[-1].append(len(Fig))
						PosNum[-1].append(1)
						FigPl[-1].append(Fig[-1])
		
			for iF in range(len(Labels)):
				for i in range(len(Labels[iF])):
					AxPl[iF][i].set_aspect(1.)
					AxPl[iF][i].autoscale(enable=True, axis='both', tight=True)
					AxPl[iF][i].set_xlabel("$R\ (m)$")
					AxPl[iF][i].set_ylabel("$Z\ (m)$")
	#				if(i == 0): AxPl[iF][i].set_title(path_label+" @ t={:.3f} s".format(tempus))
					if(i == 0): AxPl[iF][i].set_title(path_label)
					if(diff == 0):
						if(evolution[0] == 0):
							if(i == 1): AxPl[iF][i].set_title(Plasmas[iPlasma][0].ion)
						else:
							if(i == 1): AxPl[iF][i].set_title("Evol.={:d}".format(evolution[0]))
							if(i == 2): AxPl[iF][i].set_title(Plasmas[iPlasma][0].ion)
					else:
						if(i == 1): AxPl[iF][i].set_title("Ev_{:d}-Ev_{:d}".format(evolution[1],evolution[0]))
						if(i == 2): AxPl[iF][i].set_title(Plasmas[iPlasma][0].ion)


	#	Plot parameters

		for iF in range(len(Labels)):
			for i in range(len(Labels[iF])):
#				try:
					iPar = Plasmas[iPlasma][0].Triangles.VNames.index(Pars[iF][i])	

#					=========== Print section =============

#					Print total and average on whole plasma

					if(diff == 1):	TValues =  (Plasmas[iPlasma][0].Triangles.Values[iPar][TriKnots[:,0]] + Plasmas[iPlasma][0].Triangles.Values[iPar][TriKnots[:,1]] + Plasmas[iPlasma][0].Triangles.Values[iPar][TriKnots[:,2]])/3.
					else:				TValues =  Plasmas[iPlasma][0].Triangles.Values[iPar]

					if(PrintAllVolSumAve[iF][i] != 0):
						SumPar = np.sum(TriVols*TValues)
						if((PrintAllVolSumAve[iF][i] == 1) or (PrintAllVolSumAve[iF][i] == 1)):					
							Value   = SumPar
							PreName = "\tTotal        "
								
						if(PrintAllVolSumAve[iF][i] > 1):
							Value   = SumPar/AllVolume
							PreName = "\tAverage      "

						print_write(fid_stat2d, PreName+Pars[iF][i]+blanks[:15-len(Pars[iF][i])]+"= {:10.3e}".format(Value))
						if(save == "stat"):
							StatHeader += PreName[1:]+Pars[iF][i]+","
							StatValues.append(Value)
							StatFormats.append('{:17.4e},')
							if(save_norm):
								StatValues[-1] *=Facts[iF][i]
								StatHeader = StatHeader[:-1] + " (*" + FactsUnits[iF][i]+ "),"

#					Print total and average on plasma inside separatrix

					if(PrintInSepSumAve[iF][i] != 0):
						SumPar = np.sum(TriVols[InSepTri]*TValues[InSepTri])
						if((PrintInSepSumAve[iF][i] == 1) or (PrintInSepSumAve[iF][i] == 3)):
							Value   = SumPar
							PreName = "\tIn Sep. Tot  "

						if(PrintInSepSumAve[iF][i] > 1):
							Value   = SumPar/InSepVolume
							PreName = "\tIn Sep. Ave  "

						print_write(fid_stat2d, PreName+Pars[iF][i]+blanks[:15-len(Pars[iF][i])]+"= {:10.3e}".format(Value))

						if(save == "stat"):
							StatHeader += PreName[1:]+Pars[iF][i]+","
							StatValues.append(Value)
							StatFormats.append('{:17.4e},')
							if(save_norm):
								StatValues[-1] *=Facts[iF][i]
								StatHeader = StatHeader[:-1] + " (*" + FactsUnits[iF][i]+ "),"

#					Print total and average on plasma outside separatrix

					if(PrintOutSepSumAve[iF][i] != 0):
						SumPar = np.sum(TriVols*TValues) - np.sum(TriVols[InSepTri]*TValues[InSepTri])
						if((PrintOutSepSumAve[iF][i] == 1) or (PrintOutSepSumAve[iF][i] == 3)):
							Value   = SumPar
							PreName = "\tOut Sep. Tot "
						if(PrintOutSepSumAve[iF][i] > 1):
							Value   = SumPar/OutSepVolume
							PreName = "\tOut Sep. Ave "
						
						print_write(fid_stat2d, PreName+Pars[iF][i]+blanks[:15-len(Pars[iF][i])]+"= {:10.3e}".format(Value))
						if(save == "stat"):
							StatHeader += PreName[1:]+Pars[iF][i]+","
							StatValues.append(Value)
							StatFormats.append('{:17.4e},')
							if(save_norm):
								StatValues[-1] *=Facts[iF][i]
								StatHeader = StatHeader[:-1] + " (*" + FactsUnits[iF][i]+ "),"

#					Print total and average on plasma close to separatrix (but inside)

					if(PrintOnSepSumAve[iF][i] != 0):
						SumPar = np.sum(TriVols[OnSepTri]*TValues[OnSepTri])
						if((PrintOnSepSumAve[iF][i] == 1) or (PrintOnSepSumAve[iF][i] == 3)):
							Value   = SumPar
							PreName = "\tOn Sep. Tot  "
						if(PrintOnSepSumAve[iF][i] > 1):
							Value   = SumPar/OnSepVolume
							PreName = "\tOn Sep. Ave  "

						print_write(fid_stat2d, PreName+Pars[iF][i]+blanks[:15-len(Pars[iF][i])]+"= {:10.3e}".format(SumPar/OnSepVolume))
						if(save == "stat"):
							StatHeader += PreName[1:]+Pars[iF][i]+","
							StatValues.append(Value)
							StatFormats.append('{:17.4e},')
							if(save_norm):
								StatValues[-1] *=Facts[iF][i]
								StatHeader = StatHeader[:-1] + " (*" + FactsUnits[iF][i]+ "),"

#					=========== Plot section =============

					if(no_plot == 0):
						Vmin = Plasmas[iPlasma][0].Triangles.Values[iPar].min()
						Vmax = Plasmas[iPlasma][0].Triangles.Values[iPar].max()
						if((MaxRange[iF][i] != 0.) and (diff == 0)):
							MinValue = Vmax/MaxRange[iF][i]
							Plasmas[iPlasma][0].Triangles.Values[iPar] = np.where(Plasmas[iPlasma][0].Triangles.Values[iPar] > MinValue, Plasmas[iPlasma][0].Triangles.Values[iPar], MinValue)
						else:
							MinValue = Vmin

						MinSet,MinValue = search_min_max(FigNum[iF][i], PosNum[iF][i], MinMax[4])
						MaxSet,MaxValue = search_min_max( FigNum[iF][i], PosNum[iF][i], MinMax[5])
						if(not MinSet): MinValue = MinValue
						if(not MaxSet): MaxValue = Vmax

						if(diff == 1):
							Im.append(AxPl[iF][i].tricontourf(R, Z, TriKnots, Plasmas[iPlasma][0].Triangles.Values[iPar]*Facts[iF][i], cmap="jet"))
							if(Vmax > Vmin): AxPl[iF][i].tricontour(TripTriang, Plasmas[iPlasma][0].Triangles.Values[iPar]*Facts[iF][i], [0.], vmin=MinValue*Facts[iF][i], vmax=MaxValue*Facts[iF][i], colors='k')
						else:
							if((log_scale == 0) or (not PosPars[iF][i])):
								Im.append(AxPl[iF][i].tripcolor(TripTriang, Plasmas[iPlasma][0].Triangles.Values[iPar]*Facts[iF][i], vmin=MinValue*Facts[iF][i], vmax=MaxValue*Facts[iF][i]))
							else:
								if(MinValue > 0.):
									Im.append(AxPl[iF][i].tripcolor(TripTriang, Plasmas[iPlasma][0].Triangles.Values[iPar]*Facts[iF][i],norm=LogNorm(vmin=MinValue*Facts[iF][i], vmax=MaxValue*Facts[iF][i])))
								else:
									iPos = np.where(Plasmas[iPlasma][0].Triangles.Values[iPar] > 0.)[0]
									if(len(iPos) > 0):
										MinValue = np.min(Plasmas[iPlasma][0].Triangles.Values[iPar][iPos])
										VPos = np.where(Plasmas[iPlasma][0].Triangles.Values[iPar] > 0., Plasmas[iPlasma][0].Triangles.Values[iPar], MinValue)

										Im.append(AxPl[iF][i].tripcolor(TripTriang, VPos*Facts[iF][i],norm=LogNorm(vmin=MinValue*Facts[iF][i], vmax=MaxValue*Facts[iF][i])))
									else:
										Im.append(AxPl[iF][i].tripcolor(TripTriang, Plasmas[iPlasma][0].Triangles.Values[iPar]*Facts[iF][i]))

						plot2d_walls(AxPl[iF][i], Config.Walls)


						Xlims = AxPl[iF][i].get_xlim()
						Ylims = AxPl[iF][i].get_ylim()
						xMinSet,xMinValue = search_min_max(FigNum[iF][i], PosNum[iF][i], MinMax[0], Xlims[0])
						xMaxSet,xMaxValue = search_min_max( FigNum[iF][i], PosNum[iF][i], MinMax[1], Xlims[1])
						yMinSet,yMinValue = search_min_max(FigNum[iF][i], PosNum[iF][i], MinMax[2], Ylims[0])
						yMaxSet,yMaxValue = search_min_max( FigNum[iF][i], PosNum[iF][i], MinMax[3], Ylims[1])
						if(xMinSet or xMaxSet or yMinSet or yMaxSet):
							AxPl[iF][i].autoscale(enable=False, tight=True)
							AxPl[iF][i].set_xlim(xMinValue, xMaxValue)
							AxPl[iF][i].set_ylim(yMinValue, yMaxValue)

						cb = FigPl[iF][i].colorbar(Im[-1], ax=AxPl[iF][i])
						cb.set_label(Labels[iF][i])

#		=========== Export section =============

		if(save == "export"):
			
			ParExt="_"+NameExt
			if(iPlasma == 0):
				ExpPars   = ["Dens", "Temp"]
				ExpNames   = ["Dens", "Temp"]
			elif(iPlasma == 1):
				ExpPars   = ["Dens", "Temp","Nn", "Tn","vxn", "vyn", "vzn","Nm", "Tm",  "Pn",  "vxm", "vym", "vzm"]
				ExpNames   = ["Dens", "Temp","Nn", "Tn","vxn", "vyn", "vzn","Nm", "Tm",  "Pn",  "vxm", "vym", "vzm"]
			elif((iPlasma > 1) and (Plasmas[iPlasma][0].charge == 1)):
				ExpPars   = ["Nn", "Tn","vxn", "vyn", "vzn"]
				ExpNames   = ["Nn", "Tn","vxn", "vyn", "vzn"]
				ParExt=ParExt[:-1]
			else:
				ExpPars=[]

			for kPar in range(len(ExpPars)):			
				try:
						ePar = Plasmas[iPlasma][0].Triangles.VNames.index(ExpPars[kPar]+NameExt)
				except:
					print("WARNING: Not found ",ExpPars[kPar]+NameExt," for iPlasma = ", iPlasma)
					print("\tAvailable names=", Plasmas[iPlasma][0].Triangles.VNames)
					exit()
			
				np.savetxt(export_path+"/"+ExpNames[kPar]+ParExt+".txt", Plasmas[iPlasma][0].Triangles.Values[ePar] , delimiter=",", fmt="%15.7e", comments="")

	if(save != "none"):
		if(save == "stat"):
			save_stat("stat2d.csv", StatHeader, StatValues, StatFormats)
		elif((save == "print") or (save == "print_append")):
			fid_stat2d.close()
		elif((save == "export") or (no_plot != 0)):
			pass
		else:
			for i in range(len(Fig)):
				i_plot_file += 1
				if(one_plot != 1): Fig[i].set_size_inches(10.05,7.44)
				if(save == "pdf"):
					print("save pdf page=",i)
					pdf.savefig(Fig[i])
				else:
					Fig[i].savefig("cvs_plot2d_on_tri_t={:.3f}_{:d}.".format(tempus,i_plot_file)+save)
#		pyp.show(block=False)
#		pyp.close(0)
	elif(no_plot == 0):
		pyp.show()

	if(save == "pdf"):	pdf.close()

	print("cvs_plot2d_on_tri: Completed")

	return

def print_write(fid, string):
	print(string)
	if(fid != 0): fid.write(string+"\n")
	return
#
# This routine set to the minimum positive value all zero or negative values
#

def set_min_positive(Values):
	index_n = np.where(Values <= 0.)
	if(len(index_n[0]) > 0):
		index_p = np.where(Values > 0.)
		min_p   = np.min(Values[index_p])
		Values[index_n] = min_p

	return Values

def  search_min_max(FigNum, PosNum, MinMax, DefVal=None):
	nMinMax = int(len(MinMax)/3)
	if((nMinMax == 0) or (nMinMax*3 != len(MinMax))): return False, 0.

	for i in range(nMinMax):
		if((MinMax[3*i] == FigNum) and (MinMax[3*i+1] == PosNum) ): 
			
			return True, MinMax[3*i+2]

	if(DefVal != None): return True, DefVal
	else:				return False, 0.