#!/usr/bin/env python
# -*- coding: utf-8 -*-


from tkinter				import Tk, Text, filedialog, messagebox, BOTH, W, N, E, S, Menu, HORIZONTAL
from tkinter				import DoubleVar, IntVar, StringVar
from tkinter.messagebox 	import Message
from tkinter.ttk			import Entry, Menubutton, Checkbutton, Frame, Label, Style, LabelFrame, Notebook

import sys
import types
import os
import matplotlib
matplotlib.use('TkAgg')

import scipy.io
import numpy								as np
import numpy.matlib							as mat
from numpy									import arange, sin, pi

from routines.reset_config					import reset_config
from routines.utils_walls					import wall_fix_eirene_walls, wall_define_eirene_walls


from mesh.generate_soledge_mesh				import generate_soledge_mesh
from mesh.find_Xpoints						import find_Xpoints
from mesh.surround_Xpoints					import surround_Xpoints
from eirene.generate_eirene_grid			import generate_eirene_grid
from eirene.set_wall_boundary_condition		import set_wall_boundary_condition


from files.load_walls_file					import load_walls_file
from files.save_soledge_mesh_file			import save_soledge_mesh_file
from files.load_soledge_mesh_file			import load_soledge_mesh_file
from files.load_matlab_mesh_file			import load_matlab_mesh_file

from files.eqdsk_routines					import *
from files.save_npco_char					import save_npco_char
from files.save_elemente					import save_elemente
from files.save_neighbors					import save_neighbors
from files.save_puffs						import save_puffs
from files.save_surfaces					import save_surfaces
from files.save_wall_back_interp			import save_wall_back_interp
from files.save_eirene_triangles			import save_eirene_triangles
from files.load_eirene_triangles			import load_eirene_triangles
from files.load_wall_material_pump_puff		import load_wall_material_pump_puff
from files.save_wall_material_pump_puff		import save_wall_material_pump_puff

from routines.globals						import *
from routines.utils_walls					import walls_define_path
from routines.plot_fields_flux				import plot_fields_flux

from interfaces.mesh_gen					import mesh_gen
from interfaces.zone_gen					import zone_gen
from interfaces.eirene_gen					import eirene_gen
from interfaces.progressbar 				import ProgressBar
from interfaces.edit_wall					import edit_wall

class GridGen(Frame):

	def __init__(self, parent):
		Frame.__init__(self, parent)

		self.parent	= parent
		self.Config	= reset_config()

		self.MaxXpoints = 4
		self.MaxBariers = 6
		nPadx = 5
		nPady = 5
		self.parent.title("Soledge grid generation")
		self.pack(fill=BOTH, expand=True)

		self.columnconfigure(1, weight=1)
#		self.columnconfigure(1, pad=7)
#		self.rowconfigure(3, weight=1)
#		self.rowconfigure(5, pad=7)


		menubar = Menu(self.master)
		filemenu = Menu(menubar, tearoff=0)
		filemenu.add_command(label="Load Mesh", command=self.cmd_load_mesh)
		filemenu.add_command(label="Save Mesh", command=self.cmd_save_mesh)
		filemenu.add_command(label="Save Mesh As..", command=self.cmd_save_mesh_as)
		filemenu.add_separator()
		filemenu.add_command(label="Save Eirene", command=self.cmd_save_eirene)
		filemenu.add_separator()
		filemenu.add_command(label="Load Fields", command=self.cmd_load_mag_file)
		filemenu.add_command(label="Load Wall Shape", command=self.cmd_load_wall_file)
		filemenu.add_separator()
		filemenu.add_command(label="Import MatLab Mesh", command=self.cmd_import_matlab_mesh)
		filemenu.add_command(label="Import Wall&Puff", command=self.cmd_import_wall_pump_puff)
		filemenu.add_separator()
		filemenu.add_command(label="Export Fields", command=self.cmd_export_mag_file)
		filemenu.add_command(label="Export Wall shape", command=self.cmd_export_wall)
		filemenu.add_command(label="Export Wall&Puff", command=self.cmd_export_wall_pump_puff)
		filemenu.add_separator()	
		filemenu.add_command(label="Exit", command=self.quit)
		menubar.add_cascade(label="File", menu=filemenu)

		toolsmenu = Menu(menubar, tearoff=0)
		toolsmenu.add_command(label="Edit wall", 		command=self.cmd_edit_wall)

		ShowFieldMenu = Menu(self.master)
		ShowFieldMenu.add_command(label="Loaded",			command=self.cmd_plot_fields_loaded)
		ShowFieldMenu.add_command(label="Smoothed",		command=self.cmd_plot_fields_smoothed)
		ShowFieldMenu.add_command(label="Difference", 	command=self.cmd_plot_fields_difference)
		toolsmenu.add_cascade(label="Show fields",		menu=ShowFieldMenu)

		menubar.add_cascade(label="Tools", menu=toolsmenu)

		OptionsMenu = Menu(menubar, tearoff=0)

		FileTypeMenu = Menu(self.master)

		self.Soledge2DCheckValue = IntVar()
		self.Soledge2DCheckValue.set(1)
		FileTypeMenu.add_checkbutton(label="Soledge2D", 		variable=self.Soledge2DCheckValue)
		self.Soledge2DCheckValue.trace("w", self.cmd_Soledge2DCheckValue)

		self.Soledge3XCheckValue = IntVar()
		self.Soledge3XCheckValue.set(0)
		FileTypeMenu.add_checkbutton(label="Soledge3X",			variable=self.Soledge3XCheckValue)
		self.Soledge3XCheckValue.trace("w", self.cmd_Soledge3XCheckValue)
		self.FileType = FILE_TYPE_SOLEDGE2D
		
		OptionsMenu.add_cascade(label="File Type",				menu=FileTypeMenu)

		ShowWallMenu = Menu(self.master)

		self.ShowPlasmaWallCheckValue = IntVar()
		self.ShowPlasmaWallCheckValue.set(1)
		ShowWallMenu.add_checkbutton(label="Soledge", 		variable=self.ShowPlasmaWallCheckValue)
		self.ShowPlasmaWallCheckValue.trace("w", self.cmd_ShowPlasmaWallCheckValue)

		self.ShowEireneWallCheckValue = IntVar()
		self.ShowEireneWallCheckValue.set(1)
		ShowWallMenu.add_checkbutton(label="Eirene", 		variable=self.ShowEireneWallCheckValue)
		self.ShowEireneWallCheckValue.trace("w", self.cmd_ShowEireneWallCheckValue)

		self.ShowExtraWallCheckValue = IntVar()
		self.ShowExtraWallCheckValue.set(0)
		ShowWallMenu.add_checkbutton(label="Extra", 		variable=self.ShowExtraWallCheckValue)
		self.ShowExtraWallCheckValue.trace("w", self.cmd_ShowExtraWallCheckValue)
	
		OptionsMenu.add_cascade(label="Show wall",				menu=ShowWallMenu)

		EireneSurfsMenu = Menu(self.master)

		self.EireneSurfsPointKnotsCheckValue = IntVar()
		self.EireneSurfsPointKnotsCheckValue.set(1)
		EireneSurfsMenu.add_checkbutton(label="Point to nodes", 	variable=self.EireneSurfsPointKnotsCheckValue)
		self.EireneSurfsPointKnotsCheckValue.trace("w", self.cmd_EireneSurfsPointKnotsCheckValue)

		self.EireneSurfsPointFreeCheckValue = IntVar()
		self.EireneSurfsPointFreeCheckValue.set(0)
		EireneSurfsMenu.add_checkbutton(label="Free points", 		variable=self.EireneSurfsPointFreeCheckValue)
		self.EireneSurfsPointFreeCheckValue.trace("w", self.cmd_EireneSurfsPointFreeCheckValue)

		OptionsMenu.add_cascade(label="Eirene surfaces",			menu=EireneSurfsMenu)

		self.SameRaiseFactorCheckValue = IntVar()
		self.SameRaiseFactorCheckValue.set(1)
		OptionsMenu.add_checkbutton(label="Same Raise Factor",			variable=self.SameRaiseFactorCheckValue)
		self.SameRaiseFactorCheckValue.trace("w", self.cmd_SameRaiseFactorCheckValue)

		self.IgnoreFieldsErrorsCheckValue = IntVar()
		self.IgnoreFieldsErrorsCheckValue.set(0)
		OptionsMenu.add_checkbutton(label="Ignore fields errors",			variable=self.IgnoreFieldsErrorsCheckValue)
#		self.IgnoreFieldsErrorsCheckValue.trace("w", self.cmd_IgnoreFieldsErrorsCheckValue)

		menubar.add_cascade(label="Options", 					menu=OptionsMenu)

		helpmenu = Menu(menubar, tearoff=0)
		helpmenu.add_command(label="About...", command=self.cmd_donothing)
		menubar.add_cascade(label="Help", menu=helpmenu)

		self.master.config(menu=menubar)

		
		self.GenNames	 = ["Zones", "Mesh", "Eirene"] 
		self.GenWidgets  = []
		
		self.GenNotebook = Notebook(self)
		self.zone_gen = zone_gen(self.GenNotebook, self)
		self.mesh_gen = mesh_gen(self.GenNotebook, self)
		self.eirene_gen = eirene_gen(self.GenNotebook, self)
		
		self.GenNotebook.add(self.zone_gen,   text="Zones")
		self.GenNotebook.add(self.mesh_gen,   text="Mesh")
		self.GenNotebook.add(self.eirene_gen, text="Eirene")

		self.GenNotebook.grid(row=0, column=0, pady=nPady, padx=nPadx)
		self.GenNotebook.tab(1, state= "disabled")
		self.GenNotebook.tab(2, state= "disabled")
		self.cmd_Soledge2DCheckValue(0, 0, 0)
	
#		Variable definitions

		self.eqdsk_file		= False
		self.MeshFile		= False
		self.resize_mode	= False
		self.sel_mode_d		= False
		self.CheckPlot		= False
		self.InNewLimit		= False
		self.InMoveCut		= False

#	Internal functions
#===============================================================================

	def cmd_plot_fields_loaded(self):
		plot_fields_flux(self.Config, plot_type=0)

	def cmd_plot_fields_smoothed(self):
		plot_fields_flux(self.Config, plot_type=1)

	def cmd_plot_fields_difference(self):
		plot_fields_flux(self.Config, plot_type=2)

	def cmd_Soledge2DCheckValue(self, a, b, c):
		if(self.Soledge2DCheckValue.get() == 1): 
			self.Soledge3XCheckValue.set(0)
			self.zone_gen.S3xSettingsFrame.state(["disabled"])
			self.zone_gen.S3xNPhiEntry.state(["disabled"])
			self.zone_gen.S3xR0Entry.state(["disabled"])
			self.zone_gen.S3xa0Entry.state(["disabled"])
			self.zone_gen.S3xZ0Entry.state(["disabled"])
			self.zone_gen.S3xBt0Entry.state(["disabled"])
			self.zone_gen.S3xResetButton.state(["disabled"])

			self.FileType = FILE_TYPE_SOLEDGE2D

	def cmd_Soledge3XCheckValue(self, a, b, c):
		if(self.Soledge3XCheckValue.get() == 1):
			self.Soledge2DCheckValue.set(0)
			self.zone_gen.S3xSettingsFrame.state(["!disabled"])
			self.zone_gen.S3xNPhiEntry.state(["!disabled"])
			self.zone_gen.S3xR0Entry.state(["!disabled"])
			self.zone_gen.S3xa0Entry.state(["!disabled"])
			self.zone_gen.S3xZ0Entry.state(["!disabled"])
			self.zone_gen.S3xBt0Entry.state(["!disabled"])
			self.zone_gen.S3xResetButton.state(["!disabled"])
			self.FileType = FILE_TYPE_SOLEDGE3X

	def cmd_EireneSurfsPointKnotsCheckValue(self, a, b, c):
		if(self.EireneSurfsPointKnotsCheckValue.get() == 1): 
			self.EireneSurfsPointFreeCheckValue.set(0)
			self.FileType = FILE_TYPE_SOLEDGE2D

	def cmd_EireneSurfsPointFreeCheckValue(self, a, b, c):
		if(self.EireneSurfsPointFreeCheckValue.get() == 1):
			self.EireneSurfsPointKnotsCheckValue.set(0)
			self.FileType = FILE_TYPE_SOLEDGE3X

	def cmd_SameRaiseFactorCheckValue(self, a, b, c):
		if(self.SameRaiseFactorCheckValue.get() == 0): 
			self.zone_gen.RaiseFacSEntry.state(["!disabled"])
			self.zone_gen.RaiseFacEEntry.state(["!disabled"])
			self.zone_gen.RaiseFacWEntry.state(["!disabled"])
		else:
			self.zone_gen.RaiseFacSEntry.state(["disabled"])
			self.zone_gen.RaiseFacEEntry.state(["disabled"])
			self.zone_gen.RaiseFacWEntry.state(["disabled"])
			self.zone_gen.RaiseFacSEntryValue.set(self.zone_gen.RaiseFacNEntryValue.get())
			self.zone_gen.RaiseFacEEntryValue.set(self.zone_gen.RaiseFacNEntryValue.get())
			self.zone_gen.RaiseFacWEntryValue.set(self.zone_gen.RaiseFacNEntryValue.get())

			self.FileType = FILE_TYPE_SOLEDGE2D

	def cmd_ShowPlasmaWallCheckValue(self, a, b, c):
		self.zone_gen.update_figure()
		self.mesh_gen.update_figure()
		self.eirene_gen.update_figure()
		return

	def cmd_ShowEireneWallCheckValue(self, a, b, c):
		self.zone_gen.update_figure()
		self.mesh_gen.update_figure()
		self.eirene_gen.update_figure()
		return

	def cmd_ShowExtraWallCheckValue(self, a, b, c):
		self.zone_gen.update_figure()
		self.mesh_gen.update_figure()
		self.eirene_gen.update_figure()
		return

	def cmd_donothing(self, empty):
		return

	def cmd_load_eqdsk(self):
		self.eqdsk_file = filedialog.askopenfilename(filetypes=(("EQDSK files", "*.eqdsk"),("EQDSK files", "*.dat"),
															("Mat files", "*.mat)"), ("Mat files", "*.m"), ("All files", "*.*") ))

		if(self.eqdsk_file):
#			try:
				eqdsk = load_eqdsk_file(self.eqdsk_file)
				self.Config.nr			= eval(self.zone_gen.NREntryValue.get())
				self.Config.nz			= eval(self.zone_gen.NZEntryValue.get())
				self.Config.OutLength2	= eval(self.zone_gen.OutLengthEntryValue.get())
				self.Config.Smooth2		= eval(self.zone_gen.SmoothInEntryValue.get())
				rz_b_psi = eqdsk_compute_fields(eqdsk, nr=self.Config.nr, nz=self.Config.nz, smooth=self.Config.Smooth2, out_length=self.Config.OutLength2, IgnoreErrors=self.IgnoreFieldsErrorsCheckValue.get())

				self.Config.Walls				= [types.SimpleNamespace()]
				self.Config.Walls[0].Rwall		= np.copy(eqdsk.rlim)
				self.Config.Walls[0].Zwall		= np.copy(eqdsk.zlim)
				self.Config.Walls[0].Type		= EXTERNAL_PLASMA_WALL
				self.Config.Walls[0].LineType	= np.array([0,0,0])

				self.Config.r2			= np.copy(rz_b_psi[:,:,0])
				self.Config.z2			= np.copy(rz_b_psi[:,:,1])
				self.Config.Br2			= np.copy(rz_b_psi[:,:,2])
				self.Config.Bz2			= np.copy(rz_b_psi[:,:,3])
				self.Config.Bphi2		= np.copy(rz_b_psi[:,:,4])
				if(eqdsk.simag > eqdsk.sibry):
					self.Config.flux2		= np.copy(-rz_b_psi[:,:,5])
				else:
					self.Config.flux2		= np.copy(rz_b_psi[:,:,5])

				eqdsk		= 0
				rz_b_psi	= 0

				self.Config.r2D   				= self.Config.r2
				self.Config.z2D					= self.Config.z2
				self.Config.Br2D				= self.Config.Br2
				self.Config.Bphi2D				= self.Config.Bphi2
				self.Config.Bz2D				= self.Config.Bz2
				self.Config.flux2D				= self.Config.flux2
				self.Config.Smooth2D			= self.Config.Smooth2
				self.Config.extrapol_val	= 0.
				self.Config.Smooth2DExt		= 0.
				self.Config.raise_value 	= 0.
				self.Config.raise_power 	= 1
				
				self.Config.in_equ_OK	= True
				self.Config.wall_OK		= True
				self.resize_mode		= True
				walls_define_path(self.Config)
				self.post_loaded_fields()

				self.Config.flux2_x_psi = np.empty((len(self.Config.X_points)), dtype='f8')
				for k in range(len(self.Config.X_points)): self.Config.flux2_x_psi[k] = self.Config.X_points[k].psi

#			except:
#				messagebox.showerror("Open Source File", "Failed to read file\n'%s'" % eqdsk_file)
		return

#	Reading cvs wall file
#---------------------------

	def cmd_load_wall_file(self):
		wall_file = filedialog.askopenfilename(filetypes=(("cvs files", "*.csv"), ("Mat files", "*.mat"), ("Mat files", "*.m"), ("All files", "*.*")))

		if(not wall_file): return

		FileExt = os.path.splitext(os.path.split(wall_file)[1])[1]

		if(FileExt == ".csv"):
			try:
				wall	 = np.loadtxt(wall_file, delimiter=',',comments='"')

				self.Config.Walls				= [types.SimpleNamespace()]
				self.Config.Walls[0].Rwall  	= wall[:,0]
				self.Config.Walls[0].Zwall 		= wall[:,1]
				self.Config.Walls[0].Type		= EXTERNAL_PLASMA_WALL
				self.Config.Walls[0].LineType	= np.array([0,0,1])
				self.Config.Walls[0].Changed	= True
				wall = None
			except:
				messagebox.showerror("Read Wall File", "Failed to read file {:s}\n".format(wall_file))
				return

		elif(FileExt == ".txt"):
			try:
				Walls = np.load_walls_file(wall_file)
				self.Config.Walls = Walls
			except:
				messagebox.showerror("Read Wall File", "Failed to read file {:s}\n".format(wall_file))
				return
	
		elif((FileExt == ".mat") or (FileExt == ".m")):
			try:
				MatData = scipy.io.loadmat(wall_file)
				try:
					self.Config.Walls				= [types.SimpleNamespace()]
					self.Config.Walls[0].Type		= EXTERNAL_PLASMA_WALL
					self.Config.Walls[0].LineType	= np.array([0,0,0])
					self.Config.Walls[0].Rwall  	= MatData['R_wall'][0]
					self.Config.Walls[0].Zwall  	= MatData['Z_wall'][0]
				except:
					self.Config.Walls[0].Rwall		= MatData['Rwall'][:,0]
					self.Config.Walls[0].Zwall		= MatData['Zwall'][:,0]
				self.Config.Walls[0].Changed = True
				MatData	= None
			except:
				messagebox.showerror("Read Wall File", "Failed to read file {:s}\n".format(wall_file))
				return
		else:
			messagebox.showerror("Read Wall File", "Unknow file type for {:s}\n".format(wall_file))

		self.update_wall_shape()


		return


# reading mag file
#---------------------------

	def cmd_load_mag_file(self):
		mag_file = filedialog.askopenfilename(filetypes=(("EQDSK files", "*.eqdsk"),("EQDSK files", "*.dat"),("cvs files", "*.csv"),
																					("Mat files", "*.mat"), ("Mat files", "*.m"),  ("All files", "*.*")))

		if(not mag_file): return

		FileExt = os.path.splitext(os.path.split(mag_file)[1])[1]

		if((FileExt == ".eqdsk") or (FileExt == ".dat")):
			self.eqdsk_file = mag_file
			try:
				eqdsk = load_eqdsk_file(self.eqdsk_file)
				self.Config.nr			= eval(self.zone_gen.NREntryValue.get())
				self.Config.nz			= eval(self.zone_gen.NZEntryValue.get())
				self.Config.OutLength2	= eval(self.zone_gen.OutLengthEntryValue.get())
				self.Config.Smooth2		= eval(self.zone_gen.SmoothInEntryValue.get())
				rz_b_psi = eqdsk_compute_fields(eqdsk, nr=self.Config.nr, nz=self.Config.nz, smooth=self.Config.Smooth2, out_length=self.Config.OutLength2, IgnoreErrors=self.IgnoreFieldsErrorsCheckValue.get())

				self.Config.r2			= np.copy(rz_b_psi[:,:,0])
				self.Config.z2			= np.copy(rz_b_psi[:,:,1])
				self.Config.Br2			= np.copy(rz_b_psi[:,:,2])
				self.Config.Bz2			= np.copy(rz_b_psi[:,:,3])
				self.Config.Bphi2		= np.copy(rz_b_psi[:,:,4])
				if(eqdsk.simag > eqdsk.sibry):
					self.Config.flux2		= np.copy(-rz_b_psi[:,:,5])
				else:
					self.Config.flux2		= np.copy(rz_b_psi[:,:,5])

				self.Config.Walls				= [types.SimpleNamespace()]
				self.Config.Walls[0].Rwall		= np.copy(eqdsk.rlim)
				self.Config.Walls[0].Zwall		= np.copy(eqdsk.zlim)
				self.Config.Walls[0].Changed	= True
				self.Config.Walls[0].Type		= EXTERNAL_PLASMA_WALL
				self.Config.Walls[0].LineType	= np.array([0,0,0])
				walls_define_path(self.Config)
				self.Config.wall_OK		= True

				eqdsk		= 0
				rz_b_psi	= 0
			except:
				messagebox.showerror("Read Magnetic File", "Failed to read file {:s}\n".format(mag_file))
				return

		elif(FileExt == ".csv"):
			try:
				rz_b_psi = np.loadtxt(mag_file, delimiter=',',comments='"')

#				Check for data with only one istance of R and Z
				
				HasBlock= False
				iBlocks = np.where(rz_b_psi[:,0]==2000.)[0]									#check for organization (Z,R) with only first values of R &Z
				if(len(iBlocks) > 0):
					HasBlock= True
					nR = np.min(iBlocks)
					nZ = int(rz_b_psi.shape[0]/nR)
					rz_b_psi[:,2] = np.copy(rz_b_psi[:,2].reshape(nZ,nR).T.reshape(-1))
					rz_b_psi[:,3] = np.copy(rz_b_psi[:,3].reshape(nZ,nR).T.reshape(-1))
					rz_b_psi[:,4] = np.copy(rz_b_psi[:,4].reshape(nZ,nR).T.reshape(-1))
					rz_b_psi[:,5] = np.copy(rz_b_psi[:,5].reshape(nZ,nR).T.reshape(-1))
				else:																		#check 
					iBlocks = np.where(rz_b_psi[:,0]==1000.)[0]								#check for organization (R,Z) with only first values of R &Z
					if(len(iBlocks) > 0):
						HasBlock= True
						nR = np.min(iBlocks)
						nZ = int(rz_b_psi.shape[0]/nR)
			
				if(HasBlock):
					r_values = np.copy(rz_b_psi[:nR,0])
					for iR in range(nR): 			rz_b_psi[iR*nZ:(iR+1)*nZ,0] = r_values[iR]
					for iR in range(1,nR):			rz_b_psi[iR*nZ:(iR+1)*nZ,1] = rz_b_psi[:nZ,1]
					

				self.Config.r2			= np.copy(rz_b_psi[:,0])
				self.Config.z2			= np.copy(rz_b_psi[:,1])
				self.Config.Br2			= np.copy(rz_b_psi[:,2])
				self.Config.Bphi2		= np.copy(rz_b_psi[:,3])
				self.Config.Bz2			= np.copy(rz_b_psi[:,4])
				self.Config.flux2		= np.copy(rz_b_psi[:,5])
				rz_b_psi	= 0

				nTot	= self.Config.r2.shape[0]
				R0		= self.Config.r2[0]
				for i in range(nTot):
					if(R0 != self.Config.r2[i]):
						nzPoints = i											#number of points on each R
						break
				nxPoints = int(nTot/nzPoints)									#Number of points on each Z
				if nxPoints*nzPoints != nTot:
					print("\tError on Magnetic file: not matrix R,Z data")
					exit()
				self.Config.OutLength2	= -1.
				self.Config.Smooth2		= 0.
				self.Config.r2   	= self.Config.r2[:].reshape(nxPoints,nzPoints).T
				self.Config.z2		= self.Config.z2[:].reshape(nxPoints,nzPoints).T
				self.Config.Br2		= self.Config.Br2[:].reshape(nxPoints,nzPoints).T
				self.Config.Bphi2	= self.Config.Bphi2[:].reshape(nxPoints,nzPoints).T
				self.Config.Bz2		= self.Config.Bz2[:].reshape(nxPoints,nzPoints).T
				self.Config.flux2	= self.Config.flux2[:].reshape(nxPoints,nzPoints).T

			except:
				messagebox.showerror("Read Magnetic File", "Failed to read file {:s}\n".format(mag_file))
				return

		elif((FileExt == ".mat") or (FileExt == ".m")):
			try:
				MatData = scipy.io.loadmat(mag_file)

				try:
					r1 = MatData['R_psi'][0]
					z1 = MatData['Z_psi'][0]
					self.Config.r2			= mat.repmat(r1, z1.shape[0],1)
					self.Config.z2			= mat.repmat(z1, r1.shape[0],1).T

					self.Config.Br2			= MatData['Br']
					self.Config.Bphi2		= MatData['Bphi']
					self.Config.Bz2			= MatData['Bz']
					self.Config.flux2		= MatData['psi']
					if(r1.shape[0] != z1.shape[0]):
						if(self.Config.Br2.shape[0] != z1.shape[0]):	self.Config.Br2		= self.Config.Br2.T
						if(self.Config.Bphi2.shape[0] != z1.shape[0]):	self.Config.Bphi2	= self.Config.Bphi2.T
						if(self.Config.Bz2.shape[0] != z1.shape[0]):	self.Config.Bz2		= self.Config.Bz2.T
						if(self.Config.flux2.shape[0] != z1.shape[0]):	self.Config.flux2	= self.Config.flux2.T
				except:
					self.Config.r2			= MatData['r2D']
					self.Config.z2			= MatData['z2D']

					self.Config.Br2			= MatData['Br2D']
					self.Config.Bphi2		= MatData['Bphi2D']
					self.Config.Bz2			= MatData['Bz2D']
					self.Config.flux2		= MatData['flux2D']

				self.Config.OutLength2	= -1.
				self.Config.Smooth2		= 0.
				try:
					try:
						self.Config.Walls				= [types.SimpleNamespace()]
						self.Config.Walls[0].Type		= EXTERNAL_PLASMA_WALL
						self.Config.Walls[0].LineType	= np.array([0,0,0])
						self.Config.Walls[0].Rwall		= MatData['R_wall'][0]
						self.Config.Walls[0].Zwall		= MatData['Z_wall'][0]
					except:
						self.Config.Walls[0].Rwall	= MatData['Rwall'][:,0]
						self.Config.Walls[0].Zwall	= MatData['Zwall'][:,0]

					walls_define_path(self.Config)
					self.Config.wall_OK		= True
				except:
					pass

				MatData	= None; r1 = None; z1 = None

			except:
				messagebox.showerror("Read Mag File", "Failed to read file {:s}\n".format(mag_file))
				return
		else:
			messagebox.showerror("Read Mag File", "Unknow file type for {:s}\n".format(mag_file))


		self.Config.in_equ_OK	= True

		self.Config.nr = self.Config.r2.shape[1]
		self.Config.nz = self.Config.r2.shape[0]

		self.zone_gen.NREntryValue.set("{:d}".format(self.Config.nr))
		self.zone_gen.NZEntryValue.set("{:d}".format(self.Config.nz))

		self.Config.r2D   		= np.copy(self.Config.r2)
		self.Config.z2D			= np.copy(self.Config.z2)
		self.Config.Br2D		= np.copy(self.Config.Br2)
		self.Config.Bphi2D		= np.copy(self.Config.Bphi2)
		self.Config.Bz2D		= np.copy(self.Config.Bz2)
		self.Config.flux2D		= np.copy(self.Config.flux2)
		self.Config.Smooth2D		= self.Config.Smooth2
		self.Config.extrapol_val	= 0.
		self.Config.Smooth2DExt		= 0.
		self.Config.raise_value		= 0.
		self.Config.raise_power 	= 1

		self.Config.in_equ_OK	= True
		self.resize_mode		= True
		self.post_loaded_fields()

		self.Config.flux2_x_psi = np.empty((len(self.Config.X_points)), dtype='f8')
		for k in range(len(self.Config.X_points)): self.Config.flux2_x_psi[k] = self.Config.X_points[k].psi

#=============================
#	load mesh previously saved
#=============================

	def cmd_load_mesh(self, File=None):

		if(File is None):
			File = filedialog.askopenfilename(title="Select Mesh file", filetypes=(("Mesh files", "*.h5"), ("All files", "*.*") ))
			if(not File): return
#		try:
		if(True):
			self.busy()

			self.Config = load_soledge_mesh_file(File)
			
			if(self.Config.Mesh_OK):
				try:
					EireneDir = self.get_eirene_dirname(File)

					self.Eirene = load_eirene_triangles(EireneDir+"triangles.h5", FileType=self.FileType)
					self.Eirene.Wall.EWalls = wall_fix_eirene_walls(self.Config, self.Eirene.Wall.EWalls)
					self.Config.Eirene_OK = True
				except:
					self.Config.Eirene_OK = False

			self.post_loaded_mesh()
			self.notbusy()
			self.MeshFile = File

		else:
#		except:
			messagebox.showerror("Reading mesh data file", "Failed to read file: {:s}".format(self.MeshFile))


#=============================
#	post_loaded_mesh
#=============================

	def post_loaded_mesh(self):
		Config = self.Config

		if(Config.extrapol_OK):	self.zone_gen.update_table_X_points()

		if(Config.wall_OK):	walls_define_path(Config)	
			
		if(Config.in_equ_OK):
			self.zone_gen.NREntryValue.set("{:d}".format(Config.nr))
			self.zone_gen.NZEntryValue.set("{:d}".format(Config.nz))
			self.zone_gen.OutLengthEntryValue.set("{:0.2f}".format(Config.OutLength2))
			self.zone_gen.SmoothInEntryValue.set("{:0.2f}".format(Config.Smooth2D))

			if(Config.extrapol_OK):
				self.zone_gen.ExtScale.set(Config.extrapol_val)
				self.zone_gen.SmoothExtEntryValue.set("{:0.2f}".format(Config.Smooth2DExt))
				for i in range(len(self.zone_gen.RaiseFacEntryValues)):
					self.zone_gen.RaiseFacEntryValues[i].set("{:0.2e}".format(Config.raise_value[i]))
				self.zone_gen.RaisePowerOptMenuValue.set("{:d}".format(Config.raise_power))
				
				if(Config.xPoints_OK):
					self.zone_gen.BordCoreEntryValue.set("{:0.5f}".format(Config.psicore))
					self.zone_gen.BordOutEntryValue.set("{:0.5f}".format(Config.psiout))
					
					indexes	=  np.argsort(np.array([Config.X_points[k].index for k in range(len(Config.X_points))]))
					for i in range(min(len(Config.X_points),self.zone_gen.MaxXpoints)):
						self.zone_gen.BordXPtEntryPsiAValue[i].set("{:.5f}".format(Config.X_points[indexes[i]].cut[0].psilim))
						self.zone_gen.BordXPtEntryPsiCValue[i].set("{:.5f}".format(Config.X_points[indexes[i]].cut[2].psilim))
				else:
					self.post_loaded_fields()

					self.Config.flux2_x_psi = np.empty((len(self.Config.X_points)), dtype='f8')
					for k in range(len(self.Config.X_points)): self.Config.flux2_x_psi[k] = self.Config.X_points[k].psi
		
#		if(self.Config.Frontiers_OK and (not self.Config.MagZones_OK)):
#			self.zone_gen.cmd_define_mag_zones()

		self.zone_gen.resize_mode = True
		self.zone_gen.update_figure()
		
		if(Config.MagZones_OK):
			Config.MagZones_OK = False									# force update
			self.set_MagZones_OK(Config, True, Propagate=False)
		if(Config.Mesh_OK):
			Config.Mesh_OK = False										# force update
			self.set_Mesh_OK(Config, True, Propagate=False)
		if(Config.Eirene_OK):											# force update
			Config.Eirene_OK = False
			self.set_Eirene_OK(Config, True, Propagate=False)
			

#============	
#	save mesh
#============

	def cmd_save_mesh(self):
		if(not self.MeshFile): return
		if(self.MeshFile):
			self.busy()

			if(self.Config.MagGrid_OK and (not self.Config.Mesh_OK)): self.mesh_gen.cmd_soledge_mesh()

			if(self.FileType == FILE_TYPE_SOLEDGE3X):
				self.Config.Soledge3x = types.SimpleNamespace()
				self.Config.Soledge3x.NPhi		= self.zone_gen.S3xNPhiEntryValue.get()
				self.Config.Soledge3x.R0		= self.zone_gen.S3xR0EntryValue.get()
				self.Config.Soledge3x.Z0		= self.zone_gen.S3xZ0EntryValue.get()
				self.Config.Soledge3x.deltaR	= self.zone_gen.S3xa0EntryValue.get()
				self.Config.Soledge3x.Bcore		= self.zone_gen.S3xBt0EntryValue.get()

			save_soledge_mesh_file(self.MeshFile, self.Config)

#			if(self.Config.Mesh_OK) and (not self.Config.Eirene_OK):	self.eirene_gen.cmd_plasma_triangles_gen()

			if(self.Config.Eirene_OK): 
				EireneDir = self.get_eirene_dirname(self.MeshFile)
				save_eirene_triangles(EireneDir, self.Eirene, FileType=self.FileType)

			self.notbusy()
						
		return

#=============================
#	get_eirene_dirname
#=============================

	def get_eirene_dirname(self, MeshFile):
		MeshDir  = os.path.dirname(MeshFile)
		if(len(MeshDir) != 0): MeshDir = MeshDir + "/"

		BaseName = os.path.splitext(os.path.split(MeshFile)[1])[0]
		ExtName  = os.path.splitext(os.path.split(MeshFile)[1])[1]
		if(BaseName == "mesh"):
			EireneDir = MeshDir
#		elif((BaseName[0:5] == "mesh_") or (BaseName[0:5] == "mesh-")):
#			EireneDir = MeshDir + "/" +  BaseName[5:] + "/"
#		elif(BaseName[0:4] == "mesh"):
#			EireneDir = MeshDir + "/" +  BaseName[4:] + "/"
		else:
			EireneDir = MeshDir +  BaseName + "/"
				
		return EireneDir
		
#=============================
#	save mesh selecting folder
#=============================

	def cmd_save_mesh_as(self):
		if(self.MeshFile):
			MeshDir  = os.path.dirname(self.MeshFile)
			MeshFile = os.path.split(self.MeshFile)[1]
			MeshFile = filedialog.asksaveasfilename(title="Select file to store mesh data", initialfile=MeshFile, initialdir=MeshDir  )
		else:
			MeshFile = filedialog.asksaveasfilename(title="Select file to store mesh data")
		if(not MeshFile): return
		
		self.MeshFile = MeshFile
		self.cmd_save_mesh()					

#=============================
#	save eirene files
#=============================

	def cmd_save_eirene(self):
		if((self.MeshFile) and (self.Config.Eirene_OK)):
			
			EireneDir = self.get_eirene_dirname(self.MeshFile)	
			try:
				os.mkdir(EireneDir)
			except OSError:
				pass
		

			set_wall_boundary_condition(self.Eirene)

			save_eirene_triangles(EireneDir, self.Eirene, FileType=self.FileType)
			save_neighbors(EireneDir, self.Eirene, FileType=self.FileType)
			save_puffs(EireneDir, self.Eirene, FileType=self.FileType)
			save_npco_char(EireneDir, self.Eirene, FileType=self.FileType)
			save_elemente(EireneDir, self.Eirene, FileType=self.FileType)
			save_surfaces(EireneDir, self.Eirene, FileType=self.FileType)
				
#			save_wall_back_interp(EireneDir, self.Eirene)
			


#=============================
#	import matlsab mesh
#=============================

	def	cmd_import_matlab_mesh(self):
		MatFile = filedialog.askopenfilename(title="Select Matlab file to import", filetypes=(("mat", "*.mat"), ("All files", "*.*")))
		self.Config, self.Eirene = load_matlab_mesh_file(MatFile)

		self.post_loaded_mesh()

#=============================
#	import  wall, pumps and puffs
#=============================

	def cmd_import_wall_pump_puff(self):
		if(self.MeshFile):
			MeshDir  = os.path.dirname(self.MeshFile)
		else:
			MeshDir	 = "./"
				
		WallFile = filedialog.askopenfilename(title="Select file to import Wall/Pump & Puff", filetypes=(("text", "*.txt"), ("All files", "*.*")), initialdir=MeshDir )
		if(WallFile):
			load_wall_material_pump_puff(WallFile, self.Eirene)
		
			self.eirene_gen.material_tree_view_fill()
			self.eirene_gen.pump_tree_view_fill()
			self.eirene_gen.puff_tree_view_fill()
			self.eirene_gen.eirene_zones_tree_view_fill()
			self.eirene_gen.update_figure()


#=============================
#	export wall, pumps and puffs
#=============================
			
	def cmd_export_wall_pump_puff(self):
		if(self.MeshFile):
			MeshDir  = os.path.dirname(self.MeshFile)
		else:
			MeshDir	 = "./"	
					
		WallFile = filedialog.asksaveasfilename(title="Select file to export Wall/Pump & Puff", filetypes=(("text", "*.txt"), ("All files", "*.*") ), initialdir=MeshDir )
		if(WallFile):
			save_wall_material_pump_puff(WallFile, self.Eirene)

#=============================
#	export magnetic file to eqdsk, matlab or csv
#=============================

	def cmd_export_mag_file(self):
		if(self.MeshFile):
			MeshDir  = os.path.dirname(self.MeshFile)
		else:
			MeshDir	 = "./"	
					
		mag_file = filedialog.asksaveasfilename(title="Select file to export mag file", filetypes=(("cvs files", "*.csv"), ("EQDSK files", "*.eqdsk"), 
																																																	("Mat files", "*.mat)"), ("Mat files", "*.m"), ("All files", "*.*") ), initialdir=MeshDir )
		if(not mag_file): return

		FileExt = os.path.splitext(os.path.split(mag_file)[1])[1]

		if(FileExt == ".eqdsk"):
			if(not self.eqdsk_file):
				self.eqdsk_file = filedialog.askopenfilename(filetypes=(("EQDSK files", "*.eqdsk"),("All files", "*.*")))
				if(not self.eqdsk_file): return

			eqdsk = load_eqdsk_file(self.eqdsk_file)

			scale_psi_twopi = False
			if(scale_psi_twopi):
				eqdsk.sibry	/= 2.*np.pi
				eqdsk.simag	/= 2.*np.pi
				eqdsk.psirz	/= 2.*np.pi

			reverse_psi = False
			if(reverse_psi): eqdsk.psirz = -eqdsk.psirz

			eqdsk.limitr	= len(self.Config.Walls[0].Rwall)
			eqdsk.rlim		= self.Config.Walls[0].Rwall
			eqdsk.zlim		= self.Config.Walls[0].Zwall
			save_eqdsk_file(mag_file, eqdsk)

			"""
		elif(FileExt == ".csv"):
			try:
				rz_b_psi = np.loadtxt(mag_file, delimiter=',',comments='"')

				self.Config.r2			= np.copy(rz_b_psi[:,0])
				self.Config.z2			= np.copy(rz_b_psi[:,1])
				self.Config.Br2			= np.copy(rz_b_psi[:,2])
				self.Config.Bphi2		= np.copy(rz_b_psi[:,3])
				self.Config.Bz2			= np.copy(rz_b_psi[:,4])
				self.Config.flux2		= np.copy(rz_b_psi[:,5])
				rz_b_psi	= 0

				nTot	= self.Config.r2.shape[0]
				R0		= self.Config.r2[0]
				for i in range(nTot):
					if R0 != self.Config.r2[i]:
						nzPoints = i											#number of points on each R
						break

				nxPoints = int(nTot/nzPoints)									#Number of points on each Z
				if nxPoints*nzPoints != nTot:
					print("\tError on Magnetic file: not matrix R,Z data")
					exit()
				print("flux2.min(),flux2.max()=",self.Config.flux2.min(),self.Config.flux2.max())
				self.Config.r2   	= self.Config.r2[:].reshape(nxPoints,nzPoints).T
				self.Config.z2		= self.Config.z2[:].reshape(nxPoints,nzPoints).T
				self.Config.Br2		= self.Config.Br2[:].reshape(nxPoints,nzPoints).T
				self.Config.Bphi2	= self.Config.Bphi2[:].reshape(nxPoints,nzPoints).T
				self.Config.Bz2		= self.Config.Bz2[:].reshape(nxPoints,nzPoints).T
				self.Config.flux2	= self.Config.flux2[:].reshape(nxPoints,nzPoints).T
				MatData = scipy.io.savemat(mag_file)

			except:
				messagebox.showerror("Read Magnetic File", "Failed to read file {:s}\n".format(mag_file))
				return
			"""
		elif((FileExt == ".mat") or (FileExt == ".m")):
#			try:
				if(False):
					MatData = {'r2D' : self.Config.r2, 'R_psi' : self.Config.r2[0,:], 'Z_psi' : self.Config.z2[:,0],  'Br' : self.Config.Br2, 'Bphi' : self.Config.Bphi2,  'Bz' : self.Config.Bz2,  'psi' : self.Config.flux2,
										'R_wall': self.Config.Walls[0].Rwall, 'Z_wall': self.Config.Walls[0].Zwall}
				else:
					MatData = {'r2D' : self.Config.r2, 'r2D' : self.Config.r2, 'z2D' : self.Config.z2,  'Br2D' : self.Config.Br2, 'Bphi2D' : self.Config.Bphi2,  'Bz2D' : self.Config.Bz2,  'flux2D' : self.Config.flux2,
										'Rwall': self.Config.Walls[0].Rwall, 'Zwall': self.Config.Walls[0].Zwall}

				scipy.io.savemat(mag_file, MatData)
#			except:
#				messagebox.showerror("Save Magnetic File", "Failed to write file {:s}\n".format(mag_file))
#				return
		else:
			messagebox.showerror("Save Magnetic File", "Unknow file type for {:s}\n".format(mag_file))

#=============================
#	export wall to csv
#=============================

	def cmd_export_wall(self):
		if(self.MeshFile):
			MeshDir  = os.path.dirname(self.MeshFile)
		else:
			MeshDir	 = "./"	
					
		WallShapeFile = filedialog.asksaveasfilename(title="Select file to export wall shape", filetypes=(("cvs files", "*.csv"), ("EQDSK files", "*.eqdsk"), 
																																																	("Mat files", "*.mat)"), ("Mat files", "*.m"), ("All files", "*.*") ), initialdir=MeshDir )


		if(not WallShapeFile): return

		FileExt = os.path.splitext(os.path.split(WallShapeFile)[1])[1]

		if(FileExt == ".csv"):
			wall = np.array([self.Config.Walls[2].Rwall,self.Config.Walls[2].Zwall]).T
			np.savetxt(WallShapeFile, wall, delimiter=',',comments='"')

		elif(FileExt == ".eqdsk"):
			if(not self.eqdsk_file):
				self.eqdsk_file = filedialog.askopenfilename(filetypes=(("EQDSK files", "*.eqdsk"),("All files", "*.*")))
				if(not self.eqdsk_file): return

			eqdsk = load_eqdsk_file(self.eqdsk_file)

			scale_psi_twopi = False
			if(scale_psi_twopi):
				eqdsk.sibry	/= 2.*np.pi
				eqdsk.simag	/= 2.*np.pi
				eqdsk.psirz	/= 2.*np.pi

			reverse_psi = False
			if(reverse_psi): eqdsk.psirz = -eqdsk.psirz

			eqdsk.limitr	= len(self.Config.Walls[0].Rwall)
			eqdsk.rlim		= self.Config.Walls[0].Rwall
			eqdsk.zlim		= self.Config.Walls[0].Zwall
			save_eqdsk_file(WallShapeFile, eqdsk)
		else:
			messagebox.showerror("Save Wall Shape File", "Unknow file type for {:s}\n".format(WallShapeFile))


#=======================================
#	tools menu
#=======================================

	def cmd_edit_wall(self):
		Changed = edit_wall(self)
		if(Changed):
			self.update_wall_shape()


	def update_wall_shape(self):
		
		Config	= self.Config
		walls_define_path(Config)

		PlasmaWallChanged = False
		EireneWallChanged = False
		for iWall in range(len(Config.Walls)):
			Wall = Config.Walls[iWall]
			if((Wall.Type == INTERNAL_PLASMA_WALL) or (Wall.Type == EXTERNAL_PLASMA_WALL)):
				Config.wall_OK	= True
				if(Wall.Changed): PlasmaWallChanged = True
			if(((Wall.Type == INTERNAL_EIRENE_WALL) or (Wall.Type == EXTERNAL_EIRENE_WALL)) and Wall.Changed):
				EireneWallChanged = True
			Wall.Changed = False

		if(PlasmaWallChanged):
			self.zone_gen.resize_mode		= True		
			self.zone_gen.update_figure()

			ChangedSegment = False
			if(Config.MagZones_OK):												#check for alligned segments
				MagZones			= Config.MagZones
				MagMegazones		= Config.MagMegazones
				MagPMegazones		= Config.MagPMegazones
				for k in range (len(MagZones)):
					if(MagZones[k].northaligned or MagZones[k].southaligned):
						if(hasattr(MagZones[k],'subSouth')):
							for n in range(2):
								MagZones[k].subSouth[n].R = np.array([])
								MagZones[k].subSouth[n].Z = np.array([])
						if(hasattr(MagZones[k],'subNorth')):
							for n in range(2):
								MagZones[k].subNorth[n].R = np.array([])
								MagZones[k].subNorth[n].Z = np.array([])

						ChangedSegment = True
						pmz = MagZones[k].pmz
						MagPMegazones[pmz].ismeshed 	= False
						MagPMegazones[pmz].isaligned	= True
						MagPMegazones[pmz].iAlignWall	= MagZones[k].iAlignWall
						MagPMegazones[pmz].meshchanged	= True
						if(hasattr(MagPMegazones[pmz], 'refpoints')):
							MagPMegazones[pmz].refpoints.R		= np.array([])
							MagPMegazones[pmz].refpoints.Z		= np.array([])
						if(hasattr(MagPMegazones[pmz], 'subrefpoints')):
							for n in range(2):
								MagPMegazones[pmz].subrefpoints[n].R	= np.array([])
								MagPMegazones[pmz].subrefpoints[n].Z	= np.array([])

			if(ChangedSegment and Config.Segments_OK):
				self.set_Segments_OK(Config, False)
				if(self.Config.MagGrid_OK):	self.mesh_gen.cmd_clear_grid()
			else:
				self.set_Eirene_OK(Config, False)

		if(EireneWallChanged and Config.Eirene_OK):
			self.Eirene.Wall.EWalls = wall_define_eirene_walls(Config)
			self.set_Eirene_OK(Config, False)

		self.zone_gen.update_figure()
		self.mesh_gen.update_figure()
		self.eirene_gen.update_figure()
		return
		
#=======================================
#	after loading fields do some actions
#=======================================

	def post_loaded_fields(self):
		
		self.Config.Frontiers		= []
		self.Config.ZoneElements	= []
		self.Config.X_points		= []
		self.Config.MagZones		= []
		self.Config.MagMegazones	= []
		self.Config.MagPMegazones	= []

		find_Xpoints(self, self.Config)
		for k in range(len(self.Config.X_points)): self.Config.X_points[k].sel = False

		self.Config.RadArroundXp = eval(self.zone_gen.RadArroundXpEntryValue.get())*1e-2
		surround_Xpoints(self, self.Config, self.CheckPlot, RadArroundXp=self.Config.RadArroundXp)

		self.Config.equ_OK			= True
		self.Config.new_equ_OK		= True
		self.set_extrapol_OK(self.Config, True)
		
		self.zone_gen.resize_mode		= True		
		self.zone_gen.update_figure()
		self.zone_gen.update_table_X_points()
		
# set status 

	def set_extrapol_OK(self, Config, OnOff, Propagate=True):
		Config.extrapol_OK		= OnOff
		if(Config.xPoints_OK and Propagate): self.set_xPoints_OK(Config, False, Propagate=Propagate)
	
	def set_xPoints_OK(self, Config, OnOff, Propagate=True):
		Config.xPoints_OK		= OnOff
		if(Config.Frontiers_OK and Propagate): self.set_Frontiers_OK(Config, False, Propagate=Propagate)
			
	def set_Frontiers_OK(self, Config, OnOff, Propagate=True):
		Config.Frontiers_OK	= OnOff
		if(Config.MagZones_OK and Propagate): self.set_MagZones_OK(Config, False, Propagate=Propagate)
			
	def set_MagZones_OK(self, Config, OnOff, Propagate=True):
		if(OnOff and (not  Config.MagZones_OK)):
			self.GenNotebook.tab(1, state= "normal")
			self.GenNotebook.select(1)

			Config.MagZones_OK	= OnOff
			self.mesh_gen.set_optmizazion_values()
			self.mesh_gen.resize_mode = True
			self.mesh_gen.update_data()

		elif((not OnOff) and Config.MagZones_OK):
			self.GenNotebook.tab(1, state= "disabled")
			self.GenNotebook.select(0)
			
		Config.MagZones_OK	= OnOff
		if(Config.Segments_OK and Propagate): self.set_Segments_OK(Config, False, Propagate=Propagate)
			
	def set_Segments_OK(self, Config, OnOff, Propagate=True):
		Config.Segments_OK	= OnOff
		if(Config.MagGrid_OK and Propagate): 
			self.set_MagGrid_OK(Config, False, Propagate=Propagate)
			if(not OnOff): self.mesh_gen.cmd_clear_grid()
			
	def set_MagGrid_OK(self, Config, OnOff, Propagate=True):
		Config.MagGrid_OK			= OnOff
		if(Propagate): 
			Config.transp_prof_OK		= False
			Config.feedback_transp_OK 	= False	
			Config.feedback_puffing_OK 	= False	
			if(Config.Mesh_OK): self.set_Mesh_OK(Config, False, Propagate=Propagate)
			
	def set_Mesh_OK(self, Config, OnOff, Propagate=True):
		if(OnOff and (not Config.Mesh_OK)):
			self.GenNotebook.tab(2, state= "normal")
			self.GenNotebook.select(2)
			self.eirene_gen.update_data()
		if((not OnOff) and Config.Mesh_OK):
			self.GenNotebook.tab(2, state= "disabled")
			self.GenNotebook.select(1)
			self.mesh_gen.update_data()
		Config.Mesh_OK				= OnOff

		if(Propagate): Config.transp_values_OK	= False
		if(Config.Eirene_OK and Propagate): self.set_Eirene_OK(Config, False, Propagate=Propagate)
	
	def set_Eirene_OK(self, Config, OnOff, Propagate=True):
		if((not OnOff) and Config.Eirene_OK):
			Config.Eirene_OK = False
			self.Eirene = types.SimpleNamespace()
		else:
			Config.Eirene_OK = OnOff
		self.eirene_gen.update_data()

	def	busy(self):
		self.config(cursor="watch")
		self.update()

	def	notbusy(self):
		self.config(cursor="arrow")
		self.update()

def main():
	
	root = Tk()
	root.geometry("+0+0")
#	root.geometry("1200x700+0+0")
#	root.geometry("1100x800+0+0")

	app = GridGen(root)

#	app.MatFile="/home/innocente/Desktop/2016-07-14_Grid_gen_511_IP/grid_ref_refined_hom-1.mat"
#	app.cmd_import_matlab_mesh()

#	app.MeshFile="/home/innocente/mesher/mesh_high-4.h5"
#	app.MeshFile="/home/innocente/mesher/mesh_jt60sa-3.h5"
#	app.MeshFile="/home/innocente/mesher/mesh_jt60sa-grid.h5"
#	app.MeshFile="/home/innocente/mesher/mesh_JT60SA_grid-3.h5"

#	app.cmd_load_eqdsk()
	if(len(sys.argv) == 2): app.cmd_load_mesh(File = sys.argv[1])

	root.mainloop()



if __name__ == '__main__':
	main()
