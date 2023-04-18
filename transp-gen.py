#!/usr/bin/env python
# -*- coding: utf-8 -*-


from tkinter 				import Tk, Text, filedialog, messagebox, BOTH, W, N, E, S, Menu, LEFT
from tkinter 				import DoubleVar, IntVar, StringVar, OptionMenu
from tkinter.ttk 			import Entry, Menubutton, Checkbutton, Frame, Button, Label, Style, LabelFrame, Notebook
from tkinter.messagebox 	import Message

import sys
import os
import types
from math 											import floor
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg				import FigureCanvasTkAgg
from matplotlib.figure								import Figure
from matplotlib.widgets 							import Cursor
try:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2TkAgg as NavigationToolbar2Tk
except:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2Tk

import numpy										as np
from numpy											import arange, sin, pi
from interfaces.generate_profile					import generate_profile
from routines.reset_config							import reset_config
from routines.intersect_contour						import intersect_2contours
from routines.utils_walls							import plot2d_walls

from mesh.find_zones_intersections					import find_zones_intersections
from mesh.generate_transport_mesh					import generate_transport_mesh
from mesh.get_rz_core_sep							import get_rz_core_sep
from mesh.compute_mesh_intersections				import compute_mesh_intersections
from files.load_soledge_mesh_file					import load_soledge_mesh_file
from files.save_soledge_mesh_file					import save_soledge_mesh_file
from files.save_profiles_as_text					import save_profiles_as_text
from files.load_profiles_as_text					import load_profiles_as_text
from files.load_feedback_file						import load_feedback_file

from routines.globals								import *

class transp_gen(Frame):

	def __init__(self, parent):
		Frame.__init__(self, parent)

		self.parent = parent
		self.initUI()

	def initUI(self):

		nPadx = 5
		nPady = 2
		self.parent.title("Transport profiles setting")
		self.pack(fill=BOTH, expand=True)

#		self.columnconfigure(1, weight=1)

		menubar = Menu(self.master)
		filemenu = Menu(menubar, tearoff=0)
		filemenu.add_command(label="Load Mesh", 	 				command=self.cmd_load_mesh)
		filemenu.add_command(label="Save Mesh", 	 				command=self.cmd_save_soledge_transport_mesh)
		filemenu.add_command(label="Save Mesh As..",			command=self.cmd_save_soledge_transport_mesh_as)
		filemenu.add_separator()
		filemenu.add_command(label="Import flux profiles",		command=self.cmd_import_flux_profiles)
		filemenu.add_command(label="Export flux profiles",		command=self.cmd_export_flux_profiles)
		filemenu.add_command(label="Import diffusion file",	command=self.cmd_import_diffusion_file)
		filemenu.add_separator()
		filemenu.add_command(label="Exit",							command=self.cmd_quit)
		menubar.add_cascade(label="File",							menu=filemenu)

		OptionsMenu = Menu(menubar, tearoff=0)

		option2dmenu = Menu(self.master)

		self.Shading2DCheckValue = IntVar()
		self.Shading2DCheckValue.set(1)
		option2dmenu.add_checkbutton(label="Shading", 				variable=self.Shading2DCheckValue)

		self.Borders2DCheckValue = IntVar()
		self.Borders2DCheckValue.set(0)
		option2dmenu.add_checkbutton(label="Zone borders",			variable=self.Borders2DCheckValue)

		self.Contour2DCheckValue = IntVar()
		self.Contour2DCheckValue.set(0)
		option2dmenu.add_checkbutton(label="Sep. contours",		variable=self.Contour2DCheckValue)
		
		OptionsMenu.add_cascade(label="2D plots",						menu=option2dmenu)

		ShowWallMenu = Menu(self.master)

		self.ShowPlasmaWallCheckValue = IntVar()
		self.ShowPlasmaWallCheckValue.set(1)
		ShowWallMenu.add_checkbutton(label="Soledge", 		variable=self.ShowPlasmaWallCheckValue)
		self.ShowPlasmaWallCheckValue.trace("w", self.cmd_ShowWallsCheckValue)

		self.ShowEireneWallCheckValue = IntVar()
		self.ShowEireneWallCheckValue.set(0)
		ShowWallMenu.add_checkbutton(label="Eirene", 		variable=self.ShowEireneWallCheckValue)
		self.ShowEireneWallCheckValue.trace("w", self.cmd_ShowWallsCheckValue)

		self.ShowExtraWallCheckValue = IntVar()
		self.ShowExtraWallCheckValue.set(0)
		ShowWallMenu.add_checkbutton(label="Extra", 		variable=self.ShowExtraWallCheckValue)
		self.ShowExtraWallCheckValue.trace("w", self.cmd_ShowWallsCheckValue)
	
		OptionsMenu.add_cascade(label="Show wall",				menu=ShowWallMenu)

		menubar.add_cascade(label="Options", 						menu=OptionsMenu)

		helpmenu = Menu(menubar, tearoff=0)
		helpmenu.add_command(label="About...",				command=self.cmd_about)
		menubar.add_cascade(label="Help", menu=helpmenu)

		self.master.config(menu=menubar)

#		Plot area

		PlotFrame = LabelFrame(self, text="Plot")
		PlotFrame.grid(row=0, column=0,  rowspan=3, pady=nPady, padx=nPadx, sticky=E+W+N)

		self.Fig = Figure(figsize=(4.8, 5.5), dpi=100)
		self.Fig.patch.set_facecolor('white')

		# a tk.DrawingArea
		canvas = FigureCanvasTkAgg(self.Fig, master=PlotFrame)
		try:
			canvas.show()										#Matplotlib  up 2.2
		except:
			canvas.draw()
#		canvas.get_tk_widget().grid(row=1, column=0,  columnspan=5,  padx=nPadx, sticky=E+W+S+N)
		canvas.get_tk_widget().grid(row=1, column=0,  columnspan=4,  padx=nPadx)

		self.Ax = self.Fig.add_subplot(111)
		self.Cb = 0														#I need this to be able to remove it later
		self.Fig.tight_layout()

		toolbar_frame = Frame(PlotFrame)
		toolbar_frame.grid(row=0, column=0, columnspan=4, padx=nPadx, sticky=E+W+S+N)
		toolbar = NavigationToolbar2Tk(canvas, toolbar_frame )
		toolbar.update()

		iRow = 2

		EqCutBut = Button(PlotFrame, text="Eq. cut", command = self.cmd_equatorial_cut, width = 10)
		EqCutBut.grid(row=iRow, column=0, pady=nPady)

		self.DrawCutButValue = StringVar()
		self.DrawCutButValue.set("Draw cut")
		DrawCutBut = Button(PlotFrame, textvariable=self.DrawCutButValue, command = self.cmd_draw_cut, width = 10)
		DrawCutBut.grid(row=iRow, column=1, pady=nPady)

		CoordLabel = Label(PlotFrame, text="Coord.:", width = 5)
		CoordLabel.grid(row=iRow, column=2,  pady=nPady, padx=nPadx)

		self.CoordEntryValue = StringVar()
		self.CoordEntry = Entry(PlotFrame, textvariable = self.CoordEntryValue, width = 14)
		self.CoordEntry.grid(row=iRow, column=3,  pady=nPady, padx=nPadx)
		self.CoordEntry.bind("<Return>", self.cmd_enter_coord_entry)
		self.CoordEntry.state(["disabled"])


		iRow += 1

		self.SelectCutBut = Button(PlotFrame, text="Select cut", command = self.cmd_select_cut, width = 10)
		self.SelectCutBut.grid(row=iRow, column=0, pady=nPady)

		self.RemoveCutBut = Button(PlotFrame, text="Remove cut", command = self.cmd_remove_cut, width = 10)
		self.RemoveCutBut.grid(row=iRow, column=1, pady=nPady)


		ShowLabel = Label(PlotFrame, text="Zone sel.", width = 7)
		ShowLabel.grid(row=iRow, column=2,  pady=nPady, padx=nPadx)

		self.ZonesOptions = ZONE_MODES
		self.ZoneModeOptMenuValue = StringVar()
		self.ZoneModeOptMenuValue.set(self.ZonesOptions[DEF_ZONE_MODE])
		self.ZoneModeOptMenu = OptionMenu(PlotFrame, self.ZoneModeOptMenuValue, *self.ZonesOptions)
		self.ZoneModeOptMenu.grid(row=iRow, column=3,  pady=nPady, padx=nPadx, sticky=E+W)

		iRow += 1

		self.GenerateBut = Button(PlotFrame, text="Generate Mesh", command = self.cmd_generate_mesh_values, width = 25)
		self.GenerateBut.grid(row=iRow, column=0, columnspan=2, pady=nPady, padx=nPadx)


		ShowLabel = Label(PlotFrame, text="Show", width = 5)
		ShowLabel.grid(row=iRow, column=2,  pady=nPady, padx=nPadx)

		self.FigureOptions = ("Borders", "Chie", "D", "Nu", "Chi","Vpinch")
		self.FigureModeOptMenuValue = StringVar()
		self.FigureModeOptMenuValue.set(self.FigureOptions[0])
		self.FigureModeOptMenu = OptionMenu(PlotFrame, self.FigureModeOptMenuValue, *self.FigureOptions)
		self.FigureModeOptMenu.grid(row=iRow, column=3,  pady=nPady, padx=nPadx, sticky=W+E)
		self.FigureModeOptMenuValue.trace("w", self.cmd_figure_option_menu)

		iRow = 0
		
#		Transport profile area

		"""
#		This manage transport profiles different for various ions (not presently possible with soledge2d 
		
		NumberOfIonsLabel = Label(self, text="Number of Ions", width = 12)
		NumberOfIonsLabel.grid(row=iRow, column=1,  pady=nPady, padx=nPadx)
		
		self.NumberOfIonsOptions = ["1", "2", "3","4","5"]
		self.NumberOfIonsOptMenuValue = StringVar()
		self.NumberOfIonsOptMenuValue.set(self.NumberOfIonsOptions[0])
		self.NumberOfIonsOptMenu = OptionMenu(self, self.NumberOfIonsOptMenuValue, self.NumberOfIonsOptions[0], self.NumberOfIonsOptions[1], self.NumberOfIonsOptions[2], 
																						self.NumberOfIonsOptions[3], self.NumberOfIonsOptions[4])
		self.NumberOfIonsOptMenu.grid(row=0, column=2,  pady=nPady, padx=nPadx,sticky=W)

		EmptysLabel = Label(self, text=" ", width = 35)
		EmptysLabel.grid(row=iRow, column=3,  pady=nPady, padx=nPadx)
		
		iRow += 1
		"""

#		Option area
		
		OptionFrame = LabelFrame(self, text="Option parameters")
		OptionFrame.grid(row=iRow, column=1,  pady=nPady, padx=nPadx, sticky=E+W+N)

		BallooningModeLabel = Label(OptionFrame, text="Ballooning mode", justify = LEFT, width=15)
		BallooningModeLabel.grid(row=0, column=0,  pady=nPady, padx=nPadx, sticky=W)

		self.BallooningModeOptMenuValues = ("No", "1/Bt", "1/B","1/Bt^2","1/B^2")
		self.BallooningModeOptMenuValue = StringVar()
		self.BallooningModeOptMenuValue.set(self.BallooningModeOptMenuValues[DEF_BALLOONING_MODE])
		BallooningModeOptMenu = OptionMenu(OptionFrame, self.BallooningModeOptMenuValue, *self.BallooningModeOptMenuValues)
		BallooningModeOptMenu.grid(row=0, column=1,  pady=nPady, padx=nPadx, sticky=W)

		iRow += 1

#		BALLOONING_TAB_NAMES	= ["Chi_e ", "D    ", "Nu   ", "Chi_i ", "Vinch"] 
#		BALLOONING_NAMES		= ["Chie", "D", "Nu", "Chi", "Vinch"] 					#Balloning names in mesh file
		
		self.TabProfileNames	= BALLOONING_TAB_NAMES
		self.BallooningNames	= BALLOONING_NAMES
		self.FluxProfileWidgets  = []
		self.FluxProfNotebook = Notebook(self)
		for k in range(len(self.BallooningNames)):
			self.FluxProfileWidgets.append(generate_profile(self.FluxProfNotebook, Name=self.BallooningNames[k], Title="Flux profile", FigSize=(7.5, 2.)))
			self.FluxProfNotebook.add(self.FluxProfileWidgets[k].Frame, text=self.TabProfileNames[k])

		self.FluxProfNotebook.grid(row=iRow, column=1, columnspan=1,  pady=nPady, padx=nPadx, sticky=E+W+N+S)

		iRow += 1

		self.ThetaProfileWidgets  = []
		self.ThetaProfNotebook = Notebook(self)
		for k in range(len(self.BallooningNames)):
			self.ThetaProfileWidgets.append(generate_profile(self.ThetaProfNotebook, Name=self.BallooningNames[k], Title="Poloidal Profile", FigSize=(7.5, 2.)))
			self.ThetaProfNotebook.add(self.ThetaProfileWidgets[k].Frame, text=self.TabProfileNames[k])

		self.ThetaProfNotebook.grid(row=iRow, column=1, columnspan=1,  pady=nPady, padx=nPadx, sticky=E+W+N+S)


		for k in range(len(self.FluxProfileWidgets)):
			flux_profile_copy_from  = []
			theta_profile_copy_from = []
			tab_profile_copy_from   = []
			for i in range(len(self.FluxProfileWidgets)):
				if(i != k):
					flux_profile_copy_from.append(self.FluxProfileWidgets[i])
					theta_profile_copy_from.append(self.ThetaProfileWidgets[i])
					tab_profile_copy_from.append(self.TabProfileNames[i])

			self.FluxProfileWidgets[k].copy_from_set(flux_profile_copy_from, tab_profile_copy_from)
			self.ThetaProfileWidgets[k].copy_from_set(theta_profile_copy_from, tab_profile_copy_from)

		self.resize_mode			= False
		self.sel_mode_d				= False
		self.cut_sel				= 0
		self.InDrawCut				= False
		self.changed				= False
		
		self.Config					= reset_config()
		self.TranspCuts				= []
		self.CutLine				= np.empty((2,2),dtype='f8')

		self.parent.protocol("WM_DELETE_WINDOW", self.cmd_quit)
		
	def cmd_about(self):
		print("To do cmd_about")
		
#=================
#	Load mesh data
#=================

	def cmd_load_mesh(self, MeshFile=None):
#	load mesh previously saved

		if(MeshFile is None):
			MeshFile = filedialog.askopenfilename(filetypes=(("Mesh files", "*.h5"), ("All files", "*.*") ))
			if(not MeshFile): return

#		try:
		if(True):
			self.busy()

			self.changed = False
			self.Config  = load_soledge_mesh_file(MeshFile)
			if(self.Config.transp_prof_OK):
				self.TranspCuts = self.Config.TranspCuts
				for i in range(len(self.Config.TranspCuts)):
					self.Config.TranspCuts[i].Flux.Lengths = compute_mesh_intersections(self.Config, self.Config.TranspCuts[i].Flux)			#Compute conversion lengths
				
			else:
				self.TranspCuts = []
				self.Config.transp_values_OK = False
				
			for k in range(len(self.Config.MagZones)):
				self.Config.MagZones[k].IsTranpProf		= self.Config.transp_prof_OK
				self.Config.MagZones[k].IsTranspValues	= self.Config.transp_values_OK
				for sz_k in self.Config.MagZones[k].list: self.Config.Zones[sz_k].IsTranspValues = self.Config.transp_values_OK
		
			self.MeshFile		= MeshFile
			self.MeshDir  		= os.path.dirname(self.MeshFile)
			self.sel_mode_d		= False
			self.cut_sel		= -1
			self.set_profile_data()
				
			self.resize_mode 	= True
			self.update_figure()

			self.notbusy()

#		except:
#			messagebox.showerror("Open Source File", "Failed to read mesh file:{:s}\n".format(MeshFile))

				
	def cmd_load_ref_mesh(self):
		print("To do cmd_load_ref_mesh")


	def cmd_save_soledge_transport_mesh(self):
		if(not self.MeshFile): return
		self.set_transp_ok()
		if(not self.Config.transp_prof_OK):
			messagebox.showwarning("Saving file", "Mesh file can be saved only when all zones are ok")
			return

		try:
			self.busy()
			if(len(self.TranspCuts) > 0):
				self.Config.transp_prof_OK = True
				self.Config.TranspCuts = self.TranspCuts
			else:
				self.Config.transp_prof_OK = False	

			save_soledge_mesh_file(self.MeshFile, self.Config)

			self.changed = False
			self.notbusy()
		except:
			self.notbusy()
			messagebox.showerror("Save Mesh file", "Failed to write mesh file:{:s}\n".format(self.MeshFile))

	def cmd_save_soledge_transport_mesh_as(self):
		self.set_transp_ok()
		if(not self.Config.transp_prof_OK):
			messagebox.showwarning("Saving file", "Mesh file can be saved only when all zones are ok")
			return
		
		MeshFile = filedialog.asksaveasfilename(filetypes=(("Mesh files", "*.h5"), ("All files", "*.*") ))

		if(not MeshFile): return
#		try:
		if(True):
			self.busy()
			if(len(self.TranspCuts) > 0):
				self.Config.transp_prof_OK = True
				self.Config.TranspCuts = self.TranspCuts
			else:
				self.Config.transp_prof_OK = False	

			save_soledge_mesh_file(MeshFile, self.Config)
			self.MeshFile = MeshFile
			self.MeshDir  = os.path.dirname(self.MeshFile)
			self.changed  = False
			self.notbusy()
#		except:
#			messagebox.showerror("Save Mesh file", "Failed to write mesh file:{:s}\n".format(self.MeshFile))

	def cmd_import_flux_profiles(self):
		if(not self.MeshFile): return
		
		if(self.cut_sel < 0):
			messagebox.showwarning("Saving file", "Profiles can be loaded only when at least one cut is selected")
			return
		ProfFile = filedialog.askopenfilename(initialdir=self.MeshDir, filetypes=(("Comma separated values", "*.csv"),  \
																					("Tab separated values",    "*.tsv"),  \
																					("Space separated values",  "*.*")))
		if(ProfFile):
#		try:
			xName = LENGTH_TYPES[0]
			self.get_profile_data()
			load_profiles_as_text(ProfFile, self.TranspCuts[self.cut_sel].Flux, xName)
			self.set_profile_data()
#		except:
#			messagebox.showerror("Save Mesh file", "Failed to import profiles file:{:s}\n".format(ProfFile))

	def cmd_quit(self):
		if(self.changed):
			answer=messagebox.askquestion("Mesh modified", "Save changes to mesh")
			if(answer=="yes"): self.cmd_save_soledge_transport_mesh()
		self.quit()

	def cmd_export_flux_profiles(self):
		if(not self.MeshFile): return
		
		if(self.cut_sel < 0):
			messagebox.showwarning("Saving file", "Profiles can be saved only when at least one cut is selected")
			return

		ProfFile = filedialog.asksaveasfilename(initialdir=self.MeshDir, filetypes=(("Comma separated values", "*.csv"),  \
																					("Tab separated values",    "*.tsv"),  \
																					("Space separated values",  "*.*")))
		if(ProfFile):
#		try:
			self.get_profile_data()
			xName = LENGTH_TYPES[0]
			save_profiles_as_text(ProfFile, self.TranspCuts[self.cut_sel].Flux, xName)
#		except:
#			messagebox.showerror("Save Mesh file", "Failed to export profiles file:{:s}\n".format(ProfFile))

	def cmd_import_diffusion_file(self):
		if(not self.MeshFile): return
		
		if(self.cut_sel < 0):
			messagebox.showwarning("Loading file", "SOLEDGE diffusion can be imported when at least one cut is selected")
			return

		DiffFile = filedialog.askopenfilename(initialdir=self.MeshDir)

		if(DiffFile):
			load_feedback_file(DiffFile, 1, self.Config, self.TranspCuts[self.cut_sel])

			self.set_profile_data()
			self.update_figure()

	def set_transp_ok(self):
		self.Config.transp_prof_OK		= True
		self.Config.transp_values_OK	= True
		for k in range(len(self.Config.MagZones)):
			self.Config.transp_prof_OK = self.Config.transp_prof_OK and self.Config.MagZones[k].IsTranpProf
			self.Config.transp_values_OK = self.Config.transp_values_OK and self.Config.MagZones[k].IsTranspValues


#######################################
#	Option  menu show walls
#######################################
		
	def cmd_ShowWallsCheckValue(self, a, b, c):
		self.update_figure()
		return

#######################################
#	Draw a cut to set transport profile
#######################################

	def cmd_equatorial_cut(self):

		Rcore, Zcore, CoreMegazone = get_rz_core_sep(self.Config, core_and_sep=False)

		iMax = np.argmax(Rcore)
		self.CutLine[0,1] = Zcore[iMax]
		self.CutLine[1,1] = Zcore[iMax]

		Ri, Zi, is1, is2  = intersect_2contours(Rcore, Zcore, np.array([Rcore.min()-1., Rcore.max()+1.]),  np.array([self.CutLine[0,1],self.CutLine[1,1]]))

		if(len(Ri) != 2):
			print("Error in finding core center")
			return

		self.CutLine[0,0] = 0.5*(Ri[0] + Ri[1])
		self.CutLine[1,0] = 2.*self.CutLine[0,0]

		self.nCutPoints   = 2
		self.post_defined_cut()

	def cmd_draw_cut(self):
		if(not self.InDrawCut):
			self.DrawCutButValue.set("Cancel")
			self.CoordEntry.state(["!disabled"])

			self.nCutPoints = 0
			self.cid_mouse_button_press_limit = self.Fig.canvas.mpl_connect('button_press_event', self.cmd_button_press_add_point)
			self.cursor = Cursor(self.Ax, useblit=False, color='red', linewidth=1)
			self.InDrawCut = True
		else:
			self.InDrawCut = False
			self.DrawCutButValue.set("Draw cut")
			self.CoordEntry.state(["disabled"])
			if(self.cid_mouse_button_press_limit != 0):
				self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press_limit)
				self.cid_mouse_button_press_limit = 0
				self.cursor=0

		self.update_figure()

	def cmd_button_press_add_point(self, event):
		self.CutLine[self.nCutPoints,0] = event.xdata
		self.CutLine[self.nCutPoints,1] = event.ydata
		self.nCutPoints                += 1
		self.Ax.plot(self.CutLine[:self.nCutPoints,0],self.CutLine[:self.nCutPoints,1],'r.-')
		self.Fig.canvas.draw()
		if(self.nCutPoints == 2):	
			self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press_limit)
			self.cid_mouse_button_press_limit = 0
			self.cursor=0
			self.InDrawCut = False
			self.DrawCutButValue.set("Draw cut")
			self.CoordEntry.state(["disabled"])
			
			self.post_defined_cut()

	def post_defined_cut(self):

			MagZones					= self.Config.MagZones
			TranspCuts				= types.SimpleNamespace()
			SkipZones				= [Zone.IsTranpProf for Zone in MagZones]
			TranspCuts.Flux			= find_zones_intersections(self.Config, self.CutLine, SkipZones)
			TranspCuts.Flux.Lengths = compute_mesh_intersections(self.Config, TranspCuts.Flux)				#Compute lengths for conversion between L, X, Rho_v, Rho_F
			TranspCuts.sel			= False


			if(len(TranspCuts.Flux.nz) != 0):
				TranspCuts.BallooningMode = self.BallooningModeOptMenuValues.index(self.BallooningModeOptMenuValue.get())
				TranspCuts.ZoneMode		  = self.ZonesOptions.index(self.ZoneModeOptMenuValue.get())


				MagMegazones	= self.Config.MagMegazones
				
#				Build list of zone affected by the cut
#				Flux direction

				nZones				  = 0
				TranspCuts.Flux.MagZones = []
				if(TranspCuts.ZoneMode == ZONE_MODE_SINGLE):				#Only cutted zone
					for k in range(len(TranspCuts.Flux.nz)):
						TranspCuts.Flux.MagZones.append(np.array([TranspCuts.Flux.nz[k]]))
						MagZones[TranspCuts.Flux.nz[k]].IsTranpProf = True
						nZones += 1


				elif(TranspCuts.ZoneMode == ZONE_MODE_CONNECTED):
					for k in range(len(TranspCuts.Flux.nz)):

						iMZone = MagZones[TranspCuts.Flux.nz[k]].mz							#add zones magnetically connected
						zones_list   = np.empty((0),dtype='i4')
						for iZone in MagMegazones[iMZone].list:
							if(not MagZones[iZone].IsTranpProf):
								zones_list = np.append(zones_list, iZone)
								MagZones[iZone].IsTranpProf = True
								nZones += 1
						
						TranspCuts.Flux.MagZones.append(zones_list)	
				elif(TranspCuts.ZoneMode == ZONE_MODE_SAME_FLUX):
					print("To be implmented flux connected")

				elif(TranspCuts.ZoneMode == ZONE_MODE_SAME_X):								#add zones witht the same x
					TranspCuts.Flux.MagZones = []

					for k in TranspCuts.Flux.nz:
						xMin = MagZones[k].x[0, 0]-0.05
						xMax = MagZones[k].x[-1,0]+0.05
						zones_list   = np.empty((0),dtype='i4')
						for iZone in range(len(MagZones)):
							if((not MagZones[iZone].IsTranpProf) and (MagZones[iZone].x[0,0] >= xMin) and (MagZones[iZone].x[-1,0] <= xMax)):
								zones_list = np.append(zones_list, iZone)
								MagZones[iZone].IsTranpProf 	= True
								nZones += 1
									
						TranspCuts.Flux.MagZones.append(zones_list)	

				TranspCuts.Flux.Profiles = []
				for k in range(len(self.BallooningNames)):
					TranspCuts.Flux.Profiles.append(types.SimpleNamespace())
					TranspCuts.Flux.Profiles[-1].Name	 = self.BallooningNames[k]
					TranspCuts.Flux.Profiles[-1].xValues = np.array([TranspCuts.Flux.d12[0,0],TranspCuts.Flux.d12[-1,1]])
					TranspCuts.Flux.Profiles[-1].Values	 = np.ones((2), dtype='f8')

#				Theta direction

				TranspCuts.Theta		= types.SimpleNamespace()
				TranspCuts.Theta.d12	= np.empty((nZones,2), dtype='f8')
				
				dZero = 0.
				if(TranspCuts.ZoneMode == ZONE_MODE_SINGLE):						#single zone
					TranspCuts.Theta.nz		= np.copy(TranspCuts.Flux.nz)
					
					for k in range(len(TranspCuts.Theta.nz)):
						iZone = TranspCuts.Theta.nz[k]
						iCenter = floor(MagZones[iZone].gridR.shape[0]/2)
						R		= MagZones[iZone].gridR[iCenter,:]
						Z		= MagZones[iZone].gridZ[iCenter,:]
						TranspCuts.Theta.d12[k,0]  = dZero
						dZero					  += np.sum(np.sqrt((R[1:]-R[:-1])**2+(Z[1:]-Z[:-1])**2))
						TranspCuts.Theta.d12[k,1]  = dZero
						dZero 					  += 0.5								#Add interval between MagMegazones

				elif((TranspCuts.ZoneMode == ZONE_MODE_CONNECTED) or (TranspCuts.ZoneMode == ZONE_MODE_SAME_X)):
					TranspCuts.Theta.nz	= np.empty((nZones),   dtype='i4')
					StartZones 			= TranspCuts.Flux.nz
					for k in range(len(TranspCuts.Flux.nz)):
						StartZones = np.append(StartZones, TranspCuts.Flux.MagZones[k])

					if(TranspCuts.ZoneMode == ZONE_MODE_SAME_X):
						nSearch = len(StartZones)
					else:
						nSearch = len(TranspCuts.Flux.nz)
											
					kZone = 0
					for i in range(nSearch):
						iZone = StartZones[i]
						if(len(np.where(TranspCuts.Theta.nz[:kZone] == iZone)[0]) == 0):
							jMZone = MagZones[iZone].mz										#add zones magnetically connected
							for jZoneWest in MagMegazones[jMZone].list:
								if(MagZones[jZoneWest].Neighbour.west == EW_BORDER_NEIGHBOUR): break			#Search West end
									
							jZone	   = jZoneWest
							nInZones   = len(MagMegazones[jMZone].list)						#Manage circular MagMegazones
							jInZones   = 0
							FoundZone  = False					
							while ((jZone >= 0) and (jInZones < nInZones)):
								jInZones += 1
								iCenter	  = floor(MagZones[jZone].gridR.shape[0]/2)
								R		  = MagZones[jZone].gridR[iCenter,:]
								Z	 	  = MagZones[jZone].gridZ[iCenter,:]
								lZone	  = np.sum(np.sqrt((R[1:]-R[:-1])**2+(Z[1:]-Z[:-1])**2))
								ii		  = np.where(StartZones == jZone)
								if(len(ii[0]) > 0):
									TranspCuts.Theta.nz[kZone] 	   = jZone
									TranspCuts.Theta.d12[kZone,0]  = dZero
									dZero						  += lZone
									TranspCuts.Theta.d12[kZone,1]  = dZero
									jZone						   = MagZones[jZone].Neighbour.east
									FoundZone					   = True
									kZone						  += 1
								elif (FoundZone):
									dZero						  += lZone
										
							dZero += 0.5							#Add interval between MagMegazones	
								
					TranspCuts.Theta.Lengths = np.empty((0,0),dtype='f8')					#Empty array for compatibility with Flux

				elif(TranspCuts.ZoneMode == ZONE_MODE_SAME_FLUX):
					print("To be implmented flux connected")
				
					
				ii		= np.where(TranspCuts.Theta.d12[1:,0]-TranspCuts.Theta.d12[:-1,1] > 0.); ii = ii[0]
				if(len(ii) > 0):
					dKnots	= np.concatenate((np.array([TranspCuts.Theta.d12[0,0]]), TranspCuts.Theta.d12[ii,1], TranspCuts.Theta.d12[ii+1,0], np.array([TranspCuts.Theta.d12[-1,1]]))) 
				else:
					dKnots	= np.array([TranspCuts.Theta.d12[0,0], TranspCuts.Theta.d12[-1,1]])
	
				TranspCuts.Theta.Profiles = []
				for k in range(len(self.BallooningNames)):
					TranspCuts.Theta.Profiles.append(types.SimpleNamespace())
					TranspCuts.Theta.Profiles[-1].Name	  = self.BallooningNames[k]
					TranspCuts.Theta.Profiles[-1].xValues = np.copy(dKnots)
					TranspCuts.Theta.Profiles[-1].Values  = np.ones(dKnots.shape, dtype='f8')
						
				if(self.sel_mode_d): self.get_profile_data() 

				self.TranspCuts.append(TranspCuts)
				self.sel_mode_d		= True
				self.cut_sel		= len(self.TranspCuts) - 1
				self.set_profile_data()

			self.update_figure()
			
########################################
#	Define cut point based on coordinate
########################################

	def cmd_enter_coord_entry(self, event):
		Coords = self.CoordEntryValue.get()
		Coords = "".join(Coords.split())
		
		IsRelative = False
		if(Coords[0] == "@"):
			if(self.nCutPoints > 0):	IsRelative = True
			Coords = Coords[1:]
			
		iAngle = Coords.find("<")
		IsAngular = False
		if(iAngle >=0):
			Coords = Coords.replace("<",",")
			IsAngular = True
	
		try:
			Values = eval(Coords)
			if(len(Values) != 2):return
		except:
			return
		
		value =	types.SimpleNamespace()
		if(IsAngular):
			value.xdata = Values[0]*np.cos(Values[1]*np.pi/180.)
			value.ydata = Values[0]*np.sin(Values[1]*np.pi/180.)
		else:
			value.xdata = Values[0]
			value.ydata = Values[1]
			
		if(IsRelative):
			value.xdata += self.CutLine[0,0]
			value.ydata += self.CutLine[0,1]
		
		self.cmd_button_press_add_point(value)

##########################
#	Select an existing cut
##########################

	def cmd_select_cut(self):
		self.cid_mouse_button_press = self.Fig.canvas.mpl_connect('button_press_event', self.cmd_mouse_button_press_cut_selected)
		self.cursor = Cursor(self.Ax, useblit=False, color='red', linewidth=1)

	def cmd_mouse_button_press_cut_selected(self, event):
		self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press)
		self.cid_mouse_button_press = 0
		self.cursor=0

#		find closest cut

		if(self.sel_mode_d): self.get_profile_data() 			

		TranspCuts	=  self.TranspCuts
		Rpt		 = event.xdata
		Zpt		 = event.ydata
		dmin     = 1e10
		for k in range(len(TranspCuts)):
			vx	 =TranspCuts[k].Flux.r12[:,1] - TranspCuts[k].Flux.r12[:,0]
			vy	= TranspCuts[k].Flux.z12[:,1] - TranspCuts[k].Flux.z12[:,0]
			vx1	= Rpt - TranspCuts[k].Flux.r12[:,0]
			vy1	= Zpt - TranspCuts[k].Flux.z12[:,0]
			vx2	= Rpt - TranspCuts[k].Flux.r12[:,1]
			vy2	= Zpt - TranspCuts[k].Flux.z12[:,1]
			d	= np.zeros(len(vx))
			p	= (vx*vx1+vy*vy1)/(vx**2+vy**2)
			d	= np.where((p > 0) & (p<1), \
								np.abs(vx*vy1-vy*vx1)/np.sqrt(vx**2+vy**2), \
								np.where(np.sqrt(vx1**2+vy1**2) < np.sqrt(vx2**2+vy2**2), \
										 np.sqrt(vx1**2+vy1**2),  np.sqrt(vx2**2+vy2**2)))
			d = np.min(d)
			if(d < dmin):
				dmin				= d
				self.cut_sel		= k
				
		self.sel_mode_d		= True
		self.set_profile_data()
				
		for i in range(len(self.TranspCuts[self.cut_sel].Flux.nz)):
			for iZone in self.TranspCuts[self.cut_sel].Flux.MagZones[i]:
				self.Config.MagZones[iZone].IsTranspValues = False
				for sz_k in self.Config.MagZones[iZone].list: self.Config.Zones[sz_k].IsTranspValues = False
				
		self.update_figure()

###########################
#	Remove and existing cut
###########################

	def cmd_remove_cut(self):
		k = self.cut_sel
		for i in range(len(self.TranspCuts[k].Flux.nz)):
			for iZone in self.TranspCuts[k].Flux.MagZones[i]:
				self.Config.MagZones[iZone].IsTranpProf	 = False
				self.Config.MagZones[iZone].IsTranspValues = False
				for sz_k in self.Config.MagZones[iZone].list: self.Config.Zones[sz_k].IsTranspValues = False

		del self.TranspCuts[k]
		self.sel_mode_d		= False
		self.cut_sel		= -1
		self.set_profile_data()
		self.update_figure()

	def cmd_generate_mesh_values(self):
		self.get_profile_data()
		generate_transport_mesh(self, self.Config, self.TranspCuts)
		self.update_figure()
									
	def cmd_figure_option_menu(self,a,b,c):
		self.update_figure()

############################################
#	Set & get profile data on profile widget
############################################

	def set_profile_data(self):
		if(self.cut_sel >= 0):
			self.changed = True
			k = self.cut_sel
			self.BallooningModeOptMenuValue.set(self.BallooningModeOptMenuValues[self.TranspCuts[k].BallooningMode])
			self.ZoneModeOptMenuValue.set(self.ZonesOptions[self.TranspCuts[k].ZoneMode])

			for i in range(len(self.FluxProfileWidgets)):
				self.FluxProfileWidgets[i].set_data(self.TranspCuts[k].Flux,   self.TranspCuts[k].Flux.Profiles[i])
				self.ThetaProfileWidgets[i].set_data(self.TranspCuts[k].Theta, self.TranspCuts[k].Theta.Profiles[i])
		else:
			self.BallooningModeOptMenuValue.set(self.BallooningModeOptMenuValues[DEF_BALLOONING_MODE])
			self.ZoneModeOptMenuValue.set(self.ZonesOptions[DEF_ZONE_MODE])
			for i in range(len(self.FluxProfileWidgets)):
				self.FluxProfileWidgets[i].reset_data()
				self.ThetaProfileWidgets[i].reset_data()
			
	def get_profile_data(self):
		if(self.cut_sel >= 0):
			k = self.cut_sel
			self.TranspCuts[k].BallooningMode = self.BallooningModeOptMenuValues.index(self.BallooningModeOptMenuValue.get())
			self.TranspCuts[k].ZoneMode		  = self.ZonesOptions.index(self.ZoneModeOptMenuValue.get())
			for i in range(len(self.FluxProfileWidgets)):
				self.TranspCuts[k].Flux.Profiles[i],  self.TranspCuts[k].Flux.d12   = self.FluxProfileWidgets[i].get_data()
				self.TranspCuts[k].Theta.Profiles[i], self.TranspCuts[k].Theta.d12 = self.ThetaProfileWidgets[i].get_data()

			for i in range(len(self.TranspCuts[k].Flux.nz)):
				for iZone in self.TranspCuts[k].Flux.MagZones[i]:
					self.Config.MagZones[iZone].IsTranpProf	 = True
					self.Config.MagZones[iZone].IsTranspValues = False
					for sz_k in self.Config.MagZones[iZone].list: self.Config.Zones[sz_k].IsTranspValues = False


	def	busy(self):
		self.config(cursor="watch")
		self.update()

	def	notbusy(self):
		self.config(cursor="arrow")
		self.update()

#==============
#	Main figure 
#==============

#update main figure

	def update_figure(self):

		Config			= self.Config
		MagZones		= self.Config.MagZones
		Zones			= self.Config.Zones
		X_points		= self.Config.X_points
		TranspCuts		= self.TranspCuts
		Fig				= self.Fig
		Ax				= self.Ax

		FigureMode = self.FigureOptions.index(self.FigureModeOptMenuValue.get())	
		
		if(FigureMode > 0):
			iBall	  = FigureMode - 1
			nZones	  = len(Zones)
			vMin	  = np.empty(nZones, dtype='f8') 
			vMax	  = np.empty(nZones, dtype='f8')
			nOkZones  = 0
			for k in range(nZones):
				if(Zones[k].IsTranspValues):
					vMin[nOkZones] = np.min(Zones[k].Ballooning[iBall])
					vMax[nOkZones] = np.max(Zones[k].Ballooning[iBall])
					nOkZones	   += 1
			if(nOkZones == 0): return

		if(self.Cb != 0):
			Fig.clf()
			self.Ax = Ax = Fig.add_subplot(111)
			self.Cb = 0
		else:
			Ax.clear()		
		
		Ax.tick_params(axis='both', which='major', labelsize=8)

		if(FigureMode == 0):
			Xlims = Ax.get_xlim()
			Ylims = Ax.get_ylim()
	
			if(Config.Mesh_OK):
				for k in range(len(MagZones)):
					EdgeLineType	='b-'
					EdgeLineWidth	= 1
					Ax.plot(MagZones[k].gridR[0,:],  MagZones[k].gridZ[0,:],  EdgeLineType, linewidth=EdgeLineWidth)		#south
					Ax.plot(MagZones[k].gridR[:,-1], MagZones[k].gridZ[:,-1], EdgeLineType, linewidth=EdgeLineWidth)		#East
					Ax.plot(MagZones[k].gridR[-1,:], MagZones[k].gridZ[-1,:], EdgeLineType, linewidth=EdgeLineWidth)		#North
					Ax.plot(MagZones[k].gridR[:,0],  MagZones[k].gridZ[:,0],  EdgeLineType, linewidth=EdgeLineWidth)		#West
	
					iCenter = floor(MagZones[k].gridR.shape[1]/2)
					R=(MagZones[k].gridR[0,iCenter]+MagZones[k].gridR[-1,iCenter])*0.5
					Z=(MagZones[k].gridZ[0,iCenter]+MagZones[k].gridZ[-1,iCenter])*0.5
					if(MagZones[k].IsTranpProf):
						Ax.text(R,Z,"{:d}".format(k+1), color='green',  bbox=dict(facecolor='green', alpha=0.5))
					else:
						Ax.text(R,Z,"{:d}".format(k+1), color='red', bbox=dict(facecolor='red', alpha=0.5))
	
			if(Config.equ_OK):
				for k in range(len(X_points)):
					Ax.contour(Config.r2D, Config.z2D, Config.flux2D, colors="r", linestyles='-', linewidths=2, levels=[X_points[k].psi])
#					Ax.text(X_points[k].R,X_points[k].Z,"{:d}".format(X_points[k].index),size=18)
	
	
			if(len(TranspCuts) > 0):
				for k in range(len(TranspCuts)):
					for i in range(len(TranspCuts[k].Flux.nz)):
						if(self.sel_mode_d and (self.cut_sel == k)):
							Ax.plot(TranspCuts[k].Flux.r12[i,:], TranspCuts[k].Flux.z12[i,:],'c.--', linewidth=2)
						else: 
							Ax.plot(TranspCuts[k].Flux.r12[i,:], TranspCuts[k].Flux.z12[i,:],'b.--', linewidth=2)
		
#						Ax.plot(frontiers[k].P1[:,0], frontiers[k].P1[:,1],'mv')
						
			if((not ((Xlims[0] == 0.) and (Xlims[1] == 1.))) and (not self.resize_mode)):
				Ax.set_xlim(Xlims[0], Xlims[1])
				Ax.set_ylim(Ylims[0], Ylims[1])
			else:
				self.resize_mode = False
			if(Config.wall_OK): plot2d_walls(Ax, Config.Walls, self.ShowPlasmaWallCheckValue.get(), self.ShowEireneWallCheckValue.get(), self.ShowExtraWallCheckValue.get())
		else:
			vMin  = np.min(vMin[:nOkZones])
			vMax  = np.max(vMax[:nOkZones])
			if(vMin == vMax):
				vMin = 0.9*vMin
				vMax = 1.1*vMax
				
			for k in range(nZones):
				if(self.Shading2DCheckValue.get() == 1):
					shade = 'gouraud'
				else:
					shade = 'flat'

				if(Zones[k].IsTranspValues):	
					Ball = np.append(Zones[k].Ballooning[iBall], Zones[k].Ballooning[iBall][:,-1].reshape(Zones[k].Ballooning[iBall].shape[0],1), axis=1)
					Ball = np.append(Ball, Ball[-1,:].reshape(1,Ball.shape[1]), axis=0)

					imBall = Ax.pcolormesh(Zones[k].gridR,  Zones[k].gridZ, Ball,  shading=shade, vmin=vMin, vmax=vMax)

				if(self.Borders2DCheckValue.get() == 1):
					for k in range(len(Zones)):
						EdgeLineType	='k--'
						EdgeLineWidth	= 1
						Ax.plot(Zones[k].gridR[0,:],  Zones[k].gridZ[0,:],  EdgeLineType, linewidth=EdgeLineWidth)		#south
						Ax.plot(Zones[k].gridR[:,-1], Zones[k].gridZ[:,-1], EdgeLineType, linewidth=EdgeLineWidth)		#East
						Ax.plot(Zones[k].gridR[-1,:], Zones[k].gridZ[-1,:], EdgeLineType, linewidth=EdgeLineWidth)		#North
						Ax.plot(Zones[k].gridR[:,0],  Zones[k].gridZ[:,0],  EdgeLineType, linewidth=EdgeLineWidth)		#West


				if(self.Contour2DCheckValue.get() == 1):
					for k in range(len(X_points)):
						Ax.contour(Config.r2D, Config.z2D, Config.flux2D, colors="r", linestyles='-', linewidths=1, levels=[X_points[k].psi])
					
			self.Cb = Fig.colorbar(imBall, ax=Ax)
			self.Cb.ax.tick_params(labelsize=8)
			self.Cb.set_label(self.FigureOptions[FigureMode], fontsize=10)
			
			if(Config.wall_OK): plot2d_walls(Ax, Config.Walls, self.ShowPlasmaWallCheckValue.get(), self.ShowEireneWallCheckValue.get(), self.ShowExtraWallCheckValue.get())

			self.resize_mode = True

		Ax.set_aspect(1.)
		Fig.canvas.draw()

		
def main():

	root = Tk()
	root.geometry("+0+0")
#	root.geometry("1102x692+0+0")
	app = transp_gen(root)

	if(len(sys.argv) == 2): app.cmd_load_mesh(MeshFile = sys.argv[1])

	root.mainloop()

if __name__ == '__main__':
	main()

