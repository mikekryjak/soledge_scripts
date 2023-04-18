#!/usr/bin/env python
# -*- coding: utf-8 -*-


from tkinter 				import Tk, Text, filedialog, messagebox, BOTH, W, N, E, S, Menu, HORIZONTAL, LEFT, CENTER, RIGHT
from tkinter 				import DoubleVar, IntVar, StringVar, OptionMenu
from tkinter.ttk 			import Entry, Menubutton, Checkbutton, Frame, Button, Label, Style, LabelFrame, Notebook
from tkinter.messagebox 	import Message

import sys
import os
import types
import matplotlib
from math 											import floor, sqrt
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg				import FigureCanvasTkAgg
from matplotlib.figure								import Figure
from matplotlib.widgets 							import Cursor
try:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2TkAgg as NavigationToolbar2Tk
except:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2Tk

import numpy										as np
import scipy.interpolate							as interp

from interfaces.generate_profile					import generate_profile
from routines.globals								import FEEDBACK_TAB_NAMES, FEEDBACK_NAMES, LENGTH_TYPES, CORE_NEIGHBOUR, DEF_BALLOONING_MODE, DEBUG
from routines.reset_config							import reset_config
from routines.intersect_contour						import intersect_contour
from routines.exp_data_routines						import get_exp_data_point
from routines.utils_walls							import plot2d_walls

from mesh.get_rz_core_sep							import get_rz_core_sep
from mesh.find_zones_intersections					import find_zones_intersections
from mesh.generate_feedback_profiles				import generate_feedback_profiles
from mesh.generate_feedback_mesh					import generate_feedback_mesh
from mesh.compute_mesh_intersections				import compute_mesh_intersections
from files.load_soledge_mesh_file					import load_soledge_mesh_file
from files.save_soledge_mesh_file					import save_soledge_mesh_file
from files.save_feedback_plasma_files				import save_feedback_plasma_files
from files.save_feedback_transp_files				import save_feedback_transp_files
from files.save_feedback_puffing_file				import save_feedback_puffing_file
from files.save_profiles_as_text					import save_profiles_as_text
from files.load_profiles_as_text					import load_profiles_as_text
from files.load_feedback_file						import load_feedback_file
from files.load_exp_descr							import load_exp_descr


class feedback_gen(Frame):

	def __init__(self, parent):
		Frame.__init__(self, parent)

		self.parent = parent
		self.initUI()


	def initUI(self):

		nPadx = 5
		nPady = 2
		self.parent.title("FeedbackTransp profiles setting")
		self.pack(fill=BOTH, expand=True)

#		self.columnconfigure(1, weight=1)

		menubar = Menu(self.master)
		filemenu = Menu(menubar, tearoff=0)
		filemenu.add_command(label="Load Mesh", 	 					command=self.cmd_load_mesh)
		filemenu.add_command(label="Save Mesh", 	 					command=self.cmd_save_soledge_mesh)
		filemenu.add_command(label="Save Mesh As..",				command=self.cmd_save_soledge_mesh_as)
		filemenu.add_separator()
		filemenu.add_command(label="Save Transport Feedback",		command=self.cmd_save_feedback_transp)
		filemenu.add_command(label="Save Puffing Feedback",		command=self.cmd_save_feedback_puffing)
		filemenu.add_command(label="Save Plasma files",			command=self.cmd_save_plasma)
		filemenu.add_separator()
		filemenu.add_command(label="Import profiles",			command=self.cmd_import_profiles)
		filemenu.add_command(label="Export profiles",			command=self.cmd_export_profiles)
		filemenu.add_separator()
		filemenu.add_command(label="Import feedbakProf file",	command=self.cmd_import_feedbakProf_file)
		filemenu.add_command(label="Import diffusion file",	command=self.cmd_import_diffusion_file)
		filemenu.add_separator()
		filemenu.add_command(label="Exit",							command=self.quit)
		menubar.add_cascade(label="File",							menu=filemenu)


		OptionsMenu = Menu(menubar, tearoff=0)

		option2dmenu = Menu(self.master)

		self.Shading2DCheckValue = IntVar()
		self.Shading2DCheckValue.set(1)
		option2dmenu.add_checkbutton(label="Shading", 			variable=self.Shading2DCheckValue)

		self.Borders2DCheckValue = IntVar()
		self.Borders2DCheckValue.set(0)
		option2dmenu.add_checkbutton(label="Zone borders",		variable=self.Borders2DCheckValue)

		self.Contour2DCheckValue = IntVar()
		self.Contour2DCheckValue.set(0)
		option2dmenu.add_checkbutton(label="Sep. contours",		variable=self.Contour2DCheckValue)
		
		OptionsMenu.add_cascade(label="2D plots",					menu=option2dmenu)


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

		"""
		OptionsMenu = Menu(menubar, tearoff=0)

		interpmenu = Menu(self.master)
		self.InterpRadioValue = IntVar()
		self.InterpRadioValue.set(0)
		interpmenu.add_radiobutton(label="linear", 			variable=self.InterpRadioValue, value=0)
		interpmenu.add_radiobutton(label="spline", 	 		variable=self.InterpRadioValue, value=1)

		optsmenu.add_cascade(label="Interpolation",			menu=interpmenu)
		
		menubar.add_cascade(label="Options", menu=OptionsMenu)
		"""

		helpmenu = Menu(menubar, tearoff=0)
		helpmenu.add_command(label="About...",				command=self.cmd_about)
		menubar.add_cascade(label="Help", menu=helpmenu)

		self.master.config(menu=menubar)

#		Plot area

		PlotFrame = LabelFrame(self, text="Plot")
		PlotFrame.grid(row=0, column=0,  rowspan=4, pady=nPady, padx=nPadx, sticky=E+W+N)

		self.Fig = Figure(figsize=(4.8, 5.5), dpi=100)
		self.Fig.patch.set_facecolor('white')

		# a tk.DrawingArea
		canvas = FigureCanvasTkAgg(self.Fig, master=PlotFrame)
		try:
			canvas.show()										#Matplotlib  up 2.2
		except:
			canvas.draw()
#		canvas.get_tk_widget().grid(row=1, column=0,  columnspan=4,  padx=nPadx, sticky=E+W+S+N)
		canvas.get_tk_widget().grid(row=1, column=0,  columnspan=4,  padx=nPadx)

		self.Ax = self.Fig.add_subplot(111)
		self.Cb = 0														#I need this to be able to remove it later
		self.Fig.tight_layout()

		toolbar_frame = Frame(PlotFrame)
		toolbar_frame.grid(row=0, column=0, columnspan=6, padx=nPadx, sticky=E+W+S+N)
		toolbar = NavigationToolbar2Tk(canvas, toolbar_frame )
		toolbar.update()
		
		iRow = 2
		self.EqCutBut = Button(PlotFrame, text="Eq. cut", command = self.cmd_equatorial_cut, width = 10)
		self.EqCutBut.grid(row=iRow, column=0, padx=nPadx, pady=nPady, sticky=W+E)

		self.DrawCutButValue = StringVar()
		self.DrawCutButValue.set("Draw cut")
		self.DrawCutBut = Button(PlotFrame, textvariable=self.DrawCutButValue, command = self.cmd_draw_cut, width = 10)
		self.DrawCutBut.grid(row=iRow, column=1, padx=nPadx, pady=nPady, sticky=W+E)

		CoordLabel = Label(PlotFrame, text="Coord.:", width = 5)
		CoordLabel.grid(row=iRow, column=2,  pady=nPady, padx=nPadx)

		self.CoordEntryValue = StringVar()
		self.CoordEntry = Entry(PlotFrame, textvariable = self.CoordEntryValue, width = 10)
		self.CoordEntry.grid(row=iRow, column=3,  pady=nPady, padx=nPadx)
		self.CoordEntry.bind("<Return>", self.cmd_enter_coord_entry)
		self.CoordEntry.state(["disabled"])

		"""
		ShowLabel = Label(PlotFrame, text="Zone sel.", width = 7)
		ShowLabel.grid(row=iRow, column=3,  pady=nPady, padx=nPadx)

		self.Config.ZonesOptions = ["Connected", "Single"]
		self.ZoneModeOptMenuValue = StringVar()
		self.ZoneModeOptMenuValue.set(self.Config.ZonesOptions[0])
		self.ZoneModeOptMenu = OptionMenu(PlotFrame, self.ZoneModeOptMenuValue, self.Config.ZonesOptions[0], self.Config.ZonesOptions[1])
		self.ZoneModeOptMenu.grid(row=iRow, column=4,  pady=nPady, padx=nPadx)
		
		iRow += 1

		self.SelectCutBut = Button(PlotFrame, text="Select cut", command = self.cmd_select_cut, width = 10)
		self.SelectCutBut.grid(row=iRow, column=0, pady=nPady)

		self.RemoveCutBut = Button(PlotFrame, text="Remove cut", command = self.cmd_remove_cut, width = 10)
		self.RemoveCutBut.grid(row=iRow, column=1, pady=nPady)
		"""

		iRow += 1

		self.GenerateBut = Button(PlotFrame, text="Generate Mesh", command = self.cmd_generate_mesh_values, width = 25)
		self.GenerateBut.grid(row=iRow, column=0, columnspan=2, pady=nPady, padx=nPadx, sticky=W+E)

		
		ShowLabel = Label(PlotFrame, text="Show", width = 4)
		ShowLabel.grid(row=iRow, column=2,  pady=nPady, padx=nPadx)
				
		self.FeedbackTabNames	= FEEDBACK_TAB_NAMES
		self.FeedbackNames		= FEEDBACK_NAMES

		self.FigureOptions		= ["Borders"]
		for k in self.FeedbackTabNames: self.FigureOptions.append(k)
		opts = ()
		for k in self.FigureOptions: opts += (k,)	
		self.FigureModeOptMenuValue = StringVar()
		self.FigureModeOptMenuValue.set(self.FigureOptions[0])
		self.FigureModeOptMenu = OptionMenu(PlotFrame, self.FigureModeOptMenuValue, *opts)
		self.FigureModeOptMenu.grid(row=iRow, column=3,  pady=nPady, padx=nPadx, sticky=W+E)
		self.FigureModeOptMenuValue.trace("w", self.cmd_figure_option_menu)

		iRow = 0

		self.FluxProfileWidgets  = []
		self.FluxProfNotebook = Notebook(self)
		for k in range(len(self.FeedbackNames)):
			self.FluxProfileWidgets.append(generate_profile(self.FluxProfNotebook, Name=self.FeedbackNames[k], Title="Flux profile"))
			self.FluxProfNotebook.add(self.FluxProfileWidgets[k].Frame, text=self.FeedbackTabNames[k])

		self.FluxProfNotebook.grid(row=iRow, column=1, columnspan=2,  pady=nPady, padx=nPadx)

		iRow += 1
		
#		FeedbackTransp parameters

		self.FeedbackNotebook = Notebook(self)
		self.FeedbackNotebook.grid(row=iRow, column=1, columnspan=2,  pady=nPady, padx=nPadx, sticky=W+E+N)

#		PLasma FeedbackTransp frame

		TranspFeedbackFrame = Frame(self.FeedbackNotebook)
		self.FeedbackNotebook.add(TranspFeedbackFrame , text="Trasport parameters")		

		lWidth = 35
		eWidth = 15
		nPady  = 1
		iFrow  = 0

		DminLabel = Label(TranspFeedbackFrame, text = "Minimal diffusivity (m^2/s):", width = lWidth)
		DminLabel.grid(row=iFrow, column=0,  pady=nPady*2, padx=nPadx, sticky=W)

		self.DminEntryValue = DoubleVar()
		DminEntry = Entry(TranspFeedbackFrame,  textvariable = self.DminEntryValue, width = eWidth)
		DminEntry.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx)

		iFrow += 1

		DmaxLabel = Label(TranspFeedbackFrame, text = "Maximal diffusivity (m^2/s):", width = lWidth)
		DmaxLabel.grid(row=iFrow, column=0,  pady=nPady, padx=nPadx, sticky=W)

		self.DmaxEntryValue = DoubleVar()
		DmaxEntry = Entry(TranspFeedbackFrame,  textvariable = self.DmaxEntryValue, width = eWidth)
		DmaxEntry.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx)

		iFrow += 1

		KeepLabel = Label(TranspFeedbackFrame, text = "Iteration window to evolve D:", width = lWidth)
		KeepLabel.grid(row=iFrow, column=0,  pady=nPady, padx=nPadx, sticky=W)

		self.KeepEntryValue = IntVar()
		KeepEntry = Entry(TranspFeedbackFrame,  textvariable = self.KeepEntryValue, width = eWidth)
		KeepEntry.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx)

		iFrow += 1

		GainLabel = Label(TranspFeedbackFrame, text = "Gain for FeedbackTransp loop (0 no FB):", width = lWidth)
		GainLabel.grid(row=iFrow, column=0,  pady=nPady, padx=nPadx, sticky=W)

		self.GainEntryValue = DoubleVar()
		GainEntry = Entry(TranspFeedbackFrame,  textvariable = self.GainEntryValue, width = eWidth)
		GainEntry.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx)

		iFrow += 1

		GainGLabel = Label(TranspFeedbackFrame, text = "Gradient Gain for FeedbackTransp loop (0 no FB):", width = lWidth)
		GainGLabel.grid(row=iFrow, column=0,  pady=nPady, padx=nPadx, sticky=W)

		self.GainGEntryValue = DoubleVar()
		GainGEntry = Entry(TranspFeedbackFrame,  textvariable = self.GainGEntryValue, width = eWidth)
		GainGEntry.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx, sticky=W)

		iFrow += 1

		LambdaLabel = Label(TranspFeedbackFrame, text = "Lambda decay length (m):", width = lWidth)
		LambdaLabel.grid(row=iFrow, column=0,  pady=nPady, padx=nPadx, sticky=W)

		self.LambdaEntryValue = DoubleVar()
		LambdaEntry = Entry(TranspFeedbackFrame,  textvariable = self.LambdaEntryValue, width = eWidth)
		LambdaEntry.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx)

		ExtrapLambdaBut = Button(TranspFeedbackFrame, text="Extrapolate by lambda", command = self.cmd_extrapolate_prof, width = 20)
		ExtrapLambdaBut.grid(row=iFrow, column=2, pady=nPady, padx=nPadx)

		iFrow += 1

		LambdaLabel = Label(TranspFeedbackFrame, text = "Use ballooning pattern:", width = lWidth)
		LambdaLabel.grid(row=iFrow, column=0,  pady=nPady, padx=nPadx, sticky=W)

		self.BallooningModeOptMenuValues = ("No", "1/Bt", "1/B")
		self.UsePatternOptMenuValue = StringVar()
		self.UsePatternOptMenuValue.set(self.BallooningModeOptMenuValues[DEF_BALLOONING_MODE])
		UsePatternOptMenu = OptionMenu(TranspFeedbackFrame, self.UsePatternOptMenuValue, *self.BallooningModeOptMenuValues)
		UsePatternOptMenu.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx, sticky=W+E)

		iFrow += 1

		DminLabel = Label(TranspFeedbackFrame, text = "Min density (m^-3):", width = lWidth)
		DminLabel.grid(row=iFrow, column=0,  pady=nPady, padx=nPadx, sticky=W)

		self.NminEntryValue = DoubleVar()
		NminEntry = Entry(TranspFeedbackFrame,  textvariable = self.NminEntryValue, width = eWidth)
		NminEntry.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx)

		iFrow += 1

		TminLabel = Label(TranspFeedbackFrame, text = "Min temperature (eV):", width = lWidth)
		TminLabel.grid(row=iFrow, column=0,  pady=nPady, padx=nPadx, sticky=W)

		self.TminEntryValue = DoubleVar()
		TminEntry = Entry(TranspFeedbackFrame,  textvariable = self.TminEntryValue, width = eWidth)
		TminEntry.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx)

	
#		Plasma parameters

		iFrow += 1
		
		PlasmaFrame = LabelFrame(TranspFeedbackFrame, text="Plasma parameters")
		PlasmaFrame.grid(row=iFrow, column=0, columnspan=2,  pady=nPady*2, padx=nPadx, sticky=W+E+N+S)

		iFrow = 0

		NRefLabel = Label(PlasmaFrame, text = "Reference density (m^-3):", width = lWidth)
		NRefLabel.grid(row=iFrow, column=0,  pady=nPady, padx=nPadx, sticky=W)

		self.NRefEntryValue = DoubleVar()
		NRefEntry = Entry(PlasmaFrame,  textvariable = self.NRefEntryValue, width = eWidth)
		NRefEntry.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx)

		iFrow += 1

		TRefLabel = Label(PlasmaFrame, text = "Reference temperature (eV):", width = lWidth)
		TRefLabel.grid(row=iFrow, column=0,  pady=nPady, padx=nPadx, sticky=W)

		self.TRefEntryValue = DoubleVar()
		TRefEntry = Entry(PlasmaFrame,  textvariable = self.TRefEntryValue, width = eWidth)
		TRefEntry.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx)


#		Gas-puffing FeedbackTransp frame

		GasPuffingFeedbackFrame = Frame(self.FeedbackNotebook)
		self.FeedbackNotebook.add(GasPuffingFeedbackFrame , text="Gas-puffing parameters")		

		nPady = 2
		iFrow = 0

		PuffDensityLabel = Label(GasPuffingFeedbackFrame, text = "Density:", width = lWidth)
		PuffDensityLabel.grid(row=iFrow, column=0,  pady=nPady*2, padx=nPadx, sticky=W)

		self.PuffDensityEntryValue = DoubleVar()
		PuffDensityEntry = Entry(GasPuffingFeedbackFrame,  textvariable = self.PuffDensityEntryValue, width = eWidth)
		PuffDensityEntry.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx)

		self.AutoDensityCkButValue = IntVar()
		self.AutoDensityCkButValue.set(1)
		AutoDensityCkBut = Checkbutton(GasPuffingFeedbackFrame, text="Auto", variable = self.AutoDensityCkButValue, width = lWidth)
		AutoDensityCkBut.grid(row=iFrow, column=2,  pady=nPady, padx=nPadx, sticky=W+E)

		"""
		opts = ("Auto", "Manual")
		self.AutoDensityOptMenuValue = StringVar()
		self.AutoDensityOptMenuValue.set(opts[0])
		AutoDensityOptMenu = OptionMenu(GasPuffingFeedbackFrame, self.AutoDensityOptMenuValue, *opts)
		AutoDensityOptMenu.grid(row=iFrow, column=2,  pady=nPady, padx=nPadx, sticky=W+E)
		"""

		iFrow += 1

		PuffGainLabel = Label(GasPuffingFeedbackFrame, text = "Gain:", width = lWidth)
		PuffGainLabel.grid(row=iFrow, column=0,  pady=nPady, padx=nPadx, sticky=W)

		self.PuffGainEntryValue = DoubleVar()
		PuffGainEntry = Entry(GasPuffingFeedbackFrame,  textvariable = self.PuffGainEntryValue, width = eWidth)
		PuffGainEntry.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx)

		iFrow += 1

		PuffTauiLabel = Label(GasPuffingFeedbackFrame, text = "Tau_i:", width = lWidth)
		PuffTauiLabel.grid(row=iFrow, column=0,  pady=nPady, padx=nPadx, sticky=W)

		self.PuffTauiEntryValue = DoubleVar()
		PuffTauiEntry = Entry(GasPuffingFeedbackFrame,  textvariable = self.PuffTauiEntryValue, width = eWidth)
		PuffTauiEntry.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx)

		iFrow += 1

		PuffTaudLabel = Label(GasPuffingFeedbackFrame, text = "Tau_d:", width = lWidth)
		PuffTaudLabel.grid(row=iFrow, column=0,  pady=nPady, padx=nPadx, sticky=W)

		self.PuffTaudEntryValue = DoubleVar()
		PuffTaudEntry = Entry(GasPuffingFeedbackFrame,  textvariable = self.PuffTaudEntryValue, width = eWidth)
		PuffTaudEntry.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx)

		iFrow += 1

		PuffMinPuffLabel = Label(GasPuffingFeedbackFrame, text = "Min Puff (particle/s):", width = lWidth)
		PuffMinPuffLabel.grid(row=iFrow, column=0,  pady=nPady, padx=nPadx, sticky=W)

		self.PuffMinPuffEntryValue = DoubleVar()
		PuffMinPuffEntry = Entry(GasPuffingFeedbackFrame,  textvariable = self.PuffMinPuffEntryValue, width = eWidth)
		PuffMinPuffEntry.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx)

		iFrow += 1

		PuffMaxPuffLabel = Label(GasPuffingFeedbackFrame, text = "Max Puff (particle/s):", width = lWidth)
		PuffMaxPuffLabel.grid(row=iFrow, column=0,  pady=nPady, padx=nPadx, sticky=W)

		self.PuffMaxPuffEntryValue = DoubleVar()
		PuffMaxPuffEntry = Entry(GasPuffingFeedbackFrame,  textvariable = self.PuffMaxPuffEntryValue, width = eWidth)
		PuffMaxPuffEntry.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx)

		iFrow += 1

		PuffPointToFeedLabel = Label(GasPuffingFeedbackFrame, text = "Point to FeedbackTransp:", width = lWidth)
		PuffPointToFeedLabel.grid(row=iFrow, column=0,  pady=nPady, padx=nPadx, sticky=W)

		self.SetPuffButValue = StringVar()
		self.SetPuffButValue.set("Set")
		self.SetPuffBut = Button(GasPuffingFeedbackFrame, textvariable = self.SetPuffButValue, command = self.cmd_set_puffing_refpoint, width = eWidth)
		self.SetPuffBut.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx)

		self.DminEntryValue.set(1.e-3)
		self.DmaxEntryValue.set(10.)
		self.KeepEntryValue.set(1000)
		self.GainEntryValue.set("{:.4e}".format(1.e-5))
		self.GainGEntryValue.set("{:.4e}".format(1.e-5))
		self.LambdaEntryValue.set(0.0015)

		self.NminEntryValue.set(0.02)
		self.TminEntryValue.set(0.02)
		self.NRefEntryValue.set("{:.4e}".format(1e20))
		self.TRefEntryValue.set(1.e3)

		self.PuffDensityEntryValue.set("{:.4e}".format(1.e20))
		self.PuffGainEntryValue.set(5000.)
		self.PuffTauiEntryValue.set("{:.4e}".format(1.e-6))
		self.PuffTaudEntryValue.set("{:.4e}".format(1.e10))
		self.PuffMinPuffEntryValue.set("{:.4e}".format(1.5e23))
		self.PuffMaxPuffEntryValue.set("{:.4e}".format(1.5e17))

		for k in range(len(self.FluxProfileWidgets)):
			flux_profile_copy_from  = []
			tab_profile_copy_from   = []
			for i in range(len(self.FluxProfileWidgets)):
				if(i != k):
					flux_profile_copy_from.append(self.FluxProfileWidgets[i])
					tab_profile_copy_from.append(self.FeedbackTabNames[i])

			self.FluxProfileWidgets[k].copy_from_set(flux_profile_copy_from, tab_profile_copy_from)

		self.resize_mode			= False
		self.InDrawCut				= False
		self.InDrawPuff				= False
		
		self.MeshFile				= False
		self.Config					= reset_config()
		self.Config.FeedbackTransp	= 0
		self.CutLine				= np.empty((2,2),dtype='f8')
		
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
			self.Config = load_soledge_mesh_file(MeshFile)
			
			if(self.Config.feedback_transp_OK or self.Config.feedback_puffing_OK): self.set_feedback_parameters()

			if(self.Config.feedback_transp_OK):
				self.Config.FeedbackTransp.Cut.Flux.Lengths = compute_mesh_intersections(self.Config, self.Config.FeedbackTransp.Cut.Flux)			#Compute conversion lengths
				self.set_cut_zones(self.Config.FeedbackTransp.Cut)
				self.set_profile_data()

				generate_feedback_profiles(self, self.Config)
				generate_feedback_mesh(self, self.Config)
			else:
				for k in range(len(self.Config.MagZones)): self.Config.MagZones[k].IsFeedbackProf	= 0
			
			for k in range(len(self.Config.MagZones)): 	self.Config.MagZones[k].IsFeedbackValues	= self.Config.feedback_transp_OK
			for k in range(len(self.Config.Zones)):		self.Config.Zones[k].IsFeedbackValues		= self.Config.feedback_transp_OK

				
			self.resize_mode 	= True
			self.update_figure()

			self.MeshFile		= MeshFile
			self.MeshDir  		= os.path.dirname(self.MeshFile)
			if((len(self.MeshDir)==0) or ((len(self.MeshDir)==1) and  (self.MeshDir == "/"))): self.MeshDir="./"
			

		else:
#		except:
			messagebox.showerror("Open Source File", "Failed to read mesh file:{:s}\n".format(File))

	def cmd_save_soledge_mesh(self):
		if(not self.MeshFile): return
		self.set_feedback_transp_ok()
		if((not self.Config.feedback_transp_OK) and (not self.Config.feedback_puffing_OK)):
			messagebox.showwarning("Saving file", "Mesh file can be saved only when zones are ok or puffing is defined")
			return

		try:
			self.Config.FeedbackTransp = self.Config.FeedbackTransp
			save_soledge_mesh_file(self.MeshFile, self.Config)
		except:
			messagebox.showerror("Save Mesh file", "Failed to write mesh file:{:s}\n".format(self.MeshFile))


	def cmd_save_soledge_mesh_as(self):
		if(not self.MeshFile): return
		
		self.set_feedback_transp_ok()
		if((not self.Config.feedback_transp_OK) and (not self.Config.feedback_puffing_OK)):
			messagebox.showwarning("Saving file", "Mesh file can be saved only when zones are ok or puffing is defined ")
			return

		MeshFile = filedialog.asksaveasfilename(initialdir=self.MeshDir, filetypes=(("Mesh files", "*.h5"), ("All files", "*.*") ))

#		try:
		if(True):
			self.Config.FeedbackTransp = self.Config.FeedbackTransp
			save_soledge_mesh_file(MeshFile, self.Config)
			self.MeshFile = MeshFile
			self.MeshDir  = os.path.dirname(self.MeshFile)
#		except:
#			messagebox.showerror("Save Mesh file", "Failed to write mesh file:{:s}\n".format(self.MeshFile))

	def cmd_save_feedback_transp(self):
		if(not self.MeshFile): return
		
		self.set_feedback_transp_ok()
		if(not self.Config.feedback_transp_OK):
			messagebox.showwarning("Saving file", "FeedbackTransp files can be saved only when profiles are ok")
			return
		if(os.path.split(self.MeshFile)[1] == "mesh.h5"):
			FeedbackPath = self.MeshDir
		else:
			FeedbackPath = self.MeshDir+"/"+os.path.splitext(os.path.split(self.MeshFile)[1])[0]+"_Feedback/"

		save_feedback_transp_files(FeedbackPath, self.Config.FeedbackTransp)
		
	def cmd_save_plasma(self):
		if(not self.MeshFile): return
		
		self.set_feedback_transp_ok()
		if(not self.Config.feedback_transp_OK):
			messagebox.showwarning("Saving file", "Plasma files can be saved only when all zones are ok")
			return
		if(os.path.split(self.MeshFile)[1] == "mesh.h5"):
			PlasmaPath = self.MeshDir+"/Results/"
		else:
			PlasmaPath = self.MeshDir+"/"+os.path.splitext(os.path.split(self.MeshFile)[1])[0]+"_Results/"

		save_feedback_plasma_files(PlasmaPath, self.Config)

	def cmd_save_feedback_puffing(self):
		if(not self.MeshFile): return
		
		if(not self.Config.feedback_puffing_OK):
			messagebox.showwarning("Saving file", "Feedback puffing file can be saved only when point is selected")
			return
		if(os.path.split(self.MeshFile)[1] == "mesh.h5"):
			FeedbackPath = self.MeshDir+"/"
		else:
			FeedbackPath = self.MeshDir+"/"+os.path.splitext(os.path.split(self.MeshFile)[1])[0]+"_Feedback/"
			
#		If requested update target density

		if((self.Config.feedback_transp_OK) and (self.AutoDensityCkButValue.get() == 1)) :
			nTarget = self.Config.Zones[self.Config.FeedbackPuffing.kTarget].FeedbackValues[0][self.Config.FeedbackPuffing.iTarget,self.Config.FeedbackPuffing.jTarget]*self.Config.FeedbackTransp.Data.NRef
			self.PuffDensityEntryValue.set("{:.4e}".format(nTarget))

		self.get_feedback_parameters()
		save_feedback_puffing_file(FeedbackPath, self.Config.FeedbackPuffing)

	def cmd_import_profiles(self):
		if(not self.MeshFile): return
		
		if(self.Config.FeedbackTransp == 0):
			messagebox.showwarning("Loading file", "Profiles can be imported when a cut has been defined")
			return

		ProfFile = filedialog.askopenfilename(initialdir=self.MeshDir, filetypes=(("Comma separated values", "*.csv"),  \
																					("Tab separated values",    "*.tsv"),  \
																					("Exp description file",    "*.exp"),  \
																					("Space separated values",  "*.*")))

#		try:
		if(ProfFile):
			xName = LENGTH_TYPES[0]
		
			self.get_profile_data()

			FileExt = os.path.splitext(os.path.split(ProfFile)[1])[1]
			if(FileExt == ".exp"):
				exp_descr = load_exp_descr(ProfFile)								#load diags description
				Cut = self.Config.FeedbackTransp.Cut.Flux
				if(hasattr(Cut, "Lengths")):
					Lengths  = Cut.Lengths
				else:
					Lengths  = np.empty((0,0), dtype='f8')
				for i in range(len(Cut.Profiles)):
					if(DEBUG > 0): print("\tSearch for following data:    ", Cut.Profiles[i].Name)
					diag_data = get_exp_data_point(xName, Cut.Profiles[i].Name, exp_descr, Lengths, OutxOffset = Cut.d12[0,0], DiagType = "POINT")
					if(len(diag_data) == 1):
						Cut.Profiles[i].xValues = np.append(np.append(Cut.Profiles[i].xValues[0],diag_data[0].xValues),Cut.Profiles[i].xValues[-1])
						Cut.Profiles[i].Values  = np.append(np.append(diag_data[0].yValues[0],diag_data[0].yValues),diag_data[0].yValues[-1])
			else:
				load_profiles_as_text(ProfFile, self.Config.FeedbackTransp.Cut.Flux, xName, OutxOffset=self.Config.FeedbackTransp.Cut.Flux.d12[0,0])

			self.set_profile_data()
			self.update_figure()

#		except:
#			messagebox.showerror("Save Mesh file", "Failed to write profiles file:{:s}\n".format(ProfFile))

	def cmd_export_profiles(self):
		if(not self.MeshFile): return
		
		if(self.Config.FeedbackTransp == 0):
			messagebox.showwarning("Saving file", "Profiles can be saved only when a cut has been defined")
			return

		ProfFile = filedialog.asksaveasfilename(initialdir=self.MeshDir, filetypes=(("Comma separated values", "*.csv"),  \
																					("Tab separated values",    "*.tsv"),  \
																					("Space separated values",  "*.*")))

#		try:
		if(ProfFile):
			xName = LENGTH_TYPES[0]
	
			self.get_profile_data()
			save_profiles_as_text(ProfFile, self.Config.FeedbackTransp.Cut.Flux, xName)
#		except:
#			messagebox.showerror("Save Mesh file", "Failed to write profiles file:{:s}\n".format(ProfFile))


	def cmd_import_feedbakProf_file(self):
		if(not self.MeshFile): return
		
		if(self.Config.FeedbackTransp == 0):
			messagebox.showwarning("Loading file", "SOLEDGE feedbakProf can be imported when at least one cut is selected")
			return

		FeedFile = filedialog.askopenfilename(initialdir=self.MeshDir)

		if(FeedFile):
			load_feedback_file(FeedFile, 0, self.Config, self.Config.FeedbackTransp.Cut)

			self.set_profile_data()
			self.update_figure()

	def cmd_import_diffusion_file(self):
		if(not self.MeshFile): return
		
		if(self.Config.FeedbackTransp == 0):
			messagebox.showwarning("Loading file", "SOLEDGE diffusion can be imported when at least one cut is selected")
			return

		DiffFile = filedialog.askopenfilename(initialdir=self.MeshDir)

		if(DiffFile):
			load_feedback_file(DiffFile, 1, self.Config, self.Config.FeedbackTransp.Cut)

			self.set_profile_data()
			self.update_figure()


	def set_feedback_transp_ok(self):
		if((not self.Config.feedback_transp_OK) and (self.Config.FeedbackTransp != 0)):
			self.generate_mesh_values()
			self.Config.feedback_transp_OK	= True
			for k in range(len(self.Config.MagZones)):
				self.Config.feedback_transp_OK = self.Config.feedback_transp_OK and self.Config.MagZones[k].IsFeedbackValues
		

#######################################
#	Option  menu show walls
#######################################
		
	def cmd_ShowWallsCheckValue(self, a, b, c):
		self.update_figure()
		return

#######################################
#	Draw a cut to set profile
#######################################

	def cmd_equatorial_cut(self):

#		Search Core zone
		MagZones = self.Config.MagZones

		Rcore, Zcore, CoreMegazone = get_rz_core_sep(self.Config, core_and_sep = False)

		cCore			= types.SimpleNamespace()
		cCore.arc		= [types.SimpleNamespace()]
		cCore.arc[0].x	= Rcore
		cCore.arc[0].y	= Zcore

		cLine			= types.SimpleNamespace()
		cLine.arc		= [types.SimpleNamespace()]
		cLine.arc[0].x	= np.array([Rcore.min()-1., Rcore.max()+1.])

		iMax = np.argmax(Rcore)
		cLine.arc[0].y	= np.array([Zcore[iMax],Zcore[iMax]])
			
		X	= intersect_contour(cCore, cLine)
		if(len(X) != 2):
			print("Error in finding core center")
			return

		self.CutLine[0,0] = 0.5*(X[0].x + X[1].x)
		self.CutLine[1,0] = 2.*self.CutLine[0,0]

		self.CutLine[0,1] = cLine.arc[0].y[0]
		self.CutLine[1,1] = cLine.arc[0].y[0]
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

		Cut				 = types.SimpleNamespace()
		Cut.Flux 		 = find_zones_intersections(self.Config, self.CutLine)
		Cut.Flux.Lengths = compute_mesh_intersections(self.Config, Cut.Flux)				#Compute lengths for conversion between L, X, Rho_v, Rho_F
		Cut.sel			 = False
			
		if(len(Cut.Flux.nz) != 0):

			self.set_cut_zones(Cut)
				
			Cut.Flux.Profiles = []
			for k in range(len(self.FeedbackNames)):
				Cut.Flux.Profiles.append(types.SimpleNamespace())
				Cut.Flux.Profiles[-1].Name	  = self.FeedbackNames[k]
				Cut.Flux.Profiles[-1].xValues = np.array([Cut.Flux.d12[0,0],Cut.Flux.d12[-1,1]])
				Cut.Flux.Profiles[-1].Values  = np.ones((2), dtype='f8')
			
			if(self.Config.FeedbackTransp == 0): self.Config.FeedbackTransp = types.SimpleNamespace()
			self.Config.FeedbackTransp.Cut  = Cut
			self.Config.feedback_transp_OK	= False
			self.set_profile_data()
				
		self.update_figure()

	def set_cut_zones(self, Cut):
				
#		Build list of zone affected by the cut
#		Flux direction
		MagZones = self.Config.MagZones
		Zones	 = self.Config.Zones

		Cut.Flux.MagZones = []

		for iZone in range(len(MagZones)):
			MagZones[iZone].IsFeedbackProf		= 0
			MagZones[iZone].IsFeedbackValues	= False

		for iZone in range(len(Zones)): Zones[iZone].IsFeedbackValues	= False

		for k in Cut.Flux.nz:
			xMin = MagZones[k].x[0, 0]-0.01
			xMax = MagZones[k].x[-1,0]+0.01
			zones_list = np.empty((0),dtype='i4')
			for iZone in range(len(MagZones)):
				if((MagZones[iZone].x[0, 0] >= xMin) and (MagZones[iZone].x[-1, 0] <= xMax)):
					if(MagZones[iZone].IsFeedbackProf == 0): zones_list	= np.append(zones_list, iZone)
					MagZones[iZone].IsFeedbackProf = 2
				elif((((MagZones[iZone].x[0, 0] >= xMin) and (MagZones[iZone].x[0, 0] <= xMax)) or
					((MagZones[iZone].x[-1,0] >= xMin) and (MagZones[iZone].x[-1,0] <= xMax))) and
					(MagZones[iZone].IsFeedbackProf == 0)):
						zones_list						= np.append(zones_list, iZone)
						MagZones[iZone].IsFeedbackProf 	= 1
						
			Cut.Flux.MagZones.append(zones_list)	
	
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

	def cmd_generate_mesh_values(self):
		if(self.Config.FeedbackTransp == 0):
			messagebox.showwarning("Generating mesh", "Mesh generation impossible because cut is not defined")
			return
		self.generate_mesh_values()
		
	def generate_mesh_values(self):
		if(self.Config.FeedbackTransp == 0): return
		
		self.get_profile_data()
		self.get_feedback_parameters()
		generate_feedback_profiles(self, self.Config)
		generate_feedback_mesh(self, self.Config)
		self.update_figure()
									
	def cmd_figure_option_menu(self,a,b,c):
		self.update_figure()


############################################
#	Set puffing reference point
############################################

	def cmd_set_puffing_refpoint(self):
		if(self.FigureOptions.index(self.FigureModeOptMenuValue.get()) > 0): return
		if(not self.InDrawPuff):
			self.SetPuffButValue.set("Cancel")

			self.cid_mouse_button_press_puffing = self.Fig.canvas.mpl_connect('button_press_event', self.cmd_button_press_set_puff)
			self.cursor = Cursor(self.Ax, useblit=False, color='red', linewidth=1)
			self.InDrawPuff = True
		else:
			self.InDrawPuff = False
			self.SetPuffButValue.set("Set")
			if(self.cid_mouse_button_press_puffing != 0):
				self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press_puffing)
				self.cid_mouse_button_press_puffing = 0
				self.cursor=0

		self.update_figure()

	def cmd_button_press_set_puff(self, event):
		if(self.cid_mouse_button_press_puffing != 0):
			self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press_puffing)
			self.cid_mouse_button_press_puffing = 0
			self.cursor=0
	
		if(not self.Config.feedback_puffing_OK): self.Config.FeedbackPuffing = types.SimpleNamespace()
		self.Config.feedback_puffing_OK = True

		FeedbackPuffing = self.Config.FeedbackPuffing
		Zones			= self.Config.Zones
		Rpoint	= event.xdata
		Zpoint	= event.ydata
		dmin	= 1e20
		for k in range(len(Zones)):
			i,j    = np.unravel_index(np.argmin(np.sqrt((Zones[k].gridRc - Rpoint)**2 + (Zones[k].gridZc - Zpoint)**2)), Zones[k].gridRc.shape)
			dmin_k = sqrt((Zones[k].gridRc[i,j] - Rpoint)**2 + (Zones[k].gridZc[i,j] - Zpoint)**2)
			if(dmin_k < dmin):
				FeedbackPuffing.iTarget = i
				FeedbackPuffing.jTarget = j
				FeedbackPuffing.kTarget = k
				dmin					= dmin_k

		FeedbackPuffing.R = Zones[FeedbackPuffing.kTarget].gridRc[FeedbackPuffing.iTarget,FeedbackPuffing.jTarget]
		FeedbackPuffing.Z = Zones[FeedbackPuffing.kTarget].gridZc[FeedbackPuffing.iTarget,FeedbackPuffing.jTarget]
		
		self.set_feedback_transp_ok()
		if((self.Config.feedback_transp_OK) and (self.AutoDensityCkButValue.get() == 1)) :
			nTarget = Zones[FeedbackPuffing.kTarget].FeedbackValues[0][FeedbackPuffing.iTarget,FeedbackPuffing.jTarget]*self.Config.FeedbackTransp.Data.NRef
			self.PuffDensityEntryValue.set("{:.4e}".format(nTarget))

		self.cmd_set_puffing_refpoint()

############################################
#	Set & get profile data on profile widget
############################################

	def set_profile_data(self):
		for i in range(len(self.FluxProfileWidgets)):
			self.FluxProfileWidgets[i].set_data(self.Config.FeedbackTransp.Cut.Flux,   self.Config.FeedbackTransp.Cut.Flux.Profiles[i])
			
	def get_profile_data(self):
		for i in range(len(self.FluxProfileWidgets)):
			self.Config.FeedbackTransp.Cut.Flux.Profiles[i],  self.Config.FeedbackTransp.Cut.Flux.d12   = self.FluxProfileWidgets[i].get_data()

		for i in range(len(self.Config.FeedbackTransp.Cut.Flux.nz)):
			for iZone in self.Config.FeedbackTransp.Cut.Flux.MagZones[i]:
				self.Config.MagZones[iZone].IsFeedbackValues = False

	def get_feedback_parameters(self):
		if(self.Config.FeedbackTransp != 0):
			FeedbackTransp = self.Config.FeedbackTransp
			FeedbackTransp.Data				= types.SimpleNamespace()
			FeedbackTransp.Data.Dmin		= self.DminEntryValue.get()
			FeedbackTransp.Data.Dmax		= self.DmaxEntryValue.get()
			FeedbackTransp.Data.Keep		= self.KeepEntryValue.get()
			FeedbackTransp.Data.Gain		= self.GainEntryValue.get()
			FeedbackTransp.Data.GainG		= self.GainGEntryValue.get()
			FeedbackTransp.Data.Lambda		= self.LambdaEntryValue.get()
			FeedbackTransp.Data.UsePattern	= self.BallooningModeOptMenuValues.index(self.UsePatternOptMenuValue.get())
			
			FeedbackTransp.Data.Nmin		= self.NminEntryValue.get()
			FeedbackTransp.Data.Tmin		= self.TminEntryValue.get()
			FeedbackTransp.Data.NRef		= self.NRefEntryValue.get()
			FeedbackTransp.Data.TRef		= self.TRefEntryValue.get()

#			self.Config.FeedbackTransp.Data.InterpMode	= self.InterpRadioValue.get()
			self.Config.FeedbackTransp.Data.InterpMode	= 0

		if(self.Config.FeedbackPuffing != 0):
			FeedbackPuffing 		= self.Config.FeedbackPuffing

			FeedbackPuffing.AutoTarget	= self.AutoDensityCkButValue.get()
			FeedbackPuffing.nTarget		= self.PuffDensityEntryValue.get()
			FeedbackPuffing.Gain		= self.PuffGainEntryValue.get()
			FeedbackPuffing.Tau_i		= self.PuffTauiEntryValue.get()
			FeedbackPuffing.Tau_d		= self.PuffTaudEntryValue.get()
			FeedbackPuffing.MaxPuff		= self.PuffMinPuffEntryValue.get()
			FeedbackPuffing.MinPuff		= self.PuffMaxPuffEntryValue.get()


	def set_feedback_parameters(self):
		if(self.Config.feedback_transp_OK):
			FeedbackTransp = self.Config.FeedbackTransp

			self.DminEntryValue.set(FeedbackTransp.Data.Dmin)
			self.DmaxEntryValue.set(FeedbackTransp.Data.Dmax)
			self.KeepEntryValue.set(FeedbackTransp.Data.Keep)
			self.GainEntryValue.set("{:.4e}".format(FeedbackTransp.Data.Gain))
			self.GainGEntryValue.set("{:.4e}".format(FeedbackTransp.Data.GainG))
			self.LambdaEntryValue.set(FeedbackTransp.Data.Lambda)
			self.UsePatternOptMenuValue.set(self.BallooningModeOptMenuValues[FeedbackTransp.Data.UsePattern])

			self.NminEntryValue.set(FeedbackTransp.Data.Nmin)
			self.TminEntryValue.set(FeedbackTransp.Data.Tmin)
			self.NRefEntryValue.set("{:.4e}".format(FeedbackTransp.Data.NRef))
			self.TRefEntryValue.set(FeedbackTransp.Data.TRef)

#			self.InterpRadioValue.set(FeedbackTransp.Data.InterpMode)

		if(self.Config.feedback_puffing_OK):
			FeedbackPuffing = self.Config.FeedbackPuffing

			self.AutoDensityCkButValue.set(FeedbackPuffing.AutoTarget)
			self.PuffDensityEntryValue.set("{:.4e}".format(FeedbackPuffing.nTarget))
			self.PuffGainEntryValue.set(FeedbackPuffing.Gain)
			self.PuffTauiEntryValue.set(FeedbackPuffing.Tau_i)
			self.PuffTaudEntryValue.set(FeedbackPuffing.Tau_d)
			self.PuffMinPuffEntryValue.set("{:.4e}".FeedbackPuffing.MinPuff)
			self.PuffMaxPuffEntryValue.set("{:.4e}".FeedbackPuffing.MaxPuff)

	def cmd_extrapolate_prof(self):

		self.get_profile_data()

		FeedbackTransp = self.Config.FeedbackTransp
		iZoneSep = 0

		nInLambda = 3
		nGrid	  = 200
		TotDist   = sqrt((FeedbackTransp.Cut.Flux.r12[-1,1] - 	FeedbackTransp.Cut.Flux.r12[iZoneSep,1])**2 + (FeedbackTransp.Cut.Flux.z12[-1,1] - 	FeedbackTransp.Cut.Flux.z12[iZoneSep,1])**2)
		GridDist  = np.arange(nGrid+1, dtype='f8')/nGrid

		GridDist *=TotDist

		Lambda = self.LambdaEntryValue.get()*5./2.																#5./2 Factor from Lambda heat to lambda N & T
		for iProf in range(4):																					#Extrapolate Ne, Te, Ne, Ti only
			xValues	= np.copy(FeedbackTransp.Cut.Flux.Profiles[iProf].xValues)
			Values  = np.copy(FeedbackTransp.Cut.Flux.Profiles[iProf].Values)
			iFirst  = np.where(xValues <= FeedbackTransp.Cut.Flux.d12[iZoneSep,1])[0]
			nFirst  = len(iFirst)
			if(nFirst < 2): continue

			iFirst  = np.sort(iFirst)																			#Extrapolate to separatrix
			if(xValues[iFirst[-1]] <= FeedbackTransp.Cut.Flux.d12[iZoneSep,1]):
				iSepVal = Values[iFirst[-1]] + (Values[iFirst[-1]] - Values[iFirst[-2]])/(xValues[iFirst[-1]] - xValues[iFirst[-2]])*(FeedbackTransp.Cut.Flux.d12[iZoneSep,1]-xValues[iFirst[-1]])
				nFirst += 1

			nNew	= min(11,max(5,floor(TotDist/Lambda*nInLambda)))

			FeedbackTransp.Cut.Flux.Profiles[iProf].xValues = np.empty((nNew + nFirst), dtype='f8')
			FeedbackTransp.Cut.Flux.Profiles[iProf].Values  = np.empty((nNew + nFirst), dtype='f8')

			FeedbackTransp.Cut.Flux.Profiles[iProf].xValues[:len(iFirst)] = xValues[:len(iFirst)]
			FeedbackTransp.Cut.Flux.Profiles[iProf].Values[:len(iFirst)]  = Values[:len(iFirst)]

			if(nFirst > len(iFirst)):
				FeedbackTransp.Cut.Flux.Profiles[iProf].xValues[len(iFirst)] = FeedbackTransp.Cut.Flux.d12[1,0]
				FeedbackTransp.Cut.Flux.Profiles[iProf].Values[len(iFirst)]  = iSepVal
			else:
				iSepVal	= Values[len(iFirst)-1]

			if(nNew > floor(TotDist/Lambda*nInLambda)):
				dNew	= TotDist*np.arange(1,nNew+1)/nNew
			else:
				dNew	= np.append(np.arange(1,nNew)*Lambda/nInLambda, TotDist)

			VNew	= iSepVal*np.exp(-dNew/Lambda)
			FeedbackTransp.Cut.Flux.Profiles[iProf].Values[nFirst:] = VNew
				
			FeedbackTransp.Cut.Flux.Profiles[iProf].xValues[nFirst:] = dNew + FeedbackTransp.Cut.Flux.d12[iZoneSep,1]
			FeedbackTransp.Cut.Flux.Profiles[iProf].xValues[-1]	   = xValues[-1]

		self.set_profile_data()

#==============
#	Main figure 
#==============

#update main figure

	def update_figure(self):

		Config			= self.Config
		MagZones		= self.Config.MagZones
		Zones			= self.Config.Zones
		X_points		= self.Config.X_points
		FeedbackTransp	= self.Config.FeedbackTransp
		FeedbackPuffing	= self.Config.FeedbackPuffing
		Fig				= self.Fig
		Ax				= self.Ax

		FigureMode = self.FigureOptions.index(self.FigureModeOptMenuValue.get())	
		nMagZones  = len(MagZones)
		nZones	   = len(Zones)
		
		if(FigureMode > 0):
			iFead	  = FigureMode - 1
			vMin	  = np.empty(nZones, dtype='f8') 
			vMax	  = np.empty(nZones, dtype='f8')
			nOkZones  = 0
			for k in range(nZones):
				if(Zones[k].IsFeedbackValues):
					vMin[nOkZones] = np.min(Zones[k].FeedbackValues[iFead])
					vMax[nOkZones] = np.max(Zones[k].FeedbackValues[iFead])
					nOkZones	   += 1
			if(nOkZones == 0): return

		Xlims = Ax.get_xlim()
		Ylims = Ax.get_ylim()

		if(self.Cb != 0):
			Fig.clf()
			self.Ax = Ax = Fig.add_subplot(111)
			self.Cb = 0
		else:
			Ax.clear()		
		
		Ax.tick_params(axis='both', which='major', labelsize=8)
		
		if(FigureMode == 0):
	
			if(Config.Mesh_OK):
				for k in range(nMagZones):
					EdgeLineType	='b-'
					EdgeLineWidth	= 1
					MagZones[k].gridR
					Ax.plot(MagZones[k].gridR[0,:],  MagZones[k].gridZ[0,:],  EdgeLineType, linewidth=EdgeLineWidth)		#south
					Ax.plot(MagZones[k].gridR[:,-1], MagZones[k].gridZ[:,-1], EdgeLineType, linewidth=EdgeLineWidth)		#East
					Ax.plot(MagZones[k].gridR[-1,:], MagZones[k].gridZ[-1,:], EdgeLineType, linewidth=EdgeLineWidth)		#North
					Ax.plot(MagZones[k].gridR[:,0],  MagZones[k].gridZ[:,0],  EdgeLineType, linewidth=EdgeLineWidth)		#West
	
					iCenter = floor(MagZones[k].gridR.shape[1]/2)
					R=(MagZones[k].gridR[0,iCenter]+MagZones[k].gridR[-1,iCenter])*0.5
					Z=(MagZones[k].gridZ[0,iCenter]+MagZones[k].gridZ[-1,iCenter])*0.5
					if(MagZones[k].IsFeedbackProf == 2):
						Ax.text(R,Z,"{:d}".format(k+1), color='blue',  bbox=dict(facecolor='blue', alpha=0.5))
					elif(MagZones[k].IsFeedbackProf == 1):
						Ax.text(R,Z,"{:d}".format(k+1), color='blue',  bbox=dict(facecolor='cyan', alpha=0.5))
					else:
						Ax.text(R,Z,"{:d}".format(k+1), color='red', bbox=dict(facecolor='red', alpha=0.5))
	
			if(Config.equ_OK):
				for k in range(len(X_points)):
					Ax.contour(Config.r2D, Config.z2D, Config.flux2D, colors="r", linestyles='-', levels=[X_points[k].psi])
	#				Ax.text(X_points[k].R,X_points[k].Z,"{:d}".format(X_points[k].index),size=18)
	
			if(Config.wall_OK): plot2d_walls(Ax, Config.Walls, self.ShowPlasmaWallCheckValue.get(), self.ShowEireneWallCheckValue.get(), self.ShowExtraWallCheckValue.get())
	
			if(FeedbackTransp != 0):
				for i in range(len(FeedbackTransp.Cut.Flux.nz)):
					Ax.plot(FeedbackTransp.Cut.Flux.r12[i,:], FeedbackTransp.Cut.Flux.z12[i,:],'b.--', linewidth=2)
				
			if(FeedbackPuffing != 0):
				Ax.plot(FeedbackPuffing.R, FeedbackPuffing.Z,'m.',markersize=8, linewidth=3)
				Ax.text(FeedbackPuffing.R, FeedbackPuffing.Z, "PR", color='m', fontsize=18)

			if((not ((Xlims[0] == 0.) and (Xlims[0] == 1.))) and (not self.resize_mode)):
				Ax.set_xlim(Xlims[0], Xlims[1])
				Ax.set_ylim(Ylims[0], Ylims[1])
			else:
				self.resize_mode = False
		else:
			vMin  = np.min(vMin[:nOkZones])*self.Config.FeedbackTransp.Data.RefMin[iFead, 0]
			vMax  = np.max(vMax[:nOkZones])*self.Config.FeedbackTransp.Data.RefMin[iFead, 0]
			if(vMin == vMax):
				vMin = 0.9*vMin
				vMax = 1.1*vMax

			if(self.Shading2DCheckValue.get() == 1):
				shade = 'gouraud'
			else:
				shade = 'flat'
				
			for k in range(nZones):
				if(Zones[k].IsFeedbackValues):	
					zc = np.append(Zones[k].FeedbackValues[iFead], Zones[k].FeedbackValues[iFead][:,-1].reshape(Zones[k].FeedbackValues[iFead].shape[0],1), axis=1)
					zc = np.append(zc, zc[-1,:].reshape(1,zc.shape[1]), axis=0)*FeedbackTransp.Data.RefMin[iFead, 0]
					imBall = Ax.pcolormesh(Zones[k].gridR,  Zones[k].gridZ, zc,  shading=shade, vmin=vMin, vmax=vMax)
					

			if(self.Borders2DCheckValue.get() == 1):
				for k in range(nMagZones):
					EdgeLineType	='k--'
					EdgeLineWidth	= 1
					Ax.plot(MagZones[k].gridR[0,:],  MagZones[k].gridZ[0,:],  EdgeLineType, linewidth=EdgeLineWidth)		#south
					Ax.plot(MagZones[k].gridR[:,-1], MagZones[k].gridZ[:,-1], EdgeLineType, linewidth=EdgeLineWidth)		#East
					Ax.plot(MagZones[k].gridR[-1,:], MagZones[k].gridZ[-1,:], EdgeLineType, linewidth=EdgeLineWidth)		#North
					Ax.plot(MagZones[k].gridR[:,0],  MagZones[k].gridZ[:,0],  EdgeLineType, linewidth=EdgeLineWidth)		#West


			if(self.Contour2DCheckValue.get() == 1):
				for k in range(len(X_points)):
					Ax.contour(Config.r2D, Config.z2D, Config.flux2D, colors="r", linestyles='-', levels=[X_points[k].psi])

			self.Cb = Fig.colorbar(imBall, ax=Ax)
			self.Cb.ax.tick_params(labelsize=8)
			self.Cb.set_label(self.FigureOptions[FigureMode], fontsize=10)
			
			if(Config.wall_OK): plot2d_walls(Ax, Config.Walls, self.ShowPlasmaWallCheckValue.get(), self.ShowEireneWallCheckValue.get(), self.ShowExtraWallCheckValue.get())

				

			self.resize_mode = True

		Ax.set_aspect(1.)
		Fig.canvas.draw()

		
def main():

	root = Tk()
#	root.geometry("1155x675+0+0")
	app = feedback_gen(root)

	if(len(sys.argv) == 2): app.cmd_load_mesh(MeshFile = sys.argv[1])

	root.mainloop()


if __name__ == '__main__':
	main()

