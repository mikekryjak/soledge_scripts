#!/usr/bin/env python
# -*- coding: utf-8 -*-


from tkinter 				import Tk, Text, filedialog, messagebox, BOTH, W, N, E, S, Menu, HORIZONTAL, LEFT, CENTER, RIGHT
from tkinter 				import DoubleVar, IntVar, StringVar, OptionMenu
from tkinter.ttk 			import Entry, Menubutton, Checkbutton, Frame, Button, Label, Style, LabelFrame, Notebook, Treeview, Separator
from tkinter.messagebox 	import Message

try:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2TkAgg as NavigationToolbar2Tk
except:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2Tk

import sys
import os
import types
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot							as pyp
from math 											import floor, sqrt

import numpy										as np
import scipy.interpolate							as interp

from matplotlib.backends.backend_tkagg				import FigureCanvasTkAgg
from matplotlib.figure								import Figure
from matplotlib.widgets 							import Cursor
from matplotlib.backends.backend_pdf				import PdfPages

from interfaces.ask_name							import ask_name
from routines.globals								import DEBUG, MIN_AUTO_VALUE, MAX_AUTO_VALUE
from routines.reset_config							import reset_config
from routines.utils_walls							import plot2d_walls

from files.load_soledge_mesh_file					import load_soledge_mesh_file
from files.save_soledge_mesh_file					import save_soledge_mesh_file
from files.save_custom_save_file					import save_custom_save_file
from files.load_plasma_files						import load_plasma_files
from files.load_custom_save_file					import load_custom_save_file

from mesh.mag_to_zones_pos							import mag_to_zones_pos



class custom_plots(Frame):

	def __init__(self, parent):
		Frame.__init__(self, parent)

		self.parent = parent
		self.initUI()

	def initUI(self):

		nPadx = 5
		nPady = 2
		self.parent.title("Custom plots setting")
		self.pack(fill=BOTH, expand=True)

#		self.columnconfigure(1, weight=1)

		menubar = Menu(self.master)
		filemenu = Menu(menubar, tearoff=0)
		filemenu.add_command(label="Load Mesh", 	 					command=self.cmd_load_mesh)
		filemenu.add_command(label="Save Mesh", 	 					command=self.cmd_save_soledge_mesh)
		filemenu.add_command(label="Save Mesh As..",				command=self.cmd_save_soledge_mesh_as)
		filemenu.add_separator()
		filemenu.add_command(label="Load custom_save",				command=self.cmd_load_custom_save)
		filemenu.add_command(label="Save custom_save",				command=self.cmd_save_custom_save)
		filemenu.add_separator()
		filemenu.add_command(label="Exit",								command=self.quit)
		menubar.add_cascade(label="File",								menu=filemenu)

		OptionsMenu = Menu(menubar, tearoff=0)

		ShowWallMenu = Menu(self.master)

		self.ShowPlasmaWallCheckValue = IntVar()
		self.ShowPlasmaWallCheckValue.set(1)
		ShowWallMenu.add_checkbutton(label="Soledge", 		variable=self.ShowPlasmaWallCheckValue)
		self.ShowPlasmaWallCheckValue.trace("w", self.figure_option_menu_event)

		self.ShowEireneWallCheckValue = IntVar()
		self.ShowEireneWallCheckValue.set(1)
		ShowWallMenu.add_checkbutton(label="Eirene", 		variable=self.ShowEireneWallCheckValue)
		self.ShowEireneWallCheckValue.trace("w", self.figure_option_menu_event)

		self.ShowExtraWallCheckValue = IntVar()
		self.ShowExtraWallCheckValue.set(0)
		ShowWallMenu.add_checkbutton(label="Extra", 		variable=self.ShowExtraWallCheckValue)
		self.ShowExtraWallCheckValue.trace("w", self.figure_option_menu_event)

		OptionsMenu.add_cascade(label="Show wall",				menu=ShowWallMenu)

		menubar.add_cascade(label="Options", 					menu=OptionsMenu)

		helpmenu = Menu(menubar, tearoff=0)
		helpmenu.add_command(label="About...",				command=self.cmd_about)
		menubar.add_cascade(label="Help", menu=helpmenu)

		self.master.config(menu=menubar)

		self.PlotParallel	= []
		self.PlotPerp	 	= []
		self.PlotTemporal	= []
		self.CustomPlots		= [self.PlotParallel, self.PlotPerp, self.PlotTemporal]

		nPadx = 5
		nPady = 3
		lWidth = 6
		eWidth = 10
		xWidth = 10
		cWidth = 15
		
#		Plot area

		PlotFrame = LabelFrame(self, text="Plot")
		PlotFrame.grid(row=0, column=0,  rowspan=4, pady=nPady, padx=nPadx, sticky=E+W+N)

		self.Fig = Figure(figsize=(4.8, 5.8), dpi=100)
		self.Fig.patch.set_facecolor('white')

		# a tk.DrawingArea
		canvas = FigureCanvasTkAgg(self.Fig, master=PlotFrame)
		try:
			canvas.show()										#Matplotlib  up 2.2
		except:
			canvas.draw()
		canvas.get_tk_widget().grid(row=1, column=0, columnspan=4,  padx=nPadx, sticky=E+W+S+N)

		self.Ax = self.Fig.add_subplot(111)
		self.Fig.tight_layout()

		toolbar_frame = Frame(PlotFrame)
		toolbar_frame.grid(row=0, column=0, columnspan=4, padx=nPadx, sticky=E+W+S+N)

		toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
		toolbar.update()

		iRow = 2
		
		ShowLabel = Label(PlotFrame, text="Show", width = lWidth)
		ShowLabel.grid(row=iRow, column=0,  pady=nPady, padx=nPadx)

		self.FigureOptions = ["Mag. Zones", "Split Zones", "Grid", "Cent. Mesh"]
		self.FigureModeOptMenuValue = StringVar()
		self.FigureModeOptMenuValue.set(self.FigureOptions[0])
		self.FigureModeOptMenu = OptionMenu(PlotFrame, self.FigureModeOptMenuValue, *self.FigureOptions)
		self.FigureModeOptMenu.grid(row=iRow, column=1,  pady=nPady, padx=nPadx, sticky=W+E)
		self.FigureModeOptMenuValue.trace("w", self.figure_option_menu_event)


		ShowLabel = Label(PlotFrame, text="Zone sel.", width = lWidth)
		ShowLabel.grid(row=iRow, column=2,  pady=nPady, padx=nPadx)

		self.ZonesOptions = ["Single","Connected"]
		self.ZoneModeOptMenuValue = StringVar()
		self.ZoneModeOptMenuValue.set(self.ZonesOptions[0])
		self.ZoneModeOptMenu = OptionMenu(PlotFrame, self.ZoneModeOptMenuValue, *self.ZonesOptions)
		self.ZoneModeOptMenu.grid(row=iRow, column=3,  pady=nPady, padx=nPadx, sticky=E+W)

#		In/out definition
		
		InOutFrame = LabelFrame(self, text="In & out")
		InOutFrame.grid(row=0, column=1,  pady=nPady, padx=nPadx, sticky=E+W+N)

		iRow = 0

		InLabel = Label(InOutFrame, text="input", width = lWidth)
		InLabel.grid(row=iRow, column=0, columnspan=1,  pady=nPady, padx=nPadx, sticky=W)

		self.InTypeOptions = ["Results", "Evolution", "File"]
		self.InTypeOptMenuValue = StringVar()
		self.InTypeOptMenuValue.set(self.InTypeOptions[0])
		self.InTypeOptMenu = OptionMenu(InOutFrame, self.InTypeOptMenuValue, *self.InTypeOptions)
		self.InTypeOptMenu.grid(row=iRow, column=1, columnspan=1, pady=nPady, padx=nPadx, sticky=E+W)
		self.InTypeOptMenuValue.trace("w", self.in_type_option_menu_event)

		self.EvolutionsEntryValue = StringVar()
		self.EvolutionsEntryValue = StringVar()
		self.EvolutionsEntry = Entry(InOutFrame,  textvariable = self.EvolutionsEntryValue, width = eWidth)
		self.EvolutionsEntry.grid(row=iRow, column=2, pady=nPady, padx=nPadx, sticky=W)
		self.EvolutionsEntry.state(["disabled"])

		iRow += 1

		OutLabel = Label(InOutFrame, text="output", width = lWidth)
		OutLabel.grid(row=iRow, column=0, columnspan=1,  pady=nPady, padx=nPadx, sticky=W)

		self.OutTypeOptions = ["window", "png", "ps", "eps", "pdf"]
		self.OutTypeOptMenuValue = StringVar()
		self.OutTypeOptMenuValue.set(self.OutTypeOptions[0])
		self.OutTypeOptMenu = OptionMenu(InOutFrame, self.OutTypeOptMenuValue, *self.OutTypeOptions)
		self.OutTypeOptMenu.grid(row=iRow, column=1, columnspan=1, pady=nPady, padx=nPadx, sticky=E+W)
#		self.OutTypeOptMenuValue.trace("w", self.plot_type_option_menu_event)

#		Plots definition
		
		PlotSetFrame = LabelFrame(self, text="Plots")
		PlotSetFrame.grid(row=1, column=1,  pady=nPady, padx=nPadx, sticky=E+W+N)

		self.PlotTreeIds = []

		iRow = 0

		PlotTypeLabel = Label(PlotSetFrame, text="Plot type", width = lWidth)
		PlotTypeLabel.grid(row=iRow, column=0, columnspan=2,  pady=nPady, padx=nPadx, sticky=W)

		self.PlotTypeOptions = ["Parallel", "Perpendicular", "Temporal"]
		self.PlotTypeOptMenuValue = StringVar()
		self.PlotTypeOptMenuValue.set(self.PlotTypeOptions[0])
		self.PlotTypeOptMenu = OptionMenu(PlotSetFrame, self.PlotTypeOptMenuValue, *self.PlotTypeOptions)
		self.PlotTypeOptMenu.grid(row=iRow, column=1, columnspan=2, pady=nPady, padx=nPadx, sticky=E+W)
		self.PlotTypeOptMenuValue.trace("w", self.plot_type_option_menu_event)
		self.iPlotType = 0

		iRow += 1

		self.PlotTree = Treeview(PlotSetFrame, height=5)
		self.PlotTree.grid(row=iRow, column=0, columnspan=3,  pady=nPady, padx=nPadx, sticky=E+W+S+N)
		self.tree_view_define(self.PlotTree, ["Label", "Color"], [CENTER, CENTER], [cWidth,int(cWidth/2)], [W, W])
		self.PlotTree.bind("<<TreeviewSelect>>", self.plot_tree_selection_event)

		iRow += 1
		
		AddPlotBut = Button(PlotSetFrame, text="Add", command = self.cmd_add_plot, width = xWidth)
		AddPlotBut.grid(row=iRow, column=0, pady=nPady, padx=nPadx)
		
		RenaPlotBut = Button(PlotSetFrame, text="Rename", command = self.cmd_rename_plot, width = xWidth)
		RenaPlotBut.grid(row=iRow, column=2, pady=nPady, padx=nPadx)

		iRow += 1
		
		AddPlotBut = Button(PlotSetFrame, text="Copy", command = self.cmd_copy_plot, width = xWidth)
		AddPlotBut.grid(row=iRow, column=0, pady=nPady, padx=nPadx)
		
		RenaPlotBut = Button(PlotSetFrame, text="Paste", command = self.cmd_paste_plot, width = xWidth)
		RenaPlotBut.grid(row=iRow, column=2, pady=nPady, padx=nPadx)
		
		DropPlotBut = Button(PlotSetFrame, text="Drop", command = self.cmd_drop_plot, width = xWidth)
		DropPlotBut.grid(row=iRow, column=1, pady=nPady, padx=nPadx)

		iRow += 1
		
		PlotAllPlotBut = Button(PlotSetFrame, text="Plot All", command = self.cmd_plot_all, width = xWidth)
		PlotAllPlotBut.grid(row=iRow, column=0, pady=nPady, padx=nPadx)
		
		PlotPlotBut = Button(PlotSetFrame, text="Plot", command = self.cmd_plot_one, width = xWidth)
		PlotPlotBut.grid(row=iRow, column=1, pady=nPady, padx=nPadx)
		
		WritePlotBut = Button(PlotSetFrame, text="Write", command = self.cmd_write_one, width = xWidth)
		WritePlotBut.grid(row=iRow, column=2, pady=nPady, padx=nPadx)


		nPadx = 5
		nPady = 3
		lWidth = 8
		eWidth = 10
		xWidth = 8
		cWidth = 10

#		Plots setting
		
		PlotSetFrame = LabelFrame(self, text="Plot setting")
		PlotSetFrame.grid(row=0, column=2, rowspan=4, pady=nPady, padx=nPadx, sticky=E+W+N)
		self.PlotSetFrame = PlotSetFrame

		iRow = 0

		PageSettingLabel = Label(PlotSetFrame, text = "Page", width = lWidth)
		PageSettingLabel.grid(row=iRow, column=0,  pady=nPady, padx=nPadx)

		PageSettingSeparator = Separator(PlotSetFrame)
		PageSettingSeparator.grid(row=iRow, column=1,  columnspan=3,  pady=nPady, padx=0, sticky=E+W)

		iRow += 1

		nRowsLabel = Label(PlotSetFrame, text="n Rows", width = lWidth)
		nRowsLabel.grid(row=iRow, column=0,  pady=nPady, padx=nPadx, sticky=W)

		self.nRowsOptions = ["1", "2", "3", "4"]
		self.nRowsOptMenuValue = StringVar()
		self.nRowsOptMenuValue.set(self.nRowsOptions[0])
		self.nRowsOptMenu = OptionMenu(PlotSetFrame, self.nRowsOptMenuValue, *self.nRowsOptions)
		self.nRowsOptMenu.grid(row=iRow, column=1, pady=nPady, padx=nPadx, sticky=E+W)
		self.nRowsOptMenuValue.trace("w", self.nrows_option_menu_event)


		nColsLabel = Label(PlotSetFrame, text="n Cols", width = lWidth)
		nColsLabel.grid(row=iRow, column=2,  pady=nPady, padx=nPadx, sticky=W)

		self.nColsOptions = ["1", "2", "3", "4"]
		self.nColsOptMenuValue = StringVar()
		self.nColsOptMenuValue.set(self.nColsOptions[0])
		self.nColsOptMenu = OptionMenu(PlotSetFrame, self.nColsOptMenuValue, *self.nColsOptions)
		self.nColsOptMenu.grid(row=iRow, column=3, pady=nPady, padx=nPadx, sticky=E+W)
		self.nColsOptMenuValue.trace("w", self.ncols_option_menu_event)

		iRow += 1

		SameXLabel = Label(PlotSetFrame, text = "Same x:", width = lWidth)
		SameXLabel.grid(row=iRow, column=0, pady=nPady, padx=nPadx, sticky=W)

		self.SameXCkButValue = IntVar()
		self.SameXCkButValue.set(1)
		SameXCkBut = Checkbutton(PlotSetFrame, text="", variable = self.SameXCkButValue, width = eWidth)
		SameXCkBut.grid(row=iRow, column=1, padx=nPadx,  pady=nPady, sticky=W)
		self.SameXCkButValue.trace("w", self.samex_option_menu_event)

		iRow += 1

		PlotSettingLabel = Label(PlotSetFrame, text = "Sub plot", width = lWidth)
		PlotSettingLabel.grid(row=iRow, column=0,  pady=nPady, padx=nPadx)

		PlotSettingSeparator = Separator(PlotSetFrame)
		PlotSettingSeparator.grid(row=iRow, column=1,  columnspan=3,  pady=nPady, padx=0, sticky=E+W)

		
		iRow += 1
		
		ParamLabel = Label(PlotSetFrame, text="e & ions", width = lWidth)
		ParamLabel.grid(row=iRow, column=0,  pady=nPady, padx=nPadx, sticky=W)


		self.IonsOptions = ["all single", "e & first single", "e & first togeher", "same A togeher", "all togeher"]
		self.IonsOptMenuValue = StringVar()
		self.IonsOptMenuValue.set(self.IonsOptions[0])
		self.IonsOptMenu = OptionMenu(PlotSetFrame, self.IonsOptMenuValue, *self.IonsOptions)
		self.IonsOptMenu.grid(row=iRow, column=1, columnspan=3, pady=nPady, padx=nPadx, sticky=E+W)

		iRow += 1

		ParLabelsLabel = Label(PlotSetFrame, text = "Par. labels:", width = lWidth)
		ParLabelsLabel.grid(row=iRow, column=0, pady=nPady, padx=nPadx, sticky=W)

		self.ParLabelsCkButValue = IntVar()
		self.ParLabelsCkButValue.set(0)
		ParLabelsCkBut = Checkbutton(PlotSetFrame, text="", variable = self.ParLabelsCkButValue, width = eWidth)
		ParLabelsCkBut.grid(row=iRow, column=1, padx=nPadx,  pady=nPady, sticky=W)

		iRow += 1

		xLabel = Label(PlotSetFrame, text = "x Label:", width = lWidth)
		xLabel.grid(row=iRow, column=0, pady=nPady, padx=nPadx, sticky=W)

		self.xLabelEntryValue = StringVar()
		self.xLabelEntry = Entry(PlotSetFrame,  textvariable = self.xLabelEntryValue, width = eWidth)
		self.xLabelEntry.grid(row=iRow, column=1, pady=nPady, padx=nPadx, sticky=W)


		yLabelLabel = Label(PlotSetFrame, text = "y Label:", width = lWidth)
		yLabelLabel.grid(row=iRow, column=2, pady=nPady, padx=nPadx, sticky=W)

		self.yLabelEntryValue = StringVar()
		self.yLabelEntry = Entry(PlotSetFrame,  textvariable = self.yLabelEntryValue, width = eWidth)
		self.yLabelEntry.grid(row=iRow, column=3,  pady=nPady, padx=nPadx, sticky=W)

		iRow += 1

		xMinLabel = Label(PlotSetFrame, text = "x Min:", width = lWidth)
		xMinLabel.grid(row=iRow, column=0, pady=nPady, padx=nPadx, sticky=W)

		self.xMinEntryValue = StringVar()
		self.xMinEntry = Entry(PlotSetFrame,  textvariable = self.xMinEntryValue, width = eWidth)
		self.xMinEntry.grid(row=iRow, column=1, pady=nPady, padx=nPadx, sticky=W)


		xMaxLabel = Label(PlotSetFrame, text = "x Max:", width = lWidth)
		xMaxLabel.grid(row=iRow, column=2, pady=nPady, padx=nPadx, sticky=W)

		self.xMaxEntryValue = StringVar()
		self.xMaxEntry = Entry(PlotSetFrame,  textvariable = self.xMaxEntryValue, width = eWidth)
		self.xMaxEntry.grid(row=iRow, column=3,  pady=nPady, padx=nPadx, sticky=W)

		iRow += 1

		yMinLabel = Label(PlotSetFrame, text = "y Min:", width = lWidth)
		yMinLabel.grid(row=iRow, column=0, pady=nPady, padx=nPadx, sticky=W)

		self.yMinEntryValue = StringVar()
		self.yMinEntry = Entry(PlotSetFrame,  textvariable = self.yMinEntryValue, width = eWidth)
		self.yMinEntry.grid(row=iRow, column=1, pady=nPady, padx=nPadx, sticky=W)


		yMaxLabel = Label(PlotSetFrame, text = "y Max:", width = lWidth)
		yMaxLabel.grid(row=iRow, column=2, pady=nPady, padx=nPadx, sticky=W)

		self.yMaxEntryValue = StringVar()
		self.yMaxEntry = Entry(PlotSetFrame,  textvariable = self.yMaxEntryValue, width = eWidth)
		self.yMaxEntry.grid(row=iRow, column=3,  pady=nPady, padx=nPadx, sticky=W)

		iRow += 1

		xLogLabel = Label(PlotSetFrame, text = "Log x", width = lWidth)
		xLogLabel.grid(row=iRow, column=0, pady=nPady, padx=nPadx, sticky=W)

		self.xLogCkButValue = IntVar()
		self.xLogCkButValue.set(0)
		self.xLogCkBut = Checkbutton(PlotSetFrame, text="", variable = self.xLogCkButValue, width = eWidth)
		self.xLogCkBut.grid(row=iRow, column=1, padx=nPadx,  pady=nPady, sticky=W)

		xLogLabel = Label(PlotSetFrame, text = "Log y", width = lWidth)
		xLogLabel.grid(row=iRow, column=2, pady=nPady, padx=nPadx, sticky=W)

		self.yLogCkButValue = IntVar()
		self.yLogCkButValue.set(0)
		self.yLogCkBut = Checkbutton(PlotSetFrame, text="", variable = self.yLogCkButValue, width = eWidth)
		self.yLogCkBut.grid(row=iRow, column=3, padx=nPadx,  pady=nPady)

#		SubPlot definitions

		iRow += 1

		ParSettingLabel = Label(PlotSetFrame, text = "SubPlot", width = lWidth)
		ParSettingLabel.grid(row=iRow, column=0,  pady=nPady, padx=nPadx)

		ParSettingSeparator = Separator(PlotSetFrame)
		ParSettingSeparator.grid(row=iRow, column=1,  columnspan=3,  pady=nPady, padx=0, sticky=E+W)
		
		iRow += 1
		
		ParamLabel = Label(PlotSetFrame, text="Param.", width = lWidth)
		ParamLabel.grid(row=iRow, column=0,  pady=nPady, padx=nPadx, sticky=W)


		self.ParamOptions = ["density", "temperature", "velocity", "radiation", "to. rad.", "Zeff",
							"SE", "Sn",
							"Nn", "Nm", "Nti", "Tn", "Tm",
							"Vr Neut.", "Vz Neut.", "Vphi Neut." ,
							"Vr Mol.", "Vz Mol.", "Vphi Mol."]
		self.ParamSignals = ["Dens", "Temp", "velocity", "IRad", "TRad", "Zeff",
							"SE", "Sn",
							"Nn", "Nm", "Nti", "Tn", "Tm",
							"vxn", "vzn", "vyn",
							"vxm", "vzm", "vym"]
		self.ParamLabels = ["$n\ (*10^{19}\ m^{-3})$", "$T\ (eV)$", "$V\ (m/s)$", "$Rad_{Ion}\ (kW/m^3)$", "$Rad_{tot}\ (kW/m^3)$", "Zeff",
							"$S_E\ (kW*m^{-3}s^{-1})$", "$S_n\ (*10^{19}\ m^{-3}s^{-1})$",
							"$n_{neutrals}\ (*10^{19}\ m^{-3})$", "$n_{mols}\ (*10^{19}\ m^{-3})$", "Nti\ (*10^{19}\ m^{-3})", "$T_{neutrals}\ (eV)$", "$T_{mols}\ (eV)$",
							"$V_R^{neutrals}\ (*10^3\ m/s)$", "$V_Z^{neutrals}\ (*10^3\ m/s)$", "$V_{phi}^{neutrals}\ (*10^3\ m/s)$", 
							"$V_R^{mols}\ (*10^3\ m/s)$", "$V_Z^{mols}\ (*10^3\ m/s)$", "$V_{phi}^{mols}\ (*10^3\ m/s)$"]
		self.ParamFact	 = [1e-19, 1., 1., 1e-3, 1e-3, 1.,
							1e-3, 1e-19,
							1e-19, 1e-19, 1., 1., 1.,
							1e-3, 1e-3, 1e-3,
							1e-3, 1e-3, 1e-3]
		self.ParamOptMenuValue = StringVar()
		self.ParamOptMenuValue.set(self.ParamOptions[0])
		self.ParamOptMenu = OptionMenu(PlotSetFrame, self.ParamOptMenuValue, *self.ParamOptions)
		self.ParamOptMenu.grid(row=iRow, column=1, columnspan=3, pady=nPady, padx=nPadx, sticky=E+W)
#		self.ParamOptMenuValue.trace("w", self.param_option_menu_event)

		iRow += 1

		LineLabel = Label(PlotSetFrame, text="line", width = lWidth)
		LineLabel.grid(row=iRow, column=0,  pady=nPady, padx=nPadx, sticky=W)

		self.LineOptions = ["solid", "dashed", "dashdot", "dotted", "none"]
		self.LineTypeOptMenuValue = StringVar()
		self.LineTypeOptMenuValue.set(self.LineOptions[0])
		self.LineTypeOptMenu = OptionMenu(PlotSetFrame, self.LineTypeOptMenuValue, *self.LineOptions)
		self.LineTypeOptMenu.grid(row=iRow, column=1, pady=nPady, padx=nPadx, sticky=E+W)

		ColorLabel = Label(PlotSetFrame, text="color", width = lWidth)
		ColorLabel.grid(row=iRow, column=2,  pady=nPady, padx=nPadx, sticky=W)

		self.ColorOptions = ["auto", "blue", "green", "cyan", "magenta", "yellow", "black"]
		self.LineColorOptMenuValue = StringVar()
		self.LineColorOptMenuValue.set(self.ColorOptions[0])
		self.LineColorOptMenu = OptionMenu(PlotSetFrame, self.LineColorOptMenuValue, *self.ColorOptions)
		self.LineColorOptMenu.grid(row=iRow, column=3, pady=nPady, padx=nPadx, sticky=E+W)

		iRow += 1

		ColSpanTree = 4
		self.ParametersTree = Treeview(PlotSetFrame, height=10)
		self.ParametersTree.grid(row=iRow, column=0, columnspan=ColSpanTree,  pady=nPady, padx=nPadx, sticky=E+W+S+N)
		self.tree_view_define(self.ParametersTree, ["SubPlot"], [CENTER], [cWidth,], [W])
		self.ParametersTree.bind("<<TreeviewSelect>>", self.parameters_tree_selection_event)

		iRow += 1
		
		UpParBut = Button(PlotSetFrame, text="Up", command = self.cmd_up_parameter, width = xWidth)
		UpParBut.grid(row=iRow, column=2, pady=nPady, padx=nPadx)
		
		DownParBut = Button(PlotSetFrame, text="Down", command = self.cmd_down_parameter, width = xWidth)
		DownParBut.grid(row=iRow, column=3, pady=nPady, padx=nPadx)

		iRow += 1
			
		AddParBut = Button(PlotSetFrame, text="Add Sub.", command = self.cmd_add_subplot, width = xWidth)
		AddParBut.grid(row=iRow, column=0, pady=nPady, padx=nPadx)
			
		AddParBut = Button(PlotSetFrame, text="Add Par.", command = self.cmd_add_parameter, width = xWidth)
		AddParBut.grid(row=iRow, column=1, pady=nPady, padx=nPadx)
		
		DropParBut = Button(PlotSetFrame, text="Drop", command = self.cmd_drop_parameter, width = xWidth)
		DropParBut.grid(row=iRow, column=2, pady=nPady, padx=nPadx)
		
		DropParBut = Button(PlotSetFrame, text="Update", command = self.cmd_update_parameter, width = xWidth)
		DropParBut.grid(row=iRow, column=3, pady=nPady, padx=nPadx)

		"""
		iRow += 1
		
		SetBut = Button(PlotSetFrame, text="Set", command = self.cmd_aset_parameters, width = xWidth*2+nPadx)
		SetBut.grid(row=iRow, column=0, columnspan=2, pady=nPady, padx=nPadx)
		"""

		self.yMinEntryValue.set("")
		self.yMaxEntryValue.set("")
		self.ParametersTreeIds	= []
		self.iPlotNum 		= -1
		self.AddType		= -1
		self.old_plot_type	= -1
		self.CopyiPlotType	= -1
		self.CopyiSub		= -1	

		self.skip_tree_selection_event = False
#		self.param_option_menu_event = False
		self.MarkersList = ["o","+","x"]
		self.ColorsList	 = ['Red','Magenta','Yellow','Orange','Pink','Blue','Green','Cyan','Grey','Purple']

		
	def cmd_about(self):
		print("To do cmd_about")
		
#=================
#	Load mesh data
#=================

	def cmd_load_mesh(self, File=None):
#	load mesh previously saved

		if(File is None):
			File = filedialog.askopenfilename(filetypes=(("Mesh files", "*.h5"), ("All files", "*.*") ))
			if(not File): return

#		try:
		if(True):
			self.Config = load_soledge_mesh_file(File)

			if(self.Config.CustomPlots_OK):
				self.CustomPlots = self.Config.CustomPlots
				self.PlotTreeIds = self.tree_view_fill_plots(self.PlotTree, self.PlotTreeIds, self.CustomPlots[self.iPlotType])
			
			self.resize_mode 	= True
			self.update_figure()

			self.MeshFile		= File
			self.MeshDir  		= os.path.dirname(self.MeshFile)
			if((len(self.MeshDir)==0) or ((len(self.MeshDir)==1) and  (self.MeshDir == "/"))): self.MeshDir="./"
		else:
#		except:
			messagebox.showerror("Open Source File", "Failed to read mesh file:{:s}\n".format(File))

	def cmd_save_soledge_mesh(self):
		if(not self.MeshFile): return

		nPlots = 0
		for iType in range(len(self.CustomPlots)): nPlots += len(self.CustomPlots[iType])

		if(nPlots > 0):
			self.Config.CustomPlots = self.CustomPlots
			self.Config.CustomPlots_OK = True
		else:
			self.Config.CustomPlots = []
			self.Config.CustomPlots_OK = False

#		try:
		if(True):
			save_soledge_mesh_file(self.MeshFile, self.Config)
#		except:
#			messagebox.showerror("Save Mesh file", "Failed to write mesh file:{:s}\n".format(self.MeshFile))


	def cmd_save_soledge_mesh_as(self):
		if(not self.MeshFile): return	
		MeshFile = filedialog.asksaveasfilename(initialdir=self.MeshDir, filetypes=(("Mesh files", "*.h5"), ("All files", "*.*") ))

		if(len(MeshFile) > 0):
			try:
				OldMeshFile	  = self.MeshFile
				self.MeshFile = MeshFile
				self.cmd_save_soledge_mesh()
				self.MeshDir  = os.path.dirname(self.MeshFile)
			except:
				self.MeshFile = OldMeshFile
				messagebox.showerror("Save Mesh file", "Failed to write mesh file:{:s}\n".format(self.MeshFile))

	def cmd_load_custom_save(self):
		if(not self.MeshFile): return

		if(os.path.split(self.MeshFile)[1] == "mesh.h5"):
			CustomPath = self.MeshDir+"/"
		else:
			CustomPath = self.MeshDir+"/"+os.path.splitext(os.path.split(self.MeshFile)[1])[0]+"/"
			if(not os.path.isdir(CustomPath[:-1])): CustomPath = self.MeshDir+"/"

		CustomPlots=load_custom_save_file(CustomPath)
		if(len(CustomPlots) == 3):
			Zones=self.Config.Zones
			for iType in range(3):
				for iPlot in range(len(CustomPlots[iType])):
					if(iType != 2):
						kMagZones=[Zones[CustomPlots[iType][iPlot].kZones[0]].magz[0]]
						for kZone  in CustomPlots[iType][iPlot].kZones:
							if(kMagZones[-1] != Zones[kZone].magz[0]): kMagZones.append(Zones[kZone].magz[0])
						CustomPlots[iType][iPlot].kMagZones = np.array(kMagZones)
						CustomPlots[iType][iPlot].MagCoords  =np.zeros(3, dtype='i4')
						CustomPlots[iType][iPlot].MagCoords[0] = Zones[CustomPlots[iType][iPlot].kZones[0]].magz[0]
						CustomPlots[iType][iPlot].MagCoords[iType+1] = CustomPlots[iType][iPlot].Coords[iType+1] + Zones[CustomPlots[iType][iPlot].kZones[0]].magz[iType+1]
					else:
						k = CustomPlots[iType][iPlot].Coords[0]
						i  = CustomPlots[iType][iPlot].Coords[1]
						j  = CustomPlots[iType][iPlot].Coords[2]
						CustomPlots[iType][iPlot].MagCoords=np.array([Zones[k].magz[0],Zones[k].magz[1]+i,Zones[k].magz[2]+j])
						CustomPlots[iType][iPlot].kMagZones=np.array([Zones[k].magz[0]])

			self.CustomPlots=CustomPlots

			self.Config.CustomPlots = self.CustomPlots
			self.Config.CustomPlots_OK = True
			self.Config.CustomPlots = self.CustomPlots
			self.PlotTreeIds = self.tree_view_fill_plots(self.PlotTree, self.PlotTreeIds, self.CustomPlots[self.iPlotType])

			self.resize_mode 	= True
			self.update_figure()
		else:
			messagebox.showerror("Open Source File", "Failed to read custom_save.txt\n")


	def cmd_save_custom_save(self):
		if(not self.MeshFile): return

		nPlots = self.set_zones_coords()
		if(os.path.split(self.MeshFile)[1] == "mesh.h5"):
			CustomPath = self.MeshDir+"/"
		else:
			CustomPath = self.MeshDir+"/"+os.path.splitext(os.path.split(self.MeshFile)[1])[0]+"/"
			if(not os.path.isdir(CustomPath[:-1])): CustomPath = self.MeshDir+"/"

		save_custom_save_file(CustomPath, self.CustomPlots)
		
	

# ######################################
#	Plot frame
# ######################################

	def figure_option_menu_event(self,a,b,c):
		self.update_figure()


# ######################################
#	In/Out frame
# ######################################

		
	def in_type_option_menu_event(self,a,b,c):
		if(self.InTypeOptMenuValue.get() == self.InTypeOptions[1]): self.EvolutionsEntry.state(["!disabled"])
		else:														self.EvolutionsEntry.state(["disabled"])


# ######################################
#	Plot type frame
# ######################################
		
	def plot_type_option_menu_event(self,a,b,c):
		for i in range(len(self.PlotTypeOptions)):
			if(self.PlotTypeOptMenuValue.get() == self.PlotTypeOptions[i]):
				iPlotType = i
				break
		if(self.iPlotType != iPlotType):
			self.PlotTreeIds = self.tree_view_fill_plots(self.PlotTree, self.PlotTreeIds, self.CustomPlots[iPlotType])
			self.iPlotType = iPlotType
			self.iPlotNum  = -1
			self.init_page_data()
			

	def cmd_add_plot(self):
		self.AddType = self.iPlotType
		self.cid_mouse_button_press_add_point = self.Fig.canvas.mpl_connect('button_press_event', self.button_press_add_point)
		self.cursor = Cursor(self.Ax, useblit=False, color='red', linewidth=1)
		return

	def cmd_rename_plot(self):
		iPlotType = self.iPlotType
		Id = self.PlotTree.focus()
		if(len(Id)>0):
			ind = self.PlotTreeIds.index(Id)
			ParName=ask_name(self, "Label", self.CustomPlots[iPlotType][ind].Name)
			if(len(ParName) > 0):
				self.CustomPlots[iPlotType][ind].Name = ParName
			self.PlotTreeIds = self.tree_view_fill_plots(self.PlotTree, self.PlotTreeIds,  self.CustomPlots[iPlotType])
		return

	def cmd_copy_plot(self):
		Id = self.PlotTree.focus()
		if(len(Id)>0):
			self.CopyiSub = self.PlotTreeIds.index(Id)
			self.CopyiPlotType = self.iPlotType
		else:
			self.CopyiPlotType	= -1
			self.CopyiSub		= -1	
		return

	def cmd_paste_plot(self):
		if((self.CopyiPlotType < 0) or (self.CopyiSub < 0)): return

		iPlotType = self.iPlotType
		Id = self.PlotTree.focus()
		if(len(Id)>0):
			ind = self.PlotTreeIds.index(Id)
			self.CustomPlots[iPlotType][ind].nRows = self.CustomPlots[self.CopyiPlotType][self.CopyiSub].nRows
			self.CustomPlots[iPlotType][ind].nCols = self.CustomPlots[self.CopyiPlotType][self.CopyiSub].nCols
			self.CustomPlots[iPlotType][ind].SameX = self.CustomPlots[self.CopyiPlotType][self.CopyiSub].SameX

			PasteSubPlots	= self.CustomPlots[self.CopyiPlotType][self.CopyiSub].SubPlots
			SubPlots		= []
			for i  in range(len(PasteSubPlots)):
				SubPlots.append(self.copy_subplot_data(PasteSubPlots[i]))
				SubPlots[i].Params = []
				for k in range(len(PasteSubPlots[i].Params)): SubPlots[i].Params.append(self.copy_parameter_data(PasteSubPlots[i].Params[k]))
			self.CustomPlots[iPlotType][ind].SubPlots =SubPlots

			self.PlotTreeIds = self.tree_view_fill_plots(self.PlotTree, self.PlotTreeIds,  self.CustomPlots[iPlotType])
			self.update_figure()
		return

	def cmd_add_subplot(self):
		if((self.iPlotType == -1) or (self.iPlotNum == -1)): return

		SubPlots   = self.CustomPlots[self.iPlotType][self.iPlotNum].SubPlots
		SubPlots.append(types.SimpleNamespace())
		SubPlots[-1].Params = [types.SimpleNamespace()]
		self.get_subplot_data(SubPlots[-1])
		self.get_parameter_data(SubPlots[-1].Params[0])

		self.CustomPlots[self.iPlotType][self.iPlotNum].IdSubPlotPar = self.tree_view_fill_parameters(SubPlots)
		return


	def cmd_drop_plot(self):
		iPlotType = self.iPlotType
		Id = self.PlotTree.focus()
		if(len(Id)>0):
			ind = self.PlotTreeIds.index(Id)
			
			del self.CustomPlots[iPlotType][ind]
			self.PlotTreeIds = self.tree_view_fill_plots(self.PlotTree, self.PlotTreeIds,  self.CustomPlots[iPlotType])
			self.update_figure()
		return

	def cmd_plot_all(self):
		iPlotType = self.iPlotType
		Plasmas = load_plasma_files("/")
		for iPlotNum in range(len(self.CustomPlots[iPlotType])): self.plot_parameters(self.CustomPlots[iPlotType][iPlotNum], Plasmas)
		return

	def cmd_plot_one(self):
		iPlotType = self.iPlotType
		Id = self.PlotTree.focus()
		if(len(Id) > 0):
			iPlotNum  = self.PlotTreeIds.index(Id)
			self.plot_parameters(self.CustomPlots[iPlotType][iPlotNum], [])
		return

	def cmd_write_one(self):
		iPlotType = self.iPlotType
		Id = self.PlotTree.focus()
		if(len(Id) > 0):
			iPlotNum  = self.PlotTreeIds.index(Id)
			Plasmas = load_plasma_files("/")
			self.write_parameters(self.CustomPlots[iPlotType][iPlotNum], Plasmas)
		return

	def plot_tree_selection_event(self, event):
		iPlotType = self.iPlotType
		Id = self.PlotTree.focus()
		if(len(Id) > 0):
			self.set_frame_state(self.PlotSetFrame, "!disabled")
			self.iPlotType = iPlotType
			self.iPlotNum  = self.PlotTreeIds.index(Id)
			self.init_page_data()
		else:
			self.iPlotNum  = -1
		return


# ######################################
#	Parameters actions commands
# ######################################
		
	def nrows_option_menu_event(self,a,b,c):
		if((self.iPlotType > -1) and (self.iPlotNum > -1)):
			self.CustomPlots[self.iPlotType][self.iPlotNum].nRows = eval(self.nRowsOptMenuValue.get())

	def ncols_option_menu_event(self,a,b,c):
		if((self.iPlotType > -1) and (self.iPlotNum > -1)):
			self.CustomPlots[self.iPlotType][self.iPlotNum].nCols = eval(self.nColsOptMenuValue.get())

	def samex_option_menu_event(self,a,b,c):
		if((self.iPlotType > -1) and (self.iPlotNum > -1)):
			self.CustomPlots[self.iPlotType][self.iPlotNum].SameX =self.SameXCkButValue.get()

	"""
	def param_option_menu_event(self,a,b,c):
		if(self.param_option_menu_event):
			self.param_option_menu_event = False
			return
		self.clear_parameter()
		self.skip_tree_selection_event = True
	"""

	def cmd_add_subplot(self):
		if((self.iPlotType == -1) or (self.iPlotNum == -1)): return

		SubPlots   = self.CustomPlots[self.iPlotType][self.iPlotNum].SubPlots
		SubPlots.append(types.SimpleNamespace())
		SubPlots[-1].Params = [types.SimpleNamespace()]
		self.get_subplot_data(SubPlots[-1])
		self.get_parameter_data(SubPlots[-1].Params[0])

		self.CustomPlots[self.iPlotType][self.iPlotNum].IdSubPlotPar = self.tree_view_fill_parameters(SubPlots)
		return

	def cmd_add_parameter(self):
		if((self.iPlotType == -1) or (self.iPlotNum == -1)): return

		SubPlots = self.CustomPlots[self.iPlotType][self.iPlotNum].SubPlots
		SubPlots[-1].Params.append(types.SimpleNamespace())
		self.get_parameter_data(SubPlots[-1].Params[-1])

		self.CustomPlots[self.iPlotType][self.iPlotNum].IdSubPlotPar = self.tree_view_fill_parameters(SubPlots)
		return

		
	def cmd_drop_parameter(self):
		if((self.iPlotType == -1) or (self.iPlotNum == -1)): return

		Id = self.ParametersTree.focus()
		if(len(Id)>0):
			iRow = self.ParametersTreeIds.index(Id)
			iParNum	   = self.CustomPlots[self.iPlotType][self.iPlotNum].IdSubPlotPar[iRow][0]
			iSubParNum = self.CustomPlots[self.iPlotType][self.iPlotNum].IdSubPlotPar[iRow][1]
			if(iSubParNum == -1):	del self.CustomPlots[self.iPlotType][self.iPlotNum].SubPlots[iParNum]
			else:					del self.CustomPlots[self.iPlotType][self.iPlotNum].SubPlots[iParNum].Params[iSubParNum]
			self.CustomPlots[self.iPlotType][self.iPlotNum].IdSubPlotPar = self.tree_view_fill_parameters(self.CustomPlots[self.iPlotType][self.iPlotNum].SubPlots)
		return
		
	def cmd_update_parameter(self):
		if((self.iPlotType == -1) or (self.iPlotNum == -1)): return

		Id = self.ParametersTree.focus()
		if(len(Id)>0):
			iRow = self.ParametersTreeIds.index(Id)
			iParNum	= self.CustomPlots[self.iPlotType][self.iPlotNum].IdSubPlotPar[iRow][0]
			iSubParNum = self.CustomPlots[self.iPlotType][self.iPlotNum].IdSubPlotPar[iRow][1]
			self.get_subplot_data(self.CustomPlots[self.iPlotType][self.iPlotNum].SubPlots[iParNum])
			if(iSubParNum == -1): self.get_parameter_data(self.CustomPlots[self.iPlotType][self.iPlotNum].SubPlots[iParNum])
			else:				  self.get_parameter_data(self.CustomPlots[self.iPlotType][self.iPlotNum].SubPlots[iParNum].Params[iSubParNum])
			self.CustomPlots[self.iPlotType][self.iPlotNum].IdSubPlotPar = self.tree_view_fill_parameters(self.CustomPlots[self.iPlotType][self.iPlotNum].SubPlots)
		return

	def cmd_up_parameter(self):
		if((self.iPlotType == -1) or (self.iPlotNum == -1)): return

		Id = self.ParametersTree.focus()
		if(len(Id)>0):
			iRow = self.ParametersTreeIds.index(Id)
			iParNum	   = self.CustomPlots[self.iPlotType][self.iPlotNum].IdSubPlotPar[iRow][0]
			iSubParNum = self.CustomPlots[self.iPlotType][self.iPlotNum].IdSubPlotPar[iRow][1]
			if(iParNum == -1):	
				SubPlots =  self.CustomPlots[self.iPlotType][self.iPlotNum].SubPlots
			else:
				SubPlots =  self.CustomPlots[self.iPlotType][self.iPlotNum].SubPlots[iParNum].Params
				iParNum  = iSubParNum

			if(iParNum > 0):
				SubPlot = SubPlots[iParNum-1]
				SubPlots[iParNum-1] = SubPlots[iParNum]
				SubPlots[iParNum]   = SubPlot
				self.CustomPlots[self.iPlotType][self.iPlotNum].IdSubPlotPar = self.tree_view_fill_parameters(SubPlots)
		return

	def cmd_down_parameter(self):
		if((self.iPlotType == -1) or (self.iPlotNum == -1)): return

		Id = self.ParametersTree.focus()
		if(len(Id)>0):
			iRow = self.ParametersTreeIds.index(Id)
			iParNum	   = self.CustomPlots[self.iPlotType][self.iPlotNum].IdSubPlotPar[iRow][0]
			iSubParNum = self.CustomPlots[self.iPlotType][self.iPlotNum].IdSubPlotPar[iRow][1]
			if(iParNum == -1):	
				SubPlots	 =  self.CustomPlots[self.iPlotType][self.iPlotNum].SubPlots
			else:
				SubPlots	 =  self.CustomPlots[self.iPlotType][self.iPlotNum].SubPlots[iParNum].Params
				iParNum  = iSubParNum

			if(iParNum < len(SubPlots)-1):
				SubPlot			= SubPlots[iParNum+1]
				SubPlots[iParNum+1]	= SubPlots[iParNum]
				SubPlots[iParNum]		= SubPlot
				self.CustomPlots[self.iPlotType][self.iPlotNum].IdSubPlotPar = self.tree_view_fill_parameters(SubPlots)
		return

	def parameters_tree_selection_event(self, event):
		if(self.skip_tree_selection_event):
			self.skip_tree_selection_event = False
			return

		Id = self.ParametersTree.focus()
		if(len(Id)>0):
			iRow = self.ParametersTreeIds.index(Id)
			iParNum	= self.CustomPlots[self.iPlotType][self.iPlotNum].IdSubPlotPar[iRow][0]
			iSubParNum = self.CustomPlots[self.iPlotType][self.iPlotNum].IdSubPlotPar[iRow][1]
			self.set_subplot_data(self.CustomPlots[self.iPlotType][self.iPlotNum].SubPlots[iParNum])
			if(iSubParNum != -1): self.set_parameter_data(self.CustomPlots[self.iPlotType][self.iPlotNum].SubPlots[iParNum].Params[iSubParNum])
		else:
			self.iParNum = -1
		return

	def init_page_data(self):
		iPlotType = self.iPlotType
		iPlotNum  = self.iPlotNum
		if(iPlotNum > -1):
			self.nRowsOptMenuValue.set(self.nRowsOptions[self.CustomPlots[iPlotType][iPlotNum].nRows-1])
			self.nColsOptMenuValue.set(self.nColsOptions[self.CustomPlots[iPlotType][iPlotNum].nCols-1])
			self.SameXCkButValue.set(self.CustomPlots[iPlotType][iPlotNum].SameX)
			self.CustomPlots[iPlotType][iPlotNum].IdSubPlotPar = self.tree_view_fill_parameters(self.CustomPlots[iPlotType][iPlotNum].SubPlots)
		else:
			self.nRowsOptMenuValue.set(self.nRowsOptions[0])
			self.nColsOptMenuValue.set(self.nColsOptions[0])
			self.SameXCkButValue.set(1)
			IdSubPlotPar = self.tree_view_fill_parameters([])
		self.init_subplot_data()


	def init_subplot_data(self):
		self.IonsOptMenuValue.set(self.IonsOptions[0])

		self.ParLabelsCkButValue.set(0)

		self.xLabelEntryValue.set("")
		self.yLabelEntryValue.set("")

		self.xMinEntryValue.set("")
		self.xMaxEntryValue.set("")

		self.yMinEntryValue.set("")
		self.yMaxEntryValue.set("")

		self.xLogCkButValue.set(0)
		self.yLogCkButValue.set(0)
		self.init_parameter_data()
		return

	def init_parameter_data(self):
		self.ParamOptMenuValue.set(self.ParamOptions[0])

		self.LineTypeOptMenuValue.set(self.LineOptions[0])
		self.LineColorOptMenuValue.set(self.ColorOptions[0])
		return

	"""
	def clear_sub_plot_data(self):
		Id = self.ParametersTree.focus()
		if(len(Id) > 0):
			self.yMinEntryValue.set("")
			self.yMaxEntryValue.set("")
			self.xLogCkButValue.set(0)
			self.yLogCkButValue.set(0)
	"""

	def set_subplot_data(self, SubPlot):
#		self.param_option_menu_event = True

		self.IonsOptMenuValue.set(SubPlot.Ions)

		self.ParLabelsCkButValue.set(SubPlot.ParLabels)

		self.xMinEntryValue.set(SubPlot.xLabel)
		self.xMinEntryValue.set(SubPlot.yLabel)

		self.xMinEntryValue.set(SubPlot.xMin)
		self.xMaxEntryValue.set(SubPlot.xMax)

		self.yMinEntryValue.set(SubPlot.yMin)
		self.yMaxEntryValue.set(SubPlot.yMax)

		self.xLogCkButValue.set(SubPlot.xLog)
		self.yLogCkButValue.set(SubPlot.yLog)


	def set_parameter_data(self, Parameter):
		self.ParamOptMenuValue.set(Parameter.Name)
		self.LineTypeOptMenuValue.set(Parameter.LineType)
		self.LineColorOptMenuValue.set(Parameter.LineColor)
		return

	def get_subplot_data(self, SubPlot):
		SubPlot.Ions = self.IonsOptMenuValue.get()

		SubPlot.ParLabels = self.ParLabelsCkButValue.get()

		SubPlot.xLabel = self.xLabelEntryValue.get()
		SubPlot.yLabel = self.yLabelEntryValue.get()

		SubPlot.xMin = self.xMinEntryValue.get()
		SubPlot.xMax = self.xMaxEntryValue.get()

		SubPlot.yMin = self.yMinEntryValue.get()
		SubPlot.yMax = self.yMaxEntryValue.get()

		SubPlot.xLog = self.xLogCkButValue.get()
		SubPlot.yLog = self.yLogCkButValue.get()

	def copy_subplot_data(self, OldSubPlot):
		SubPlot = types.SimpleNamespace()

		SubPlot.Ions		= OldSubPlot.Ions
		SubPlot.ParLabels 	= OldSubPlot.ParLabels
		SubPlot.xLabel		= OldSubPlot.xLabel
		SubPlot.yLabel		= OldSubPlot.yLabel
		SubPlot.xMin		= OldSubPlot.xMin
		SubPlot.xMax		= OldSubPlot.xMax
		SubPlot.yMin		= OldSubPlot.yMin
		SubPlot.yMax		= OldSubPlot.yMax
		SubPlot.xLog		= OldSubPlot.xLog
		SubPlot.yLog		= OldSubPlot.yLog
		return SubPlot

	def get_parameter_data(self, Parameter):
		Parameter.Name = self.ParamOptMenuValue.get()
		Parameter.LineType  = self.LineTypeOptMenuValue.get()
		Parameter.LineColor = self.LineColorOptMenuValue.get()
		return

	def copy_parameter_data(self, OldParameter):
		Parameter = types.SimpleNamespace()

		Parameter.Name		= OldParameter.Name
		Parameter.LineType  = OldParameter.LineType
		Parameter.LineColor = OldParameter.LineColor
		return Parameter


	def tree_view_fill_parameters(self, SubPlots):
		ParNames = []
		IdSubPlotPar = []
		for k in range(len(SubPlots)):
			ParNames.append("======= SubPlot # {:d} ==========".format(k+1)[:26])
			IdSubPlotPar.append([k, -1])
			for ko in range(len(SubPlots[k].Params)):
				IdSubPlotPar.append([k, ko])
				ParNames.append("    " + SubPlots[k].Params[ko].Name)
		self.ParametersTreeIds = self.tree_view_fill_params(self.ParametersTree, self.ParametersTreeIds, ParNames)
		return 	IdSubPlotPar
	
	def safe_get_entry_value(self, String, DefValue):
		if(len(String) > 0):
			try:
				return eval(String)
			except:
				return DefValue
		else:
			return DefValue
	
	def safe_set_entry_value(self, Value, DefValue):
		if(Value != DefValue):
			return "{:g}".format(Value)
		else:
			return ""

# ######################################
#	general routines
# ######################################

	def tree_view_define(self, Tree, Headings, HeadAnchors, xWidths, TextAnchors):
		col_nums=()
		for k in range(len(Headings)):  col_nums += ("{:d}".format(k),)
		Tree["columns"] = col_nums
		Tree['show'] = 'headings'
		for k in range(len(Headings)): 
			Tree.column(col_nums[k],  width=xWidths[k], anchor=TextAnchors[k])
			Tree.heading(col_nums[k], text = Headings[k], anchor=HeadAnchors[k])

	def tree_view_fill_plots(self, Tree, TreeIds, Values):
		if(len(TreeIds) > 0):
			Tree.delete(*TreeIds)
			TreeIds = []
			
		for k in range(len(Values)):
			TreeValues = (Values[k].Name,)	
			TreeValues += (self.ColorsList[k],)	
			TreeIds.append(Tree.insert("" , "end",    text="", values=TreeValues))

		self.set_frame_state(self.PlotSetFrame, "disabled")
		return TreeIds


	def tree_view_fill_params(self, Tree, TreeIds, Values):
		if(len(TreeIds) > 0):
			Tree.delete(*TreeIds)
			TreeIds = []
			
		for k in range(len(Values)):
			TreeValues = (Values[k],)	
			TreeIds.append(Tree.insert("" , "end",    text="", values=TreeValues))
		return TreeIds

	def button_press_add_point(self, event):
		if(self.cid_mouse_button_press_add_point != 0):
			self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press_add_point)
			self.cid_mouse_button_press_add_point = 0
			self.cursor=0

			Config	= self.Config																#Find nearest center mesh
			MagZones	= self.Config.MagZones
			nMagZones	= len(MagZones)
			kMin	= -1
			iMIn	= -1
			jMin	= -1
			dMin	= 1e30
			Rp		= event.xdata
			Zp		= event.ydata
			for k in range(nMagZones):
				d = np.sqrt((Rp - MagZones[k].gridRc)**2 + (Zp - MagZones[k].gridZc)**2)
				(iiMin, jjMin) = np.unravel_index(np.argmin(d), MagZones[k].gridRc.shape)
				ddMin = d[iiMin, jjMin]
				if(ddMin < dMin):
					dMin = ddMin
					kMin = k
					iMin = iiMin
					jMin = jjMin

			MagCoords = np.array([kMin, iMin, jMin])

		iPlotType = self.iPlotType
		PlotLabels = [ "parallel_{:d}", "perpendicular_{:d}", "temporal_{:d}"]

		self.CustomPlots[iPlotType].append(types.SimpleNamespace())
		self.CustomPlots[iPlotType][-1].Name			= PlotLabels[iPlotType].format(len(self.CustomPlots[iPlotType]))
		self.CustomPlots[iPlotType][-1].SubPlots		= [] 
		self.CustomPlots[iPlotType][-1].IdSubPlotPar	= [] 
		self.CustomPlots[iPlotType][-1].nRows			= 1 
		self.CustomPlots[iPlotType][-1].nCols			= 1 
		self.CustomPlots[iPlotType][-1].SameX			= 1 
		self.CustomPlots[iPlotType][-1].MagCoords 	= MagCoords

		if((self.ZoneModeOptMenuValue.get() != self.ZonesOptions[0]) and (iPlotType != 2)): 
			if(iPlotType == 0): self.CustomPlots[iPlotType][-1].kMagZones = self.get_west_east_zones(MagCoords[0], self.Config.MagZones) 
			else:				self.CustomPlots[iPlotType][-1].kMagZones = self.get_south_north_zones(MagCoords[0], self.Config.MagZones)
		else:
			self.CustomPlots[iPlotType][-1].kMagZones = np.array([MagCoords[0]])

		self.PlotTreeIds = self.tree_view_fill_plots(self.PlotTree, self.PlotTreeIds, self.CustomPlots[iPlotType])
		self.iPlotNum  = -1
		self.init_page_data()

		self.update_figure()


	def get_south_north_zones(self, k, Zones, kMagZones=[]):
		MagZones = self.Config.MagZones

		kZones = np.array([k])
		South = Zones[k].Neighbour.south
		while (South >= 0):
			if((len(kMagZones) == 1) and (len(np.where(South == MagZones[kMagZones[0]].list)[0]) == 0)): break		
			kZones = np.append(South, kZones)
			South = Zones[South].Neighbour.south

		North = Zones[k].Neighbour.north
		while (North >= 0):
			if((len(kMagZones) == 1) and (len(np.where(North == MagZones[kMagZones[0]].list)[0]) == 0)): break		
			kZones = np.append(kZones, North)
			North	  = Zones[North].Neighbour.north

		return kZones


	def get_west_east_zones(self, k, Zones, kMagZones=[]):
		MagZones = self.Config.MagZones

		kZones = np.array([k])
		West = Zones[k].Neighbour.west
		while (West >= 0):
			if((len(kMagZones) == 1) and (len(np.where(West == MagZones[kMagZones[0]].list)[0]) == 0)): break
			kZones = np.append(West, kZones)
			West = Zones[West].Neighbour.west
			if(West == k): West = -100

		if(West > -100):
			East = Zones[k].Neighbour.east
			while (East >= 0):
				if((len(kMagZones) == 1) and (len(np.where(East == MagZones[kMagZones[0]].list)[0]) == 0)): break
				kZones = np.append(kZones, East)
				East = Zones[East].Neighbour.east
		return kZones


	def set_zones_coords(self):
		nPlots = 0
		for iType in range(len(self.CustomPlots)):
			for iPlot in range(len(self.CustomPlots[iType])):
				nPlots += 1
				self.set_zones_coords_plot(self.CustomPlots[iType][iPlot], iType)

		return nPlots

	def set_zones_coords_plot(self, CustomPlot, iType):
		CustomPlot.Coords = mag_to_zones_pos(self.Config, CustomPlot.MagCoords)
		if(iType == 0):   CustomPlot.kZones = self.get_west_east_zones(CustomPlot.Coords[0],   self.Config.Zones, kMagZones=CustomPlot.kMagZones)
		elif(iType == 1): CustomPlot.kZones = self.get_south_north_zones(CustomPlot.Coords[0], self.Config.Zones, kMagZones=CustomPlot.kMagZones)
		else:			  CustomPlot.kZones = np.array([CustomPlot.Coords])
			

	def set_frame_state(self, frame, state):
		frame.state(["disabled"])
		for child in frame.winfo_children(): 
			try:
				child.state([state])
			except:
				if(state == "disabled"):	child.configure(state="disabled")
				else:						child.configure(state="active")


# ######################################
#	Plot & write routine
# ######################################
		

#	Plot routine
#	=======

	def plot_parameters(self, CustomPlot, Plasmas):
		iPlotType = self.iPlotType
		iPlotNum  = self.iPlotNum

		save = self.OutTypeOptMenuValue.get()

		if(save == "pdf"):	pdf = PdfPages("custom_plots.".format(save))   #pdf in one file only

		if(len(Plasmas) == 0): Plasmas = load_plasma_files("/")

		SubPlots  = CustomPlot.SubPlots
		iPlotVar, nPlots = self.get_plasma_var(CustomPlot.SubPlots, Plasmas)			#get parameters index
		Dist, Rc, Zc, Chi, ic, jc = self.compute_mesh_dist(CustomPlot)					#compute mesh distance and coordinate
		iIn = np.where(Chi == 0)[0];
		Rc = 0; Zc = 0; Chi = 0

		if(iIn[0] == 0):
			Dist = Dist[iIn]- Dist[iIn[0]]
		else:
			iIn  = iIn[::-1]
			Dist = -(Dist[iIn] - Dist[iIn[0]])

		nRows		= CustomPlot.nRows
		nCols		= CustomPlot.nCols
		PlotPerFig	= nRows*nCols
		nFigs		= int(nPlots/PlotPerFig)
		if(nFigs*PlotPerFig < nPlots): nFigs += 1

		Fig = []
		Ax  = []
		for i in range(nFigs):
			Fig.append(pyp.figure())
			for k in range(min(PlotPerFig,nPlots-i*PlotPerFig)):
				if(((i == 0) and (k == 0)) or (CustomPlot.SameX == 0)):	Ax.append(Fig[-1].add_subplot(nRows,nCols,k+1))
				else:													Ax.append(Fig[-1].add_subplot(nRows,nCols,k+1, sharex = Ax[0]))

				Ax[-1].locator_params(axis='x',nbins=4)

			Fig[-1].tight_layout(pad=2., w_pad=3., h_pad=3.)

		for figure in Fig:  figure.patch.set_facecolor('white')


		nLinesOnAxes   = np.zeros(nPlots, dtype='i4')

		OrderByIon = 1
		if(OrderByIon == 0):
			print("TO DOOOOOOOOOOOOOOOO")
		else:
			for iSub in range(len(iPlotVar)):
				for iVal in range(len(iPlotVar[iSub])):
					iPar		= iPlotVar[iSub][iVal][0]
					iParOpt		= iPlotVar[iSub][iVal][1]
					iPlasma		= iPlotVar[iSub][iVal][2]
					iValue		= iPlotVar[iSub][iVal][3]
					iPlot		= iPlotVar[iSub][iVal][4]-1
					FirstLine	= iPlotVar[iSub][iVal][5]
					nLines		= iPlotVar[iSub][iVal][6]
					nSamePlasma	= iPlotVar[iSub][iVal][7]
					nSamePar	= iPlotVar[iSub][iVal][8]
					Values = np.empty(0, dtype = 'f8')

					for iz in range(len(CustomPlot.kZones)):
						kz = CustomPlot.kZones[iz]
						if(Plasmas[iPlasma][kz].Nx == Plasmas[iPlasma][kz].Values[iValue].shape[0]):
							if(iPlotType == 0):	Values = np.append(Values, Plasmas[iPlasma][kz].Values[iValue][ic,:])
							else:				Values = np.append(Values, Plasmas[iPlasma][kz].Values[iValue][:,jc])
						else:
							if(iPlotType == 0):	Values = np.append(Values, Plasmas[iPlasma][kz].Values[iValue][ic,1:-1])
							else:				Values = np.append(Values, Plasmas[iPlasma][kz].Values[iValue][1:-1,jc])

					if(nLines > 1):
						if(nSamePlasma == nLines):
							LineLabel = self.ParamLabels[iParOpt]
						else:
							if(nSamePar == nLines): LineLabel = Plasmas[iPlasma][0].ion
							else:				    LineLabel = self.ParamLabels[iParOpt] + " " + Plasmas[iPlasma][0].ion

					if(SubPlots[iSub].Params[iPar].LineColor == "auto"):
						if(SubPlots[iSub].ParLabels):	Ax[iPlot].plot(Dist, Values[iIn]*self.ParamFact[iParOpt], linestyle = SubPlots[iSub].Params[iPar].LineType, color=self.ColorsList[int(min(len(self.ColorsList)-1,nLinesOnAxes[iPlot-1]))], label=LineLabel)
						else:							Ax[iPlot].plot(Dist, Values[iIn]*self.ParamFact[iParOpt], linestyle = SubPlots[iSub].Params[iPar].LineType, color=self.ColorsList[int(min(len(self.ColorsList)-1,nLinesOnAxes[iPlot-1]))])
					else:
						if(SubPlots[iSub].ParLabels):	Ax[iPlot].plot(Dist, Values[iIn]*self.ParamFact[iParOpt], linestyle = SubPlots[iSub].Params[iPar].LineType, color=SubPlots[iSub].Params[iPar].LineColor, label=LineLabel)
						else:							Ax[iPlot].plot(Dist, Values[iIn]*self.ParamFact[iParOpt], linestyle = SubPlots[iSub].Params[iPar].LineType, color=SubPlots[iSub].Params[iPar].LineColor)

					nLinesOnAxes[iPlot-1] += 1

					if(nLinesOnAxes[iPlot-1] == nLines):
						SamePlasma 	= iPlotVar[iSub][iVal][6]
						Ax[iPlot].autoscale(enable=True, axis='both', tight=True)

						if((SubPlots[iSub].ParLabels) and (nLines > 1)):	Ax[iPlot].legend(fontsize='small', loc='upper right')
						if(nSamePar == nLines):	Ax[iPlot].set_title(Plasmas[iPlasma][0].ion)

						Ax[iPlot].set_xlabel("d (m)")
						if(nSamePar == nLines): 		Ax[iPlot].set_ylabel(self.ParamLabels[iParOpt])

						if(SubPlots[iSub].xLog == 0):	Ax[iPlot].set_yscale('linear')
						else:							Ax[iPlot].set_yscale('log')
						if(SubPlots[iSub].yLog == 0):	Ax[iPlot].set_yscale('linear')
						else:							Ax[iPlot].set_yscale('log')

						if(len(SubPlots[iSub].xMin) != 0): Ax[iPlot].set_xlim(left = eval(SubPlots[iSub].xMin))
						if(len(SubPlots[iSub].xMax) != 0): Ax[iPlot].set_xlim(right = eval(SubPlots[iSub].xMax))

						if(len(SubPlots[iSub].yMin) != 0): Ax[iPlot].set_ylim(bottom = eval(SubPlots[iSub].yMin))
						if(len(SubPlots[iSub].yMax) != 0): Ax[iPlot].set_ylim(top = eval(SubPlots[iSub].yMax))

		if(save != "window"):
			i_plot_file	= 0
			for i in range(len(Fig)):
				i_plot_file += 1
				Fig[i].set_size_inches(10.05,7.44)
				if(save == "pdf"):
					pdf.savefig(Fig[i])
				else:
					Fig[i].savefig("custom_plots_{:d}.{:s}".format(i_plot_file, save))

			pyp.show(block=False)
			pyp.close()
		else:
			pyp.show()

		if(save == "pdf"):	pdf.close()


#	Write routine
#	=======

	def write_parameters(self, CustomPlot, Plasmas):
		iPlotType = self.iPlotType
		iPlotNum  = self.iPlotNum

		if(len(Plasmas) == 0): Plasmas = load_plasma_files("/")

		SubPlots  = CustomPlot.SubPlots
		iPlotVar, nPlots = self.get_plasma_var(CustomPlot.SubPlots, Plasmas)			#get parameters index
		Dist, Rc, Zc, Chi, ic, jc = self.compute_mesh_dist(CustomPlot)					#compute mesh distance and coordinate
		iIn = np.where(Chi == 0)[0]; Chi = 0
		if(iIn[0] == 0):
			Dist = Dist[iIn]- Dist[iIn[0]]
		else:
			iIn  = iIn[::-1]
			Dist = -(Dist[iIn] - Dist[iIn[0]])

		Rc   = Rc[iIn]
		Zc	 = Zc[iIn]

		nWrite = 0
		for iSub in range(len(iPlotVar)):
			for iVal in range(len(iPlotVar[iSub])): nWrite += 1

		save_cvs	  = np.empty((len(Dist), nWrite+3), dtype = 'f8')
		save_cvs[:,0] = Dist
		save_cvs[:,1] = Rc
		save_cvs[:,2] = Zc
		Header		  = "dist, Rc, Zc"


		nWrite = 2
		OrderByIon = 1
		if(OrderByIon == 0):
			print("TO DOOOOOOOOOOOOOOOO")
		else:
			for iSub in range(len(iPlotVar)):
				for iVal in range(len(iPlotVar[iSub])):
					nWrite += 1
					iPar	= iPlotVar[iSub][iVal][0]
					iParOpt = iPlotVar[iSub][iVal][1]
					iPlasma = iPlotVar[iSub][iVal][2]
					iValue  = iPlotVar[iSub][iVal][3]
					iPlot  	= iPlotVar[iSub][iVal][4]-1
					Values = np.empty(0, dtype = 'f8')
					for iz in range(len(CustomPlot.kZones)):
						kz = CustomPlot.kZones[iz]
						if(Plasmas[iPlasma][kz].Nx == Plasmas[iPlasma][kz].Values[iValue].shape[0]):
							if(iPlotType == 0):	Values = np.append(Values, Plasmas[iPlasma][kz].Values[iValue][ic,:])
							else:				Values = np.append(Values, Plasmas[iPlasma][kz].Values[iValue][:,jc])
						else:
							if(iPlotType == 0):	Values = np.append(Values, Plasmas[iPlasma][kz].Values[iValue][ic,1:-1])
							else:				Values = np.append(Values, Plasmas[iPlasma][kz].Values[iValue][1:-1,jc])

					save_cvs[:,nWrite] = Values[iIn]
#					Header		  	  += "," + Plasmas[iPlasma][0].VNames[iValue]
					if(self.ParamOptions[iParOpt][-1] != "."): ParName = self.ParamOptions[iParOpt]
					else:										ParName = self.ParamOptions[iParOpt][:-1]

					if(iPlasma == 0): 	ParName += "_e"
					elif(iPlasma == 1):	ParName += "_i"
					else:				ParName += SigName + Plasmas[iPlasma][0].ion[:-1] 

					Header		  	  += "," + ParName
					
		np.savetxt(CustomPlot.Name+".csv", save_cvs, header=Header, delimiter=",", fmt="%15.7e", comments="")
				
#=====================
#	Get plasma Indexes
#=====================

	def get_plasma_var(self, SubPlots, Plasmas):

		nPlots = 0
		iPlotVar   = []
		iPlotPlasmas = np.zeros(len(Plasmas), dtype = 'i4')
		for iSub in range(len(SubPlots)):
			iPlotVar.append([])
			iPlotPlasmas[:] = 0
			for iPar in range(len(SubPlots[iSub].Params)):
				Name = SubPlots[iSub].Params[iPar].Name
				index = self.indexes_upper(self.ParamOptions, Name)
				if(len(index) == 1):
					iParOpt = index[0]
					SigName = self.ParamSignals[iParOpt]
				else:
					print("Error in code: zero or more than one signal for one parameter")
					print("Parameters[iPar].Name=",Name)
					print("self.ParamSignal=",self.ParamSignal)
					print("index=",index)
					exit()

				Found_in_0		= False
				for iPlasma in range(len(Plasmas)):
					if(iPlasma == 0): 	SigName_ion = SigName + "e"
					elif(iPlasma == 1):	SigName_ion = SigName + "i"
					else:				SigName_ion = SigName + Plasmas[iPlasma][0].ion[:-1] 

					index = self.indexes_upper(Plasmas[iPlasma][0].VNames, SigName_ion)

					NewVar  = False
					NewPlot = 0
					if(len(index) == 1):
						iValue = index[0]

						if(SubPlots[iSub].Ions == self.IonsOptions[0]):				#all single
							NewVar  = True
							NewPlot = 1

						elif(SubPlots[iSub].Ions == self.IonsOptions[1]):			#e & first single
							if(iPlasma < 2):
								NewVar  = True
								NewPlot = 1

						elif(SubPlots[iSub].Ions == self.IonsOptions[2]):			#e & first togeher
							if(iPlasma < 2):
								NewVar  = True
								if((iPlasma == 0) or (not Found_in_0)):
									NewPlot = 1
									Found_in_0 = True

						elif(SubPlots[iSub].Ions == self.IonsOptions[3]):			#same A togeher
							NewVar  = True
							if((iPlasma < 2) or (Plasmas[iPlasma][0].mass != Plasmas[iPlasma-1][0].mass)): NewPlot = 1

						elif(SubPlots[iSub].Ions == self.IonsOptions[4]):			#all togeher
							NewVar  = True
							if(iPlasma == 0): NewPlot = 1

						if(NewVar):
							if(iPlotPlasmas[iPlasma] != 0):
								iPlot	= iPlotPlasmas[iPlasma]
								NewPlot = 0
								SamePlasma = 1
							else:
								nPlots = nPlots + NewPlot
								iPlotPlasmas[iPlasma] = nPlots
								SamePlasma = NewPlot

							iPlotVar[-1].append([iPar, iParOpt, iPlasma, iValue, iPlotPlasmas[iPlasma], NewPlot, 0, SamePlasma, 0])

		nLines 	    = np.zeros(nPlots, dtype = 'i4')
		nSamePlasma = np.zeros(nPlots, dtype = 'i4')
		nSamePar	= np.zeros(nPlots, dtype = 'i4')
		iParOptOld	= -np.ones(nPlots, dtype = 'i4')
		for iSub in range(len(iPlotVar)):
			for iVal in range(len(iPlotVar[iSub])):
				iParOpt		= iPlotVar[iSub][iVal][1]
				iPlot  		= iPlotVar[iSub][iVal][4] -1
				FirstLine	= iPlotVar[iSub][iVal][5]
				SamePlasma  = iPlotVar[iSub][iVal][7]
				nLines[iPlot]	   += 1 
				nSamePlasma[iPlot] += SamePlasma
				if(FirstLine == 1): 
					nSamePar[iPlot]	  = 1
					iParOptOld[iPlot] = iParOpt
				elif(iParOptOld[iPlot] == iParOpt): nSamePar[iPlot] += 1
				else:								iParOptOld[iPlot] = -1

		for iSub in range(len(iPlotVar)):
			for iVal in range(len(iPlotVar[iSub])):
				iPlot  		= iPlotVar[iSub][iVal][4] -1
				iPlotVar[iSub][iVal][6] = nLines[iPlot]
				iPlotVar[iSub][iVal][7] = nSamePlasma[iPlot]
				iPlotVar[iSub][iVal][8] = nSamePar[iPlot]


#		iPlotVar[plot_index][[plasma_index, variable_index]]
#		iPlotPos[plot_index][plasma_index]=number of plot

		return iPlotVar, nPlots

	def indexes_upper(self, list_strs, str):
		return  [i for i, j in enumerate(list_strs) if j.upper() == str.upper()]
			
#=====================
#	Get mesh distance
#=====================

	def compute_mesh_dist(self, CustomPlot):
		iPlotType = self.iPlotType

		self.set_zones_coords_plot(CustomPlot, iPlotType)			#define zones parameters

#		Compute distance

		Zones = self.Config.Zones

		Rc  = np.empty(0, dtype = 'f8')
		Zc  = np.empty(0, dtype = 'f8')
		Chi = np.empty(0, dtype = 'i4')
		ic = CustomPlot.Coords[1]
		jc = CustomPlot.Coords[2]
		for iz in range(len(CustomPlot.kZones)):
			kz = CustomPlot.kZones[iz]
			if(iPlotType == 0):
				Rc  = np.append(Rc,  Zones[kz].gridRc[ic,:])
				Zc  = np.append(Zc,  Zones[kz].gridZc[ic,:])
				Chi = np.append(Chi, Zones[kz].Chi[ic,:])
			elif(iPlotType == 1):
				Rc  = np.append(Rc,  Zones[kz].gridRc[:,jc])
				Zc  = np.append(Zc,  Zones[kz].gridZc[:,jc])
				Chi = np.append(Chi, Zones[kz].Chi[:,jc])

		Dist = np.append(0.,np.cumsum(np.sqrt((Rc[1:]-Rc[:-1])**2 + (Zc[1:]-Zc[:-1])**2)))

		return Dist, Rc, Zc, Chi, ic, jc
		

			
#==============
#	Main figure 
#==============

#update main figure

	def update_figure(self):

		Config			= self.Config
		MagZones		= self.Config.MagZones
		Zones			= self.Config.Zones
		Fig				= self.Fig
		Ax				= self.Ax

		nMagZones		= len(MagZones)
		nZones			= len(Zones)
		
		Ax.clear()
	
		Xlims = Ax.get_xlim()
		Ylims = Ax.get_ylim()
		
		Ax.tick_params(axis='both', which='major', labelsize=8)
		Ax.autoscale(axis='both', tight=True)
		Ax.axis('on')

		PlotWall = True
		Ax.set_aspect(1.)
		
		if(self.FigureModeOptMenuValue.get() == self.FigureOptions[0]):				#Mag zone
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
				Ax.text(R,Z,"{:d}".format(k+1), color='green', size = "x-small")

		elif(self.FigureModeOptMenuValue.get() == self.FigureOptions[1]):											#Splitted sizes
			for k in range(len(Zones)):
				EdgeLineType	='b-'
				EdgeLineWidth	= 1
				Ax.plot(Zones[k].gridR[0,:],  Zones[k].gridZ[0,:],  EdgeLineType, linewidth=EdgeLineWidth)		#south
				Ax.plot(Zones[k].gridR[:,-1], Zones[k].gridZ[:,-1], EdgeLineType, linewidth=EdgeLineWidth)		#East
				Ax.plot(Zones[k].gridR[-1,:], Zones[k].gridZ[-1,:], EdgeLineType, linewidth=EdgeLineWidth)		#North
				Ax.plot(Zones[k].gridR[:,0],  Zones[k].gridZ[:,0],  EdgeLineType, linewidth=EdgeLineWidth)		#West
	
				iCenter = floor(Zones[k].gridR.shape[1]/2)
				R=(Zones[k].gridR[0,iCenter]+Zones[k].gridR[-1,iCenter])*0.5
				Z=(Zones[k].gridZ[0,iCenter]+Zones[k].gridZ[-1,iCenter])*0.5
				Ax.text(R,Z,"{:d}".format(k+1), color='green', size = "x-small")
	
		elif(self.FigureModeOptMenuValue.get() == self.FigureOptions[2]):				#Show grid
			EdgeLineType	='b-'
			EdgeLineWidth	= 1
			for k in range(nMagZones):
				if(MagZones[k].meshortho):
					Ax.plot(MagZones[k].gridR,   MagZones[k].gridZ,   'g-')
					Ax.plot(MagZones[k].gridR.T, MagZones[k].gridZ.T, 'g-')
				else:
					Ax.plot(MagZones[k].gridR,   MagZones[k].gridZ,   '-', color=[1.,0.8,0.2])
					Ax.plot(MagZones[k].gridR.T, MagZones[k].gridZ.T, '-', color=[1.,0.8,0.2])

				Ax.plot(MagZones[k].gridR[0,:],  MagZones[k].gridZ[0,:],  EdgeLineType, linewidth=EdgeLineWidth)		#south
				Ax.plot(MagZones[k].gridR[:,-1], MagZones[k].gridZ[:,-1], EdgeLineType, linewidth=EdgeLineWidth)		#East
				Ax.plot(MagZones[k].gridR[-1,:], MagZones[k].gridZ[-1,:], EdgeLineType, linewidth=EdgeLineWidth)		#North
				Ax.plot(MagZones[k].gridR[:,0],  MagZones[k].gridZ[:,0],  EdgeLineType, linewidth=EdgeLineWidth)		#West
			
		elif(self.FigureModeOptMenuValue.get() == self.FigureOptions[3]):			#Cent. Mesh
			for k in range(nZones):
				Ax.plot(Zones[k].gridR,   Zones[k].gridZ,   'b-')
				Ax.plot(Zones[k].gridR.T, Zones[k].gridZ.T, 'b-')
				iIn, jIn = np.where(Zones[k].Chi == 0)
				Ax.plot(Zones[k].gridRc[iIn, jIn],  Zones[k].gridZc[iIn, jIn],  'r.')
				iOut, jOut = np.where(Zones[k].Chi == 1)
				Ax.plot(Zones[k].gridRc[iOut, jOut],  Zones[k].gridZc[iOut, jOut],  'y.')
				iIn = 0; jIn = 0; iOut = 0; jOut = 0

		else:
			print("ERROR invalid option check code")

#			if(self.align_psimode):
#				psimin = eval(self.PsiMinEntryValue.get())
#				psimax = eval(self.PsiMaxEntryValue.get())

#				Ax.contour(Config.r2D, Config.z2D, Config.flux2D, levels = [psimin, psimax], colors ='m', linestyles='dashed')
#				Ax.contour(Config.r2D, Config.z2D, Config.flux2D, levels = [psimin, psimax], colors ='m', linestyles='dashed')
		if(Config.wall_OK and PlotWall): plot2d_walls(Ax, Config.Walls, self.ShowPlasmaWallCheckValue.get(), self.ShowEireneWallCheckValue.get(), self.ShowExtraWallCheckValue.get())

		self.markers = ["o","+","x"]
		for iType in range(len(self.CustomPlots)):
			for iPlot in range(len(self.CustomPlots[iType])):
				k = self.CustomPlots[iType][iPlot].MagCoords[0]
				i = self.CustomPlots[iType][iPlot].MagCoords[1]
				j = self.CustomPlots[iType][iPlot].MagCoords[2]
				for iz in range(len(self.CustomPlots[iType][iPlot].kMagZones)):
					kz = self.CustomPlots[iType][iPlot].kMagZones[iz]
					if(iType == 0):
						if(iz == 0):
							Ax.plot(MagZones[kz].gridRc[i,:], MagZones[kz].gridZc[i,:], "-", color=self.ColorsList[iPlot])
						else:
							Ax.plot(np.append(Re, MagZones[kz].gridRc[i,:]), np.append(Ze, MagZones[kz].gridZc[i,:]), "-", color=self.ColorsList[iPlot])
						Re = MagZones[kz].gridRc[i,-1]
						Ze = MagZones[kz].gridZc[i,-1]

					elif(iType == 1):
						if(iz == 0):
							Ax.plot(MagZones[kz].gridRc[:,j], MagZones[kz].gridZc[:,j], "-", color=self.ColorsList[iPlot])
						else:
							Ax.plot(np.append(Re, MagZones[kz].gridRc[:,j]), np.append(Ze, MagZones[kz].gridZc[:,j]), "-", color=self.ColorsList[iPlot])
						Re = MagZones[kz].gridRc[-1,j]
						Ze = MagZones[kz].gridZc[-1,j]


				Ax.plot(MagZones[k].gridRc[i,j], MagZones[k].gridZc[i,j], self.MarkersList[iType], color=self.ColorsList[iPlot])
				

		if((not ((Xlims[0] == 0.) and (Xlims[1] == 1.))) and (not self.resize_mode)):
			Ax.set_xlim(Xlims[0], Xlims[1])
			Ax.set_ylim(Ylims[0], Ylims[1])
		else:
			self.resize_mode = False
	
		Fig.canvas.draw()

		
def main():

	root = Tk()
#	root.geometry("1200x675+0+0")
	app = custom_plots(root)

	if(len(sys.argv) == 2): app.cmd_load_mesh(File = sys.argv[1])

	root.mainloop()


if __name__ == '__main__':
	main()

