#!/usr/bin/env python
# -*- coding: utf-8 -*-


from tkinter				import Tk, Text, OptionMenu, filedialog, messagebox, BOTH, W, N, E, S, Menu, HORIZONTAL
from tkinter				import DoubleVar, IntVar, StringVar
from tkinter.messagebox 	import Message
from tkinter.ttk			import Entry, Menubutton, Checkbutton, Frame, Button, Label, Style, LabelFrame, Scale, Treeview

import types
import os

import numpy								as np
from numpy									import arange, sin, pi
from math									import floor

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg		import FigureCanvasTkAgg
from matplotlib.figure						import Figure
from matplotlib.widgets 					import Cursor
try:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2TkAgg as NavigationToolbar2Tk
except:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2Tk

from routines.globals						import *
from routines.smooth_fields					import smooth_fields 
from routines.bezier						import bezier
from routines.part_contour					import part_contour
from routines.contour_better				import contour_better
from routines.intersect_contour				import intersect_contour
from routines.find_closest_segment			import find_closest_segment
from routines.reset_config					import reset_config
from routines.utils_walls					import plot2d_walls
from routines.find_core						import find_core


from mesh.plot_ortho_lines					import plot_ortho_lines
from mesh.getArcLim							import getArcLim
from mesh.extrapolate_psi_out				import extrapolate_psi_out
from mesh.define_frontier					import define_frontier

from mesh.find_zone_corners					import find_zone_corners
from mesh.find_segment_theta				import find_segment_theta
from mesh.find_cote_psi						import find_cote_psi
from mesh.remove_double_cote				import remove_double_cote
from mesh.find_zones						import find_zones
from mesh.analyze_topology					import analyze_topology
from mesh.find_megazones_outside_plasma		import find_megazones_outside_plasma
from mesh.get_magzone_of_point				import get_magzone_of_point
from mesh.delete_magzone					import delete_magzone
from mesh.define_OMP_segment				import define_OMP_segment
from mesh.surround_Xpoints					import surround_Xpoints

from interfaces.mesh_gen					import mesh_gen
from interfaces.ask_values					import ask_values

class zone_gen(Frame):

	def __init__(self, parent, root):
		Frame.__init__(self, parent)

		self.root = root

		self.MaxXpoints = 6
		self.MaxBariers = 6
		nPadx = 5
		nPady = 3

		self.pack(fill=BOTH, expand=True)
		self.columnconfigure(1, weight=1)

#		Plot area

		PlotFrame = Frame(self)
		PlotFrame.grid(row=1, column=0, rowspan=5,  pady=nPady, padx=nPadx, sticky=E+W+N)

		self.Fig = Figure(figsize=(4.8, 5.8), dpi=100)
		self.Fig.patch.set_facecolor('white')

		# a tk.DrawingArea
		canvas = FigureCanvasTkAgg(self.Fig, master=PlotFrame)
		try:
			canvas.show()										#Matplotlib  up 2.2
		except:
			canvas.draw()
		canvas.get_tk_widget().grid(row=1, column=0,  padx=nPadx, sticky=E+W+S+N)

		self.Ax = self.Fig.add_subplot(111)
		self.Fig.tight_layout()


		toolbar_frame = Frame(PlotFrame)
		toolbar_frame.grid(row=0, column=0, padx=nPadx, sticky=E+W+S+N)
		toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
		toolbar.update()

#	Extrapolation
#	=======

		xWidth = 12

		FieldsFrame = LabelFrame(self, text="Fields & Flux")

		FieldsFrame.grid(row=1, column=1, columnspan=2, pady=nPady, padx=nPadx, sticky=E+W+N)
		iFrow = 0
		NRLabel = Label(FieldsFrame, text="n. points R:", width = xWidth)
		NRLabel.grid(row=iFrow,column=0,  pady=nPady, padx=nPadx)

		self.NREntryValue = StringVar()
		NREntry = Entry(FieldsFrame,  textvariable = self.NREntryValue, width = xWidth)
		NREntry.grid(row=iFrow,column=1,  pady=nPady, padx=nPadx)

		NZLabel = Label(FieldsFrame, text="n. points Z:", width = xWidth)
		NZLabel.grid(row=iFrow,column=2,  pady=nPady, padx=nPadx)

		self.NZEntryValue = StringVar()
		NZEntry = Entry(FieldsFrame,  textvariable = self.NZEntryValue, width = xWidth)
		NZEntry.grid(row=iFrow,column=3,  pady=nPady, padx=nPadx)

		OutlengthLabel = Label(FieldsFrame, text="Out length (m):", width = xWidth)
		OutlengthLabel.grid(row=iFrow,column=4,  pady=nPady, padx=nPadx)

		self.OutLengthEntryValue = StringVar()
		OutlengthEntry = Entry(FieldsFrame,  textvariable = self.OutLengthEntryValue, width = xWidth)
		OutlengthEntry.grid(row=iFrow,column=5,  pady=nPady, padx=nPadx)

		iFrow += 1

		SmoothLabel = Label(FieldsFrame, text="Smooth inp (%):", width = xWidth)
		SmoothLabel.grid(row=iFrow,column=0,  pady=nPady, padx=nPadx)

		self.SmoothInEntryValue = StringVar()
		SmoothInEntry = Entry(FieldsFrame,  textvariable = self.SmoothInEntryValue, width = xWidth)
		SmoothInEntry.grid(row=iFrow,column=1,  pady=nPady, padx=nPadx)

		SmoothLabel = Label(FieldsFrame, text="Smooth ext (%):", width = xWidth)
		SmoothLabel.grid(row=iFrow,column=2, pady=nPady, padx=nPadx)

		self.SmoothExtEntryValue = StringVar()
		SmoothEntry = Entry(FieldsFrame,  textvariable = self.SmoothExtEntryValue, width = xWidth)
		SmoothEntry.grid(row=iFrow,column=3,  pady=nPady, padx=nPadx)

		self.ExtFacLabelValue = StringVar()
		self.ExtFacLabelValue.set("Ext. (0.00):")
		ExtLabel = Label(FieldsFrame, textvariable = self.ExtFacLabelValue, width=xWidth)
		ExtLabel.grid(row=iFrow,column=4,  pady=nPady, padx=nPadx)

		self.ExtScale = Scale(FieldsFrame, from_ = 0., to = 1., command = self.cmd_move_ext_scale)
		self.ExtScale.grid(row=iFrow,column=5, pady=nPady, padx=nPadx)

		iFrow += 1

		RaiseFactLabel = Label(FieldsFrame, text="Raise factor:", width = xWidth)
		RaiseFactLabel.grid(row=iFrow,column=0,  pady=nPady, padx=nPadx)

		RaiseNSLabel = Label(FieldsFrame, text="N S", width = xWidth)
		RaiseNSLabel.grid(row=iFrow,column=1,  pady=nPady, padx=nPadx)

		self.RaiseFacNEntryValue = StringVar()
		self.RaiseFacNEntry = Entry(FieldsFrame,  textvariable = self.RaiseFacNEntryValue, width = xWidth)
		self.RaiseFacNEntry.grid(row=iFrow,column=2,  pady=nPady, padx=nPadx)

		self.RaiseFacSEntryValue = StringVar()
		self.RaiseFacSEntry = Entry(FieldsFrame,  textvariable = self.RaiseFacSEntryValue, width = xWidth)
		self.RaiseFacSEntry.grid(row=iFrow,column=3,  pady=nPady, padx=nPadx)
		self.RaiseFacSEntry.state(["disabled"])

		RaisePowerLabel = Label(FieldsFrame, text="Raise power:", width = xWidth)
		RaisePowerLabel.grid(row=iFrow,column=4,  pady=nPady, padx=nPadx)

		self.RaisePowerOptions = ["1", "2", "3","4","5"]
		self.RaisePowerOptMenuValue = StringVar()
		self.RaisePowerOptMenuValue.set(self.RaisePowerOptions[0])
		self.RaisePowerOptMenu = OptionMenu(FieldsFrame, self.RaisePowerOptMenuValue, *self.RaisePowerOptions)
		self.RaisePowerOptMenu.grid(row=iFrow, column=5,  pady=0, sticky=W+E)

		iFrow += 1
		RaiseEWLabel = Label(FieldsFrame, text="E W", width = xWidth)
		RaiseEWLabel.grid(row=iFrow,column=1,  pady=nPady, padx=nPadx)

		self.RaiseFacEEntryValue = StringVar()
		self.RaiseFacEEntry = Entry(FieldsFrame,  textvariable = self.RaiseFacEEntryValue, width = xWidth)
		self.RaiseFacEEntry.grid(row=iFrow,column=2,  pady=nPady, padx=nPadx)
		self.RaiseFacEEntry.state(["disabled"])

		self.RaiseFacWEntryValue = StringVar()
		self.RaiseFacWEntry = Entry(FieldsFrame,  textvariable = self.RaiseFacWEntryValue, width = xWidth)
		self.RaiseFacWEntry.grid(row=iFrow,column=3,  pady=nPady, padx=nPadx)
		self.RaiseFacWEntry.state(["disabled"])


		self.RaiseFacEntryValues = [self.RaiseFacNEntryValue, self.RaiseFacSEntryValue, self.RaiseFacEEntryValue, self.RaiseFacWEntryValue]

		iFrow += 1
		self.SmoothExtBut = Button(FieldsFrame, text="Smooth & Extrapolate", command = self.cmd_smooth_extrapolate, width = 2*xWidth)
		self.SmoothExtBut.grid(row=iFrow, column=0,  columnspan=2, pady=nPady)

#		Show area
#		=======


		ShowFrame = LabelFrame(self, text="Psi")
		ShowFrame.grid(row=2, column=1, columnspan=2,  pady=nPady, padx=nPadx, sticky=E+W+N)

		iFrow = 0

		self.FluxCkButValue = IntVar()
		self.FluxCkBut = Checkbutton(ShowFrame, text="Show surf.", variable = self.FluxCkButValue, command = self.cmd_flux_ck_button, width = xWidth)
		self.FluxCkBut.grid(row=iFrow,column=0,  pady=nPady, padx=nPadx, sticky=E+W+N)

		self.FluxScaleValue = DoubleVar()
		self.FluxScale = Scale(ShowFrame, orient= HORIZONTAL, from_ = 0., to_ = 100., \
								variable = self.FluxScaleValue, command = self.cmd_move_flux_scale)
		self.FluxScale.grid(row=iFrow,column=1, pady=nPady, padx=nPadx, sticky=E+W+N)

		self.FluxEntryValue = StringVar()
		self.FluxEntry = Entry(ShowFrame,  textvariable = self.FluxEntryValue, width = xWidth)
		self.FluxEntry.grid(row=iFrow,column=2,  pady=nPady, padx=nPadx, sticky=E+W+N)
		self.FluxEntry.bind("<Return>", self.cmd_enter_flux_entry)

		self.BordXptDetCkButValue = IntVar()
		self.BordXptDetCkBut = Checkbutton(ShowFrame, text="Xp details", variable =self.BordXptDetCkButValue, command=self.cmd_ck_bord_xpt_det, width = xWidth)
		self.BordXptDetCkBut.grid(row=iFrow, column=3,  pady=nPady, padx=nPadx)

		self.RemoveXptOptMenu = Menubutton(ShowFrame, text='Remove Xpt.', width = xWidth)
		self.RemoveXptOptMenu.grid(row=iFrow, column=4,  pady=nPady, padx=nPadx, sticky=E+W)
		self.RemoveXptOptMenu.menu = Menu(self.RemoveXptOptMenu, tearoff=0)
		self.RemoveXptOptMenu['menu'] = self.RemoveXptOptMenu.menu


#	Borders
#	=====

		nPady  = 4
		xWidth = 12

		BordFrame = LabelFrame(self, text="Borders")
		BordFrame.grid(row=3, column=1, pady=nPady, padx=nPadx, sticky=E+W+N)

		iRow = 0
		self.BordCoreLabel = Label(BordFrame, text="Psi core:", width = xWidth)
		self.BordCoreLabel.grid(row=iRow, column=0,  pady=nPady, padx=nPadx)

		self.BordCoreEntryValue = StringVar()
		self.BordCoreEntry = Entry(BordFrame, textvariable = self.BordCoreEntryValue, width = xWidth)
		self.BordCoreEntry.grid(row=iRow,column=1,  pady=nPady, padx=nPadx)

		self.BordOutLabel = Label(BordFrame, text="Psi out:", width = xWidth)
		self.BordOutLabel.grid(row=iRow, column=2,  pady=nPady, padx=nPadx)

		self.BordOutEntryValue = StringVar()
		self.BordOutEntry = Entry(BordFrame, textvariable = self.BordOutEntryValue, width = xWidth)
		self.BordOutEntry.grid(row=iRow, column=3,  pady=nPady, padx=nPadx)
	
		iRow += 1
		BordXPtLabelH	= []
		HeaderXpt 			= ["X points:", "Psi Value", "Psi border A", "Psi border C"]
		HeaderWidth		= [xWidth, xWidth, xWidth, xWidth]
		for i in range(len(HeaderXpt)):
			BordXPtLabelH.append( Label(BordFrame, text=HeaderXpt[i],width=HeaderWidth[i]))
			BordXPtLabelH[i].grid(row=iRow,column=i,  pady=nPady, padx=nPadx)

		self.BordXPtLabelNValue = []
		self.BordXPtLabelN		= []
		self.BordXPtLabelPsi	= []
		self.BordXPtEntryPsiA	= []
		self.BordXPtEntryPsiC	= []
		self.BordXPtLabelPsiValue = []
		self.BordXPtEntryPsiAValue = []
		self.BordXPtEntryPsiCValue = []
		for i in range(self.MaxXpoints):
			iRow += 1
			self.BordXPtLabelNValue.append(StringVar())
			self.BordXPtLabelNValue[i].set("")
			self.BordXPtLabelN.append( Label(BordFrame, textvariable = self.BordXPtLabelNValue[i], background="white", width = xWidth))
			self.BordXPtLabelN[i].grid(row=iRow,column=0,  pady=1, padx=nPadx)

			self.BordXPtLabelPsiValue.append(StringVar())
			self.BordXPtLabelPsiValue[i].set("")
			self.BordXPtLabelPsi.append( Label(BordFrame, textvariable = self.BordXPtLabelPsiValue[i], background="white", width = xWidth))
			self.BordXPtLabelPsi[i].grid(row=iRow,column=1,  pady=1, padx=nPadx)

			self.BordXPtEntryPsiAValue.append(StringVar())
			self.BordXPtEntryPsiAValue[i].set("")
			self.BordXPtEntryPsiA.append( Entry(BordFrame, textvariable = self.BordXPtEntryPsiAValue[i], width = xWidth))
			self.BordXPtEntryPsiA[i].grid(row=iRow,column=2,  pady=1, padx=nPadx)

			self.BordXPtEntryPsiCValue.append(StringVar())
			self.BordXPtEntryPsiCValue[i].set("")
			self.BordXPtEntryPsiC.append( Entry(BordFrame, textvariable = self.BordXPtEntryPsiCValue[i], width = xWidth))
			self.BordXPtEntryPsiC[i].grid(row=iRow,column=3,  pady=1, padx=nPadx)

		iRow += 1

		self.RadArroundXpLabel = Label(BordFrame, text="Surr. Xp. (cm):", width = xWidth)
		self.RadArroundXpLabel.grid(row=iRow, column=0,  pady=nPady, padx=nPadx, sticky=E+W)

		self.RadArroundXpEntryValue = StringVar()
		self.RadArroundXpEntry = Entry(BordFrame, textvariable = self.RadArroundXpEntryValue, width = xWidth)
		self.RadArroundXpEntry.grid(row=iRow,column=1,  pady=nPady, padx=nPadx)

		self.DescentStepLabel = Label(BordFrame, text="Des. step (mm):", width = xWidth)
		self.DescentStepLabel.grid(row=iRow, column=2,  pady=nPady, padx=nPadx, sticky=E+W)

		self.DescentStepEntryValue = StringVar()
		self.DescentStepEntry = Entry(BordFrame, textvariable = self.DescentStepEntryValue, width = xWidth)
		self.DescentStepEntry.grid(row=iRow, column=3,  pady=nPady, padx=nPadx)

		iRow += 1
		self.SplitMagZonesButton = Button(BordFrame, text="Split Zones", command= self.cmd_split_mag_zones, width = xWidth)
		self.SplitMagZonesButton.grid(row=iRow, column=0, pady=nPady, padx=nPadx)

		self.SplitMagZonesButton = Button(BordFrame, text="Conn. Xp.", command= self.cmd_connect_x_points, width = xWidth)
		self.SplitMagZonesButton.grid(row=iRow, column=1, pady=nPady, padx=nPadx)


#	Limits
#	====

		xWidth = 12

		LimitsFrame = LabelFrame(self, text="Limits")
		LimitsFrame.grid(row=4, column=1, columnspan=2, pady=nPady, padx=nPadx, sticky=E+W+N)

		iRow = 0
		self.NewLimitButtonValue = StringVar()
		self.NewLimitButtonValue.set("New limit")
		self.NewLimitButton = Button(LimitsFrame, textvariable=self.NewLimitButtonValue, command = self.cmd_new_limit, width = xWidth)
		self.NewLimitButton.grid(row=iRow, column=0, pady=2, padx=nPadx)

		self.MovePointLimitButtonValue = StringVar()
		self.MovePointLimitButtonValue.set("Move point")
		self.MovePointLimitButton = Button(LimitsFrame, textvariable=self.MovePointLimitButtonValue, command = self.cmd_move_point_limit, width = xWidth)
		self.MovePointLimitButton.grid(row=iRow, column=1, pady=2, padx=nPadx)
		self.MovePointLimitButton.state(["disabled"])

		self.RemovePointLimitButtonValue = StringVar()
		self.RemovePointLimitButtonValue.set("Remove pt.")
		self.RemovePointLimitButton = Button(LimitsFrame, textvariable=self.RemovePointLimitButtonValue, command = self.cmd_remove_point_limit, width = xWidth)
		self.RemovePointLimitButton.grid(row=iRow, column=2, pady=2, padx=nPadx)
		self.RemovePointLimitButton.state(["disabled"])

		self.RemoveLimitButtonValue = StringVar()
		self.RemoveLimitButtonValue.set("Remove limit")
		self.RemoveLimitButton = Button(LimitsFrame, textvariable=self.RemoveLimitButtonValue, command = self.cmd_remove_limit, width = xWidth)
		self.RemoveLimitButton.grid(row=iRow, column=3, pady=2, padx=nPadx)
		self.RemoveLimitButton.state(["!disabled"])

		iRow += 1
		self.MoveCutButtonValue = StringVar()
		self.MoveCutButtonValue.set("Move X Cut")
		self.MoveCutButton = Button(LimitsFrame, textvariable=self.MoveCutButtonValue, command=self.cmd_move_x_cut, width = xWidth)
		self.MoveCutButton.grid(row=iRow, column=0, pady=2, padx=nPadx)

		self.MoveCutAddButtonValue = StringVar()
		self.MoveCutAddButtonValue.set("Add points")
		self.MoveCutAddButton = Button(LimitsFrame, textvariable=self.MoveCutAddButtonValue, command=self.cmd_move_x_cut_add_point, width = xWidth)
		self.MoveCutAddButton.grid(row=iRow, column=1, pady=2, padx=nPadx)
		self.MoveCutAddButton.state(["disabled"])

		self.MoveCutMoveButtonValue = StringVar()
		self.MoveCutMoveButtonValue.set("Move point")
		self.MoveCutMoveButton = Button(LimitsFrame, textvariable=self.MoveCutMoveButtonValue, command=self.cmd_move_x_cut_move_point, width = xWidth)
		self.MoveCutMoveButton.grid(row=iRow, column=2, pady=2, padx=nPadx)
		self.MoveCutMoveButton.state(["disabled"])

		self.MoveCutRemoveButtonValue = StringVar()
		self.MoveCutRemoveButtonValue.set("Remove pt.")
		self.MoveCutRemoveButton = Button(LimitsFrame, textvariable=self.MoveCutRemoveButtonValue, command=self.cmd_move_x_cut_remove_point, width = xWidth)
		self.MoveCutRemoveButton.grid(row=iRow, column=3, pady=2, padx=nPadx)
		self.MoveCutRemoveButton.state(["disabled"])

		iRow += 1

		self.DefineMagZonesButton = Button(LimitsFrame, text="Define Zones", command= self.cmd_define_mag_zones, width = xWidth)
		self.DefineMagZonesButton.grid(row=iRow, column=0, pady=nPady, padx=nPadx)

		self.DeleteMagZonesButton = Button(LimitsFrame, text="Delete Zone", command= self.cmd_delete_mag_zones, width = xWidth)
		self.DeleteMagZonesButton.grid(row=iRow, column=1, pady=nPady, padx=nPadx)

#	S3X setting
#	=======

		self.S3xSettingsFrame = LabelFrame(self, text="S3X setting")
		self.S3xSettingsFrame.grid(row=3, column=2, pady=nPady, padx=nPadx, sticky=E+W+N+S)

		iRow = 0
		self.S3xNPhiLabel = Label(self.S3xSettingsFrame, text="Nphi:", width = xWidth)
		self.S3xNPhiLabel.grid(row=iRow, column=0,  pady=nPady, padx=nPadx)

		self.S3xNPhiEntryValue = StringVar()
		self.S3xNPhiEntry = Entry(self.S3xSettingsFrame, textvariable = self.S3xNPhiEntryValue, width = xWidth)
		self.S3xNPhiEntry.grid(row=iRow,column=1,  pady=nPady, padx=nPadx)

		iRow += 1
		self.S3xR0Label = Label(self.S3xSettingsFrame, text="Ro (m):", width = xWidth)
		self.S3xR0Label.grid(row=iRow, column=0,  pady=nPady, padx=nPadx)

		self.S3xR0EntryValue = StringVar()
		self.S3xR0Entry = Entry(self.S3xSettingsFrame, textvariable = self.S3xR0EntryValue, width = xWidth)
		self.S3xR0Entry.grid(row=iRow,column=1,  pady=nPady, padx=nPadx)

		iRow += 1
		S3xa0Label = Label(self.S3xSettingsFrame, text="ao (m):", width = xWidth)
		S3xa0Label.grid(row=iRow, column=0,  pady=nPady, padx=nPadx)

		self.S3xa0EntryValue = StringVar()
		self.S3xa0Entry = Entry(self.S3xSettingsFrame, textvariable = self.S3xa0EntryValue, width = xWidth)
		self.S3xa0Entry.grid(row=iRow,column=1,  pady=nPady, padx=nPadx)

		iRow += 1
		S3xZ0Label = Label(self.S3xSettingsFrame, text="Zo (m):", width = xWidth)
		S3xZ0Label.grid(row=iRow, column=0,  pady=nPady, padx=nPadx)

		self.S3xZ0EntryValue = StringVar()
		self.S3xZ0Entry = Entry(self.S3xSettingsFrame, textvariable = self.S3xZ0EntryValue, width = xWidth)
		self.S3xZ0Entry.grid(row=iRow,column=1,  pady=nPady, padx=nPadx)

		iRow += 1
		S3xBt0Label = Label(self.S3xSettingsFrame, text="Bto (T):", width = xWidth)
		S3xBt0Label.grid(row=iRow, column=0,  pady=nPady, padx=nPadx)

		self.S3xBt0EntryValue = StringVar()
		self.S3xBt0Entry = Entry(self.S3xSettingsFrame, textvariable = self.S3xBt0EntryValue, width = xWidth)
		self.S3xBt0Entry.grid(row=iRow,column=1,  pady=nPady, padx=nPadx)

		iRow += 1
		self.S3xResetButton = Button(self.S3xSettingsFrame, text="Reset", command = self.cmd_reset_S3Xsettings, width = xWidth)
		self.S3xResetButton.grid(row=iRow, column=1, pady=2, padx=nPadx)

#		self.S3xSettingsFrame.state(["disabled"])

#		Variable definitions

		self.ZoneElements				= []
		self.Cotes						= []

		self.resize_mode				= False
		self.sel_mode_d					= False
		self.CheckPlot					= False
		self.InNewLimit					= False
		self.InRemoveLimit				= False
		self.InMoveCut					= False
		self.zone_corners_OK			= False

#		Starting default values

		self.NREntryValue.set(0)
		self.NZEntryValue.set(0)
		self.OutLengthEntryValue.set(-1.)
		self.SmoothInEntryValue.set(0.)
		self.SmoothExtEntryValue.set(0.)		
		self.RadArroundXpEntryValue.set(CIRCLE_RADIUS*1.e2)
		self.DescentStepEntryValue.set(STEP_FOLLOW_GRAD*1.e3)
		self.ExtScale.set(0.1)
		for RaiseFacEntryValue in self.RaiseFacEntryValues: RaiseFacEntryValue.set("{:0.2e}".format(1.e-4))

		self.FluxEntryValue.set("{:f}".format(0.))
		self.FluxScaleValue.set(0.)

		self.FluxCkButValue.set(0)
		self.BordXptDetCkButValue.set(0)

		self.update_table_X_points()

#		S3X Settings

		self.S3xNPhiEntryValue.set(1)
		self.S3xR0EntryValue.set(0.)
		self.S3xa0EntryValue.set(0.)
		self.S3xZ0EntryValue.set(0.)
		self.S3xBt0EntryValue.set(0.)

#		Temporaneo

#		self.BordCoreEntryValue.set("3.43")
#		self.BordOutEntryValue.set("6.8")
#		self.BordXPtEntryPsiAValue[0].set("4.05")
#		self.BordXPtEntryPsiCValue[0].set("3.43")

#	Internal functions
#===============================================================================

	def update_table_X_points(self):
		X_points = self.root.Config.X_points
		
		indexes	=  np.argsort(np.array([X_points[k].index for k in range(len(X_points))]))
		for i in range(min(len(X_points),self.MaxXpoints)):
			self.BordXPtLabelN[i].state(["!disabled"])
			self.BordXPtLabelPsi[i].state(["!disabled"])
			self.BordXPtEntryPsiA[i].state(["!disabled"])
			self.BordXPtEntryPsiC[i].state(["!disabled"])

			self.BordXPtLabelNValue[i].set("{:d}".format(i+1))
			self.BordXPtLabelPsiValue[i].set("{:.5f}".format(X_points[indexes[i]].psi))
#			RemoveXptValues.append("Xpt N. {:d}".format(i+1))

		for i in range(len(self.root.Config.X_points), self.MaxXpoints):
			self.BordXPtLabelN[i].state(["disabled"])
			self.BordXPtLabelPsi[i].state(["disabled"])
			self.BordXPtEntryPsiA[i].state(["disabled"])
			self.BordXPtEntryPsiC[i].state(["disabled"])
			self.BordXPtLabelNValue[i].set("")
			self.BordXPtLabelPsiValue[i].set("")
			self.BordXPtEntryPsiAValue[i].set("")
			self.BordXPtEntryPsiCValue[i].set("")

		iLast = self.RemoveXptOptMenu.menu.index('end')
		if(iLast != None): self.RemoveXptOptMenu.menu.delete(0,iLast)
		if(len(X_points) > 0): self.RemoveXptOptMenu.menu.add_command(label="Remove N .1", command=lambda: self.cmd_remove_xpt(1))
		if(len(X_points) > 1): self.RemoveXptOptMenu.menu.add_command(label="Remove N .2", command=lambda: self.cmd_remove_xpt(2))
		if(len(X_points) > 2): self.RemoveXptOptMenu.menu.add_command(label="Remove N .3", command=lambda: self.cmd_remove_xpt(3))
		if(len(X_points) > 3): self.RemoveXptOptMenu.menu.add_command(label="Remove N .4", command=lambda: self.cmd_remove_xpt(4))
		if(len(X_points) > 4): self.RemoveXptOptMenu.menu.add_command(label="Remove N .5", command=lambda: self.cmd_remove_xpt(5))
		if(len(X_points) > 5): self.RemoveXptOptMenu.menu.add_command(label="Remove N .6", command=lambda: self.cmd_remove_xpt(6))

	"""
		for i in range(len(self.Frontiers)):
			self.BordPsiLabelLimN[i].state(["!disabled"])
			self.BordPsiLabelLimMin[i].state(["!disabled"])
			self.BordPsiLabelLimMax[i].state(["!disabled"])
			self.BordPsiLabelLimNValue[i].set("{:d}".format(i+1))
			self.BordPsiLabelLimMinValue[i].set("{:3.2e}".format(Frontiers[i].psimin))
			self.BordPsiLabelLimMaxValue[i].set("{:3.2e}".format(Frontiers[i].psimax))

		for i in range(len(self.root.Config.X_points), self.MaxBariers):
			self.BordPsiLabelLimN[i].state(["disabled"])
			self.BordPsiLabelLimMin[i].state(["disabled"])
			self.BordPsiLabelLimMax[i].state(["disabled"])
			self.BordPsiLabelLimNValue[i].set("")
			self.BordPsiLabelLimMinValue[i].set("")
			self.BordPsiLabelLimMaxValue[i].set("")
	"""
#	move flux scale

	def cmd_move_flux_scale(self, slider):
		if(not hasattr(self.root.Config,'r2D')):	return

		val					= self.FluxScaleValue.get()
		flux_min			= np.min(self.root.Config.flux2D)
		flux_max			= np.max(self.root.Config.flux2D)
		manual_psi_val		= flux_min + (flux_max - flux_min)*val/100.
		self.FluxEntryValue.set("{:f}".format(manual_psi_val))

		if(self.FluxCkButValue.get() == 1): self.update_figure()

#	Flux entry changed

	def cmd_enter_flux_entry(self, event):
		if((not hasattr(self.root.Config,'r2D')) or (self.FluxCkButValue.get() != 1)): return
		self.update_figure()

	def cmd_flux_ck_button(self):
		if(not hasattr(self.root.Config,'r2D')): return
		self.update_figure()

#	Smooth2Dng and extrapolation region
#-------------------------------------

#	move Extrapolating scale

	def cmd_move_ext_scale(self, slider):
		self.root.Config.extrapol_val	= self.ExtScale.get()
		self.ExtFacLabelValue.set("Ext. ({:4.2f}):".format(self.root.Config.extrapol_val))

	def cmd_ck_bord_xpt_det(self):
		if(not hasattr(self.root.Config,'r2D')):	return
		self.update_figure()

#	plot smooth and extrapolation

	def cmd_smooth_extrapolate(self):
		self.root.Config.nr				= eval(self.NREntryValue.get())
		self.root.Config.nz				= eval(self.NZEntryValue.get())
		self.root.Config.OutLength2		= eval(self.OutLengthEntryValue.get())
		self.root.Config.Smooth2D		= eval(self.SmoothInEntryValue.get())
		self.root.Config.Smooth2DExt	= eval(self.SmoothExtEntryValue.get())
		self.root.Config.extrapol_val	= self.ExtScale.get()
		if(self.root.SameRaiseFactorCheckValue.get() == 1): 
			self.root.Config.raise_value 	= np.array([eval(self.RaiseFacEntryValues[0].get()), eval(self.RaiseFacEntryValues[0].get()),
														eval(self.RaiseFacEntryValues[0].get()), eval(self.RaiseFacEntryValues[0].get())])
		else:
			self.root.Config.raise_value 	= np.array([eval(self.RaiseFacEntryValues[0].get()), eval(self.RaiseFacEntryValues[1].get()),
														eval(self.RaiseFacEntryValues[2].get()), eval(self.RaiseFacEntryValues[3].get())])
		self.root.Config.raise_power 	= eval(self.RaisePowerOptMenuValue.get())
		
		smooth_fields(self.root.Config)
		if(self.root.Config.extrapol_val > 0.): extrapolate_psi_out(self.root.Config)
		self.root.post_loaded_fields()
		
#	Frontiers
#=====================

#	Manage new limit
#	==========

	def cmd_new_limit(self):
		if(not self.InNewLimit):
			self.NewLimitButtonValue.set("Close")
			self.RemovePointLimitButton.state(["!disabled"])

			self.list_res = np.empty((0,2),dtype='f8')
			self.cid_mouse_button_press_limit = 0
			self.InNewLimit = True
			self.cid_mouse_button_press_limit = self.Fig.canvas.mpl_connect('button_press_event', self.cmd_button_press_limit)
			self.cursor = Cursor(self.Ax, useblit=False, color='red', linewidth=1)
		else:
			self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press_limit)
			self.cid_mouse_button_press_limit = 0
			self.cursor=0

			self.InNewLimit = False
			self.NewLimitButtonValue.set("New limit")
			self.RemovePointLimitButton.state(["disabled"])
			self.InLimitAddPoints = False

#			look for intersections with X_point radial limits (should be two)
			
			self.root.Config.Frontiers.append(define_frontier(self.root.Config, self.list_res))

			if(len(self.root.Config.Frontiers) >= 2*len(self.root.Config.X_points)): self.root.set_Frontiers_OK(self.root.Config, True)
			self.update_figure()


#	Limits move one point
#	==============

	def cmd_move_point_limit(self):
		print("To do move point")


#	Limits remove one point
#	==============

	def cmd_remove_point_limit(self):
		print("To do move point")

	def cmd_remove_limit(self):
		if(not self.InRemoveLimit):
			self.RemoveLimitButtonValue.set("Cancel")
			self.InRemoveLimit = True
			self.cid_mouse_button_press_limit = self.Fig.canvas.mpl_connect('button_press_event', self.cmd_button_press_limit)
			self.cursor = Cursor(self.Ax, useblit=False, color='red', linewidth=1)
		else:
			self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press_limit)
			self.cid_mouse_button_press_limit = 0
			self.cursor=0

			self.InRemoveLimit = False
			self.RemoveLimitButtonValue.set("Remove limit")
			self.update_figure()
	
	def cmd_button_press_limit(self, event):
		Rpt		  = event.xdata
		Zpt		  = event.ydata

		if(self.InNewLimit):
			self.list_res = np.append(self.list_res, [[Rpt, Zpt]], axis=0)
			self.Ax.plot(self.list_res[:,0],self.list_res[:,1],'b.-')
			self.Fig.canvas.draw()
		elif(self.InRemoveLimit):
			Frontiers = self.root.Config.Frontiers
			dMin = 1e40
			iMin = -1
			for i in range(len(Frontiers)):
				dF = np.min((Frontiers[i].R - Rpt)**2 + (Frontiers[i].Z - Zpt)**2)
				if(dF < dMin):
					dMin = dF
					iMin = i

			if(len(self.root.Config.Frontiers) == 2*len(self.root.Config.X_points)): self.root.set_Frontiers_OK(self.root.Config, False)
			del self.root.Config.Frontiers[iMin]									#Remove Frontier
			self.cmd_remove_limit()													#Reset button and update figure

#	========================================0
#	Remove Xpt

	def cmd_remove_xpt(self, Xpt_index):
		ret_status = messagebox.askokcancel("Remove X point", "Do you want remove X point n. {:d}".format(Xpt_index), icon=messagebox.QUESTION, default=messagebox.OK)
		X_points = self.root.Config.X_points
		if(ret_status):
			for k in range(len(self.root.Config.X_points)):
				if(X_points[k].index == Xpt_index-1):
					del X_points[k]
					break
			for k in range(len(X_points)):
				if(X_points[k].index >= Xpt_index): X_points[k].index -= 1		
				
			self.update_table_X_points()
			self.update_figure()

#================
#	Move X branch
#================

	def cmd_move_x_cut(self):
		if(not self.InMoveCut):
			ret_status = messagebox.askokcancel("Move X branch", "1) choose branch to move\n2) set end point\n3) add bezier points:\n" + \
												"    a) select segmeent\n    b) set control point", icon=messagebox.INFO, default=messagebox.OK)
			if(not ret_status): return
			
			self.MoveCutButtonValue.set("Close")
			self.InAddPoint		= False
			self.InMovePoint	= False
			self.InRemovePoint	= False
			self.InMoveCut		= self.set_move_x_button_press_event()
			self.X_point_sel	= -1
		else:
			self.MoveCutButtonValue.set("Move X Cut")
			self.MoveCutAddButton.state(["disabled"])
			self.MoveCutMoveButton.state(["disabled"])
			self.MoveCutRemoveButton.state(["disabled"])
			self.InMoveCut = self.reset_move_x_button_press_event()

			Config		= self.root.Config
			X_points	= self.root.Config.X_points
			c1 	   		= types.SimpleNamespace()
			c1.arc 		= [types.SimpleNamespace()]
			c1.arc[0].x, c1.arc[0].y = bezier(self.rNewBranch, self.zNewBranch)
			c2			= contour_better(Config.r2D, Config.z2D, Config.flux2D, [X_points[self.X_point_sel].cut[self.cut_sel].psilim])
			X			= intersect_contour(c1,c2)

			cin 	    = types.SimpleNamespace()
			cin.x		= c1.arc[0].x
			cin.y		= c1.arc[0].y
			p1 	   		= types.SimpleNamespace()
			p1.x		= X[0].x
			p1.y		= X[0].y
			cout		= part_contour(cin,p1)
			self.Ax.plot(cout.arc[0].x, cout.arc[0].y, 'g-', linewidth=3)
			X_points[self.X_point_sel].cut[self.cut_sel].R		= cout.arc[0].x
			X_points[self.X_point_sel].cut[self.cut_sel].Z		= cout.arc[0].y
			X_points[self.X_point_sel].cut[self.cut_sel].type	= 3
			self.sel_mode_d										= False
			self.update_figure()
				
	def cmd_move_x_cut_add_point(self):
		if(not self.InAddPoint):
			self.MoveCutAddButtonValue.set("Cancel")
			self.InAddPoint = self.set_move_x_button_press_event()
			self.iPtSeg    = -1
		else:
			self.MoveCutAddButtonValue.set("Add point")
			self.InAddPoint = self.reset_move_x_button_press_event()

	def cmd_move_x_cut_move_point(self):
		if(not self.InMovePoint):
			self.MoveCutMoveButtonValue.set("Cancel")
			self.InMovePoint = self.set_move_x_button_press_event()
			self.iPtSeg    = -1
		else:
			self.MoveCutMoveButtonValue.set("Move point")
			self.InMovePoint = self.reset_move_x_button_press_event()

	def cmd_move_x_cut_remove_point(self):
		if(not self.InRemovePoint):
			self.MoveCutRemoveButtonValue.set("Cancel")
			self.InRemovePoint = self.set_move_x_button_press_event()
		else:
			self.MoveCutRemoveButtonValue.set("Remove point")
			self.InRemovePoint = self.reset_move_x_button_press_event()
		

	def set_move_x_button_press_event(self):
		self.cid_mouse_button_press = self.Fig.canvas.mpl_connect('button_press_event', self.cmd_move_x_button_press)
		self.cursor = Cursor(self.Ax, useblit=False, color='red', linewidth=1)
		return True

	def reset_move_x_button_press_event(self):
		if(self.cid_mouse_button_press != 0):
			self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press)
			self.cid_mouse_button_press = 0
			self.cursor=0
		return False

	
#	Mouses events for branch movement

	def cmd_move_x_button_press(self, event):
		Rpt		 = event.xdata
		Zpt		 = event.ydata
		reset	 = self.reset_move_x_button_press_event()
		X_points = self.root.Config.X_points


		if(self.InAddPoint):
			if(self.iPtSeg == -1):
				d, self.iPtSeg = find_closest_segment(Rpt, Zpt, self.rNewBranch, self.zNewBranch)
				self.Ax.plot(self.rNewBranch[self.iPtSeg:self.iPtSeg+2], self.zNewBranch[self.iPtSeg:self.iPtSeg+2],'r-')
				self.Fig.canvas.draw()
				self.InAddPoint = self.set_move_x_button_press_event()
			else:
				self.cmd_move_x_cut_add_point()
				self.rNewBranch = np.concatenate((self.rNewBranch[:self.iPtSeg+1], np.array([Rpt]), self.rNewBranch[self.iPtSeg+1:]))
				self.zNewBranch = np.concatenate((self.zNewBranch[:self.iPtSeg+1], np.array([Zpt]), self.zNewBranch[self.iPtSeg+1:]))
				rBezier, zBezier = bezier(self.rNewBranch, self.zNewBranch)
				self.update_figure()
				self.Ax.plot(self.rNewBranch, self.zNewBranch,'k.-')
				self.Ax.plot(rBezier, zBezier,'b-')
				self.Fig.canvas.draw()
				self.MoveCutMoveButton.state(["!disabled"])
				self.MoveCutRemoveButton.state(["!disabled"])


		elif(self.InMovePoint):
			if(self.iPtSeg == -1):
				d = np.sqrt((self.rNewBranch[1:-1]-Rpt)**2 + (self.zNewBranch[1:-1]-Zpt)**2)
				self.iPtSeg = np.argmin(d)+1
				self.Ax.plot(self.rNewBranch[self.iPtSeg], self.zNewBranch[self.iPtSeg],'ro')
				self.Fig.canvas.draw()
				self.InMovePoint = self.set_move_x_button_press_event()
			else:
				self.cmd_move_x_cut_move_point()
				self.InMovePoint  = self.reset_move_x_button_press_event()
				self.rNewBranch[self.iPtSeg] = Rpt
				self.zNewBranch[self.iPtSeg] = Zpt
				rBezier, zBezier = bezier(self.rNewBranch, self.zNewBranch)
				self.update_figure()
				self.Ax.plot(self.rNewBranch, self.zNewBranch,'k.-')
				self.Ax.plot(rBezier, zBezier,'b-')
				self.Fig.canvas.draw()

		elif(self.InRemovePoint):
			self.cmd_move_x_cut_remove_point()
			d = np.sqrt((self.rNewBranch[1:-1]-Rpt)**2 + (self.zNewBranch[1:-1]-Zpt)**2)
			self.iPtSeg = np.argmin(d)+1
			self.rNewBranch = np.append(self.rNewBranch[:self.iPtSeg], self.rNewBranch[self.iPtSeg+1:])
			self.zNewBranch = np.append(self.zNewBranch[:self.iPtSeg], self.zNewBranch[self.iPtSeg+1:])

			rBezier, zBezier = bezier(self.rNewBranch, self.zNewBranch)
			self.update_figure()
			self.Ax.plot(self.rNewBranch, self.zNewBranch,'k.-')
			self.Ax.plot(rBezier, zBezier,'b-')
			self.Fig.canvas.draw()
			if(len(self.rNewBranch) < 3):
				self.MoveCutMoveButton.state(["disabled"])
				self.MoveCutRemoveButton.state(["disabled"])
				
		elif(self.InMoveCut):										# Move Branch
			if(self.X_point_sel == -1):
				dmin = 1e10											# find closest segment
				for k in range(len(X_points)):
					for n in range(4):
						d, i = find_closest_segment(Rpt, Zpt, X_points[k].cut[n].R, X_points[k].cut[n].Z)
						if(d < dmin):
							dmin				= d
							self.X_point_sel	= k
							self.cut_sel		= n
							self.sel_mode_d		= True
				
				self.update_figure()		
#				ret_status = messagebox.showinfo("Move X branch", "Choose chose new end point and then add control points")
				self.InMoveCut = self.set_move_x_button_press_event()
			else:
				InMoveCut		= self.reset_move_x_button_press_event()
						
				self.rNewBranch = np.array([X_points[self.X_point_sel].R, Rpt])
				self.zNewBranch = np.array([X_points[self.X_point_sel].Z, Zpt])
				vec	= np.array([[X_points[self.X_point_sel].R, X_points[self.X_point_sel].Z],[Rpt, Zpt]])
				rBezier, zBezier = bezier(self.rNewBranch, self.zNewBranch)
				self.update_figure()
				self.Ax.plot(self.rNewBranch, self.zNewBranch,'k.-')
				self.Ax.plot(rBezier, zBezier,'b-')
				self.Fig.canvas.draw()
				
				self.MoveCutAddButton.state(["!disabled"])


#=========================
#	Split MagZones button
#=========================

	def cmd_split_mag_zones(self):
		Config	 = self.root.Config
		X_points = self.root.Config.X_points

		Config.psicore		= eval(self.BordCoreEntryValue.get())
		Config.psiout		= eval(self.BordOutEntryValue.get())
		Config.DescentStep	= eval(self.DescentStepEntryValue.get())*1e-3
		
		indexes	=  np.argsort(np.array([X_points[k].index for k in range(len(X_points))]))
		for k in range(len(indexes)):
			X_points[indexes[k]].cut[0].psilim = eval(self.BordXPtEntryPsiAValue[k].get())
			X_points[indexes[k]].cut[1].psilim = Config.psiout
			X_points[indexes[k]].cut[2].psilim = eval(self.BordXPtEntryPsiCValue[k].get())
			X_points[indexes[k]].cut[3].psilim = Config.psiout

#		FoundOrtho = plot_ortho_lines(self.root, Config, step=DescentStep)
		FoundOrtho = plot_ortho_lines(self.root, Config, step=STEP_FOLLOW_GRAD)
		if(FoundOrtho): FoundArcs  = getArcLim(self.root, Config)
		else:			FoundArcs  = False
		if(Config.Frontiers_OK): Config.Frontiers = []
		self.root.set_xPoints_OK(Config, FoundArcs)

		self.update_figure()


	def cmd_connect_x_points(self):
		Values = ask_values(self, ["Xp_1","Xp_2"], title="Connecting X-points")					#Ask X-points to connect
		if(len(Values) != 2): return

		Config	 = self.root.Config
		X_points = self.root.Config.X_points
		indexes	=  np.argsort(np.array([X_points[k].index for k in range(len(X_points))]))
		iX 		= [eval(Values[0])-1,eval(Values[1])-1]													#set X-points

		psi_ave = 0.5*(X_points[indexes[iX[0]]].psi + X_points[indexes[iX[1]]].psi)						#Average psi
		X_points[indexes[iX[0]]].psi = psi_ave
		X_points[indexes[iX[1]]].psi = psi_ave

		cpsi = contour_better(Config.r2D,Config.z2D,Config.flux2D, [psi_ave])							#Find new position for X-points on the new flux surface
		for kX in range(2):
			d_min = np.empty(len(cpsi.arc), dtype='f8')
			x_min = np.empty(len(cpsi.arc), dtype='f8')
			y_min = np.empty(len(cpsi.arc), dtype='f8')
			for iarc in range(len(cpsi.arc)):
				dist  = (cpsi.arc[iarc].x - X_points[indexes[iX[kX]]].R)**2 + (cpsi.arc[iarc].y - X_points[indexes[iX[kX]]].Z)**2
				i_min = np.argmin(dist)
				d_min[iarc] = dist[i_min]
				x_min[iarc] = cpsi.arc[iarc].x[i_min]
				y_min[iarc] = cpsi.arc[iarc].y[i_min]
			i_min = np.argmin(d_min)
			X_points[indexes[iX[kX]]].R = x_min[i_min]
			X_points[indexes[iX[kX]]].Z = y_min[i_min]

		surround_Xpoints(self.root, Config, self.root.CheckPlot, RadArroundXp=Config.RadArroundXp)		#surround new X-points
		self.update_table_X_points()																	#Update table

#=========================			
#	Define MagZones button
#=========================

	def cmd_define_mag_zones(self):
		Config = self.root.Config
		if(not Config.Frontiers_OK):
			messagebox.showwarning("Grid-Gen", "At least {:d} limits must be defines".format(2*len(self.root.Config.X_points)))
			return
		
		self.zone_corners_OK = find_zone_corners(self.root, Config)
		if(self.zone_corners_OK):
			IsOk = find_segment_theta(self, Config)
			if(IsOk): IsOk = find_cote_psi(self, Config)
			if(IsOk): IsOk = remove_double_cote(Config)		
			if(IsOk): IsOk = find_zones(self,Config)
			if(IsOk): IsOk = analyze_topology(self, Config)
			if(IsOk): IsOk = find_megazones_outside_plasma(self, Config)
			if(IsOk): IsOk = define_OMP_segment(Config)
			if(IsOk):
				Config.Split = []
				self.root.set_MagZones_OK(Config, True)
				self.update_figure()

	def cmd_delete_mag_zones(self):
		if(not self.root.Config.MagZones_OK): return

		ret_status = messagebox.askokcancel("Delete Zone", "Remove a zone selecting it with the cursor", icon=messagebox.INFO, default=messagebox.OK)
		if(not ret_status): return

		self.cid_mouse_button_press_delete_zone = self.Fig.canvas.mpl_connect('button_press_event', self.cmd_button_press_delete_zone)
		self.cursor = Cursor(self.Ax, useblit=False, color='red', linewidth=1)

	
	def cmd_button_press_delete_zone(self, event):
		self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press_delete_zone)
		self.cid_mouse_button_press_delete_zone = 0
		self.cursor = 0
		
		kMagZone = get_magzone_of_point(self.root.Config, event.xdata, event.ydata)
		if(kMagZone == -1):
			ret_status = messagebox.showinfo("Delete Zone", "the selected point is not in a zone. Skipping...")
			return
		else:
			ret_status = messagebox.askokcancel("Delete Zone", "Do you want remove Zone n. {:d}".format(kMagZone+1), icon=messagebox.QUESTION, default=messagebox.OK)
			if(not ret_status): return

			self.root.set_MagZones_OK(self.root.Config, False)
			delete_magzone(self.root.Config, kMagZone)

			self.root.set_MagZones_OK(self.root.Config, True)
			self.update_figure()


#=========================			
#	S3X settings
#=========================

	def cmd_reset_S3Xsettings(self):					#like matlab routine reset_S3Xsettings

		Config		= self.root.Config
		X_points	= Config.X_points
		Rcore, Zcore, Bcore = find_core(Config)							# Locate center of plasma chamber
	
#		Determines a by looking at intersection between the separatrix and the horizontal line going through the core center

		psisep = 1.e10
		for k in range(len(Config.X_points)):
			if(X_points[k].psi < psisep): psisep = X_points[k].psi

		csep = contour_better(Config.r2D, Config.z2D, Config.flux2D, [psisep])

		ccore 	   		= types.SimpleNamespace()
		ccore.arc 		= [types.SimpleNamespace()]
		ccore.arc[0].x	= np.array([Rcore, Rcore+100.])
		ccore.arc[0].y	= np.array([Zcore, Zcore])
		X = intersect_contour(csep,ccore)

		if(len(X) > 0): deltaR = np.max(np.array([x.x for x in X])) - Rcore
		else:			deltaR = X[0].x - Rcore

		self.S3xR0EntryValue.set("{:.3f}".format(Rcore))				# Sets the values in the box
		self.S3xa0EntryValue.set("{:.3f}".format(deltaR))
		self.S3xZ0EntryValue.set("{:.3f}".format(Zcore))
		self.S3xBt0EntryValue.set("{:.3f}".format(Bcore))

#==============
#	Main figure 
#==============

#update main figure

	def update_figure(self):

		Config		= self.root.Config
		Frontiers	= Config.Frontiers
		MagZones	= self.root.Config.MagZones
		X_points	= self.root.Config.X_points
		Fig			= self.Fig
		Ax			= self.Ax

		Xlims = Ax.get_xlim()
		Ylims = Ax.get_ylim()

		Ax.clear()

		Ax.tick_params(axis='both', which='major', labelsize=8)
		Ax.autoscale(axis='both', tight=True)
		if(Config.equ_OK):
			for k in range(len(X_points)):
				Ax.contour(Config.r2D, Config.z2D, Config.flux2D, colors="r", levels=[X_points[k].psi])
				Ax.text(X_points[k].R,X_points[k].Z,"{:d}".format(X_points[k].index+1),size=18)

		if(Config.wall_OK): plot2d_walls(Ax, Config.Walls, self.root.ShowPlasmaWallCheckValue.get(), self.root.ShowEireneWallCheckValue.get(), self.root.ShowExtraWallCheckValue.get())

		if(self.FluxCkButValue.get() == 1):
			Config.manual_psi_val	 = eval(self.FluxEntryValue.get())
			Ax.contour(Config.r2D, Config.z2D, Config.flux2D, colors="c", levels=[Config.manual_psi_val])

		if(Config.xPoints_OK):
			for k in range(len(X_points)):
				for n in range(4):
					Ax.plot(X_points[k].cut[n].R,    X_points[k].cut[n].Z,    'b-', linewidth=2)
					Ax.plot(X_points[k].cut[n].psiR, X_points[k].cut[n].psiZ, 'r--')


		if(self.BordXptDetCkButValue.get() == 1):
			for k in range(len(X_points)):
#				Ax.text(X_points[k].branch[0].R, X_points[k].branch[0].Z, '1')
#				Ax.text(X_points[k].branch[1].R, X_points[k].branch[0].Z, '2')
#				Ax.text(X_points[k].branch[2].R, X_points[k].branch[0].Z, '3')
#				Ax.text(X_points[k].branch[3].R, X_points[k].branch[0].Z, '4')

				labels = ['A', 'B', 'C', 'D']
				for n in range(4): Ax.text(X_points[k].cut[n].Rs, X_points[k].cut[n].Zs, labels[n])

		for k in range(len(Frontiers)):
			if(Frontiers[k].sel):
				Ax.plot(Frontiers[k].R, Frontiers[k].Z,'c.--', linewidth=2)
			else:
				Ax.plot(Frontiers[k].R, Frontiers[k].Z,'b.--', linewidth=2)

			Ax.plot(Frontiers[k].P1[:,0], Frontiers[k].P1[:,1],'mv')

		if(Config.MagZones_OK):
			for k in range(len(MagZones)):
				sCenter = floor(len(MagZones[k].south.R)/2)
				nCenter = floor(len(MagZones[k].north.R)/2)
				R=(MagZones[k].south.R[sCenter]+MagZones[k].north.R[nCenter])*0.5
				Z=(MagZones[k].south.Z[sCenter]+MagZones[k].north.Z[nCenter])*0.5

				Ax.plot(MagZones[k].south.R, MagZones[k].south.Z,'r-',	linewidth=2)
				Ax.plot(MagZones[k].north.R,MagZones[k].north.Z,	'b--',	linewidth=2)
				Ax.plot(MagZones[k].east.R,MagZones[k].east.Z,	'k-', 	linewidth=2)
				Ax.plot(MagZones[k].west.R,MagZones[k].west.Z,	'g--',	linewidth=2)
				Ax.text(R,Z,"{:d}".format(k+1))
		elif(self.zone_corners_OK):
			for k in range(len(MagZones)):
				if(hasattr(MagZones[k].south,'R') and hasattr(MagZones[k].north,'R')):
					sCenter = floor(len(MagZones[k].south.R)/2)
					nCenter = floor(len(MagZones[k].north.R)/2)
					R=(MagZones[k].south.R[sCenter]+MagZones[k].north.R[nCenter])*0.5
					Z=(MagZones[k].south.Z[sCenter]+MagZones[k].north.Z[nCenter])*0.5
					Ax.text(R,Z,"{:d}".format(k+1))
				elif(hasattr(MagZones[k].west,'R') and hasattr(MagZones[k].east,'R')):
					wCenter = floor(len(MagZones[k].west.R)/2)
					eCenter = floor(len(MagZones[k].east.R)/2)
					R=(MagZones[k].west.R[wCenter]+MagZones[k].east.R[eCenter])*0.5
					Z=(MagZones[k].west.Z[wCenter]+MagZones[k].east.Z[eCenter])*0.5
					Ax.text(R,Z,"{:d}".format(k+1))
			

		if(self.sel_mode_d):
			Ax.plot(X_points[self.X_point_sel].cut[self.cut_sel].R, \
				X_points[self.X_point_sel].cut[self.cut_sel].Z,				'm-',linewidth=3)

		Ax.set_aspect(1.)
		if((not ((Xlims[0] == 0.) and (Xlims[1] == 1.))) and (not self.resize_mode)):
			Ax.set_xlim(Xlims[0], Xlims[1])
			Ax.set_ylim(Ylims[0], Ylims[1])
		else:
			self.resize_mode = False

		Fig.canvas.draw()

