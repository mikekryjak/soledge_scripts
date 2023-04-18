#!/usr/bin/env python
# -*- coding: utf-8 -*-

from tkinter import *
import tkinter		as tk
from tkinter		import Tk, Text, filedialog, messagebox, BOTH, W, N, E, S, Menu, HORIZONTAL
from tkinter		import DoubleVar, IntVar, StringVar, OptionMenu
from tkinter.ttk	import Entry, Menubutton, Checkbutton, Frame, Button, Label, Style, LabelFrame, Scale, Separator

import types
import matplotlib
from tkinter.messagebox import Message
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg		import FigureCanvasTkAgg
from matplotlib.figure						import Figure
from matplotlib.widgets 					import Cursor
try:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2TkAgg as NavigationToolbar2Tk
except:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2Tk

from math									import floor, sqrt, asin, atan2, sin, cos
import os
import numpy								as np
import scipy.interpolate					as interp
from numpy									import arange, sin, pi

from interfaces.set_properties				import set_properties

from files.load_text_data					import load_text_data
from files.load_walls_file					import load_walls_file
from files.save_walls_file					import save_walls_file
from routines.find_closest_segment			import find_closest_segment
from routines.utils_walls					import walls_copy, wall_copy, plot2d_walls, walls_find_closest_segment, walls_find_closest_point
from routines.globals						import LENGTH_TYPES, DRAWING_WALL

class EditWall(Toplevel):

	def __init__(self, parent):
		Toplevel.__init__(self, parent, Title = None)

		self.result	= None
		self.transient(parent)
		self.parent	= parent
		self.title("Edit Wall")
		self.protocol("WM_DELETE_WINDOW", self.cmd_cancel)

		menubar = Menu(self.master)
		filemenu = Menu(menubar, tearoff=0)
		filemenu.add_command(label="Load Walls", command=self.cmd_load_walls)
		filemenu.add_command(label="Add Walls", command=self.cmd_add_walls)
		filemenu.add_command(label="Save Single Wall As..", command=self.cmd_save_wall)
		filemenu.add_command(label="Save All Walls As..", command=self.cmd_save_walls)
		menubar.add_cascade(label="File", menu=filemenu)

		OptionsMenu = Menu(menubar, tearoff=0)

		ShowWallMenu = Menu(self.master)

		self.ShowPlasmaWallCheckValue = IntVar()
		self.ShowPlasmaWallCheckValue.set(parent.ShowPlasmaWallCheckValue.get())
		ShowWallMenu.add_checkbutton(label="Soledge", 		variable=self.ShowPlasmaWallCheckValue)
		self.ShowPlasmaWallCheckValue.trace("w", self.cmd_ShowPlasmaWallCheckValue)

		self.ShowEireneWallCheckValue = IntVar()
		self.ShowEireneWallCheckValue.set(parent.ShowEireneWallCheckValue.get())
		ShowWallMenu.add_checkbutton(label="Eirene", 		variable=self.ShowEireneWallCheckValue)
		self.ShowEireneWallCheckValue.trace("w", self.cmd_ShowEireneWallCheckValue)

		self.ShowExtraWallCheckValue = IntVar()
		self.ShowExtraWallCheckValue.set(parent.ShowExtraWallCheckValue.get())
		ShowWallMenu.add_checkbutton(label="Extra", 		variable=self.ShowExtraWallCheckValue)
		self.ShowExtraWallCheckValue.trace("w", self.cmd_ShowExtraWallCheckValue)
	
		OptionsMenu.add_cascade(label="Show wall",				menu=ShowWallMenu)

		menubar.add_cascade(label="Options", 					menu=OptionsMenu)

		self.config(menu=menubar)

		nPadx  = 5
		nPady  = 5
		xWidth = 12

		self.Frame = Frame(self)
		self.Frame.pack(fill=BOTH, expand=True)
		self.Frame.columnconfigure(1, weight=1)
		self.Frame.grid(row=0, column=0, rowspan=4,  pady=nPady, padx=nPadx, sticky=E+W+N)

		iRow = 0
		"""
		self.FileOptions = ["Files", "Load Contour", "Save Contour", "Save Contour As..","Save EQDSK"]
		self.FileOptionsValue = StringVar()
		self.FileOptionsValue.set(self.FileOptions[0])
		opts = ()
		for k in self.FileOptions: opts += (k,)	
		self.FileOptionsMenu = OptionMenu(self.Frame, self.FileOptionsValue, *opts)
		self.FileOptionsMenu.grid(row=0, column=0,  pady=nPady, sticky=W+E)
		self.FileOptionsValue.trace("w", self.cmd_file)
		iRow += 1
		"""

#		Plot area
		self.Fig = Figure(figsize=(6.0, 5.2), dpi=100)
		self.Fig.patch.set_facecolor('white')

		# a tk.DrawingArea
		canvas = FigureCanvasTkAgg(self.Fig, master=self.Frame)
		try:
			canvas.show()										#Matplotlib  up 2.2
		except:
			canvas.draw()
		canvas.get_tk_widget().grid(row=iRow+1, column=0, columnspan=6,  padx=nPadx, sticky=E+W+S+N)

		self.Ax = self.Fig.add_subplot(111)
		self.Ax.tick_params(axis='both', which='major', labelsize=8)
		self.Fig.tight_layout()

		toolbar_frame = Frame(self.Frame, width=50)
		toolbar_frame.grid(row=iRow, column=0, columnspan=4, padx=nPadx, sticky=W+S+N)
		toolbar = NavigationToolbar2Tk(canvas, toolbar_frame )

		toolbar.update()
		iRow += 2

		EditLabel = Label(self.Frame, text = "Edit", width = xWidth)
		EditLabel.grid(row=iRow, column=0,  pady=nPady, padx=nPadx)

		EditSeparator = Separator(self.Frame)
		EditSeparator.grid(row=iRow, column=1,  columnspan=5,  pady=nPady, padx=0, sticky=E+W)
		iRow += 1

		self.AddPointButValue = StringVar()
		self.AddPointButValue.set("Add point")
		self.AddPointBut = Button(self.Frame, textvariable=self.AddPointButValue, command = self.cmd_add_point, width = xWidth)
		self.AddPointBut.grid(row=iRow, column=0, pady=nPady, padx=nPadx)

		self.RemovePointButValue = StringVar()
		self.RemovePointButValue.set("Remove point")
		self.RemovePointBut = Button(self.Frame, textvariable=self.RemovePointButValue, command = self.cmd_remove_point, width = xWidth)
		self.RemovePointBut.grid(row=iRow, column=1, pady=nPady, padx=nPadx)

		self.MovePointButValue = StringVar()
		self.MovePointButValue.set("Move point")
		self.MovePointBut = Button(self.Frame, textvariable=self.MovePointButValue, command = self.cmd_move_point, width = xWidth)
		self.MovePointBut.grid(row=iRow, column=2, pady=nPady, padx=nPadx)

		TranspValueLabel = Label(self.Frame, text="Value:", width = xWidth)
		TranspValueLabel.grid(row=iRow, column=3,  pady=nPady, padx=nPadx)

		self.CoordEntryValue = StringVar()
		self.CoordEntry = Entry(self.Frame,  textvariable = self.CoordEntryValue, width = xWidth)
		self.CoordEntry.grid(row=iRow, column=4,  pady=nPady, padx=nPadx)
		self.CoordEntry.bind("<Return>", self.cmd_enter_coord_entry)
		self.CoordEntry.state(["disabled"])
		iRow += 1

		self.OffsetPointsButValue = StringVar()
		self.OffsetPointsButValue.set("Offset points")
		self.OffsetPointsBut = Button(self.Frame, textvariable=self.OffsetPointsButValue, command = self.cmd_offset_points, width = xWidth)
		self.OffsetPointsBut.grid(row=iRow, column=0, pady=nPady, padx=nPadx)

		self.RotatePointsButValue = StringVar()
		self.RotatePointsButValue.set("Rotate points")
		self.RotatePointsBut = Button(self.Frame, textvariable=self.RotatePointsButValue, command = self.cmd_rotate_points, width = xWidth)
		self.RotatePointsBut.grid(row=iRow, column=1, pady=nPady, padx=nPadx)
		iRow += 1

		self.NewWallButValue = StringVar()
		self.NewWallButValue.set("New wall")
		self.NewWallBut = Button(self.Frame, textvariable=self.NewWallButValue, command = self.cmd_new_wall, width = xWidth)
		self.NewWallBut.grid(row=iRow, column=0, pady=nPady, padx=nPadx)

		self.DropWallButValue = StringVar()
		self.DropWallButValue.set("Drop wall")
		self.DropWallBut = Button(self.Frame, textvariable=self.DropWallButValue, command = self.cmd_drop_wall, width = xWidth)
		self.DropWallBut.grid(row=iRow, column=1, pady=nPady, padx=nPadx)

		self.SplitWallButValue = StringVar()
		self.SplitWallButValue.set("Split wall")
		self.SplitWallBut = Button(self.Frame, textvariable=self.SplitWallButValue, command = self.cmd_split_wall, width = xWidth)
		self.SplitWallBut.grid(row=iRow, column=2, pady=nPady, padx=nPadx)

		self.DoubleWallButValue = StringVar()
		self.DoubleWallButValue.set("Duplicate wall")
		self.DoubleWallBut = Button(self.Frame, textvariable=self.DoubleWallButValue, command = self.cmd_double_wall, width = xWidth)
		self.DoubleWallBut.grid(row=iRow, column=3, pady=nPady, padx=nPadx)

		self.ConnectWallsButValue = StringVar()
		self.ConnectWallsButValue.set("Connect walls")
		self.ConnectWallsBut = Button(self.Frame, textvariable=self.ConnectWallsButValue, command = self.cmd_connect_walls, width = xWidth)
		self.ConnectWallsBut.grid(row=iRow, column=4, pady=nPady, padx=nPadx)

		iRow += 1
		self.PropertiesButValue = StringVar()
		self.PropertiesButValue.set("Properties")
		self.PropertiesBut = Button(self.Frame, textvariable=self.PropertiesButValue, command = self.cmd_properties_wall, width = xWidth)
		self.PropertiesBut.grid(row=iRow, column=0, pady=nPady, padx=nPadx)

		iRow += 1
		InfoLabel = Label(self.Frame, text = "Info", width = xWidth)
		InfoLabel.grid(row=iRow, column=0,  pady=nPady, padx=nPadx)

		InfoSeparator = Separator(self.Frame)
		InfoSeparator.grid(row=iRow, column=1,  columnspan=5,  pady=nPady, padx=0, sticky=E+W)
		iRow += 1

		self.AngleButValue = StringVar()
		self.AngleButValue.set("Get Angle")
		self.AngleBut = Button(self.Frame, textvariable=self.AngleButValue, command = self.cmd_info_angle, width = xWidth)
		self.AngleBut.grid(row=iRow, column=0, pady=nPady)

		self.DistanceButValue = StringVar()
		self.DistanceButValue.set("Get Dist")
		self.DistanceBut = Button(self.Frame, textvariable=self.DistanceButValue, command = self.cmd_info_distance, width = xWidth)
		self.DistanceBut.grid(row=iRow, column=1, pady=nPady)
		iRow += 1

		ExitSeparator = Separator(self.Frame)
		ExitSeparator.grid(row=iRow, column=0,  columnspan=6,  pady=nPady, padx=0, sticky=E+W)
		iRow += 1

		self.OKBut = Button(self.Frame, text="OK", command = self.cmd_exit, width = xWidth)
		self.OKBut.grid(row=iRow, column=0, pady=nPady)

		self.CancelBut = Button(self.Frame, text="Cancel", command = self.cmd_cancel, width = xWidth)
		self.CancelBut.grid(row=iRow, column=1, pady=nPady)

#		Starting default values
		self.SaveFile		= False
		self.resize_mode	= True
		self.InButton		= False 
		self.InSaveWall		= False
		self.InAddPoint		= False
		self.InRemovePoint	= False
		self.InMovePoint	= False
		self.InOffsetPoints	= False
		self.InRotatePoints	= False
		self.InProperties	= False
		self.InSplitWall	= False
		self.InConnectWalls	= False
		self.InDoubleWall	= False
		self.InNewWall		= False
		self.InDropWall		= False
		self.InCloseWall	= False
		self.InInfoAngle	= False
		self.InInfoDistance	= False

		self.Walls			= walls_copy(self.parent.Config.Walls)
		self.update_figure()

		self.wait_visibility()
		self.grab_set()
		self.wait_window(self)

#	Internal functions
#===============================================================================

	def cmd_cancel(self):
		self.result	= False
		self.parent.focus_set()
		self.destroy()

	def cmd_exit(self):
		self.result	= True
		self.parent.Config.Walls = walls_copy(self.Walls)
		self.parent.focus_set()
		self.destroy()

	def cmd_donothing(self, empty):
		return

#	Load walls data from file

	def cmd_load_walls(self):
		wall_file = filedialog.askopenfilename(filetypes=(("wall files", "*.txt"), ("All files", "*.*")))
		if(not wall_file): return
		try:
			Walls = load_walls_file(wall_file)
			self.Walls = Walls
			self.update_figure()
		except:
			messagebox.showerror("Read Walls File", "Failed to read file {:s}\n".format(wall_file))
			return


	def cmd_add_walls(self):
		wall_file = filedialog.askopenfilename(filetypes=(("wall files", "*.txt"), ("wall files", "*.csv"), ("All files", "*.*")))
		if(not wall_file): return
		try:
			Walls = load_walls_file(wall_file)
			for Wall in Walls: self.Walls.append(Wall)
			self.update_figure()
		except:
			messagebox.showerror("Read Wall File", "Failed to read file {:s}\n".format(wall_file))
			return

	def cmd_save_walls(self):
		wall_file = filedialog.asksaveasfilename(filetypes=(("wall files", "*.txt"), ("All files", "*.*")))
		if(not wall_file): return
		try:
			save_walls_file(wall_file,self.Walls)
		except:
			messagebox.showerror("Save Walls File", "Failed to save file {:s}\n".format(wall_file))
			return


	def cmd_save_wall(self):
		if(not self.InSaveWall):
			if(self.InButton): return
			self.InButton = self.InSaveWall = self.set_move_button_press_event()
		else:
			self.InButton = self.InSaveWall = self.reset_move_button_press_event()

			wall_file = filedialog.asksaveasfilename(filetypes=(("wall files", "*.txt"), ("All files", "*.*")))
			if(not wall_file): return
			try:
				save_walls_file(wall_file,[self.Walls[self.iWall]])
			except:
				messagebox.showerror("Save Wall File", "Failed to save file {:s}\n".format(wall_file))
				return

	def cmd_ShowPlasmaWallCheckValue(self, a, b, c):
		self.update_figure()

	def cmd_ShowEireneWallCheckValue(self, a, b, c):
		self.update_figure()

	def cmd_ShowExtraWallCheckValue(self, a, b, c):
		self.update_figure()

#	Edit ==============================================	
			
	def cmd_add_point(self):
		if(not self.InAddPoint):
			if(self.InButton): return
			print("\n\t1) Click above the segment on which to add point")
			print(   "\t2) Click on new position or write new coordinate and press Enter")
			print(   "\t3) To end press Cancel\n")
			self.iPtSeg     = -1
			self.AddPointButValue.set("end")
			self.CoordEntry.state(["!disabled"])
			self.InButton = self.InAddPoint = self.set_move_button_press_event()
		else:
			self.AddPointButValue.set("Add point")
			self.CoordEntry.state(["disabled"])
			self.InButton = self.InAddPoint = self.reset_move_button_press_event()
			self.update_figure()

	def cmd_move_point(self):
		if(not self.InMovePoint):
			if(self.InButton): return
			print("\n\t1) Click above the point to select it")
			print(   "\t2) Click on new position or write new coordinate and press Enter\n")
			print(   "\t3) To end press Cancel\n")
			self.iPtSeg     = -1
			self.MovePointButValue.set("End")
			self.CoordEntry.state(["!disabled"])
			self.InButton = self.InMovePoint = self.set_move_button_press_event()
		else:
			self.MovePointButValue.set("Move point")
			self.CoordEntry.state(["disabled"])
			self.InButton = self.InMovePoint = self.reset_move_button_press_event()
			self.update_figure()

	def cmd_remove_point(self):
		if(not self.InRemovePoint):
			if(self.InButton): return
			print("\n\t1) Click above the point to remove it")
			print(   "\t2) To end press Cancel\n")
			self.iPtSeg     = -1
			self.RemovePointButValue.set("End")
			self.InButton = self.InRemovePoint = self.set_move_button_press_event()
		else:
			self.RemovePointButValue.set("Remove point")
			self.InButton = self.InRemovePoint = self.reset_move_button_press_event()
			self.update_figure()
	
	def cmd_offset_points(self):
		if(not self.InOffsetPoints):
			if(self.InButton): return
			print("\n\t1) Select the wall to offset a sub part")
			print(   "\t2) Define rectangolar zone with two points")
			print(   "\t3) Enter offset (DeltaR,DeltaZ)in entry and press return\n")
			print(   "\t4) To end press Cancel\n")
			self.iWall	= -1
			self.iPtSeg	= -1
			self.OffsetPointsButValue.set("Cancel")
			self.CoordEntry.state(["!disabled"])
			self.InButton = self.InOffsetPoints = self.set_move_button_press_event()
		else:
			self.OffsetPointsButValue.set("Offset points")
			self.CoordEntry.state(["disabled"])
			self.InButton = self.InOffsetPoints = self.reset_move_button_press_event()
			self.update_figure()

	def cmd_rotate_points(self):
		if(not self.InRotatePoints):
			if(self.InButton): return
			print("\n\t1) Select the wall to rotate a sub part")
			print(   "\t2) Define rectangolar zone with two points")
			print(   "\t3) Set the rotation point")
			print(   "\t4) Enter angle in degree in entry and press return\n")
			print(   "\t5) To end press Cancel\n")
			self.iWall	= -1
			self.iPtSeg	= -1
			self.RotatePointsButValue.set("Cancel")
			self.CoordEntry.state(["!disabled"])
			self.InButton = self.InRotatePoints = self.set_move_button_press_event()
		else:
			self.RotatePointsButValue.set("Rotate points")
			self.CoordEntry.state(["disabled"])
			self.InButton = self.InRotatePoints = self.reset_move_button_press_event()
			self.update_figure()

#	=====================================================0

	def cmd_properties_wall(self):
		if(not self.InProperties):
			if(self.InButton): return
			print("\n\tClick above the wall to select it")
			self.iWall = -1
			self.PropertiesButValue.set("End")
			self.InButton = self.InProperties = self.set_move_button_press_event()
		else:
			self.PropertiesButValue.set("Properties")
			self.InButton = self.InProperties = self.reset_move_button_press_event()
			if(self.iWall > -1):
				OldType = self.Walls[self.iWall].Type
				self.Walls[self.iWall].Type, self.Walls[self.iWall].LineType = set_properties(self.parent, self.iWall+1, self.Walls[self.iWall].Type, self.Walls[self.iWall].LineType)
				if((OldType == DRAWING_WALL) and (self.Walls[self.iWall].Type != DRAWING_WALL)):
					Wall = self.Walls[self.iWall]
					if((Wall.Rwall[0] != Wall.Rwall[-1]) or (Wall.Zwall[0] != Wall.Zwall[-1])):					#If needed close wall
						Wall.Rwall = np.append(Wall.Rwall, Wall.Rwall[0])
						Wall.Zwall = np.append(Wall.Zwall, Wall.Zwall[0])

					Area = np.sum((Wall.Zwall[1:]+Wall.Zwall[:-1])*(Wall.Rwall[1:]-Wall.Rwall[0:-1]))
					if(Area > 0.):																				#if needed order counter clockwise
						Wall.Rwall = Wall.Rwall[::-1]
						Wall.Zwall = Wall.Zwall[::-1]
					Wall.Clockwise = -1
				
				self.update_figure()

#	=====================================================0

	def cmd_new_wall(self):
		if(not self.InNewWall):
			if(self.InButton): return
			print("\n\t1) Click points of the wall")
			print(   "\t2) close pressing End ")
			self.NewWallButValue.set("End")
			self.CoordEntry.state(["!disabled"])

			Walls = self.Walls
			Walls.append(types.SimpleNamespace())
			Walls[-1].Rwall = np.empty(0, dtype='f8')
			Walls[-1].Zwall = np.empty(0, dtype='f8')
			Walls[-1].Type	= DRAWING_WALL
			Walls[-1].LineType = np.array([0,0,0], dtype='i4')

			self.InButton = self.InNewWall = self.set_move_button_press_event()
		else:
			Walls = self.Walls
			if(len(Walls[-1].Rwall) > 1):
				Walls[-1].SubEZones  = []
				Walls[-1].Clockwise = 0
				Walls[-1].Closed	 = False
				Walls[-1].Changed	 = True
				self.result = True
			else:
				del Walls[len(Walls)-1]

			self.NewWallButValue.set("New Wall")
			self.CoordEntry.state(["disabled"])
			self.InButton = self.InNewWall = self.reset_move_button_press_event()

	def cmd_drop_wall(self):
		if(not self.InDropWall):
			if(self.InButton): return
			print("\n\t1) Click above the wall to remove it")
			self.iPtSeg     = -1
			self.DropWallButValue.set("End")
			self.InButton = self.InDropWall = self.set_move_button_press_event()
		else:
			self.DropWallButValue.set("Drop Wall")
			self.InButton = self.InDropWall = self.reset_move_button_press_event()


	def cmd_split_wall(self):
		if(not self.InSplitWall):
			if(self.InButton): return
			print("\n\t1) Click above the first point of the new wall")
			print(   "\t2) Click above the last point of the new wall")
			self.iPtSet1     = -1
			self.SplitWallButValue.set("Cancel")
			self.InButton = self.InSplitWall = self.set_move_button_press_event()
		else:
			self.SplitWallButValue.set("Split Wall")
			self.InButton = self.InSplitWall = self.reset_move_button_press_event()
			self.update_figure()

	def cmd_double_wall(self):
		if(not self.InDoubleWall):
			if(self.InButton): return
			self.iPtSeg     = -1
			self.DoubleWallButValue.set("Cancel")
		else:
			self.DoubleWallButValue.set("Doplicate Wall")
			self.InButton = self.InDoubleWall = self.reset_move_button_press_event()
			self.update_figure()

	def cmd_connect_walls(self):
		if(not self.InConnectWalls):
			if(self.InButton): return
			print("\n\t1) Click above the first point/(middle segment)/last point of wall section to connect")
			print(   "\t2) Click above the first point/(middle segment)/last point of wall section to be connected")
			print(   "\t3) The second wall will be kept as it is while the first one will be changed.")
			self.iPtSet1	= []
			self.iPtSet2	= []
			self.iWall1		= -1
			self.iWall2		= -1
			self.ConnectWallsButValue.set("Cancel")
			self.CoordEntry.state(["!disabled"])
			self.InButton = self.InConnectWalls = self.set_move_button_press_event()
		else:
			self.ConnectWallsButValue.set("Connect Walls")
			self.CoordEntry.state(["disabled"])
			self.InButton = self.InConnectWalls = self.reset_move_button_press_event()
			self.update_figure()

#	=====================================================0


	def set_move_button_press_event(self):
		self.cid_mouse_button_press = self.Fig.canvas.mpl_connect('button_press_event', self.cmd_move_button_press)
		self.cursor = Cursor(self.Ax, horizOn=True, vertOn=True, useblit=False, color='red', linewidth=1)
		return True

	def reset_move_button_press_event(self):
		if(self.cid_mouse_button_press != 0):
			self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press)
			self.cid_mouse_button_press = 0
			self.cursor=0
		return False


#	######################	
#	Mouses events for branch movement
#	######################

	def cmd_move_button_press(self, event):
		Rpt		 = event.xdata
		Zpt		 = event.ydata

		self.Rpt = Rpt
		self.Zpt = Zpt
		Walls	 = self.Walls
		if(self.InAddPoint):
			if(self.iPtSeg == -1):
				self.iWall, self.iPtSeg = walls_find_closest_segment(self, Rpt, Zpt, self.ShowPlasmaWallCheckValue.get(), self.ShowEireneWallCheckValue.get(), self.ShowExtraWallCheckValue.get())

				self.Ax.plot(Walls[self.iWall].Rwall[self.iPtSeg:self.iPtSeg+2], Walls[self.iWall].Zwall[self.iPtSeg:self.iPtSeg+2],'r-')
				self.Fig.canvas.draw()
			else:
				Walls[self.iWall].Rwall = np.concatenate((Walls[self.iWall].Rwall[:self.iPtSeg+1], np.array([Rpt]), Walls[self.iWall].Rwall[self.iPtSeg+1:]))
				Walls[self.iWall].Zwall = np.concatenate((Walls[self.iWall].Zwall[:self.iPtSeg+1], np.array([Zpt]), Walls[self.iWall].Zwall[self.iPtSeg+1:]))
				if(hasattr(Walls[self.iWall],'SubEZones')): Walls[self.iWall].SubEZones=[]
				Walls[self.iWall].Changed = True
				self.iPtSeg	= -1
				self.result = True
				self.update_figure()

		elif(self.InMovePoint):
			if(self.iPtSeg == -1):
				self.iWall, self.iPtSeg = walls_find_closest_point(self, Rpt, Zpt, self.ShowPlasmaWallCheckValue.get(), self.ShowEireneWallCheckValue.get(), self.ShowExtraWallCheckValue.get())

				self.Ax.plot(Walls[self.iWall].Rwall[self.iPtSeg], Walls[self.iWall].Zwall[self.iPtSeg],'ro')
				self.Fig.canvas.draw()
			else:
				Walls[self.iWall].Rwall[self.iPtSeg] = Rpt
				Walls[self.iWall].Zwall[self.iPtSeg] = Zpt
				if(hasattr(Walls[self.iWall],'SubEZones')): Walls[self.iWall].SubEZones=[]
				Walls[self.iWall].Changed = True
				self.iPtSeg = -1
				self.result = True
				self.update_figure()

		elif(self.InRemovePoint):
			self.iWall, self.iPtSeg = walls_find_closest_point(self, Rpt, Zpt, self.ShowPlasmaWallCheckValue.get(), self.ShowEireneWallCheckValue.get(), self.ShowExtraWallCheckValue.get())
			if(self.iPtSeg == 0):
				if(Walls[self.iWall].Closed):
					Walls[self.iWall].Rwall = np.append(Walls[self.iWall].Rwall[1:-1],Walls[self.iWall].Rwall[1])
					Walls[self.iWall].Zwall = np.append(Walls[self.iWall].Zwall[1:-1],Walls[self.iWall].Zwall[1])
				else:
					Walls[self.iWall].Rwall = Walls[self.iWall].Rwall[1:]
					Walls[self.iWall].Zwall = Walls[self.iWall].Zwall[1:]
			elif (self.iPtSeg+1 == Walls[self.iWall].Rwall[:self.iPtSeg].size):
				if(Walls[self.iWall].Closed):
					Walls[self.iWall].Rwall = np.append(Walls[self.iWall].Rwall[1:-1],Walls[self.iWall].Rwall[1])
					Walls[self.iWall].Zwall = np.append(Walls[self.iWall].Zwall[1:-1],Walls[self.iWall].Zwall[1])
				else:
					Walls[self.iWall].Rwall = Walls[self.iWall].Rwall[:-1]
					Walls[self.iWall].Zwall = Walls[self.iWall].Zwall[:-1]
			else:
				Walls[self.iWall].Rwall = np.append(Walls[self.iWall].Rwall[:self.iPtSeg], Walls[self.iWall].Rwall[self.iPtSeg+1:])
				Walls[self.iWall].Zwall = np.append(Walls[self.iWall].Zwall[:self.iPtSeg], Walls[self.iWall].Zwall[self.iPtSeg+1:])

			if(hasattr(Walls[self.iWall],'SubEZones')): Walls[self.iWall].SubEZones=[]
			Walls[self.iWall].Changed = True
			self.result = True
			self.update_figure()

#	=====================================================0

		elif(self.InOffsetPoints or  self.InRotatePoints):
			if(self.iWall	== -1):
				self.iWall, iPtSeg = walls_find_closest_segment(self, Rpt, Zpt, self.ShowPlasmaWallCheckValue.get(), self.ShowEireneWallCheckValue.get(), self.ShowExtraWallCheckValue.get())
				self.Ax.plot(Walls[self.iWall].Rwall, Walls[self.iWall].Zwall,'r.-')
				self.Fig.canvas.draw()
			elif(self.iPtSeg == -1):
				self.iPtSeg = 0
				self.OldRpt = Rpt
				self.OldZpt = Zpt
				self.Ax.plot(Rpt, Zpt,'go')
				self.Fig.canvas.draw()
			elif(self.iPtSeg == 0):
				if(self.InOffsetPoints ): self.iPtSeg = 2
				else:					  self.iPtSeg = 1
				self.NewRpt = Rpt
				self.NewZpt = Zpt
				self.Ax.plot(Rpt, Zpt,'go')
				self.Ax.plot([self.OldRpt,self.NewRpt,self.NewRpt,self.OldRpt,self.OldRpt], [self.OldZpt,self.OldZpt,self.NewZpt,self.NewZpt,self.OldZpt],'g--')
				self.Fig.canvas.draw()
			elif(self.iPtSeg == 1):
				self.iPtSeg		= 2
				self.CenterRpt	= Rpt
				self.CenterZpt	= Zpt
				self.Ax.plot(Rpt, Zpt,'bo')
				self.Fig.canvas.draw()
			elif(self.iPtSeg == 2):
				Rwall = Walls[self.iWall].Rwall
				Zwall = Walls[self.iWall].Zwall
				Rs = min(self.OldRpt, self.NewRpt)
				Re = max(self.OldRpt, self.NewRpt)
				Zs = min(self.OldZpt, self.NewZpt)
				Ze = max(self.OldZpt, self.NewZpt)
				iPts = np.where((Rwall >= Rs) & (Rwall <= Re) & (Zwall >= Zs) & (Zwall <= Ze))[0]
				if(len(iPts) > 0):
					if(self.InOffsetPoints):
						Rwall[iPts] = Rwall[iPts] + Rpt
						Zwall[iPts] = Zwall[iPts] + Zpt
						self.cmd_offset_points()
					else:
						RadAngle   = event.xdata*np.pi/180.
						DeltaR	   = Rwall[iPts] - self.CenterRpt
						DeltaZ	   = Zwall[iPts] - self.CenterZpt
						Rwall[iPts]= DeltaR*cos(RadAngle) - DeltaZ*sin(RadAngle) + self.CenterRpt
						Zwall[iPts]= DeltaR*sin(RadAngle) + DeltaZ*cos(RadAngle) + self.CenterZpt
						self.cmd_rotate_points()

#	=====================================================0

		if(self.InNewWall):
			Walls[-1].Rwall = np.append(Walls[-1].Rwall, Rpt)
			Walls[-1].Zwall = np.append(Walls[-1].Zwall, Zpt)
			self.update_figure()

		elif(self.InSplitWall):
			if(self.iPtSet1 == -1):
				self.iWall, self.iPtSet1 = walls_find_closest_point(self, Rpt, Zpt, self.ShowPlasmaWallCheckValue.get(), self.ShowEireneWallCheckValue.get(), self.ShowExtraWallCheckValue.get())

				self.Ax.plot(Walls[self.iWall].Rwall[self.iPtSet1], Walls[self.iWall].Zwall[self.iPtSet1],'ro')
				self.Fig.canvas.draw()
			else:
				d, iPtSet2 = find_closest_segment(Rpt, Zpt, Walls[self.iWall].Rwall, Walls[self.iWall].Zwall)
				self.Ax.plot(Walls[self.iWall].Rwall[iPtSet2], Walls[self.iWall].Zwall[iPtSet2],'ro')
				self.Fig.canvas.draw()

				iPt1 = min(iPtSet2,self.iPtSet1)
				iPt2 = max(iPtSet2,self.iPtSet1)

				Walls.append(types.SimpleNamespace())
				dStartEnd = sqrt((Walls[self.iWall].Rwall[-1]-Walls[self.iWall].Rwall[0])**2+(Walls[self.iWall].Zwall[-1]-Walls[self.iWall].Zwall[0])**2)
				if(dStartEnd < 1e-6):																						#closed wall
					Walls[-1].Rwall = np.append(Walls[self.iWall].Rwall[iPt1:iPt2+1],Walls[self.iWall].Rwall[iPt1])
					Walls[-1].Zwall = np.append(Walls[self.iWall].Zwall[iPt1:iPt2+1],Walls[self.iWall].Zwall[iPt1])
				else:
					Walls[-1].Rwall = np.copy(Walls[self.iWall].Rwall[iPt1:iPt2+1])
					Walls[-1].Zwall = np.copy(Walls[self.iWall].Zwall[iPt1:iPt2+1])

				Walls[-1].Type		= Walls[self.iWall].Type
				Walls[-1].LineType	= np.copy(Walls[self.iWall].LineType)
				Walls[-1].Closed	= Walls[self.iWall].Closed
				Walls[-1].Changed	= True

				Walls[self.iWall].Rwall = np.append(Walls[self.iWall].Rwall[:iPt1+1],Walls[self.iWall].Rwall[iPt2:])
				Walls[self.iWall].Zwall = np.append(Walls[self.iWall].Zwall[:iPt1+1],Walls[self.iWall].Zwall[iPt2:])
				Walls[self.iWall].Changed = True

				if(hasattr(Walls[self.iWall],'SubEZones')):
					Walls[-1].SubEZones = []
					Walls[self.iWall].SubEZones = []
				self.iPtSeg	= -1

				self.result = True
				self.update_figure()
				self.cmd_split_wall()

		elif(self.InConnectWalls):
			if(len(self.iPtSet1) == 0):
				self.iWall1, iPt = walls_find_closest_point(self, Rpt, Zpt, self.ShowPlasmaWallCheckValue.get(), self.ShowEireneWallCheckValue.get(), self.ShowExtraWallCheckValue.get())
				self.iPtSet1 = [iPt]
				self.Ax.plot(Walls[self.iWall1].Rwall[iPt], Walls[self.iWall1].Zwall[iPt],'ro')
				self.Ax.text(Walls[self.iWall1].Rwall[iPt], Walls[self.iWall1].Zwall[iPt], "1", color='r', ha="left",fontsize=10)		
				self.Fig.canvas.draw()
			elif(len(self.iPtSet1) < 2):
				iWall1, iPt = walls_find_closest_segment(self, Rpt, Zpt, self.ShowPlasmaWallCheckValue.get(), self.ShowEireneWallCheckValue.get(), self.ShowExtraWallCheckValue.get())
				if(self.iWall1 != iWall1): return
				self.iPtSet1.append(iPt+0.5)
				self.Ax.plot(Walls[iWall1].Rwall[iPt:iPt+2], Walls[iWall1].Zwall[iPt:iPt+2],'r-')
				self.Ax.text(0.5*(Walls[iWall1].Rwall[iPt]+Walls[iWall1].Rwall[iPt+1]), 0.5*(Walls[iWall1].Zwall[iPt]+Walls[iWall1].Zwall[iPt+1]), "2", color='r', ha="left",fontsize=10)		
				self.Fig.canvas.draw()
			elif(len(self.iPtSet1) < 3):
				iWall1, iPt = walls_find_closest_point(self, Rpt, Zpt, self.ShowPlasmaWallCheckValue.get(), self.ShowEireneWallCheckValue.get(), self.ShowExtraWallCheckValue.get())
				if(self.iWall1 != iWall1): return
				self.iPtSet1.append(iPt)
				self.Ax.plot(Walls[iWall1].Rwall[iPt], Walls[iWall1].Zwall[iPt],'ro')
				self.Ax.text(Walls[iWall1].Rwall[iPt], Walls[iWall1].Zwall[iPt], "3", color='r', ha="left",fontsize=10)		
				self.Fig.canvas.draw()
			elif(len(self.iPtSet2) == 0):
				self.iWall2, iPt = walls_find_closest_point(self, Rpt, Zpt, self.ShowPlasmaWallCheckValue.get(), self.ShowEireneWallCheckValue.get(), self.ShowExtraWallCheckValue.get())
				if(self.iWall2 == self.iWall1): return
				self.iPtSet2 = [iPt]
				self.Ax.plot(Walls[self.iWall2].Rwall[self.iPtSet2], Walls[self.iWall2].Zwall[self.iPtSet2],'go')
				self.Ax.text(Walls[self.iWall2].Rwall[self.iPtSet2[0]], Walls[self.iWall2].Zwall[self.iPtSet1[0]], "1", color='g', ha="left",fontsize=10)		
				self.Fig.canvas.draw()
			elif(len(self.iPtSet2) < 2):
				iWall2, iPt = walls_find_closest_segment(self, Rpt, Zpt, self.ShowPlasmaWallCheckValue.get(), self.ShowEireneWallCheckValue.get(), self.ShowExtraWallCheckValue.get())
				if(self.iWall2 != iWall2): return
				self.iPtSet2.append(iPt+0.5)
				self.Ax.plot(Walls[iWall2].Rwall[iPt:iPt+2], Walls[iWall2].Zwall[iPt:iPt+2],'g-')
				self.Ax.text(0.5*(Walls[iWall2].Rwall[iPt]+Walls[iWall2].Rwall[iPt+1]), 0.5*(Walls[iWall2].Zwall[iPt]+Walls[iWall2].Zwall[iPt+1]), "2", color='g', ha="left",fontsize=10)		
				self.Fig.canvas.draw()
			elif(len(self.iPtSet2) < 3):
				iWall2, iPt = walls_find_closest_point(self, Rpt, Zpt, self.ShowPlasmaWallCheckValue.get(), self.ShowEireneWallCheckValue.get(), self.ShowExtraWallCheckValue.get())
				if(self.iWall2 != iWall2): return
				self.iPtSet2.append(iPt)
				self.Ax.plot(Walls[iWall2].Rwall[iPt], Walls[iWall2].Zwall[iPt],'go')
				self.Ax.text(Walls[iWall2].Rwall[iPt], Walls[iWall2].Zwall[iPt], "3", color='g', ha="left",fontsize=10)		
				self.Fig.canvas.draw()
				if(len(self.iPtSet2) == 3):
					iWall1 = self.iWall1
					iWall2 = self.iWall2
					Wall1	= Walls[iWall1]
					Wall2	= Walls[iWall2]
					iPtSet1 = self.iPtSet1
					iPtSet2 = self.iPtSet2

					iPtMin1 = min(iPtSet1[0], iPtSet1[2])
					iPtMax1 = max(iPtSet1[0], iPtSet1[2])
					if((iPtMin1 < iPtSet1[1]) and (iPtMax1 > iPtSet1[1])):	MinToMax1 = 1
					else:													MinToMax1 = -1
					if(iPtMin1 == iPtSet1[2]): MinToMax1 = -MinToMax1

					iPtMin2 = min(iPtSet2[0], iPtSet2[2])
					iPtMax2 = max(iPtSet2[0], iPtSet2[2])				
					if((iPtMin2 < iPtSet2[1]) and (iPtMax2 > iPtSet2[1])):	
						iPts2	  = np.arange(iPtMin2,iPtMax2+1)
						MinToMax2 = 1
					else:
						if(Wall2.Closed): 	iPts2 = np.append(np.arange(iPtMax2-len(Wall2.Rwall),-1), np.arange(0,iPtMin2+1))
						else:				iPts2 = np.arange(iPtMax2-len(Wall2.Rwall),iPtMin2+1)
						MinToMax2 = -1
					if(iPtMin2 == iPtSet2[2]): MinToMax2 = -MinToMax2

					if(MinToMax1 != MinToMax2): iPts2 = iPts2[::-1]

					if((iPtMin1 <  iPtSet1[1]) and (iPtMax1 > iPtSet1[1])):	
						Wall1.Rwall = np.concatenate((Wall1.Rwall[:iPtMin1+1],Wall2.Rwall[iPts2],Wall1.Rwall[iPtMax1:]))
						Wall1.Zwall = np.concatenate((Wall1.Zwall[:iPtMin1+1],Wall2.Zwall[iPts2],Wall1.Zwall[iPtMax1:]))
					else:
						if(Wall1.Closed):
							Wall1.Rwall = np.concatenate((Wall1.Rwall[int(0.5*(iPtMin1+iPtMax1)):iPtMax1+1], Wall2.Rwall[iPts2], Wall1.Rwall[:int(0.5*(iPtMin1+iPtMax1))+1]))								
							Wall1.Zwall = np.concatenate((Wall1.Zwall[int(0.5*(iPtMin1+iPtMax1)):iPtMax1+1], Wall2.Zwall[iPts2], Wall1.Zwall[:int(0.5*(iPtMin1+iPtMax1))+1]))								
						else:
							Wall1.Rwall = np.append(Wall1.Rwall[iPtMin1:iPtMax1+1],	Wall2.Rwall[iPts2])								
							Wall1.Zwall = np.append(Wall1.Zwall[iPtMin1:iPtMax1+1],	Wall2.Zwall[iPts2])	
												
					self.cmd_connect_walls()

		elif(self.InDoubleWall):
			self.iWall, self.iPtSeg = walls_find_closest_point(self, Rpt, Zpt, self.ShowPlasmaWallCheckValue.get(), self.ShowEireneWallCheckValue.get(), self.ShowExtraWallCheckValue.get())
			Walls.append(wall_copy(Walls[self.iWall]))
			Walls[-1].Changed  = True
			self.result = True
			self.update_figure()
			self.cmd_double_wall()

		elif(self.InDropWall):
			self.iWall, self.iPtSeg = walls_find_closest_point(self, Rpt, Zpt, self.ShowPlasmaWallCheckValue.get(), self.ShowEireneWallCheckValue.get(), self.ShowExtraWallCheckValue.get())
			del Walls[self.iWall]
			for iWall in range(self.iWall,len(Walls)): Walls[iWall].Changed = True
			self.result = True
			self.update_figure()
			self.cmd_drop_wall()

#	=====================================================0

		elif(self.InInfoAngle):
			reset	 = self.reset_move_button_press_event()

			iWall, i = walls_find_closest_segment(self, Rpt, Zpt, self.ShowPlasmaWallCheckValue.get(), self.ShowEireneWallCheckValue.get(), self.ShowExtraWallCheckValue.get())
			R1   = Walls[iWall].Rwall[i]; R2 = Walls[iWall].Rwall[i+1]; dR = R2 - R1
			Z1   = Walls[iWall].Zwall[i]; Z2 = Walls[iWall].Zwall[i+1]; dZ = Z2 - Z1
			if((dR != 0.) and (dZ != 0.)):
				ms = (Z2-Z1)/(R2-R1)
				Rc = (Rpt/ms + ms*R1 + Zpt - Z1)/(ms + 1./ms)
				Zc = ms*(Rc - R1) + Z1
			elif(dR != 0):
				Rc = Rpt
				Zc = Z1
			else:
				Rc = R1
				Zc = Zpt

			Config = self.parent.Config
			f_Br = interp.RectBivariateSpline(Config.r2D[0,:], Config.z2D[:,0], Config.Br2D.T)										#RectBivariateSpline wants [x,y] order
			Br	 = f_Br.ev(Rc, Zc); f_Br = None

			f_Bz = interp.RectBivariateSpline(Config.r2D[0,:], Config.z2D[:,0], Config.Bz2D.T)
			Bz	 = f_Bz.ev(Rc, Zc); f_Bz = None

			f_Bphi = interp.RectBivariateSpline(Config.r2D[0,:], Config.z2D[:,0], Config.Bphi2D.T)
			Bphi   = f_Bphi.ev(Rc, Zc); f_Bphi = None

			GrazAngle = abs(asin((-dZ*Br + dR*Bz)/(sqrt(dR**2+dZ**2)*sqrt(Br**2+Bz**2+Bphi**2))))*180./np.pi
			SegAngle  = atan2(dZ, dR)*180./np.pi

			print("\tRc  = ",Rc)
			print("\tZc  = ",Zc)
			print("\tBr  = ",Br)
			print("\tBz  = ",Bz)
			print("\tBphi= ",Bphi)

			messagebox.showinfo("Angles", "Grazing angle = {:f}\n".format(GrazAngle)+
											"Segment angle = {:f}\n".format(SegAngle)+
											"Br            = {:f}\n".format(float(Br))+
											"Bz            = {:f}\n".format(float(Bz))+
											"Bphi          = {:f}".format(float(Bphi)))

			self.cmd_info_angle()

		elif(self.InInfoDistance):
			if(self.iPtSeg == -1):
				self.iPtSeg = 1
				self.OldRpt = Rpt
				self.OldZpt = Zpt
				self.Ax.plot(Rpt, Zpt,'ro')
				self.Fig.canvas.draw()
			else:
				self.Ax.plot([self.OldRpt,Rpt], [self.OldZpt,Zpt],'ro-')
				self.Fig.canvas.draw()
				self.iPtSeg = -1
				print("\tDelta_R  = ",Rpt-self.OldRpt)
				print("\tDelta_Z  = ",Zpt-self.OldZpt)
				print("\tDelat_d  = ",sqrt((Rpt-self.OldRpt)**2+(Zpt-self.OldZpt)**2))
				messagebox.showinfo("Angles", "Delta_d = {:f}\n".format(sqrt((Rpt-self.OldRpt)**2+(Zpt-self.OldZpt)**2))+
												"Delta_R = {:f}\n".format(Rpt-self.OldRpt)+
												"Delta_Z = {:f}\n".format(Zpt-self.OldZpt))

				self.cmd_info_distance()

		elif(self.InSaveWall or self.InProperties):
			self.iWall, i = walls_find_closest_segment(self, Rpt, Zpt, self.ShowPlasmaWallCheckValue.get(), self.ShowEireneWallCheckValue.get(), self.ShowExtraWallCheckValue.get())
			if(self.InSaveWall): self.cmd_save_wall()
			else:				 self.cmd_properties_wall()		 

	
# #######################################
#	Define cut point based on coordinate
# #######################################

	def cmd_enter_coord_entry(self, event):
		Coords = self.CoordEntryValue.get()
		Coords = "".join(Coords.split())
		
		IsRelative = False
		if(Coords[0] == "@"):
			if(self.iPtSeg != -1):	IsRelative = True
			Coords = Coords[1:]
			
		iAngle = Coords.find("<")
		if(iAngle >=0):
			Coords = Coords.replace("<",",")
			IsAngular = True
		else: IsAngular = False
	
		try:
			Values = eval(Coords)
		except:
			return

		try:
			nValues = len(Values)
		except:
			nValues = 1
		
		value =	types.SimpleNamespace()
		if(nValues == 2):
			if(IsAngular):
				value.xdata = Values[0]*cos(Values[1]*np.pi/180.)
				value.ydata = Values[0]*sin(Values[1]*np.pi/180.)
			else:
				value.xdata = Values[0]
				value.ydata = Values[1]
				
			if(IsRelative):
				value.xdata += self.Rpt
				value.ydata += self.Zpt
		elif((nValues == 1) and self.InRotatePoints and (self.iPtSeg == 2)):
			value.xdata = value.ydata = Values
		else:
			return

		self.cmd_move_button_press(value)

	def cmd_info_angle(self):
		if(not self.InInfoAngle):
			if(self.InButton): return
			self.iPtSeg     = -1
			self.AngleButValue.set("Cancel")
			self.InInfoAngle = self.InButton = self.set_move_button_press_event()
		else:
			self.AngleButValue.set("Get Angle")
			self.InInfoAngle = self.InButton = self.reset_move_button_press_event()
		return

	def cmd_info_distance(self):
		if(not self.InInfoDistance):
			if(self.InButton): return
			self.iPtSeg     = -1
			self.DistanceButValue.set("Cancel")
			self.InInfoDistance = self.InButton = self.set_move_button_press_event()
		else:
			self.DistanceButValue.set("Get Distance")
			self.InInfoDistance = self.InButton = self.reset_move_button_press_event()
			self.update_figure()
		return

#=================================================================

	def update_figure(self):
		Fig			= self.Fig
		Ax			= self.Ax
		Config		= self.parent.Config
		X_points	= Config.X_points
		
		Xlims = Ax.get_xlim()
		Ylims = Ax.get_ylim()

		Ax.clear()
		Ax.tick_params(axis='both', which='major', labelsize=8)
		Ax.autoscale(axis='both', tight=True)

		if(Config.equ_OK):
			for k in range(len(X_points)):
				Ax.contour(Config.r2D, Config.z2D, Config.flux2D, colors="r", levels=[X_points[k].psi])
				Ax.text(X_points[k].R,X_points[k].Z,"{:d}".format(X_points[k].index+1),size=18)

		plot2d_walls(Ax, self.Walls, self.ShowPlasmaWallCheckValue.get(), self.ShowEireneWallCheckValue.get(), self.ShowExtraWallCheckValue.get(), marker=".", end_marker=True)

		if((not ((Xlims[0] == 0.) and (Xlims[1] == 1.))) and (not self.resize_mode)):
			Ax.set_xlim(Xlims[0], Xlims[1])
			Ax.set_ylim(Ylims[0], Ylims[1])
		else:
			self.resize_mode = False

		Ax.set_aspect(1.)
		Fig.canvas.draw()


def edit_wall(MainWindow):
	d = EditWall(MainWindow)
	return d.result
