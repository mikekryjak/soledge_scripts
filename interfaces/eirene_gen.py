#!/usr/bin/env python
# -*- coding: utf-8 -*-


from tkinter		import Toplevel, Tk, Text, OptionMenu, filedialog, messagebox, BOTH, W, N, E, S, HORIZONTAL, CENTER
from tkinter		import DoubleVar, IntVar, StringVar
from tkinter.ttk	import Entry, Menubutton, Checkbutton, Frame, Button, Label, Style, LabelFrame, Scale, Separator, Treeview

import types
from tkinter.messagebox import Message

import os
import types

from math									import sqrt
import numpy								as np
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
from routines.intersect_contour				import intersect_2contours

from routines.utils_walls					import plot2d_walls, wall_esubzone_copy, wall_define_eirene_walls

from eirene.generate_vacuum_triangles		import generate_vacuum_triangles
from eirene.generate_eirene_grid			import generate_eirene_grid
from eirene.set_wall_material_pump			import set_wall_material_pump
from eirene.set_wall_puff					import set_wall_puff

from interfaces.ask_name					import ask_name
from interfaces.ask_values					import ask_values


ColorlistMat	= ['Black','Blue','Green','Cyan','Grey','Purple']
ColorlistSurf	= ['Green','Cyan','Grey','Purple','Magenta','Yellow','Orange','Pink']
ColorlistPump	= ['Red','Magenta','Yellow','Orange','Pink','Blue','Green','Cyan','Grey','Purple','Brown','Greenyellow','Crimson','Navy','Magenta','Firebrick']
ColorlistPuff	= ['Red','Magenta','Yellow','Orange','Pink','Blue','Green','Cyan','Grey','Purple','Brown','Greenyellow','Crimson','Navy','Magenta','Firebrick']
#Colors		  = [[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 1, 0], [0.3, 0.3, 0.3], [0.5, 0, 0.5]]
#Colorsp		  = [[1, 0, 0], [1, 0, 1], [1, 0.5, 0], [1, 0.75, 0.8]]

class eirene_gen(Frame):

	def __init__(self, parent, root):
		
		Frame.__init__(self, parent)

		self.root		= root

		self.EireneZonesTreeIds = []
		self.EireneSurfsTreeIds = []
		self.MaterialTreeIds 	= []
		self.PumpTreeIds	 	= []
		self.PuffTreeIds	 	= []
		
		nPadx = 5
		nPady = 3
		
		self.win_frame = Frame(self)
		self.win_frame.pack(fill=BOTH, expand=True)

		self.win_frame.columnconfigure(1, weight=1)
		
#		Plot Eirene area

		EirenePlotFrame = Frame(self.win_frame)
		EirenePlotFrame.grid(row=0, column=0, rowspan=3,  pady=nPady, padx=nPadx, sticky=E+W+N)

		self.Fig = Figure(figsize=(4.8, 5.8), dpi=100)
		self.Fig.patch.set_facecolor('white')

		# a tk.DrawingArea
		canvas = FigureCanvasTkAgg(self.Fig, master=EirenePlotFrame)
		try:
			canvas.show()										#Matplotlib  up 2.2
		except:
			canvas.draw()
		canvas.get_tk_widget().grid(row=1, column=0, columnspan=2,  padx=nPadx, sticky=E+W+S+N)

		self.Ax = self.Fig.add_subplot(111)
		self.Cb = 0														#I need this to be able to remove it later
		self.Fig.tight_layout()

		toolbar_frame = Frame(EirenePlotFrame)
		toolbar_frame.grid(row=0, column=0, columnspan=2, padx=nPadx, sticky=E+W+S+N)
		toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
		toolbar.update()

		ShowEireneLabel = Label(EirenePlotFrame, text = "Show")
		ShowEireneLabel.grid(row=2, column=0,  pady=nPady, padx=nPadx)

		self.ShowFigureOptions = ["Wall", "Wall Triangles", "Triangles","Triangles & Zone"]
		self.ShowFigureOptMenuValue = StringVar()
		self.ShowFigureOptMenuValue.set(self.ShowFigureOptions[0])
		opts = ()
		for k in self.ShowFigureOptions: opts += (k,)	
		self.ShowFigureOptMenu = OptionMenu(EirenePlotFrame, self.ShowFigureOptMenuValue, *opts)
		self.ShowFigureOptMenu.grid(row=2, column=1,  pady=nPady, padx=nPadx, sticky=W+E)
		self.ShowFigureOptMenuValue.trace("w", self.cmd_show_in_figure)


#		Grid 0,1
#		=====

		
#		Eirene vacuum parameters
#		==============
		
		EireneFrame = LabelFrame(self.win_frame, text="Eirene triangles")
		EireneFrame.grid(row=0, column=1,  pady=nPady, padx=nPadx, sticky=E+W+N)

		iFrow = 0
		xWidth = 11
		cWidth = 14

		self.EireneZonesTree = Treeview(EireneFrame, height=4)
		self.EireneZonesTree.grid(row=iFrow, column=0, columnspan=3,  pady=nPady, padx=nPadx, sticky=E+W+S+N)
		self.tree_view_define(self.EireneZonesTree, ["N.", "dR knots", "dZ knots"], [CENTER, CENTER, CENTER], [cWidth,cWidth,cWidth], [CENTER, CENTER, CENTER])

		iFrow += 1
		
		self.EireneSplitZoneButValue = StringVar()
		self.EireneSplitZoneButValue.set("Split Zone")
		EireneSplitZoneBut = Button(EireneFrame, textvariable=self.EireneSplitZoneButValue, command = self.cmd_eirene_split_zone, width = xWidth)
		EireneSplitZoneBut.grid(row=iFrow, column=0, pady=nPady, padx=nPadx)
		
		self.EireneDropSplitZoneButValue = StringVar()
		self.EireneDropSplitZoneButValue.set("Drop Split")
		EireneDropSplitZoneBut = Button(EireneFrame, textvariable=self.EireneDropSplitZoneButValue, command = self.cmd_eirene_drop_split_zone, width = xWidth)
		EireneDropSplitZoneBut.grid(row=iFrow, column=1, pady=nPady, padx=nPadx)
		
		self.EireneSetZoneBut = Button(EireneFrame, text="Set", command = self.cmd_eirene_set_zone, width = xWidth)
		self.EireneSetZoneBut.grid(row=iFrow, column=2, pady=nPady, padx=nPadx)

		iFrow += 1
		
		self.ToCenterCkButValue = IntVar()
		self.ToCenterCkButValue.set(0)
		ToCenterCkBut = Checkbutton(EireneFrame, text="To center", variable = self.ToCenterCkButValue)
		ToCenterCkBut.grid(row=iFrow,column=0, padx=nPadx,  pady=nPady)
		
		GenerateBut = Button(EireneFrame, text="Generate triangles", command = self.cmd_generate_triangles, width = 2*xWidth+4)
		GenerateBut.grid(row=iFrow, column=1, columnspan=2, pady=nPady, padx=nPadx)


#		Grid 0,2
#		=====

#		Eirene surfaces frame
		
		EireneSurfsFrame = LabelFrame(self.win_frame, text="Eirene surfaces")
		EireneSurfsFrame.grid(row=0, column=2,  pady=nPady, padx=nPadx, sticky=E+W+N)

		iFrow = 0
		xWidth = 11
		cWidth = 14

		self.EireneSurfsTree = Treeview(EireneSurfsFrame, height=4)
		self.EireneSurfsTree.grid(row=iFrow, column=0, columnspan=3,  pady=nPady, padx=nPadx, sticky=E+W+S+N)
		self.tree_view_define(self.EireneSurfsTree, ["Name", "R. Coef.", "Color"], [CENTER, CENTER, CENTER], [cWidth, cWidth,cWidth], [CENTER, CENTER, CENTER])

		iFrow += 1
		
		AddEireneSurfsBut = Button(EireneSurfsFrame, text="Add", command = self.cmd_add_surface, width = xWidth)
		AddEireneSurfsBut.grid(row=iFrow, column=0, pady=nPady, padx=nPadx)
		
		DropEireneSurfsBut = Button(EireneSurfsFrame, text="Drop", command = self.cmd_drop_surface, width = xWidth)
		DropEireneSurfsBut.grid(row=iFrow, column=1, pady=nPady, padx=nPadx)
		
		RenaEireneSurfsBut = Button(EireneSurfsFrame, text="Rename", command = self.cmd_rename_surface, width = xWidth)
		RenaEireneSurfsBut.grid(row=iFrow, column=2, pady=nPady, padx=nPadx)

		iFrow += 1

		self.SetEireneSurfsButValue = StringVar()
		self.SetEireneSurfsButValue.set("Set")
		SetEireneSurfsBut = Button(EireneSurfsFrame, textvariable=self.SetEireneSurfsButValue, command = self.cmd_set_surface, width = xWidth)
		SetEireneSurfsBut.grid(row=iFrow, column=0, pady=nPady, padx=nPadx)

#		Grid 1,1
#		=====
		
#		Wall frame
		
		MaterialsFrame = LabelFrame(self.win_frame, text="Plasma wall materials")
		MaterialsFrame.grid(row=1, column=1,  pady=nPady, padx=nPadx, sticky=E+W+N)

		iFrow = 0
		xWidth = 11
		cWidth = 20

		self.MaterialTree = Treeview(MaterialsFrame, height=4)
		self.MaterialTree.grid(row=iFrow, column=0, columnspan=3,  pady=nPady, padx=nPadx, sticky=E+W+S+N)
		self.tree_view_define(self.MaterialTree, ["Material", "Color"], [CENTER, CENTER], [cWidth,cWidth], [CENTER, CENTER])

		iFrow += 1
		
		AddMaterialBut = Button(MaterialsFrame, text="Add", command = self.cmd_add_material, width = xWidth)
		AddMaterialBut.grid(row=iFrow, column=0, pady=nPady, padx=nPadx)
		
		DropMaterialBut = Button(MaterialsFrame, text="Drop", command = self.cmd_drop_material, width = xWidth)
		DropMaterialBut.grid(row=iFrow, column=1, pady=nPady, padx=nPadx)

		
		RenaMaterialBut = Button(MaterialsFrame, text="Rename", command = self.cmd_rename_material, width = xWidth)
		RenaMaterialBut.grid(row=iFrow, column=2, pady=nPady, padx=nPadx)

		iFrow += 1
		
		SetMaterialBut = Button(MaterialsFrame, text="Set", command = self.cmd_set_material, width = xWidth)
		SetMaterialBut.grid(row=iFrow, column=0, pady=nPady, padx=nPadx)

#		Grid 1,2
#		=====
	
#		Pumps frame
		
		PumpsFrame = LabelFrame(self.win_frame, text="Pumps")
		PumpsFrame.grid(row=1, column=2,  pady=nPady, padx=nPadx, sticky=E+W+N)

		iFrow = 0
		xWidth = 11
		cWidth = 20

		self.PumpTree = Treeview(PumpsFrame, height=4)
		self.PumpTree.grid(row=iFrow, column=0, columnspan=3,  pady=nPady, padx=nPadx, sticky=E+W+S+N)
		self.tree_view_define(self.PumpTree, ["Pump", "Color"], [CENTER, CENTER], [cWidth,cWidth], [CENTER, CENTER])
	
		iFrow += 1
		AddPumpBut = Button(PumpsFrame, text="Add", command = self.cmd_add_pump, width = xWidth)
		AddPumpBut.grid(row=iFrow, column=0, pady=nPady, padx=nPadx)
		
		DropPumpBut = Button(PumpsFrame, text="Drop", command = self.cmd_drop_pump, width = xWidth)
		DropPumpBut.grid(row=iFrow, column=1, pady=nPady, padx=nPadx)
		
		RenaPumpBut = Button(PumpsFrame, text="Rename", command = self.cmd_rename_pump, width = xWidth)
		RenaPumpBut.grid(row=iFrow, column=2, pady=nPady, padx=nPadx)

		iFrow += 1
		
		SetPumpBut = Button(PumpsFrame, text="Set", command = self.cmd_set_pump, width = xWidth)
		SetPumpBut.grid(row=iFrow, column=0, pady=nPady, padx=nPadx)


#		Grid 2,1
#		=====
	
#		Puffs frame
		
		PuffsFrame = LabelFrame(self.win_frame, text="Puffs")
		PuffsFrame.grid(row=2, column=1,  pady=nPady, padx=nPadx, sticky=E+W+N)

		iFrow = 0
		xWidth = 11
		cWidth = 20

		self.PuffTree = Treeview(PuffsFrame, height=4)
		self.PuffTree.grid(row=iFrow, column=0, columnspan=3,  pady=nPady, padx=nPadx, sticky=E+W+S+N)
		self.tree_view_define(self.PuffTree,  ["Puff", "Color"], [CENTER, CENTER], [cWidth,cWidth], [CENTER, CENTER])

		iFrow += 1
		
		AddPuffBut = Button(PuffsFrame, text="Add", command = self.cmd_add_puff, width = xWidth)
		AddPuffBut.grid(row=iFrow, column=0, pady=nPady, padx=nPadx)
		
		DropPuffBut = Button(PuffsFrame, text="Drop", command = self.cmd_drop_puff, width = xWidth)
		DropPuffBut.grid(row=iFrow, column=1, pady=nPady, padx=nPadx)
		
		RenaPuffBut = Button(PuffsFrame, text="Rename", command = self.cmd_rename_puff, width = xWidth)
		RenaPuffBut.grid(row=iFrow, column=2, pady=nPady, padx=nPadx)

		iFrow += 1
		SetPuffBut = Button(PuffsFrame, text="Set", command = self.cmd_set_puff, width = xWidth)
		SetPuffBut.grid(row=iFrow, column=0, pady=nPady, padx=nPadx)


#	End Grid
#	=====
#		--------------------------------------------------
		self.EWalls	  		= []

		self.xPtSet	= np.zeros((3), dtype='f8')
		self.yPtSet	= np.zeros((3), dtype='f8')

		self.InSetSurface		= False
		self.InEireneSplitZone	= False

					
#	Internal functions
#===============================================================================
		
	def cmd_show_in_figure(self, a, b, c):
		self.update_figure()

#	Generate Plasma triangles command

	def cmd_generate_triangles(self):
		if(self.root.Config.Mesh_OK):
			ToTheCenter = self.ToCenterCkButValue.get() 
			self.root.Eirene = generate_eirene_grid(self.root, self.root.Config, self.EWalls, ToTheCenter)
			self.eirene_zones_tree_view_fill()
			self.root.set_Eirene_OK(self.root.Config, True)

#	Materials commands
#	============
		
	def cmd_add_material(self):
		Values=ask_values(self, ["Material"], title="Materials")
		if(len(Values) > 0):
			self.root.Eirene.Wall.Material.append(types.SimpleNamespace())
			self.root.Eirene.Wall.Material[-1].Name = Values[0]
			self.material_tree_view_fill()
		return

	def cmd_drop_material(self):
		if(len(self.root.Eirene.Wall.Material) == 1): return
		
		Id = self.MaterialTree.focus()
		if(len(Id)>0):
			ind = self.MaterialTreeIds.index(Id)
			
			del self.root.Eirene.Wall.Material[ind]
			self.material_tree_view_fill()
			
			iMat = np.where(~self.root.Eirene.Wall.IsPump & (self.root.Eirene.Wall.TypeMat >= ind))[0]
			if(len(iMat) > 0):
				self.root.Eirene.Wall.TypeMat[iMat] = np.where(self.root.Eirene.Wall.TypeMat > 0, self.root.Eirene.Wall.TypeMat-1, 0)
				self.update_figure()
		return

	def cmd_rename_material(self):
		Id = self.MaterialTree.focus()
		if(len(Id)>0):
			ind = self.MaterialTreeIds.index(Id)
			Values=ask_values(self, ["Material"], [self.root.Eirene.Wall.Material[ind].Name], title="Materials")
			if(len(Values) > 0):
				self.root.Eirene.Wall.Material[ind].Name = Values[0]
				self.material_tree_view_fill()
		return

	def cmd_set_material(self):
		Id = self.MaterialTree.focus()
		if(len(Id)>0):
			ind = self.MaterialTreeIds.index(Id)
			self.change_wall(ind, False)
		return

#	Pumps commands
#	===========

	def cmd_add_pump(self):
		Values=ask_values(self, ["Name"], title="Pumps")
		if(len(Values) > 0):
			self.root.Eirene.Wall.Pump.append(types.SimpleNamespace())
			self.root.Eirene.Wall.Pump[-1].Name = Values[0]
			self.pump_tree_view_fill()
		return

	def cmd_drop_pump(self):
		Id = self.PumpTree.focus()
		if(len(Id)>0):
			ind = self.PumpTreeIds.index(Id)
			del self.root.Eirene.Wall.Pump[ind]
			self.pump_tree_view_fill()			
		
			iPump = np.where(self.root.Eirene.Wall.IsPump & (self.root.Eirene.Wall.TypeMat == ind))[0]
			if(len(iPump) > 0):
				self.root.Eirene.Wall.IsPump[iPump]  = False
				self.root.Eirene.Wall.TypeMat[iPump] = 0

			iPump = np.where(self.root.Eirene.Wall.IsPump & (self.root.Eirene.Wall.TypeMat > ind))[0]
			if(len(iPump) > 0): self.root.Eirene.Wall.TypeMat[iPump] -= 1
			self.update_figure()
		return
	
	def cmd_rename_pump(self):
		Id = self.PumpTree.focus()
		if(len(Id)>0):
			ind = self.PumpTreeIds.index(Id)
			Values = ask_values(self, ["Name"], [self.root.Eirene.Wall.Pump[ind].Name], title="Pump")
			if(len(Values) > 0):
				self.root.Eirene.Wall.Pump[ind].Name = Values[0]
				self.pump_tree_view_fill()
		return

	def cmd_set_pump(self):
		Id = self.PumpTree.focus()
		if(len(Id)>0):
			ind = self.PumpTreeIds.index(Id)
			self.change_wall(ind, True)
		return

#	Puffs commands
#	===========

	def cmd_add_puff(self):
		Values = ask_values(self,["Name:","R:","Z:"], title="Puff_{:d}".format(len(self.root.Eirene.Wall.Puff)+1))		
		if(len(Values) == 3):
			try:
				set_wall_puff(self.root.Eirene, -1, eval(Values[1]), eval(Values[2]))
				self.PuffInd = -1														#new puff
				self.root.Eirene.Wall.Puff[-1].Name = Values[0]
				
				self.puff_tree_view_fill()
				self.update_figure()
			except:
				pass

	def cmd_drop_puff(self):
		Id = self.PuffTree.focus()
		if(len(Id)>0):
			ind = self.PuffTreeIds.index(Id)
			del self.root.Eirene.Wall.Puff[ind]
			self.puff_tree_view_fill()
			self.update_figure()
		return

	def cmd_rename_puff(self):
		Id = self.PuffTree.focus()
		if(len(Id)>0):
			ind = self.PuffTreeIds.index(Id)
			Values = ask_values(self, ["Name"], [self.root.Eirene.Wall.Puff[ind].Name], title="Puffings")
			if(len(Values) > 0):
				self.root.Eirene.Wall.Puff[ind].Name = Values[0]
				self.puff_tree_view_fill()
				self.update_figure()
		return

	def cmd_set_puff(self):
		if(self.ShowFigureOptMenuValue.get() != self.ShowFigureOptions[0]): return

		self.PuffInd = -1										#New puff
		self.cid_mouse_button_press_set_puff = self.Fig.canvas.mpl_connect('button_press_event', self.cmd_button_press_set_puff)
		self.cursor = Cursor(self.Ax, horizOn=True, vertOn=True, useblit=False, color='red', linewidth=1)


#	In plasma surfaces commands
#	==================
		
	def cmd_add_surface(self):
		if(self.root.Config.Eirene_OK):
			Values=ask_values(self, ["Name","R. Coef"], title="Surfaces")
			if(len(Values) > 0):
				self.root.Eirene.Surfaces.append(types.SimpleNamespace())
				self.root.Eirene.Surfaces[-1].Name  = Values[0]
				self.root.Eirene.Surfaces[-1].Rcoef = eval(Values[1])
				self.root.Eirene.Surfaces[-1].R	   = []
				self.root.Eirene.Surfaces[-1].Z	   = []
				self.surface_tree_view_fill()
		else:
			print("\n\tATTENTION!! You can set surfaces after grid generation")

		return

	def cmd_drop_surface(self):
		if(not self.root.Config.Eirene_OK): 		return
		if(len(self.root.Eirene.Surfaces) == 0):	return
		
		Id = self.EireneSurfsTree.focus()
		if(len(Id)>0):
			ind = self.EireneSurfsTreeIds.index(Id)
			
			del self.root.Eirene.Surfaces[ind]
			self.surface_tree_view_fill()
			self.update_figure()
		else:
			print("\n\tATTENTION!! You must select one surface")

		return

	def cmd_rename_surface(self):
		if(not self.root.Config.Eirene_OK): 		return
		if(len(self.root.Eirene.Surfaces) == 0):	return

		Id = self.EireneSurfsTree.focus()
		if(len(Id)>0):
			ind = self.EireneSurfsTreeIds.index(Id)
			Values=ask_values(self, ["Name", "R. coef."], [self.root.Eirene.Surfaces[ind].Name,"{:.3f}".format(self.root.Eirene.Surfaces[ind].Rcoef)], title="Surfaces")
			if(len(Values) > 0):
				self.root.Eirene.Surfaces[ind].Name  = Values[0]
				self.root.Eirene.Surfaces[ind].Rcoef = eval(Values[1])
				self.surface_tree_view_fill()
		else:
			print("\n\tATTENTION!! You must select one surface")

		return

	def cmd_set_surface(self):
		if(not self.root.Config.Eirene_OK): 		return
		if(len(self.root.Eirene.Surfaces) == 0):	return

		if(self.ShowFigureOptMenuValue.get() != self.ShowFigureOptions[2]):
			print("\n\tATTENTION!! You must set plot in Triangles mode")
			return

		if(not self.InSetSurface):
			Id = self.EireneSurfsTree.focus()
			if(len(Id)>0):
				self.InSetSurface = True
				self.SetEireneSurfsButValue.set("End")
				self.IndSurf	= self.EireneSurfsTree.index(Id)
				self.xPtSurf	= []
				self.yPtSurf	= []

				self.cid_mouse_button_press_set_surf_points = self.Fig.canvas.mpl_connect('button_press_event', self.cmd_button_press_set_surf_points)
				self.cursor = Cursor(self.Ax, horizOn=True, vertOn=True, useblit=False, color='red', linewidth=1)
			else:
				print("\n\tATTENTION!! You must select one surface")
		else:
			self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press_set_surf_points)
			self.cid_mouse_button_press_set_surf_points = 0
			self.cursor=0
			self.root.Eirene.Surfaces[self.IndSurf].R = np.array(self.xPtSurf)
			self.root.Eirene.Surfaces[self.IndSurf].Z = np.array(self.yPtSurf)
			self.InSetSurface = False
			self.SetEireneSurfsButValue.set("Set")
			self.update_figure()


	def cmd_button_press_set_surf_points(self, event):

		if(self.root.EireneSurfsPointKnotsCheckValue.get() == 1):
			iKnot = np.argmin((self.root.Eirene.RKnots - event.xdata)**2 + (self.root.Eirene.ZKnots - event.ydata)**2)
			self.xPtSurf.append(self.root.Eirene.RKnots[iKnot])
			self.yPtSurf.append(self.root.Eirene.ZKnots[iKnot])
		else:
			self.xPtSurf.append(event.xdata)
			self.yPtSurf.append(event.ydata)

		self.update_figure()
		self.Ax.plot(self.xPtSurf,self.yPtSurf,'r.-')
		self.Fig.canvas.draw()


#	eirene zones commands
#	==================
		

	def cmd_eirene_split_zone(self):
		if(not self.InEireneSplitZone):
			Id = self.EireneZonesTree.focus()
			if(len(Id)>0):
				print("\n\t1) Click outside the zone border")
				print(   "\t2) draw cut left clicking on points")
				print(   "\t3) After the last point outside zone border push button Close")
				print(   "\t   To remove last point click with right button")
				print(   "\t   you can start and/or end with a zone border point clicking with central button")
				self.InEireneSplitZone = True
				self.EireneSplitZoneButValue.set("Close")
				self.EireneDropSplitZoneButValue.set("Cancel")
				ind = self.EireneZonesTreeIds.index(Id)
				self.iEWallSplit	=  self.EireneZonesIndexes[ind][0]
				self.iESubZoneSplit =  self.EireneZonesIndexes[ind][1]
				self.xPtSplitZone	= []
				self.yPtSplitZone	= []

				self.cid_mouse_button_press_set_eirene_split_zone_points = self.Fig.canvas.mpl_connect('button_press_event', self.cmd_button_press_set_eirene_split_zone_points)
				self.cursor = Cursor(self.Ax, horizOn=True, vertOn=True, useblit=False, color='red', linewidth=1)
			else:
				print("\n\tATTENTION!! You must select the zone to split")

		else:
			self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press_set_eirene_split_zone_points)
			self.cid_mouse_button_press_set_eirene_split_zone_points = 0
			self.cursor=0
			self.InEireneSplitZone  = False
			self.EireneSplitZoneButValue.set("Split Zone")
			self.EireneDropSplitZoneButValue.set("Drop Splits")

			EWall 	 = self.EWalls[self.iEWallSplit]
			ESubZone = EWall.ESubZones[self.iESubZoneSplit]
			if(len(self.yPtSplitZone) > 1):
				self.xPtSplitZone = np.array(self.xPtSplitZone)
				self.yPtSplitZone = np.array(self.yPtSplitZone)
				if(self.LastCentral):
					d2 = (ESubZone.Rbord - self.xPtSplitZone[0])**2 + (ESubZone.Zbord - self.yPtSplitZone[0])**2
					is_end = np.argmin(d2)
					xi_end = -1
				else:
					xi_end, yi_end, is1, is_end	= intersect_2contours(self.xPtSplitZone[-2:], self.yPtSplitZone[-2:], ESubZone.Rbord, ESubZone.Zbord)
					if(len(xi_end) == 1):
						xi_end = xi_end[0]
						yi_end = yi_end[0]
						is_end = is_end[0]
					elif((len(xi_end) == 2) and (len(self.yPtSplitZone) == 2)):
						i_end = np.argmax((xi_end- self.xPtSplitZone[0])**2 + (yi_end - self.yPtSplitZone[0])**2)
						xi_end = xi_end[i_end]
						yi_end = yi_end[i_end]
						is_end = is_end[i_end]
					elif(len(xi_end) != 1):
						messagebox.showerror("Split zone", "Invalid starting intersection with zone")
						self.update_figure()
						return
			else:
				self.update_figure()
				self.Fig.canvas.draw()
				return

			se_point = [0, 0]
			se_is	 = [self.is_start, is_end]
			if(self.xi_start != -1):
				self.xPtSplitZone[0]  = self.xi_start
				self.yPtSplitZone[0]  = self.yi_start
			else:
				se_point[0]		  	  = -1
				self.xPtSplitZone  = self.xPtSplitZone[1:]
				self.yPtSplitZone  = self.yPtSplitZone[1:]

			if(xi_end != -1):
				self.xPtSplitZone[-1]  = xi_end
				self.yPtSplitZone[-1]  = yi_end
			else:
				se_point[1]		  	   = True
				self.xPtSplitZone  = self.xPtSplitZone[:-1]
				self.yPtSplitZone  = self.yPtSplitZone[:-1]
	

			if(self.is_start > is_end):
				self.xPtSplitZone = self.xPtSplitZone[::-1]
				self.yPtSplitZone = self.yPtSplitZone[::-1]
				se_point = se_point[::-1]
				se_is	 = se_is[::-1]
			

			NewESubZone = wall_esubzone_copy(ESubZone)

			NewESubZone.Rbord = np.append(ESubZone.Rbord[se_is[0]+1+se_point[0]:se_is[1]+1],	self.xPtSplitZone[::-1])
			NewESubZone.Zbord = np.append(ESubZone.Zbord[se_is[0]+1+se_point[0]:se_is[1]+1],	self.yPtSplitZone[::-1])
			NewESubZone.Rbord = np.append(NewESubZone.Rbord, NewESubZone.Rbord[0])
			NewESubZone.Zbord = np.append(NewESubZone.Zbord, NewESubZone.Zbord[0])

			ESubZone.Rbord = np.concatenate((ESubZone.Rbord[:se_is[0]+1],	self.xPtSplitZone,	 ESubZone.Rbord[se_is[1]+1+se_point[1]:]))
			ESubZone.Zbord = np.concatenate((ESubZone.Zbord[:se_is[0]+1],	self.yPtSplitZone,	 ESubZone.Zbord[se_is[1]+1+se_point[1]:]))
			EWall.ESubZones.append(NewESubZone)
			self.eirene_zones_tree_view_fill()

		self.update_figure()

	def cmd_button_press_set_eirene_split_zone_points(self, event):
		if((event.button == 3) and (len(self.xPtSplitZone) > 0)):
			del self.xPtSplitZone[-1]
			del self.yPtSplitZone[-1]
			self.update_figure()
			if(len(self.yPtSplitZone) > 0):	self.Ax.plot(self.xPtSplitZone,self.yPtSplitZone,'r.-')
			self.Fig.canvas.draw()
			return
				
		self.xPtSplitZone.append(event.xdata)
		self.yPtSplitZone.append(event.ydata)
		self.LastCentral = False

		self.update_figure()
		self.Ax.plot(self.xPtSplitZone,self.yPtSplitZone,'r.-')
		self.Fig.canvas.draw()

		if(len(self.yPtSplitZone) == 2):
			ESubZone = self.EWalls[self.iEWallSplit].ESubZones[self.iESubZoneSplit]
			if(self.LastCentral):
				d2 = (ESubZone.Rbord - self.xPtSplitZone[0])**2 + (ESubZone.Zbord - self.yPtSplitZone[0])**2
				self.is_start = np.argmin(d2)
				self.xi_start = -1
			else:
				xi_start, yi_start, is1, is_start = intersect_2contours(np.array(self.xPtSplitZone[:2]), np.array(self.yPtSplitZone[:2]), ESubZone.Rbord, ESubZone.Zbord)
				if(len(xi_start) == 1):
					self.xi_start = xi_start[0]
					self.yi_start = yi_start[0]
					self.is_start = is_start[0]
				elif(len(xi_start) == 2):
					i_start = np.argmin((xi_start - self.xPtSplitZone[0])**2 + (yi_start- self.yPtSplitZone[0])**2)
					self.xi_start = xi_start[i_start]
					self.yi_start = yi_start[i_start]
					self.is_start = is_start[i_start]
				else:
					messagebox.showerror("Split sone", "Invalid starting intersection with zone")
					self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press_set_eirene_split_zone_points)
					self.cid_mouse_button_press_set_eirene_split_zone_points = 0
					self.cursor=0
					self.update_figure()
					self.Fig.canvas.draw()
					self.InEireneSplitZone  = False
					self.EireneSplitZoneButValue.set("Split Zone")
					self.EireneDropSplitZoneButValue.set("Drop Splits")
					return

		if(event.button == 2): self.LastCentral = True

	def cmd_eirene_drop_split_zone(self):
		if(not self.InEireneSplitZone):
			Id = self.EireneZonesTree.focus()
			if(len(Id)>0):
				ind = self.EireneZonesTreeIds.index(Id)
				self.iEWallSplit   =  self.EireneZonesIndexes[ind][0]
				self.iESubZoneSplit =  self.EireneZonesIndexes[ind][1]
				EWall = self.EWalls[self.iEWallSplit]
				CWall  = EWall.Wall															#Config wall
				EWall.ESubZones = []
				EWall.ESubZones.append(types.SimpleNamespace())
				EWall.ESubZones[-1].dr_knots = (CWall.Rwall.max()-CWall.Rwall.min())/20
				EWall.ESubZones[-1].dz_knots = (CWall.Zwall.max()-CWall.Zwall.min())/20
				EWall.ESubZones[-1].Rbord	 = np.copy(CWall.Rwall) 
				EWall.ESubZones[-1].Zbord	 = np.copy(CWall.Zwall)
				self.eirene_zones_tree_view_fill()
			else:
				print("\n\tATTENTION!! You must select one of the zones to drop")

		else:
			self.InEireneSplitZone  = False
			self.EireneSplitZoneButValue.set("Split Zone")
			self.EireneDropSplitZoneButValue.set("Drop Splits")
			self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press_set_eirene_split_zone_points)
			self.cid_mouse_button_press_set_eirene_split_zone_points = 0
			self.cursor=0

		self.update_figure()
		return

	def cmd_eirene_set_zone(self):
		Id = self.EireneZonesTree.focus()
		if(len(Id)>0):
			EWalls = self.EWalls
			ind = self.EireneZonesTreeIds.index(Id)
			iEWall	  =  self.EireneZonesIndexes[ind][0]
			iESubZone =  self.EireneZonesIndexes[ind][1]
			Values=ask_values(self, ["dr knots","dz knots"], ["{:.4f}".format(EWalls[iEWall].ESubZones[iESubZone].dr_knots),"{:.4f}".format(EWalls[iEWall].ESubZones[iESubZone].dz_knots)], title="Eirene zone")
			if(len(Values) > 0):
				EWalls[iEWall].ESubZones[iESubZone].dr_knots = eval(Values[0])
				EWalls[iEWall].ESubZones[iESubZone].dz_knots = eval(Values[1])
				self.eirene_zones_tree_view_fill()
		else:
			print("\n\tATTENTION!! You must select the zones to set")
		return
				
#	Tree functions

	def tree_view_define(self, Tree, Headings, HeadAnchors, xWidths, TextAnchors):
		Tmp = ()
		for k in range(len(Headings)): Tmp += ("{:d}".format(k),)	
		Tree["columns"] = Tmp

		Tree['show'] = 'headings'
		for k in range(len(Headings)): 
			Tree.column("{:d}".format(k),  width=xWidths[k], anchor=TextAnchors[k])
			Tree.heading("{:d}".format(k), text = Headings[k], anchor=HeadAnchors[k])

	def tree_view_fill(self, Tree, TreeIds, Values):
		if(len(TreeIds) > 0):
			Tree.delete(*TreeIds)
			TreeIds = []
			
		for k in range(len(Values[0])):
			TreeValues = ()
			for i in range(len(Values)): TreeValues += (Values[i][k],)	
			TreeIds.append(Tree.insert("" , "end",    text="", values=TreeValues))
		return TreeIds
																	
		
	def material_tree_view_fill(self):
		if(len(self.root.Eirene.Wall.Material) > 0): 
			MatNames = [self.root.Eirene.Wall.Material[k].Name for k in range(len(self.root.Eirene.Wall.Material))]
			self.MaterialTreeIds = self.tree_view_fill(self.MaterialTree,  self.MaterialTreeIds, (MatNames, ColorlistMat[:len(MatNames)]))
		else:
			self.MaterialTreeIds = self.tree_view_fill(self.MaterialTree,  self.MaterialTreeIds, ([], []))
		
	def eirene_zones_tree_view_fill(self):
		EWalls = self.EWalls
		if(len(EWalls) > 0):
			self.EireneZonesIndexes = []
			SubZoneNames  = []
			dr_knots	  = []
			dz_knots	  = []
			for iEWall in range(len(EWalls)):
				ESubZones = EWalls[iEWall].ESubZones
				for iESubZone in range(len(ESubZones)):
					self.EireneZonesIndexes.append([iEWall,iESubZone])
					SubZoneNames.append("{:d},{:d}".format(EWalls[iEWall].iWall+1,iESubZone+1))
					dr_knots.append("{:0.4f}".format(ESubZones[iESubZone].dr_knots))
					dz_knots.append("{:0.4f}".format(ESubZones[iESubZone].dz_knots))
			self.EireneZonesTreeIds = self.tree_view_fill(self.EireneZonesTree,  self.EireneZonesTreeIds, (SubZoneNames, dr_knots, dz_knots))
		else:
			self.EireneZonesTreeIds = self.tree_view_fill(self.EireneZonesTree,  self.EireneZonesTreeIds, ([], [], []))

		
	def surface_tree_view_fill(self):
		if(len(self.root.Eirene.Surfaces) > 0): 
			SurfNames  = [self.root.Eirene.Surfaces[k].Name  for k in range(len(self.root.Eirene.Surfaces))]
			SurfRcoefs = [self.root.Eirene.Surfaces[k].Rcoef for k in range(len(self.root.Eirene.Surfaces))]
			self.EireneSurfsTreeIds = self.tree_view_fill(self.EireneSurfsTree,  self.EireneSurfsTreeIds, (SurfNames, SurfRcoefs, ColorlistSurf[:len(SurfNames)]))
		else:
			self.EireneSurfsTreeIds = self.tree_view_fill(self.EireneSurfsTree,  self.EireneSurfsTreeIds, ([], [], []))
		
	def pump_tree_view_fill(self):
		if(len(self.root.Eirene.Wall.Material) > 0): 
			PumpNames = [self.root.Eirene.Wall.Pump[k].Name for k in range(len(self.root.Eirene.Wall.Pump))]
			self.PumpTreeIds = self.tree_view_fill(self.PumpTree, self.PumpTreeIds, (PumpNames, ColorlistPump[:len(PumpNames)]))
		else:
			self.PumpTreeIds = self.tree_view_fill(self.PumpTree, self.PumpTreeIds, ([], []))

	def puff_tree_view_fill(self):
		if(len(self.root.Eirene.Wall.Puff) > 0): 
			PuffNames = [self.root.Eirene.Wall.Puff[k].Name for k in range(len(self.root.Eirene.Wall.Puff))]
			self.PuffTreeIds = self.tree_view_fill(self.PuffTree, self.PuffTreeIds, (PuffNames, ColorlistPuff[:len(PuffNames)]))
		else:
			self.PuffTreeIds = self.tree_view_fill(self.PuffTree, self.PuffTreeIds, ([], []))


	def change_wall(self, IndMatPump, IsPump):

		if(self.ShowFigureOptMenuValue.get() != self.ShowFigureOptions[0]): return

		self.IndMatPump	= IndMatPump
		self.IsPump		= IsPump
		self.nPtSet		= 0

		self.cid_mouse_button_press_set_wall_points = self.Fig.canvas.mpl_connect('button_press_event', self.cmd_button_press_set_wall_points)
		self.cursor = Cursor(self.Ax, horizOn=True, vertOn=True, useblit=False, color='red', linewidth=1)


	def cmd_button_press_set_wall_points(self, event):
		self.xPtSet[self.nPtSet]	= event.xdata
		self.yPtSet[self.nPtSet]	= event.ydata
		self.nPtSet += 1
		if(self.nPtSet == 3):
			self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press_set_wall_points)
			self.cid_mouse_button_press_set_wall_points = 0
			self.cursor=0

			set_wall_material_pump(self.root.Eirene, self.IsPump, self.IndMatPump, self.xPtSet, self.yPtSet)
			self.update_figure()


	def cmd_button_press_set_puff(self, event):
		if(self.cid_mouse_button_press_set_puff != 0):
			self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press_set_puff)
			self.cid_mouse_button_press_set_puff = 0
			self.cursor=0
		
		set_wall_puff(self.root.Eirene, self.PuffInd, event.xdata, event.ydata)	
		if(self.PuffInd < 0): self.root.Eirene.Wall.Puff[-1].Name = "{:d}".format(len(self.root.Eirene.Wall.Puff)+1)
	
		self.puff_tree_view_fill()
		self.update_figure()

#===============
#	Update data 
#===============

	def update_data(self):
		if(self.root.Config.Eirene_OK):
			self.ToCenterCkButValue.set(self.root.Eirene.ToTheCenter) 
			self.EWalls = self.root.Eirene.Wall.EWalls
			self.surface_tree_view_fill()
			self.material_tree_view_fill()
			self.pump_tree_view_fill()
			self.puff_tree_view_fill()
		else:			
			self.EireneSurfsTreeIds = self.tree_view_fill(self.EireneSurfsTree,  self.EireneSurfsTreeIds, ([], [], []))
			self.MaterialTreeIds	= self.tree_view_fill(self.MaterialTree,  self.MaterialTreeIds, ([], []))		
			self.PumpTreeIds	 	= self.tree_view_fill(self.PumpTree, self.PumpTreeIds, ([], []))
			self.PuffTreeIds	 	= self.tree_view_fill(self.PuffTree, self.PuffTreeIds, ([], []))
			self.EWalls			 	= wall_define_eirene_walls(self.root.Config)
		self.eirene_zones_tree_view_fill()
		self.update_figure()
		
#================
#	Eirene Figure 
#================

	def update_figure(self):
		Config		= self.root.Config

		Fig			= self.Fig
		Ax			= self.Ax
	
		Xlims = Ax.get_xlim()
		Ylims = Ax.get_ylim()
		
		if(self.Cb != 0):
			Fig.clf()
			self.Ax = Ax = Fig.add_subplot(111)
			self.Cb = 0
		else:
			Ax.clear()

		Ax.tick_params(axis='both', which='major', labelsize=8)
		Ax.autoscale(axis='both', tight=True)
		
		if(Config.wall_OK):
			plot2d_walls(Ax, Config.Walls, self.root.ShowPlasmaWallCheckValue.get(), self.root.ShowEireneWallCheckValue.get(), self.root.ShowExtraWallCheckValue.get())
			if(self.ShowFigureOptMenuValue.get() == self.ShowFigureOptions[0]):
				EWalls = self.EWalls
				for iEWall in range(len(EWalls)):
					EWall = EWalls[iEWall]
					for iESubZone in range(len(EWall.ESubZones)):
						ESubZone = EWall.ESubZones[iESubZone]
						if(len(EWall.ESubZones) > 0):
							Ax.plot(ESubZone.Rbord, ESubZone.Zbord, color="gray", linestyle="-", linewidth=1, marker=".", markersize = 6)
						Rc	  = 0.5*(np.min(ESubZone.Rbord)+np.max(ESubZone.Rbord))
						Zc	  = 0.5*(np.min(ESubZone.Zbord)+np.max(ESubZone.Zbord))
						Ax.text(Rc, Zc, "({:d},{:d})".format(EWall.iWall+1,iESubZone+1), color='m', ha="center",fontsize=10)		
					if(Config.Eirene_OK and hasattr(EWall,"RWbord")):
						Ax.plot(EWall.RWbord, EWall.ZWbord, color="blue", linestyle="-", linewidth=1, marker=".", markersize = 6)


		if(Config.Eirene_OK):
			X_points	= self.root.Config.X_points
			Eirene		= self.root.Eirene

			if(self.ShowFigureOptMenuValue.get() == self.ShowFigureOptions[1]):																					#Plot Wall Triangles
				Triangles = np.array([Eirene.WallTriangles.p1, Eirene.WallTriangles.p2, Eirene.WallTriangles.p3]).T
				Ax.triplot(Eirene.RKnots, Eirene.ZKnots, Triangles, 'b-')
				Rc = (Eirene.RKnots[Eirene.WallTriangles.p1] + Eirene.RKnots[Eirene.WallTriangles.p2] + Eirene.RKnots[Eirene.WallTriangles.p3])/3.
				Zc = (Eirene.ZKnots[Eirene.WallTriangles.p1] + Eirene.ZKnots[Eirene.WallTriangles.p2] + Eirene.ZKnots[Eirene.WallTriangles.p3])/3.
				Ax.plot(Rc,  Zc,  'm.')
				for n in range(len(Rc)):
					Ax.text(Rc[n], Zc[n], "{:d}".format(n+1), color='m', fontsize=10)

			if((self.ShowFigureOptMenuValue.get() == self.ShowFigureOptions[2]) or (self.ShowFigureOptMenuValue.get() == self.ShowFigureOptions[3])):			#Plot Triangles

				Triangles = np.array([Eirene.Triangles.p1, Eirene.Triangles.p2, Eirene.Triangles.p3]).T
				Ax.triplot(Eirene.RKnots, Eirene.ZKnots, Triangles, 'b-')

			if(self.ShowFigureOptMenuValue.get() == self.ShowFigureOptions[3]):																					#Plot Triangles & Zone

				Rc  = (Eirene.RKnots[Eirene.Triangles.p1] + Eirene.RKnots[Eirene.Triangles.p2] + Eirene.RKnots[Eirene.Triangles.p3])/3.
				Zc  = (Eirene.ZKnots[Eirene.Triangles.p1] + Eirene.ZKnots[Eirene.Triangles.p2] + Eirene.ZKnots[Eirene.Triangles.p3])/3.
				In  = np.where(Eirene.Triangles.PlasmaVacuum == 0)[0]
				Out = np.where(Eirene.Triangles.PlasmaVacuum == 1)[0]
				Ax.plot(Rc[In],  Zc[In],  'r.')
				Ax.plot(Rc[Out], Zc[Out], 'y.')


			if(self.ShowFigureOptMenuValue.get() != self.ShowFigureOptions[1]):																					#Plot Wall point
				for iWall in range(len(Eirene.Wall.R12)):
					IsPump  = Eirene.Wall.IsPump[iWall]
					TypeMat = Eirene.Wall.TypeMat[iWall]

					iPumps  = np.where(IsPump)[0]
					iMats   = np.where(~IsPump)[0]
					if((len(iPumps) > 0) or (len(iMats) > 0)):
						R12 = Eirene.Wall.R12[iWall]
						Z12 = Eirene.Wall.Z12[iWall]
						WallRc = 0.5*(R12[:,0]+R12[:,1])
						WallZc = 0.5*(Z12[:,0]+Z12[:,1])

					if(len(iPumps) > 0):
						for k in range(len(Eirene.Wall.Pump)):
							iPump = np.where(TypeMat[iPumps] == k)[0]
							Ax.plot(R12[iPumps[iPump],:].T, Z12[iPumps[iPump],:].T,'-', linewidth=2, color=ColorlistPump[k])
							Ax.plot(WallRc[iPumps[iPump]], WallZc[iPumps[iPump]],'.', markersize = 6 , color=ColorlistPump[k])

					if(len(iMats)):
						for k in range(len(Eirene.Wall.Material)):
							iMat = np.where(TypeMat[iMats] == k)[0]
							Ax.plot(R12[iMats[iMat],:].T, Z12[iMats[iMat],:].T,'-', linewidth=2, color=ColorlistMat[k])
							Ax.plot(WallRc[iMats[iMat]], WallZc[iMats[iMat]],'.', markersize=6, color=ColorlistMat[k])
						
				for k in range(len(X_points)):
					Ax.contour(Config.r2D, Config.z2D, Config.flux2D, levels = [X_points[k].psi], colors ='r', linestyles='dashed')

			for n in range(len(Eirene.Wall.Puff)):
				Ax.plot(Eirene.Wall.Puff[n].R,Eirene.Wall.Puff[n].Z,'.',markersize=8, linewidth=3, color=ColorlistPuff[n])
				Ax.text(Eirene.Wall.Puff[n].R,Eirene.Wall.Puff[n].Z, Eirene.Wall.Puff[n].Name, color=ColorlistPuff[n], fontsize=12)


			for k in range(len(Eirene.Surfaces)):
				if((len(Eirene.Surfaces[k].R) == 0) or (self.InSetSurface and (k == self.IndSurf))): continue
				Ax.plot(Eirene.Surfaces[k].R, Eirene.Surfaces[k].Z, '-', linewidth=2, color=ColorlistSurf[k])

		else:
			if(not Config.wall_OK): Ax.clear()

		Ax.set_aspect(1.)
		if(not ((Xlims[0] == 0.) and (Xlims[1] == 1.))):
			Ax.set_xlim(Xlims[0], Xlims[1])
			Ax.set_ylim(Ylims[0], Ylims[1])
	
		Fig.canvas.draw()
