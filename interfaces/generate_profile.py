#!/usr/bin/env python
# -*- coding: utf-8 -*-

import tkinter		as tk
from tkinter		import Tk, Text, filedialog, messagebox, BOTH, W, N, E, S, Menu, HORIZONTAL
from tkinter		import DoubleVar, IntVar, StringVar, OptionMenu
from tkinter.ttk	import Entry, Menubutton, Checkbutton, Frame, Button, Label, Style, LabelFrame, Scale

import types
from tkinter.messagebox import Message
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg		import FigureCanvasTkAgg
from matplotlib.figure						import Figure
from matplotlib.widgets 					import Cursor
try:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2TkAgg as NavigationToolbar2Tk
except:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2Tk

from math									import floor
import os
import numpy								as np
from numpy									import arange, sin, pi

from files.load_text_data					import load_text_data
from routines.set_profile_from_filedata		import set_cut_profile_from_filedata
from routines.globals						import LENGTH_TYPES

class generate_profile(Frame):

	def __init__(self, parent, Name="D", Title="Radial profile", FigSize=(7.5, 2.3)):

		Frame.__init__(self, parent)
		
		nPadx	= 5
		nPady	= 2
		xWidth  = 11
		self.parent = parent

		self.Name = Name
		self.Frame = LabelFrame(parent, text=Title)

		self.CopyFromBut = Button(self.Frame, text="Load from..", command = self.cmd_load_from, width = xWidth)
		self.CopyFromBut.grid(row=0, column=4, pady=nPady)

		self.CopyFromOptMenu = Menubutton(self.Frame, text='Copy from..', width = xWidth)
		self.CopyFromOptMenu.grid(row=0, column=5,  pady=nPady, padx=nPadx)
		self.CopyFromOptMenu.menu = Menu(self.CopyFromOptMenu, tearoff=0)
		self.CopyFromOptMenu['menu'] = self.CopyFromOptMenu.menu

#		Plot area

		self.Fig = Figure(figsize=FigSize, dpi=100)
		self.Fig.patch.set_facecolor('white')

		# a tk.DrawingArea
		canvas = FigureCanvasTkAgg(self.Fig, master=self.Frame)
		try:
			canvas.show()										#Matplotlib  up 2.2
		except:
			canvas.draw()
		canvas.get_tk_widget().grid(row=1, column=0, columnspan=6,  padx=nPadx, sticky=E+W+S+N)

		self.Ax = self.Fig.add_subplot(111)
		self.Ax.tick_params(axis='both', which='major', labelsize=8)
		self.Fig.tight_layout()

		toolbar_frame = Frame(self.Frame, width=50)
		toolbar_frame.grid(row=0, column=0, columnspan=4, padx=nPadx, sticky=W+S+N)
		toolbar = NavigationToolbar2Tk(canvas, toolbar_frame )

		toolbar.update()

		self.AddPointButValue = StringVar()
		self.AddPointButValue.set("Add point")
		self.AddPointBut = Button(self.Frame, textvariable=self.AddPointButValue, command = self.cmd_add_point, width = xWidth)
		self.AddPointBut.grid(row=2, column=0, pady=nPady)

		self.RemovePointButValue = StringVar()
		self.RemovePointButValue.set("Remove point")
		self.RemovePointBut = Button(self.Frame, textvariable=self.RemovePointButValue, command = self.cmd_remove_point, width = xWidth)
		self.RemovePointBut.grid(row=2, column=1, pady=nPady)

		self.MovePointButValue = StringVar()
		self.MovePointButValue.set("Move point")
		self.MovePointBut = Button(self.Frame, textvariable=self.MovePointButValue, command = self.cmd_move_point, width = xWidth)
		self.MovePointBut.grid(row=2, column=2, pady=nPady)

		self.MoveZeroButValue = StringVar()
		self.MoveZeroButValue.set("Move zero")
		self.MoveZeroBut = Button(self.Frame, textvariable=self.MoveZeroButValue, command = self.cmd_move_zero, width = xWidth)
		self.MoveZeroBut.grid(row=2, column=3, pady=nPady)

		TranspValueLabel = Label(self.Frame, text="Value:", width = 5)
		TranspValueLabel.grid(row=2, column=4,  pady=nPady, padx=0, sticky=E)

		self.CoordEntryValue = StringVar()
		self.CoordEntry = Entry(self.Frame,  textvariable = self.CoordEntryValue, width = 12)
		self.CoordEntry.grid(row=2, column=5,  pady=nPady, padx=nPadx, sticky=E)
		self.CoordEntry.bind("<Return>", self.cmd_enter_entry_value)
		self.CoordEntry.state(["disabled"])

		self.InAddPoint		= False
		self.InRemovePoint	= False
		self.InMovePoint	= False
		self.InMoveZero		= False

#	Load profile data from file
		
	def cmd_load_from(self):
		LoadFile = filedialog.askopenfilename(filetypes=(("Comma separated values", "*.csv"), \
														 ("Tab separated values",    "*.tsv"), \
														 ("Space separated values",  "*.*")))
		if(not LoadFile): return

		Headers, FileData = load_text_data(LoadFile)

		if(len(Headers) == 0): return
	
		if(Headers[0] == "1"):
			Headers[0] = LENGTH_TYPES[0]								#If none assume DIst
			Headers[1] = self.Name

#		set_cut_profile_from_filedata(Headers, FileData, self, "Length", self.Lengths, xOffset=self.d12[0,0])
		xName = LENGTH_TYPES[0]
		set_cut_profile_from_filedata(Headers, FileData, self, xName, self.Lengths, OutxOffset=self.d12[0,0])
		
		self.resize_mode = True
		self.set_iKnots()
		self.update_figure()
		return
			
#	Copy profile from other one
			
	def copy_from_set(self, Profiles, TabNames):
		self.profiles = []
		for k in range(len(Profiles)):
			self.profiles.append(Profiles[k])
			if(k==0): 	self.CopyFromOptMenu.menu.add_command(label=TabNames[k], command=lambda: self.cmd_copy_from_set(0))		
			elif(k==1): self.CopyFromOptMenu.menu.add_command(label=TabNames[k], command=lambda: self.cmd_copy_from_set(1))		
			elif(k==2): self.CopyFromOptMenu.menu.add_command(label=TabNames[k], command=lambda: self.cmd_copy_from_set(2))		
			elif(k==3): self.CopyFromOptMenu.menu.add_command(label=TabNames[k], command=lambda: self.cmd_copy_from_set(3))		
			elif(k==4): self.CopyFromOptMenu.menu.add_command(label=TabNames[k], command=lambda: self.cmd_copy_from_set(4))		
			elif(k==5): self.CopyFromOptMenu.menu.add_command(label=TabNames[k], command=lambda: self.cmd_copy_from_set(5))		

#			self.CopyFromOptMenu.menu.add_command(label=profiles[k].Name, command=lambda: self.cmd_copy_from_set(profiles[k]))		

	def set_data(self, ZoneDatas, Profile):
		self.nz			= np.copy(ZoneDatas.nz)
		self.d12		= np.copy(ZoneDatas.d12)
		self.d			= np.copy(ZoneDatas.d12.reshape(ZoneDatas.d12.size))

		if(hasattr(ZoneDatas, "Lengths")): self.Lengths = np.copy(ZoneDatas.Lengths)

		ii		= np.where(ZoneDatas.d12[1:,0]-ZoneDatas.d12[:-1,1] > 0.); ii = ii[0]
		if(len(ii) > 0):
			self.dKnots	= np.concatenate((np.array([ZoneDatas.d12[0,0]]), ZoneDatas.d12[ii,1], ZoneDatas.d12[ii+1,0], np.array([ZoneDatas.d12[-1,1]]))) 
		else:
			self.dKnots	= np.array([ZoneDatas.d12[0,0], ZoneDatas.d12[-1,1]])

		self.xValues        = np.copy(Profile.xValues)
		self.Values   		= np.copy(Profile.Values)
		self.resize_mode	= True
		self.set_iKnots()
		self.update_figure()

	def set_iKnots(self):
		self.iKnots	= np.empty(self.dKnots.shape, dtype='i4')
		for k in range(len(self.iKnots)):
			ii = np.where(self.xValues == self.dKnots[k])
			self.iKnots[k] = ii[0][0]			

	def reset_data(self):
		self.nz			= 0
		self.d12		= 0
		self.d			= 0
		self.dKnots		= 0
		self.iKnots		= 0
		self.xValues	= 0
		self.Values		= 0
		self.Ax.clear()
		self.Fig.canvas.draw()

	def get_data(self):
		Profile = types.SimpleNamespace()
		Profile.Name	= self.Name
		Profile.xValues = np.copy(self.xValues)
		Profile.Values	= np.copy(self.Values)
		
		return Profile, np.copy(self.d12)

	def set_zero(self, dMove):
		self.d   			-= dMove
		self.dKnots			-= dMove
		self.d12			-= dMove
		self.xValues		-= dMove 
		self.resize_mode	 = True
		self.update_figure()
		
	def cmd_copy_from_set(self, iProfile):
		Profile, tmp_d12 = self.profiles[iProfile].get_data()
		self.xValues = Profile.xValues
		self.Values  = Profile.Values
		self.set_iKnots()
		self.resize_mode = True
		self.update_figure()
		
	def cmd_add_point(self):
		if(not self.InAddPoint):
			self.AddPointButValue.set("Cancel")
			self.CoordEntry.state(["!disabled"])

			self.cid_mouse_button_press = self.Fig.canvas.mpl_connect('button_press_event', self.cmd_button_press)
			self.cursor = Cursor(self.Ax, useblit=False, color='red', linewidth=1)
			self.InAddPoint = True
		else:
			self.InAddPoint = False
			self.AddPointButValue.set("Add point")
			self.CoordEntry.state(["disabled"])
			if(self.cid_mouse_button_press != 0):
				self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press)
				self.cid_mouse_button_press = 0
				self.cursor=0
			self.update_figure()

	def cmd_remove_point(self):
		if(not self.InRemovePoint):
			self.RemovePointButValue.set("Cancel")
			self.CoordEntry.state(["!disabled"])

			self.cid_mouse_button_press = self.Fig.canvas.mpl_connect('button_press_event', self.cmd_button_press)
			self.cursor = Cursor(self.Ax, useblit=False, color='red', linewidth=1)
			self.InRemovePoint = True
		else:
			self.InRemovePoint = False
			self.RemovePointButValue.set("Remove point")
			self.CoordEntry.state(["disabled"])
			if(self.cid_mouse_button_press != 0):
				self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press)
				self.cid_mouse_button_press = 0
				self.cursor=0
			self.update_figure()
		
	def cmd_move_point(self):
		if(not self.InMovePoint):
			self.MovePointButValue.set("Cancel")
			self.CoordEntry.state(["!disabled"])

			self.cid_mouse_button_press = self.Fig.canvas.mpl_connect('button_press_event', self.cmd_button_press)
			self.cursor = Cursor(self.Ax, useblit=False, color='red', linewidth=1)
			self.InMovePoint = True
			self.iMove		 = -1
			self.iKnot		 = -1
		else:
			self.InMovePoint = False
			self.MovePointButValue.set("Move point")
			self.CoordEntry.state(["disabled"])
			if(self.cid_mouse_button_press != 0):
				self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press)
				self.cid_mouse_button_press = 0
				self.cursor=0
			self.update_figure()

		
	def cmd_move_zero(self):
		if(not self.InMoveZero):
			self.MoveZeroButValue.set("Cancel")
			self.CoordEntry.state(["!disabled"])

			self.cid_mouse_button_press = self.Fig.canvas.mpl_connect('button_press_event', self.cmd_button_press)
			self.cursor = Cursor(self.Ax, useblit=False, color='red', linewidth=1)
			self.InMoveZero = True
			self.iMoveZero	= -1
		else:
			self.InMoveZero = False
			self.MoveZeroButValue.set("Move zero")
			self.CoordEntry.state(["disabled"])
			if(self.cid_mouse_button_press != 0):
				self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press)
				self.cid_mouse_button_press = 0
				self.cursor=0
			self.resize_mode = True
			self.update_figure()

#	Manage entry coordinate

	def cmd_enter_entry_value(self, event):
		try:
			Values = eval(self.CoordEntryValue.get())
			if(not isinstance(Values, tuple)): Values = (Values,)
			if((len(Values) != 1) and (len(Values) != 2)): return
		except:
			return
		
		value =	types.SimpleNamespace()
		if(self.InMovePoint and (self.iKnot > -1)):
			value.xdata = self.dKnots[self.iKnot]		
			if(len(Values) == 1):	
				value.ydata = Values[0]
			else:
				value.ydata = Values[1]			
			
		elif(self.InMoveZero):
			if(len(Values) != 1): return
			
			value.xdata = self.xValues[0]
			value.ydata = 0.
		else:
			if(len(Values) != 2):return
			value.xdata = Values[0]
			value.ydata = Values[1]
			
		self.cmd_button_press(value)

#	Manage Add point, move point, remove point and move zero

	def cmd_button_press(self, event):
		xPt = event.xdata
		yPt = event.ydata
		
#		Add point

		if(self.InAddPoint):
			if((xPt <= self.xValues[0]) or (xPt >= self.xValues[-1])): return
			iKnot = np.where((self.xValues[self.iKnots[:-1]] <= xPt) & (self.xValues[self.iKnots[1:]] > xPt)); iKnot = iKnot[0]
			if(floor(iKnot/2)*2 != iKnot):	return					#Points cannot be added between odd-even nodes
						
			x	= np.append(self.xValues, xPt) 
			y	= np.append(self.Values, yPt)
			a = np.argsort(x)
			self.xValues = x[a]
			self.Values  = y[a]
			self.set_iKnots()
			self.cmd_add_point()									#Reset plot and buttons
		
#		Move point

		if(self.InMovePoint):
			if((xPt < self.xValues[0]) or (xPt > self.xValues[-1])): return
			
			if(self.iMove == -1):									#find nearest point
				self.iMove = np.argmin(np.abs(self.xValues - xPt))
				ii = np.where(self.iKnots == self.iMove)			#Check if is node
				if(len(ii[0])> 0):
					self.iKnot = ii[0]					
					self.Ax.plot(self.xValues[self.iMove], self.Values[self.iMove],'ro')
				else:
					self.iKnot = -1					
					self.Ax.plot(self.xValues[self.iMove], self.Values[self.iMove],'go')
			else:
				if((xPt < self.xValues[max(self.iMove-1,0)]) or (xPt > self.xValues[min(self.iMove+1, len(self.xValues)-1)])): return
				if(self.iKnot < 0): self.xValues[self.iMove] = xPt	#x cannot be changed for nodes
				self.Values[self.iMove]  = yPt
				self.set_iKnots()
				self.cmd_move_point()								#Reset plot and buttons
		
#		Remove point

		if(self.InRemovePoint):
			if((xPt <= self.xValues[0]) or (xPt >= self.xValues[-1])): return
			
			iRemove = np.argmin(np.abs(self.xValues - xPt))
			ii = np.where(self.iKnots == iRemove)					#nodes cannot be removed
			if(len(ii[0]) == 0):
				self.xValues = np.delete(self.xValues, iRemove)
				self.Values  = np.delete(self.Values,  iRemove)
				self.set_iKnots()
				self.cmd_remove_point()								#Reset plot and buttons
		
#		Move zero

		if(self.InMoveZero):
			if((xPt < self.xValues[0]) or (xPt > self.xValues[-1])): return
			
			iMoveZero = np.argmin(np.abs(self.d - xPt))				#find nearest zone limit
			dMove 				 = self.d[iMoveZero]
			self.d				-= dMove
			self.dKnots			-= dMove
			self.d12			-= dMove
			self.xValues		-= dMove 
			self.resize_mode	 = True
			self.cmd_move_zero()									#Reset plot and buttons
			for k in range(len(self.profiles)):
				self.profiles[k].set_zero(dMove)	

	def update_figure(self):
		Fig			= self.Fig
		Ax			= self.Ax
		
		Vmax = np.max(self.Values)
		Vmin = np.min(self.Values)

		Xlims = Ax.get_xlim()
		Ylims = Ax.get_ylim()
		Ax.clear()
		
		Ax.tick_params(axis='both', which='major', labelsize=8)
		Ax.autoscale(axis='both', tight=True)
		Ax.plot(self.xValues,self.Values,'g.-')
		Ax.plot(self.xValues[self.iKnots],self.Values[self.iKnots],'r.')
		
		for k in range(len(self.nz)):
			d = 0.5*(self.d12[k,0]+self.d12[k,1])
			Ax.text(d,0.5*Vmax,"{:d}".format(self.nz[k]+1), fontsize=8, horizontalalignment='center')
			
		for k in range(len(self.nz)):
			Ax.axvline(x=self.d12[k,0], color='b', linestyle='dashed')

		if(self.resize_mode):
			Ax.set_ylim(0., Vmax*1.1)
			self.resize_mode = False
		else:
			Ax.set_xlim(Xlims[0], Xlims[1])
			Ax.set_ylim(Ylims[0], Ylims[1])
			
		Xlims = Ax.get_xlim()
		Ylims = Ax.get_ylim()
#		for k in range(len(self.nz)):
#			Ax.plot([self.d12[k,0], self.d12[k,0]],Ylims,'b--')
#			Ax.plot([self.d12[k,1], self.d12[k,1]],Ylims,'b--')
	
		Fig.canvas.draw()

