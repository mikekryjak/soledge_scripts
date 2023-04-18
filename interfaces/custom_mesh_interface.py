#!/usr/bin/env python
# -*- coding: utf-8 -*-

import tkinter				as tk
from tkinter				import Tk, Toplevel, BOTH, W, N, E, S, Menu, HORIZONTAL
from tkinter				import DoubleVar, IntVar, StringVar, OptionMenu
from tkinter.ttk			import Checkbutton, Frame, Button, Label, LabelFrame
from tkinter				import messagebox

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg		import FigureCanvasTkAgg
from matplotlib.figure						import Figure
from matplotlib.widgets 					import Cursor
try:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2TkAgg as NavigationToolbar2Tk
except:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2Tk

import types
from math									import floor
import numpy								as np
from numpy									import arange, sin, pi
from scipy									import interpolate

from routines.bezier						import bezier, bezier3
from routines.utils_walls					import plot2d_walls

class CustomMeshInterface(Toplevel):

	def __init__(self, parent, Config, X_points, pt_num, xBezier, yBezier, R, Z, CurveType = 0, d_force = np.array([],dtype='f8')):

		Toplevel.__init__(self, parent)

		nPadx	= 5
		nPady	= 5
		xWidth  = 15

		self.parent		= parent
		self.result		= types.SimpleNamespace()
		self.result.d	= np.array([])

		self.title("Custom mesh")
		self.geometry("1200x700+{:d}+{:d}".format(parent.winfo_rootx()+50,parent.winfo_rooty()+50))

		self.protocol("WM_DELETE_WINDOW", self.cancel)


		win_frame = Frame(self)
		win_frame.pack(fill=BOTH, expand=True)

		win_frame.columnconfigure(1, weight=1)

#		Plot area 1

		plot_frame = Frame(win_frame)
#		plot_frame.grid(row=0, column=0, padx=nPadx, sticky=E+W+S+N)
		plot_frame.pack(fill=BOTH, expand=True)

		self.Fig1 = Figure(figsize=(4.8, 5), dpi=100)
		self.Fig1.patch.set_facecolor('white')

		# a tk.DrawingArea
		canvas = FigureCanvasTkAgg(self.Fig1, master=plot_frame)
		try:
			canvas.show()										#Matplotlib  up 2.2
		except:
			canvas.draw()

		canvas.get_tk_widget().grid(row=1, column=0,  padx=nPadx, sticky=E+W+S+N)

		toolbar_frame = Frame(plot_frame, width=40)
		toolbar_frame.grid(row=0, column=0, columnspan=3, padx=nPadx, sticky=E+W+S+N)
		toolbar = NavigationToolbar2Tk(canvas, toolbar_frame )

		toolbar.update()
		self.Ax1 = self.Fig1.add_subplot(111)
		self.Ax1.tick_params(axis='both', which='major', labelsize=8)

#		Plot area 2
		
		self.Fig2 = Figure(figsize=(4.8, 5), dpi=100)
		self.Fig2.patch.set_facecolor('white')

		# a tk.DrawingArea
		canvas = FigureCanvasTkAgg(self.Fig2, master=plot_frame)
		try:
			canvas.show()										#Matplotlib  up 2.2
		except:
			canvas.draw()
		canvas.get_tk_widget().grid(row=1, column=1,  padx=nPadx, sticky=E+W+S+N)

		toolbar_frame = Frame(plot_frame, width=40)
		toolbar_frame.grid(row=0, column=1, padx=nPadx, sticky=E+W+S+N)
		toolbar = NavigationToolbar2Tk(canvas, toolbar_frame )

		toolbar.update()
		self.Ax2 = self.Fig2.add_subplot(111)
		self.Ax2.tick_params(axis='both', which='major', labelsize=8)


		button_frame = Frame(win_frame)
#		button_frame.grid(row=1, column=0, padx=nPadx, sticky=E+W+S+N)
		button_frame.pack(fill=BOTH, expand=True)
	
		self.AddPointButValue = StringVar()
		self.AddPointButValue.set("Add control point")
		self.AddPointBut = Button(button_frame, textvariable=self.AddPointButValue, command = self.cmd_add_point, width = xWidth)
		self.AddPointBut.grid(row=0, column=0, padx=nPadx, pady=nPady, sticky=W)

		self.MovePointButValue = StringVar()
		self.MovePointButValue.set("Move control point")
		self.MovePointBut = Button(button_frame,textvariable=self.MovePointButValue, command = self.cmd_move_point, width = xWidth)
		self.MovePointBut.grid(row=0, column=1, padx=nPadx, pady=nPady, sticky=W)

		self.RemovePointButValue = StringVar()
		self.RemovePointButValue.set("Remove control point")
		self.RemovePointBut = Button(button_frame,textvariable=self.RemovePointButValue, command = self.cmd_remove_point, width = xWidth)
		self.RemovePointBut.grid(row=0, column=2, padx=nPadx, pady=nPady, sticky=W)

		self.ApplyButValue = Button(button_frame,text="Apply", command = self.cmd_apply, width = xWidth)
		self.ApplyButValue.grid(row=0, column=3, padx=nPadx, pady=nPady, sticky=W)
		self.ApplyButValue.state(["!disabled"])

		self.PreviewCkButValue = IntVar()
		PreviewCkBut = Checkbutton(button_frame, text="Preview", variable = self.PreviewCkButValue, command = self.cmd_preview, width = xWidth)
		PreviewCkBut.grid(row=0,column=4, padx=nPadx,  pady=nPady, sticky=E)

		self.Config	   = Config
		self.X_points  = X_points
		self.npt 	   = pt_num
		self.R		   = np.copy(R)
		self.Z		   = np.copy(Z)
		self.d_force   = d_force
		self.CurveType = CurveType
		if(self.CurveType == 0):											#Bezier curve
			if((len(xBezier.shape) == 1) and (len(xBezier) > 1) and (len(xBezier) == len(yBezier))):
				self.x_ctr	  = np.copy(xBezier)
				self.y_ctr	  = np.copy(yBezier)
			else:															# Default values
				self.x_ctr	  = np.array([1., self.npt], dtype='f8')
				self.y_ctr	  = np.array([0., 1.] , dtype= 'f8')
		else:																#Cubic spline
			if((len(xBezier.shape) == 2) and (xBezier.shape[0] > 1) and (xBezier.shape[0] == yBezier.shape[0])):
				self.x_ctr	  = np.copy(xBezier)
				self.y_ctr	  = np.copy(yBezier)
			else:															# Default values
				self.x_ctr	  = np.array([[1e10, 1., (self.npt-1)/4.],	[self.npt-(self.npt-1)/4., self.npt, 1.e10]],	dtype='f8')
				self.y_ctr	  = np.array([[1e10, 0., 0.25],				[0.75, 1., 1.e10]],							dtype= 'f8')

		self.dist	  = np.zeros(len(self.R), dtype='f8')
		self.dist[1:] = np.cumsum(np.sqrt((self.R[1:]-self.R[:-1])**2+(self.Z[1:]-self.Z[:-1])**2))
		self.dist	  = self.dist/self.dist[-1]

		self.f_psi	 = interpolate.RectBivariateSpline(self.Config.r2D[0,:], self.Config.z2D[:,0], self.Config.flux2D.T, s=0)

		self.InAddPoint		= False
		self.InMovePoint	= False
		self.InRemovePoint	= False
		self.iMove			= -1

		self.compute_meshpoint()
		self.update_figure1()
		self.update_figure2()

#		self.initial_focus.focus_set()

		self.grab_set()
		self.wait_window(self)
	
	def cmd_add_point(self):
		if(not self.InAddPoint):
			self.AddPointButValue.set("Cancel")

			self.cid_mouse_button_press = self.Fig1.canvas.mpl_connect('button_press_event', self.cmd_button_press)
			self.cursor = Cursor(self.Ax1, useblit=False, color='red', linewidth=1)
			self.InAddPoint = True
		else:
			self.InAddPoint = False
			self.AddPointButValue.set("Add control point")
			if(self.cid_mouse_button_press != 0):
				self.Fig1.canvas.mpl_disconnect(self.cid_mouse_button_press)
				self.cid_mouse_button_press = 0
				self.cursor=0
			self.update_figure1()
			self.update_figure2()

	def cmd_move_point(self):
		if(not self.InMovePoint):
			self.MovePointButValue.set("Cancel")

			self.cid_mouse_button_press = self.Fig1.canvas.mpl_connect('button_press_event', self.cmd_button_press)
			self.cursor = Cursor(self.Ax1, useblit=False, color='red', linewidth=1)
			self.InMovePoint = True
			self.iMove		 = -1
		else:
			self.InMovePoint = False
			self.MovePointButValue.set("Move control point")
			if(self.cid_mouse_button_press != 0):
				self.Fig1.canvas.mpl_disconnect(self.cid_mouse_button_press)
				self.cid_mouse_button_press = 0
				self.cursor=0
			self.update_figure1()
			self.update_figure2()


	def cmd_remove_point(self):
		if(not self.InRemovePoint):
			self.RemovePointButValue.set("Cancel")

			self.cid_mouse_button_press = self.Fig1.canvas.mpl_connect('button_press_event', self.cmd_button_press)
			self.cursor = Cursor(self.Ax1, useblit=False, color='red', linewidth=1)
			self.InRemovePoint = True
		else:
			self.InRemovePoint = False
			self.RemovePointButValue.set("Remove control point")
			if(self.cid_mouse_button_press != 0):
				self.Fig1.canvas.mpl_disconnect(self.cid_mouse_button_press)
				self.cid_mouse_button_press = 0
				self.cursor=0
			self.update_figure1()
			self.update_figure2()

	
	def cmd_button_press(self, event):
		x = event.xdata
		y = event.ydata

#		Add point
#		=====

		if(self.InAddPoint):
			if(self.CurveType == 0):											#Bezier curve
				self.x_ctr = np.append(self.x_ctr, x)
				self.y_ctr = np.append(self.y_ctr, y)
				b		   = np.argsort(self.x_ctr)
				self.x_ctr = self.x_ctr[b]
				self.y_ctr = self.y_ctr[b]
			else:																# Cubic spline
				inewpt = np.max(np.where(self.x_ctr[:,1] < x)[0])+1				# Find position of x with respect to existing control triplets

				x_ctr_new = np.zeros((self.x_ctr.shape[0]+1,3), dtype='f8')	# Generate extended and sorted array of control points triplets
				y_ctr_new = np.zeros((self.x_ctr.shape[0]+1,3), dtype='f8')
				x_ctr_new[:inewpt,:] = self.x_ctr[:inewpt,:]
				x_ctr_new[inewpt,1] = x
				x_ctr_new[inewpt+1:,:] = self.x_ctr[inewpt:,:]
				y_ctr_new[:inewpt,:] = self.y_ctr[:inewpt,:]
				y_ctr_new[inewpt,1] = y
				y_ctr_new[inewpt+1:,:] = self.y_ctr[inewpt:,:]

#				We set the slope of the derivative to be  equal to the slope between the 2 surrounding triplets
				xprev	= x_ctr_new[inewpt-1,1]
				yprev	= y_ctr_new[inewpt-1,1]
				xnext	= x_ctr_new[inewpt+1,1]
				ynext	= y_ctr_new[inewpt+1,1]
				slope	= (ynext-yprev)/(xnext-xprev)
				dx		= (xnext-xprev)/6
				dy		= dx*slope
				x_ctr_new[inewpt,0] = x - dx
				y_ctr_new[inewpt,0] = y - dy
				x_ctr_new[inewpt,2] = x + dx
				y_ctr_new[inewpt,2] = y + dy

				self.x_ctr = x_ctr_new											# Storing final result in x_ctr
				self.y_ctr = y_ctr_new

			self.compute_meshpoint()
			self.cmd_add_point()

#		Move point
#		======

		if(self.InMovePoint):
			if(self.iMove == -1):
				if(self.CurveType == 0):											#Bezier curve
					d = np.sqrt((self.x_ctr-x)**2+(self.y_ctr-y)**2)				#find nearest point
					a = np.argmin(d)
					if((a != 0) and (a != self.npt-1)):
						self.Ax1.plot(self.x_ctr[a], self.y_ctr[a], 'ms', linewidth=3)
						self.iMove = a
				else:																#Cubic spline curve
#					Compute distance of click position to all control points in triplets (not only central)
					d = np.sqrt((self.x_ctr-x)**2+(self.y_ctr-y)**2)

					ictrmin,itypemin = np.unravel_index(np.argmin(d),d.shape)		# Find the closest one
					if((ictrmin == 0) and (itypemin != 2)):							#  Trying to move left or central point in first triplet => not allowed
						messagebox.showerror("Custom segment", "Selected point cannot be moved.")
					elif((ictrmin == self.x_ctr.shape[0]-1) and (itypemin !=0)):	# Trying to move right or central point in last triplet => not allowed
						messagebox.showerror("Custom segment", "Selected point cannot be moved.")
					else:															# Highlight selected point
						self.Ax1.plot(self.x_ctr[ictrmin, itypemin], self.y_ctr[ictrmin, itypemin], 'ms', linewidth=3)
						self.iMove = (ictrmin,itypemin)

			else:
				if(self.CurveType == 0):											#Bezier curve
					self.x_ctr[self.iMove] = x
					self.y_ctr[self.iMove] = y
					b			  = np.argsort(self.x_ctr)							# Sort the list of control points
					self.x_ctr 	  = self.x_ctr[b]
					self.y_ctr	  = self.y_ctr[b]
				else:																#Cubic spline curve
#					Update point value but keep a trace of the slope
					ictrmin = self.iMove[0]
					itypemin = self.iMove[1]
					dx = self.x_ctr[ictrmin, 2] - self.x_ctr[ictrmin, 1]
					dy = self.y_ctr[ictrmin, 2] - self.y_ctr[ictrmin, 1]
					self.x_ctr[ictrmin, itypemin] = x;
					self.y_ctr[ictrmin, itypemin] = y;
					if(itypemin != 1):						# If the point is a derivative type, we need to update the one on the other side to preserve C1 continuity
						dx = x - self.x_ctr[ictrmin, 1]
						dy = y - self.y_ctr[ictrmin, 1]
						itype2 = (itypemin + 2)%4
						self.x_ctr[ictrmin,itype2] = self.x_ctr[ictrmin,1] - dx
						self.y_ctr[ictrmin,itype2] = self.y_ctr[ictrmin,1] - dy
					else:									# If the point is a central type, we need to move with it the derivative ones
						self.x_ctr[ictrmin, 0] = x - dx
						self.y_ctr[ictrmin, 0] = y - dy
						self.x_ctr[ictrmin, 2] = x + dx
						self.y_ctr[ictrmin, 2] = y + dy

						b = np.argsort(self.x_ctr[:, 1])	# We also resort the list of control points
						self.x_ctr = self.x_ctr[b,:]
						self.y_ctr = self.y_ctr[b,:]

				self.compute_meshpoint()
				self.cmd_move_point()
		
#		Remove point
#		========

		if(self.InRemovePoint):
			if(self.CurveType == 0):											#Bezier curve
				d = np.sqrt((self.x_ctr-x)**2+(self.y_ctr-y)**2)				#find nearest point
				a = np.argmin(d)
				if((a != 0) and (a != self.npt-1)):
					self.Ax1.plot(self.x_ctr[a], self.y_ctr[a], 'ms')
					answer=messagebox.askquestion("Custom segment", "Highlighted control point will be suppressed\nAre you sure you want to proceed?")
					if(answer=="yes"):
						self.x_ctr = np.append(self.x_ctr[:a], self.x_ctr[a+1:])
						self.y_ctr = np.append(self.y_ctr[:a], self.y_ctr[a+1:])
						self.compute_meshpoint()
						self.cmd_remove_point()
			else:																#Cubic spline curve
				d = np.sqrt((self.x_ctr[:,1]-x)**2+(self.y_ctr[:,1]-y)**2)		#find nearest central point
				a = np.argmin(d)
				if((a != 0) and (a != self.npt-1)):
					self.Ax1.plot(self.x_ctr[a,:], self.y_ctr[a,:], 'ms-', linewidth=3)
					answer=messagebox.askquestion("Custom segment", "Highlighted control point will be suppressed\nAre you sure you want to proceed?")
					if(answer=="yes"):
						self.x_ctr = np.append(self.x_ctr[:a,:], self.x_ctr[a+1:,:], axis=1)
						self.y_ctr = np.append(self.y_ctr[:a,:], self.y_ctr[a+1:,:], axis=1)
						self.compute_meshpoint()
						self.cmd_remove_point()

	def cmd_preview(self):
		self.update_figure2()

	def cancel(self, event=None):
		self.parent.focus_set()
		self.destroy()

	def cmd_apply(self):
		self.result.d		= np.copy(self.d)
		self.result.xBezier = np.copy(self.x_ctr)
		self.result.yBezier = np.copy(self.y_ctr)
		self.destroy()

#	Internal routines
#	##########

	def compute_meshpoint(self):
		if(self.CurveType == 0):	self.x_bex, self.y_bex = bezier(self.x_ctr, self.y_ctr)		#Bezier curve
		else:						self.x_bex, self.y_bex = bezier3(self.x_ctr, self.y_ctr)

		d = np.zeros(self.npt, dtype='f8')														#Interpolate position of mesh points on that curve
		d[0]	= 0
		d[-1]	= 1
		d[1:-1]	= np.interp(np.arange(2,self.npt), self.x_bex, self.y_bex)

#		Find d-MI point the closest to the d_force targets and move them

		id_forced = np.empty(0,dtype='i4')
		for iforce in range(len(self.d_force)):
			iclosest = np.argmin(np.abs(d-self.d_force[iforce]))		# Find closest mesh point
			iclosest = min(max(iclosest,1),len(d)-2)					#Correct if we found the first one or the last one (forbidden)
			d[iclosest] = self.d_force[iforce]							# Move
			id_forced = np.append(id_forced, iclosest)					# Mark as forced

#		Interpolate R,Z position of mesh points along segment to mesh

		self.d			= np.copy(d)
		self.id_forced	= np.copy(id_forced)
		self.Rp			= np.zeros(len(self.d), dtype='f8')
		self.Zp			= np.zeros(len(self.d), dtype='f8')
		self.Rp[0]		= self.R[0]
		self.Rp[-1]		= self.R[-1]
		self.Zp[0]		= self.Z[0]
		self.Zp[-1]		= self.Z[-1]
		self.Rp[1:-1]	= np.interp(self.d[1:-1], self.dist, self.R)
		self.Zp[1:-1]	= np.interp(self.d[1:-1], self.dist, self.Z)

#		Disable the "Apply" button if the d_MI function is not increasing

		if(np.min(np.diff(d)) < 0.):
			messagebox.showerror("Custom segment", "Current point distribution is invalid:\n\tpoints must be ordered in a monotonic way along both x and y axes.")
			self.ApplyButValue.state(["disabled"])
		else:
			self.ApplyButValue.state(["!disabled"])

	def update_figure1(self):
		Ax1			= self.Ax1

		Ax1.clear()
		Ax1.plot(self.x_bex, self.y_bex,'k-')
		Ax1.plot(np.arange(1,self.npt+1), self.d, 'ro')
		if(self.CurveType == 0):								#Bezier curve
			Ax1.plot(self.x_ctr, self.y_ctr,'bs-')
		else:
			Ax1.plot(self.x_ctr[:,1],self.y_ctr[:,1],'sb')
			for ictrl in range(self.x_ctr.shape[0]):
				Ax1.plot(self.x_ctr[ictrl,:],self.y_ctr[ictrl,:],'-+b')
		if(len(self.id_forced,) > 0): 
			Ax1.plot(self.id_forced, self.d[self.id_forced],'bo', label="Forced")
			Ax1.legend(loc="lower right")

		Ax1.set_xlim(1., self.npt)
		Ax1.set_ylim(0., 1.)
		self.Fig1.canvas.draw()

	def update_figure2(self):
		Config		= self.Config
		X_points	= self.X_points
		Ax2			= self.Ax2
		
		Xlims = Ax2.get_xlim()
		Ylims = Ax2.get_ylim()
		Ax2.clear()
		plot2d_walls(Ax2, Config.Walls)
		for k in range(len(X_points)):
			Ax2.contour(Config.r2D, Config.z2D, Config.flux2D, levels=[X_points[k].psi],  colors ='r', linestyles='solid')

		Ax2.plot(self.R,self.Z,'g-')
		Ax2.text(self.R[0],self.Z[0],'0',fontsize=16)
		Ax2.plot(self.Rp,self.Zp,'k.')
		Ax2.text(self.R[-1],self.Z[-1],'1',fontsize=16)
		if(self.PreviewCkButValue.get()):
			psi =self.f_psi.ev(self.Rp, self.Zp)
			Ax2.contour(Config.r2D, Config.z2D, Config.flux2D, levels = psi,  colors ='b', linestyles='solid')

		Ax2.set_aspect(1.)
		if(not ((Xlims[0] == 0.) and (Xlims[1] == 1.))):
			Ax2.set_xlim(Xlims[0], Xlims[1])
			Ax2.set_ylim(Ylims[0], Ylims[1])
		self.Fig2.canvas.draw()


def custom_mesh_interface(root, Config, X_points, pt_num, xBezier, yBezier, R, Z, CurveType = 0, d_force = np.array([],dtype='f8')):
	d = CustomMeshInterface(root, Config, X_points, pt_num, xBezier, yBezier, R, Z, CurveType = CurveType, d_force=d_force)

	return d.result
