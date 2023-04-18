from tkinter import *
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg		import FigureCanvasTkAgg
try:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2TkAgg as NavigationToolbar2Tk
except:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2Tk

from matplotlib.figure						import Figure
from matplotlib.widgets 					import Cursor
from routines.utils_walls					import plot2d_walls

class AskForXpoint(Toplevel):

	def __init__(self, parent, MeshGenerator, level, point, title = None):

		Toplevel.__init__(self, parent)
		self.transient(parent)

		if title:
			self.title(title)

		self.parent = parent

		self.result = None

		body = Frame(self)
		self.initial_focus = self.body(body, MeshGenerator, level, point)
		body.pack(padx=5, pady=5)

		self.buttonbox()

		self.grab_set()

		if not self.initial_focus:
			self.initial_focus = self

		self.protocol("WM_DELETE_WINDOW", self.cancel)

		self.geometry("+%d+%d" % (parent.winfo_rootx()+50,
								  parent.winfo_rooty()+50))

		self.initial_focus.focus_set()

		self.grab_set()
		self.wait_window(self)


	def body(self, master, MeshGenerator, level, point):

		Fig = Figure(figsize=(4.8, 5.8), dpi=100)
		self.Fig.patch.set_facecolor('white')

		# a tk.DrawingArea

		canvas = FigureCanvasTkAgg(Fig, master=master)
		try:
			canvas.show()										#Matplotlib  up 2.2
		except:
			canvas.draw()

		Ax = Fig.add_subplot(111)
		self.Fig.tight_layout()

		toolbar_frame = Frame(master)
		toolbar = NavigationToolbar2Tk(canvas, toolbar_frame )

		toolbar.update()

		toolbar_frame.pack(fill=X)
		canvas.get_tk_widget().pack(fill=X)

		Ax.contour(MeshGenerator.r2D, MeshGenerator.z2D, MeshGenerator.flux2D, levels =[level])
		Ax.plot(point.x,point.y,'ks')
		plot2d_walls(Ax, MeshGenerator.Config.Walls)
		Ax.set_aspect(1.)
		Fig.canvas.draw()

	def buttonbox(self):
		# add standard button box. override if you don't want the
		# standard buttons

		box = Frame(self)

		w = Button(box, text="X Point", width=10, command=self.X_point, default=ACTIVE)
		w.pack(side=LEFT, padx=5, pady=5)
		w = Button(box, text="O Point", width=10, command=self.O_point)
		w.pack(side=LEFT, padx=5, pady=5)

		box.pack()

	def cancel(self, event=None):
		# put focus back to the parent window
		self.parent.focus_set()
		self.destroy()

	def X_point(self, event=None):
		self.result = "X"
		self.destroy()

	def O_point(self, event=None):
		self.result = "O"
		self.destroy()


def ask_for_Xpoint(root, MeshGenerator, level, point, title='X point or O point'):
	d = AskForXpoint(root, MeshGenerator, level, point)

	return d.result