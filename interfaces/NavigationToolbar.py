import numpy as np
import tkinter as tk
import matplotlib
matplotlib.use('TkAgg')
import matplotlib as mpl
from matplotlib.patches import Rectangle
try:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2TkAgg as NavigationToolbar2Tk
except:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2Tk

# custom toolbar 

class NavigationToolbar(NavigationToolbar2Tk):
	def __init__(self,canvas_,parent_):
		self.toolitems = (
						    ('Save', 'Save the figure', 'filesave', 'save_figure'),
						    ('Home', 'Reset original view', 'home', 'home'),
						    ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'),
						    ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'))
		NavigationToolbar2Tk.__init__(self,canvas_,parent_)

	def pan(self):
		NavigationToolbar2Tk.pan(self)
		self.mode = "pan" 											#this replace "pan/zoom"
		self.set_message(self.mode)

	def zoom(self):
		NavigationToolbar2Tk.zoom(self)
		self.mode = "zoom" 											#this replace "zoom rect"
		self.set_message(self.mode)

"""
						    ('Home', 'Reset original view', 'home', 'home'),
						    ('Back', 'Back to  previous view', 'back', 'back'),
						    ('Forward', 'Forward to next view', 'forward', 'forward'),
						    (None, None, None, None),
						    ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'),
						    ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'))
#						    (None, None, None, None),
#						    ('Subplots', 'Configure subplots', 'subplots', 'configure_subplots'),
#						    ('Save', 'Save the figure', 'filesave', 'save_figure'),
#						    )
"""
