import tkinter as tk
from tkinter import ttk

class show_zones_neighbours(tk.Tk):

	def __init__(self, parent, Zones):
#		self.parent = parent
		self.window	= tk.Toplevel(parent)
		self.window.wm_title("Zone Neighbour")

		self.tree = ttk.Treeview(self.window,)
		self.tree["columns"]=("0","1","2","3","4")
		self.tree['show'] = 'headings'
		Headings = ["Zone", "North", "South", "East", "West" ]
		xWidths  = 50
		for k in range(5): 
			self.tree.column("{:d}".format(k),  width=xWidths, anchor=tk.CENTER)
			self.tree.heading("{:d}".format(k), text = Headings[k], anchor=tk.CENTER)

		for k in range(len(Zones)):
			Values = ("{:d}".format(k+1), 
					  "{:d}".format(Zones[k].Neighbour.north+1), "{:d}".format(Zones[k].Neighbour.south+1),
				 	  "{:d}".format(Zones[k].Neighbour.east+1), "{:d}".format(Zones[k].Neighbour.west+1))
			self.tree.insert("" , "end",    text="", values=Values)																	

		self.tree.pack()
		self.window.update()
		
	def __del__(self):
		self.window.destroy()

