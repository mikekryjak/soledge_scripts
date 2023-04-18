import tkinter as tk
from tkinter import ttk

class ProgressBar(tk.Tk):

	def __init__(self, parent, Title="Progress bar", Label="Waiting...", Maximum=1., Value=0):
#		self.parent = parent
		self.window	= tk.Toplevel(parent)
		self.window.wm_title(Title)

		self.progress = ttk.Progressbar(self.window, orient="horizontal", length=200, mode="determinate")
		self.progress.pack()
		self.label = ttk.Label(self.window, text=Label)
		self.label.pack()

		self.progress["maximum"] = Maximum
		self.progress["value"] = Value
#		self.parent.update()							#it works  either parent or window
		self.window.update()
		
	def __del__(self):
		self.window.destroy()

	def Update(self, Value, Label=""):
		self.progress["value"] = Value
		if(len(Label)> 0): self.label.text = Label
#		self.parent.update()							#it works  either parent or window
		self.window.update()
