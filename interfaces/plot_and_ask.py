from tkinter import *
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg		import FigureCanvasTkAgg
from matplotlib.figure						import Figure
from matplotlib.widgets 					import Cursor
try:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2TkAgg as NavigationToolbar2Tk
except:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2Tk

class PlotAndAsk(Toplevel):

	def __init__(self, parent, LinesData=[], TextsData=[], AskTexts=[], title = None):

		Toplevel.__init__(self, parent)
		self.transient(parent)

		if title:
			self.title(title)

		self.parent = parent

		self.EntryValues = []
		self.result 	 = None

		body = Frame(self)
		self.initial_focus = self.body(body, LinesData, TextsData)
		body.pack(padx=5, pady=5)

		if(len(AskTexts) > 0): self.entrybox(AskTexts)

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


	def body(self, master, LinesData, TextsData):

		Fig = Figure(figsize=(4.8, 5.8), dpi=100)
		Fig.patch.set_facecolor('white')

		# a tk.DrawingArea
		canvas = FigureCanvasTkAgg(Fig, master=master)
		try:
			canvas.show()										#Matplotlib  up 2.2
		except:
			canvas.draw()

		Ax = Fig.add_subplot(111)
		Fig.tight_layout()

		toolbar_frame = Frame(master)
		toolbar = NavigationToolbar2Tk(canvas, toolbar_frame )
		toolbar.update()

		toolbar_frame.pack(fill=X)
		canvas.get_tk_widget().pack(fill=X)

		if(len(LinesData) > 1):
			Xarrs		= LinesData[0]
			Yarrs		= LinesData[1]
			if(len(LinesData) > 2): LineStyles = LinesData[2]
			else:					LineStyles = ['k-' for k in range(len(Xarrs))]

			if(len(LinesData) > 3): 
				LineLabels	= LinesData[3]
				for k in range(len(Xarrs)):
					if((len(Xarrs[k].shape) == 2) and  (len(Yarrs[k].shape) == 2) and (Yarrs[k].shape[1] == 3)):
						Ax.triplot(Xarrs[k][0,:], Xarrs[k][1,:], Yarrs[k], LineStyles[k], label=LineLabels[k])
					else:
						Ax.plot(Xarrs[k], Yarrs[k], LineStyles[k], label=LineLabels[k])

				Ax.legend(fontsize='small', loc='lower left')
			else:
				for k in range(len(Xarrs)):
					if((len(Xarrs[k].shape) == 2) and  (len(Yarrs[k].shape) == 2) and (Yarrs[k].shape[1] == 3)):
						Ax.triplot(Xarrs[k][0,:], Xarrs[k][1,:], Yarrs[k], LineStyles[k])					
					else:
						Ax.plot(Xarrs[k], Yarrs[k], LineStyles[k])


		if(len(TextsData) == 5):
			for k in range(len(TextsData[0])):
				Ax.text(TextsData[0][k], TextsData[1][k], TextsData[2][k], color = TextsData[3][k], fontsize=TextsData[4][k])
		elif(len(TextsData) == 3):
			for k in range(len(TextsData[0])):
				Ax.text(TextsData[0][k], TextsData[1][k], TextsData[2][k], color = 'k', fontsize=18)

		Ax.set_aspect(1.)
		Fig.canvas.draw()

	def entrybox(self, AskTexts):

		MaxLen = 0
		for k in range(len(AskTexts)): MaxLen = max(MaxLen, len(AskTexts[k]))

		for k in range(len(AskTexts)):
			box = Frame(self)

			w=Label(box, text=AskTexts[k], width=MaxLen)
			w.pack(side=LEFT, padx=5, pady=5)
			self.EntryValues.append(StringVar())
			w = Entry(box, textvariable = self.EntryValues[-1], width=4)
			w.pack(side=LEFT, padx=5, pady=5)

			box.pack()

	def buttonbox(self):
		box = Frame(self)

		w = Button(box, text="OK", width=10, command=self.OK, default=ACTIVE)
		w.pack(side=LEFT, padx=5, pady=5)
		w = Button(box, text="Cancel", width=10, command=self.Stop)
		w.pack(side=LEFT, padx=5, pady=5)

		box.pack()

	def cancel(self, event=None):
		# put focus back to the parent window
		self.parent.focus_set(0)
		self.result = (0,)
		self.destroy()

	def OK(self, event=None):
		self.result = (1,)
		for E in self.EntryValues: self.result += (E.get(),)
		self.destroy()

	def Stop(self, event=None):
		self.result = (0,)
		self.destroy()


def plot_and_ask(root, LinesData=[], TextsData=[], AskTexts=[], title="Plot"):
	d = PlotAndAsk(root, LinesData, TextsData, AskTexts, title=title)

	return d.result