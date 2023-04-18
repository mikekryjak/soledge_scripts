from tkinter 				import *
import numpy				as np
from routines.globals		import *

class SetProperties(Toplevel):

	def __init__(self, parent, WallNum, Type, LineType, title = None):

		Toplevel.__init__(self, parent)
		self.result = []

		self.transient(parent)
		self.parent = parent
		self.title("Wall # {:d} properties".format(WallNum))
		self.protocol("WM_DELETE_WINDOW", self.cancel)

		body = Frame(self)
		body.pack(padx=5, pady=5)
		self.initial_focus = self.body(body, Type, LineType)

		self.buttonbox()

		self.wait_visibility()
#		self.initial_focus.focus_set()
		self.grab_set()
		self.wait_window()


	def body(self, master, Type, LineType):

		nPadx  = 5
		nPady  = 5
		xWidth = 10

		iRow = 0

		w = Label(master, text = "Type", width = xWidth)
		w.grid(row=iRow,column=0,  pady=nPady, padx=nPadx)

		self.WallTypeOptions = WALL_TYPE_NAMES
		self.WallTypeOptMenuValue = StringVar()
		self.WallTypeOptMenuValue.set(self.WallTypeOptions[WALL_TYPE_INDEX[Type]])
		w = OptionMenu(master, self.WallTypeOptMenuValue, *self.WallTypeOptions)
		w.grid(row=iRow, column=1, pady=nPady, padx=nPadx, sticky=E+W)

		iRow += 1

		w = Label(master, text = "Line color", width = xWidth)
		w.grid(row=iRow,column=0,  pady=nPady, padx=nPadx)

		self.WallColorOptions = LINE_COLORS
		self.WallColorOptMenuValue = StringVar()
		self.WallColorOptMenuValue.set(self.WallColorOptions[LineType[0]])
		w = OptionMenu(master, self.WallColorOptMenuValue, *self.WallColorOptions)
		w.grid(row=iRow, column=1, pady=nPady, padx=nPadx, sticky=E+W)

		iRow += 1

		w = Label(master, text="Line type", width=xWidth)
		w.grid(row=iRow,column=0, pady=nPady, padx=nPadx)

#		https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/linestyles.html?highlight=linestyle

		self.WallLineOptions = LINE_STYLES
		self.WallLineOptMenuValue = StringVar()
		self.WallLineOptMenuValue.set(self.WallLineOptions[LineType[1]])
		w = OptionMenu(master, self.WallLineOptMenuValue, *self.WallLineOptions)
		w.grid(row=iRow, column=1,  pady=nPady, padx=nPadx, sticky=E+W)

		iRow += 1

		w = Label(master, text="Line width", width=xWidth)
		w.grid(row=iRow,column=0, pady=nPady, padx=nPadx)

#		https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/linestyles.html?highlight=linestyle

		self.WallWidthOptions = LINE_WIDTH
		self.WallWidthOptMenuValue = StringVar()
		self.WallWidthOptMenuValue.set(self.WallWidthOptions[LineType[2]])
		w = OptionMenu(master, self.WallWidthOptMenuValue, *self.WallWidthOptions)
		w.grid(row=iRow, column=1, pady=nPady, padx=nPadx, sticky=E+W)

	def buttonbox(self):
		# add standard button box. override if you don't want the
		# standard buttons

		box = Frame(self)

		w = Button(box, text="Ok", width=10, command=self.OkButton, default=ACTIVE)
		w.pack(side=LEFT, padx=5, pady=5)
		w = Button(box, text="Cancel", width=10, command=self.CancelButton)
		w.pack(side=LEFT, padx=5, pady=5)

		box.pack()

	def cancel(self, event=None):
		# put focus back to the parent window
		self.parent.focus_set()
		self.result = []
		self.destroy()

	def OkButton(self, event=None):
		Type	  = self.get_option_menu_index(self.WallTypeOptMenuValue.get(), self.WallTypeOptions)
		LineColor = self.get_option_menu_index(self.WallColorOptMenuValue.get(), self.WallColorOptions)
		LineStyle = self.get_option_menu_index(self.WallLineOptMenuValue.get(), self.WallLineOptions)
		LineWidth = self.get_option_menu_index(self.WallWidthOptMenuValue.get(), self.WallWidthOptions)
		self.result = [WALL_TYPE_VALUES[Type], LineColor, LineStyle, LineWidth]
		self.destroy()

	def CancelButton(self, event=None):
		self.result = []
		self.destroy()

	def get_option_menu_index(self, Value, ArrValues):
		for k in range(len(ArrValues)):
			if(Value == ArrValues[k]):
				return k
		return -1

def set_properties(root, WallNum, Type, LineType, title='Get values'):
	d = SetProperties(root, WallNum, Type, LineType)

	if(len(d.result) == 0):	return Type, LineType
	else:
		return d.result[0], np.array(d.result[1:], dtype='i4')
