from tkinter import *
from tkinter.messagebox 	import Message

class AskToSplitZones(Toplevel):

	def __init__(self, parent, title = None):

		Toplevel.__init__(self, parent)
		self.transient(parent)

		if(title): self.title(title)

		self.parent = parent

		self.result = ()

		nPadx	= 5
		nPady	= 5
		RowBase = 0
		body_frame = Frame(self)
		body_frame.pack(fill=BOTH, expand=True)
#		self.initial_focus = self.body(body, RowBase, nPadx, nPady)
		RowBase = self.body(body_frame, RowBase, nPadx, nPady)

		RowBase += 1 
		self.buttonbox(body_frame, RowBase, nPadx, nPady)

		self.grab_set()

		self.initial_focus = body_frame

		self.protocol("WM_DELETE_WINDOW", self.cancel)

		self.geometry("+%d+%d" % (parent.winfo_rootx()+200,
								  parent.winfo_rooty()+200))

		self.initial_focus.focus_set()

		self.grab_set()
		self.wait_window(self)


	def body(self, master, RowBase, nPadx, nPady):

		w = Label(master, text="Num. zone to split", width=20)
		w.grid(row=RowBase, column=0, padx=nPadx, pady=nPady, sticky=W)
		
		self.NumZoneEntryValue = StringVar()
		self.NumZoneEntryValue.set("")
		w = Entry(master,  textvariable = self.NumZoneEntryValue, width = 10)
		w.grid(row=RowBase, column=1, padx=nPadx, pady=nPady, sticky=E)

		RowBase += 1
		w = Label(master, text="Num. sub zones", width=20)
		w.grid(row=RowBase, column=0, padx=nPadx, pady=nPady, sticky=W)
		
		self.NumSubZonesEntryValue = StringVar()
		self.NumSubZonesEntryValue.set("")
		w = Entry(master,  textvariable = self.NumSubZonesEntryValue, width = 10)
		w.grid(row=RowBase, column=1, padx=nPadx, pady=nPady, sticky=E)
		return RowBase

	def buttonbox(self, master, RowBase, nPadx, nPady):
		
		# add standard button box. override if you don't want the
		# standard buttons
		w = Button(master, text="Ok", width=10, command=self.OkButton, default=ACTIVE)
		w.grid(row=RowBase, column=0, padx=nPadx, pady=nPady)
		
		w = Button(master, text="Cancel", width=10, command=self.CancelButton)
		w.grid(row=RowBase, column=1, padx=nPadx, pady=nPady)
		
		return RowBase


	def cancel(self, event=None):		# put focus back to the parent window
		self.parent.focus_set()
		self.destroy()

	def OkButton(self, event=None):
		try:
			self.result = (eval(self.NumZoneEntryValue.get()),eval(self.NumSubZonesEntryValue.get()))
		except:
			messagebox.showwarning("Invalid datas")
			return
		
		self.destroy()

	def CancelButton(self, event=None):
		self.destroy()


def ask_to_split_zones(root, title='Get split data'):
	d = AskToSplitZones(root)

	return d.result