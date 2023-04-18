from tkinter import *

class AskName(Toplevel):

	def __init__(self, parent, labels, values=[], title = None):

		Toplevel.__init__(self, parent)
		self.transient(parent)

		if(title):
			self.title(title)

		self.parent = parent

		self.result = None

		body = Frame(self)
		body.pack(padx=5, pady=5)
		self.initial_focus = self.body(body, labels, values)

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


	def body(self, master, labels, values):

		nPadx = 5
		nPady = 5
		self.EntryValues = []
		for i in range(len(labels)):
			w = Label(master, text=labels[i], width=20)
			w.grid(row=i, column=0,  pady=nPady, padx=nPadx, sticky=E+W+N)
			
			self.EntryValues.append(StringVar())
			if(len(values) > i): self.EntryValues[-1].set(values[i])
			else:				 self.EntryValues[-1].set("")

			w = Entry(master,  textvariable = self.EntryValues[-1], width = 20)
			w.grid(row=i, column=1,  pady=nPady, padx=nPadx, sticky=E+W+N)

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
		self.destroy()

	def OkButton(self, event=None):
		self.result = []
		for EntryValue in self.EntryValues: self.result.append(EntryValue.get())
		self.destroy()

	def CancelButton(self, event=None):
		self.result = []
		self.destroy()


def ask_values(root, labels, values=[], title='Get values'):
	d = AskName(root, labels, values, title=title)

	return d.result