from tkinter import *

class AskName(Toplevel):

	def __init__(self, parent, label, value="", title = None):

		Toplevel.__init__(self, parent)
		self.transient(parent)

		if(title):
			self.title(title)

		self.parent = parent

		self.result = None

		body = Frame(self)
		self.initial_focus = self.body(body, label, value)
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


	def body(self, master, label, value):

		w = Label(master, text=label, width=10)
		w.pack(side=LEFT, padx=5, pady=5)
		
		self.NameEntryValue = StringVar()
		self.NameEntryValue.set(value)
		w = Entry(master,  textvariable = self.NameEntryValue, width = 10)
		w.pack(side=LEFT, padx=5, pady=5)

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
		self.result = self.NameEntryValue.get()
		self.destroy()

	def CancelButton(self, event=None):
		self.result = ""
		self.destroy()


def ask_name(root, label, value="", title='Get name'):
	d = AskName(root, label, value)

	return d.result