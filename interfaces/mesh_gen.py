
from tkinter		import Toplevel, Tk, Text, OptionMenu, filedialog, messagebox, BOTH, W, N, E, S, HORIZONTAL, LEFT
from tkinter		import DoubleVar, IntVar, StringVar
from tkinter.ttk	import Entry, Menubutton, Checkbutton, Frame, Button, Label, Style, LabelFrame, Scale, Separator

import types
from tkinter.messagebox import Message

from math 									import floor, sqrt
import numpy								as np
from numpy									import arange, sin, pi

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg		import FigureCanvasTkAgg
from matplotlib.figure						import Figure
from matplotlib.widgets 					import Cursor
from matplotlib.path	 					import Path
import matplotlib.pyplot					as plt
try:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2TkAgg as NavigationToolbar2Tk
except:
	from matplotlib.backends.backend_tkagg		import NavigationToolbar2Tk

from mesh.mesh_all_zones					import mesh_all_zones
from mesh.conform_mesh						import conform_mesh
from mesh.tweak_mesh						import tweak_mesh
from mesh.subdivide_selected_segment		import subdivide_selected_segment
from mesh.mesh_segment						import mesh_segment
from mesh.regularize_mesh_size				import regularize_mesh_size
from mesh.generate_soledge_mesh				import generate_soledge_mesh
from mesh.vertical_split					import vertical_split
from mesh.horizontal_split					import horizontal_split
from mesh.recompute_megazones				import recompute_megazones
from mesh.get_magzone_of_point				import get_zone_of_point

from routines.find_closest_segment			import find_closest_segment
from routines.utils_routines				import safe_eval
from routines.utils_walls					import plot2d_walls

from interfaces.show_zones_neighbours		import show_zones_neighbours
from interfaces.eirene_gen					import eirene_gen
from interfaces.ask_to_split_zones			import ask_to_split_zones

class mesh_gen(Frame):

	def __init__(self, parent, root):
		Frame.__init__(self, parent)

		self.SkipEvents = True

		self.root = root

		nPadx  = 5
		nPady  = 5
		xWidth = 12

		self.win_frame = Frame(self)
		self.win_frame.pack(fill=BOTH, expand=True)

		self.win_frame.columnconfigure(1, weight=1)

		PlotFrame = Frame(self.win_frame)
		PlotFrame.grid(row=0, column=0, rowspan=4,  pady=nPady, padx=nPadx, sticky=E+W+N)

#		Plot area

		self.Fig = Figure(figsize=(4.8, 5.8), dpi=100)
		self.Fig.patch.set_facecolor('white')

		# a tk.DrawingArea
		canvas = FigureCanvasTkAgg(self.Fig, master=PlotFrame)
		try:
			canvas.show()										#Matplotlib  up 2.2
		except:
			canvas.draw()
		canvas.get_tk_widget().grid(row=1, column=0,  columnspan=3,  padx=nPadx, sticky=E+W+S+N)

		self.Ax = self.Fig.add_subplot(111)
		self.Cb = 0														#I need this to be able to remove it later
		self.Fig.tight_layout()

		toolbar_frame = Frame(PlotFrame)
		toolbar_frame.grid(row=0, column=0, columnspan=3,  padx=nPadx, sticky=E+W+S+N)
		toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
		toolbar.update()

		self.InfoLabelValue = StringVar()
		InfoLabel = Label(PlotFrame, textvariable = self.InfoLabelValue, justify = LEFT, width = 30)
		InfoLabel.grid(row=2, column=0,  pady=nPady, padx=nPadx)
		self.InfoLabelValue.set("")

		ShowLabel = Label(PlotFrame, text = "Show", anchor=E, width = xWidth)
		ShowLabel.grid(row=2, column=1,  pady=nPady, padx=nPadx)

		self.ShowInFigureOptions = ["Segments", "Grid", "Quality", "Size", "Cent. Mesh", "Chi", "Br", "Bz", "Bt", "x mesh", "z mesh","Zones Map","Zones Bar", "Zones","Mesh Index."]
		self.ShowInFigureOptMenuValue = StringVar()
		self.ShowInFigureOptMenuValue.set(self.ShowInFigureOptions[0])
		opts = ()
		for k in self.ShowInFigureOptions: opts += (k,)	
		self.ShowInFigureOptMenu = OptionMenu(PlotFrame, self.ShowInFigureOptMenuValue, *opts)
		self.ShowInFigureOptMenu.grid(row=2, column=2,  pady=nPady, sticky=W+E)
		self.ShowInFigureOptMenuValue.trace("w", self.cmd_show_in_figure)

		
#		Tools area
		
		ToolsFrame = LabelFrame(self.win_frame, text="Tools")
		ToolsFrame.grid(row=0, column=1,  pady=nPady, padx=nPadx, sticky=E+W+N)

		iFrow = 0
		
		ClearSegMeshBut = Button(ToolsFrame, text="Clear Segments", command = self.cmd_clear_segments, width = xWidth)
		ClearSegMeshBut.grid(row=iFrow, column=0, padx=nPadx, pady=nPady)
		
		ClearGridMeshBut = Button(ToolsFrame, text="Clear Grid", command = self.cmd_clear_grid, width = xWidth)
		ClearGridMeshBut.grid(row=iFrow, column=1, padx=nPadx, pady=nPady)

		ShowNeighLabel = Label(ToolsFrame, text = "Show Neighbors")
		ShowNeighLabel.grid(row=iFrow,column=2,  pady=nPady, padx=nPadx)

		self.ShowNeighCkButValue = IntVar()
		self.ShowNeighCkButValue.set(0)
		ShowNeighCkBut = Checkbutton(ToolsFrame, text="", variable = self.ShowNeighCkButValue, command = self.cmd_show_neighbors, width = 2*xWidth)
		ShowNeighCkBut.grid(row=iFrow,column=3, padx=nPadx,  pady=nPady)

#		Segment MeshGen area
		
		SegmentMeshGenFrame = LabelFrame(self.win_frame, text="Segments")
		SegmentMeshGenFrame.grid(row=1, column=1,  pady=nPady, padx=nPadx, sticky=E+W+N)

		iFrow = 0
		
		SelectSegmentBut = Button(SegmentMeshGenFrame, text="Select segment", command = self.cmd_select_segment, width = xWidth)
		SelectSegmentBut.grid(row=iFrow, column=0, pady=nPady)
		
		SelectSegmentBut = Button(SegmentMeshGenFrame, text="Deselect all", command = self.cmd_deselect_all, width = xWidth)
		SelectSegmentBut.grid(row=iFrow, column=1, pady=nPady)
		
		SelectSegmentBut = Button(SegmentMeshGenFrame, text="Build segment", command = self.cmd_mesh_segment, width = xWidth)
		SelectSegmentBut.grid(row=iFrow, column=2, pady=nPady)

		iFrow += 1

		NumPointsLabel = Label(SegmentMeshGenFrame, text = "Number of points", width = 1.5*xWidth)
		NumPointsLabel.grid(row=iFrow, column=0, columnspan=1,  pady=nPady, padx=nPadx)

		self.NumPointsEntryValue = IntVar()
		NumPointsEntry = Entry(SegmentMeshGenFrame,  textvariable = self.NumPointsEntryValue, width = xWidth)
		NumPointsEntry.grid(row=iFrow, column=1, columnspan=1,  pady=nPady, padx=nPadx)

		iFrow += 1

		RefinementLabel = Label(SegmentMeshGenFrame, text = "Refinement", width = xWidth)
		RefinementLabel.grid(row=iFrow, column=0,  pady=nPady, padx=nPadx)


		RefinementSeparator = Separator(SegmentMeshGenFrame)
		RefinementSeparator.grid(row=iFrow, column=1,  columnspan=4,  pady=nPady, padx=0, sticky=E+W)

		iFrow += 1

		RefinementTypeLabel = Label(SegmentMeshGenFrame, text = "Type", width = xWidth)
		RefinementTypeLabel.grid(row=iFrow,column=0,  pady=nPady, padx=nPadx)

		self.RefinementTypeOptions = ["None", "Linear", "Bezier", "Spline"]
		self.RefinementTypeOptMenuValue = StringVar()
		self.RefinementTypeOptMenuValue.set(self.RefinementTypeOptions[0])
		opts = ()
		for k in self.RefinementTypeOptions: opts += (k,)	
		self.RefinementTypeOptMenu = OptionMenu(SegmentMeshGenFrame, self.RefinementTypeOptMenuValue, *opts)
		self.RefinementTypeOptMenu.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx, sticky=E+W)
		self.RefinementTypeOptMenuValue.trace("w", self.cmd_refinement_type)

		iFrow += 1
		
		RefinementSideLabel = Label(SegmentMeshGenFrame, text = "Side", width = xWidth)
		RefinementSideLabel.grid(row=iFrow,column=0,  pady=nPady, padx=nPadx)

		self.RefinementSideOptions = ["Left", "Right", "Both"]
		self.RefinementSideOptMenuValue = StringVar()
		self.RefinementSideOptMenuValue.set(self.RefinementSideOptions[0])
		opts = ()
		for k in self.RefinementSideOptions: opts += (k,)	
		self.RefinementSideOptMenu = OptionMenu(SegmentMeshGenFrame, self.RefinementSideOptMenuValue, *opts)
		self.RefinementSideOptMenu.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx, sticky=E+W)
#		self.RefinementSideOptMenu.state(["disabled"])

		RefinementModeLabel = Label(SegmentMeshGenFrame, text = "Mode", width = xWidth)
		RefinementModeLabel.grid(row=iFrow,column=2,  pady=nPady, padx=nPadx)

		self.RefinementModeOptions = ["Relative", "Absolute"]
		self.RefinementModeOptMenuValue = StringVar()
		self.RefinementModeOptMenuValue.set(self.RefinementModeOptions[0])
		opts = ()
		for k in self.RefinementModeOptions: opts += (k,)	
		self.RefinementModeOptMenu = OptionMenu(SegmentMeshGenFrame, self.RefinementModeOptMenuValue, *opts)
		self.RefinementModeOptMenu.grid(row=iFrow, column=3,  pady=nPady, padx=nPadx, sticky=E+W)
		self.RefinementModeOptMenuValue.trace("w", self.cmd_refinement_mode)
#		self.RefinementModeOptMenu.state(["disabled"])

		iFrow += 1

		RefinementSizeLabel = Label(SegmentMeshGenFrame, text = "Refinement size", width = xWidth)
		RefinementSizeLabel.grid(row=iFrow,column=0,  pady=nPady, padx=nPadx)

		self.RefinementSizeLeftLabelValue = StringVar()
		self.RefinementSizeLeftLabelValue.set("Left")
		RefinementSizeLeftLabel = Label(SegmentMeshGenFrame, textvariable = self.RefinementSizeLeftLabelValue, width = xWidth)
		RefinementSizeLeftLabel.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx)

		self.RefinementSizeLeftEntryValue = DoubleVar()
		self.RefinementSizeLeftEntry = Entry(SegmentMeshGenFrame,  textvariable = self.RefinementSizeLeftEntryValue, width = xWidth)
		self.RefinementSizeLeftEntry.grid(row=iFrow, column=2,  pady=nPady, padx=nPadx)
		self.RefinementSizeLeftEntry.state(["disabled"])

		self.RefinementSizeRightLabelValue = StringVar()
		self.RefinementSizeRightLabelValue.set("Right")
		RefinementSizeRightLabel = Label(SegmentMeshGenFrame, textvariable = self.RefinementSizeRightLabelValue, width = xWidth)
		RefinementSizeRightLabel.grid(row=iFrow, column=3,  pady=nPady, padx=nPadx)

		self.RefinementSizeRightEntryValue = DoubleVar()
		self.RefinementSizeRightEntry = Entry(SegmentMeshGenFrame,  textvariable = self.RefinementSizeRightEntryValue, width = xWidth)
		self.RefinementSizeRightEntry.grid(row=iFrow, column=4,  pady=nPady, padx=nPadx)
		self.RefinementSizeRightEntry.state(["disabled"])

		iFrow += 1

		AlignementsLabel = Label(SegmentMeshGenFrame, text = "Alignements", width = xWidth)
		AlignementsLabel.grid(row=iFrow,column=0,  pady=nPady, padx=nPadx)

		AlignementsSeparator = Separator(SegmentMeshGenFrame)
		AlignementsSeparator.grid(row=iFrow,column=1,  columnspan=4,  pady=nPady, padx=0, sticky=W+E)

		iFrow += 1

		AlignementModeLabel = Label(SegmentMeshGenFrame, text = "Align with:", width = xWidth)
		AlignementModeLabel.grid(row=iFrow, column=0,  pady=nPady, padx=nPadx)

		self.AlignementTypeOptions = ["Flux", "Wall"]
		self.AlignementTypeOptMenuValue = StringVar()
		self.AlignementTypeOptMenuValue.set(self.AlignementTypeOptions[0])
		self.AlignementTypeOptMenu = OptionMenu(SegmentMeshGenFrame, self.AlignementTypeOptMenuValue, self.AlignementTypeOptions[0], self.AlignementTypeOptions[1])
		self.AlignementTypeOptMenu.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx, sticky=E+W)
		self.AlignementTypeOptMenuValue.trace("w", self.cmd_align_type)

		ShowForceLabel = Label(SegmentMeshGenFrame, text = "Force Ortho")
		ShowForceLabel.grid(row=iFrow,column=2,  pady=nPady, padx=nPadx)

		self.ForceOrthoCkButValue = IntVar()
		self.ForceOrthoCkButValue.set(0)
		ForceOrthoCkBut = Checkbutton(SegmentMeshGenFrame, text="", variable = self.ForceOrthoCkButValue, width = 2*xWidth)
		ForceOrthoCkBut.grid(row=iFrow,column=3, padx=nPadx,  pady=nPady)

		iFrow += 1

		PsiMinLabel = Label(SegmentMeshGenFrame, text = "Psi min:", width = xWidth)
		PsiMinLabel.grid(row=iFrow, column=0,  pady=nPady, padx=nPadx)

		self.PsiMinEntryValue = StringVar()
		self.PsiMinEntry = Entry(SegmentMeshGenFrame,  textvariable = self.PsiMinEntryValue, width = xWidth)
		self.PsiMinEntry.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx)
		self.PsiMinEntry.state(["disabled"])

		self.PsiMinScaleValue = DoubleVar()
		self.PsiMinScale = Scale(SegmentMeshGenFrame, orient= HORIZONTAL, from_ = 0., to_ = 1., \
								variable = self.PsiMinScaleValue, command = self.cmd_move_flux_min_scale)
		self.PsiMinScale.grid(row=iFrow,column=2,  columnspan=2, pady=nPady, padx=nPadx, sticky=E+W)
		self.PsiMinScale.state(["disabled"])

		iFrow += 1

		PsiMaxLabel = Label(SegmentMeshGenFrame, text = "Psi max:", width = xWidth)
		PsiMaxLabel.grid(row=iFrow,column=0,  pady=nPady, padx=nPadx)

		self.PsiMaxEntryValue = StringVar()
		self.PsiMaxEntry = Entry(SegmentMeshGenFrame,  textvariable = self.PsiMaxEntryValue, width = xWidth)
		self.PsiMaxEntry.grid(row=iFrow,column=1,  pady=nPady, padx=nPadx)
		self.PsiMaxEntry.state(["disabled"])

		self.PsiMaxScaleValue = DoubleVar()
		self.PsiMaxScale = Scale(SegmentMeshGenFrame, orient= HORIZONTAL, from_ = 0., to_ = 1., \
								variable = self.PsiMaxScaleValue, command = self.cmd_move_flux_max_scale)
		self.PsiMaxScale.grid(row=iFrow,column=2,  columnspan=2, pady=nPady, padx=nPadx, sticky=E+W)
		self.PsiMaxScale.state(["disabled"])

#		Optimization area
		
		OptimizationFrame = LabelFrame(self.win_frame, text="Optimization parameters")
		OptimizationFrame.grid(row=2, column=1,  pady=nPady, padx=nPadx, sticky=E+W+N)

		xWidth	= 11
		iFrow 	= 0

		OrthogonalLabel = Label(OptimizationFrame, text = "Orthogonal", width = xWidth)
		OrthogonalLabel.grid(row=iFrow, column=0,  pady=nPady, padx=nPadx)

		self.OptimizationEntryValue = StringVar()
		self.OptimizationEntry = Entry(OptimizationFrame,  textvariable = self.OptimizationEntryValue, width = 8)
		self.OptimizationEntry.grid(row=iFrow, column=1,  pady=nPady, padx=nPadx)

		self.OptimizationScaleValue = DoubleVar()
		OptimizationScale = Scale(OptimizationFrame, orient= HORIZONTAL, from_ = 0., to_ = 1., \
								variable = self.OptimizationScaleValue, command = self.cmd_move_optimization_scale)
		OptimizationScale.grid(row=iFrow,column=2,  columnspan=3, pady=nPady, padx=nPadx, sticky=E+W)

		EquidistantLabel = Label(OptimizationFrame, text = "Equidistant", width = xWidth)
		EquidistantLabel.grid(row=iFrow,column=5,  pady=nPady, padx=nPadx)

		iFrow += 1

		OptdMinLabel = Label(OptimizationFrame, text = "d min (mm)", width = xWidth)
		OptdMinLabel.grid(row=iFrow,column=0,  pady=nPady, padx=nPadx)

		self.OptdMinEntryValue = DoubleVar()
		OptdMinEntry = Entry(OptimizationFrame,  textvariable = self.OptdMinEntryValue, width = 8)
		OptdMinEntry.grid(row=iFrow,column=1,  pady=nPady, padx=nPadx)

		OptLenghtLabel = Label(OptimizationFrame, text = "Lenght (1000)", width = xWidth)
		OptLenghtLabel.grid(row=iFrow,column=2,  pady=nPady, padx=nPadx)

		self.OptLenghtEntryValue = DoubleVar()
		OptLenghtEntry = Entry(OptimizationFrame,  textvariable = self.OptLenghtEntryValue, width = 8)
		OptLenghtEntry.grid(row=iFrow,column=3,  pady=nPady, padx=nPadx)


		OptTargetLabel = Label(OptimizationFrame, text = "Target (0.75)", width = xWidth)
		OptTargetLabel.grid(row=iFrow,column=4,  pady=nPady, padx=nPadx)

		self.OptTargetEntryValue = DoubleVar()
		OptTargetEntry = Entry(OptimizationFrame,  textvariable = self.OptTargetEntryValue, width = 8)
		OptTargetEntry.grid(row=iFrow,column=5,  pady=nPady, padx=nPadx)

		iFrow += 1

		OptdRangeLabel = Label(OptimizationFrame, text = "Range (0.2)", width = xWidth)
		OptdRangeLabel.grid(row=iFrow,column=0,  pady=nPady, padx=nPadx)

		self.OptRangeEntryValue = DoubleVar()
		OptRangeEntry = Entry(OptimizationFrame,  textvariable = self.OptRangeEntryValue, width = 8)
		OptRangeEntry.grid(row=iFrow,column=1,  pady=nPady, padx=nPadx)

		OptAlphaLabel = Label(OptimizationFrame, text = "Alpha (0.1)", width = xWidth)
		OptAlphaLabel.grid(row=iFrow,column=2,  pady=nPady, padx=nPadx)

		self.OptAlphaEntryValue = DoubleVar()
		OptAlphaEntry = Entry(OptimizationFrame,  textvariable = self.OptAlphaEntryValue, width = 8)
		OptAlphaEntry.grid(row=iFrow,column=3,  pady=nPady, padx=nPadx)

		OptdrMinLabel = Label(OptimizationFrame, text = "dr min (mm)", width = xWidth)
		OptdrMinLabel.grid(row=iFrow,column=4,  pady=nPady, padx=nPadx)

		self.OptdrMinEntryValue = DoubleVar()
		OptdrMinEntry = Entry(OptimizationFrame,  textvariable = self.OptdrMinEntryValue, width = 8)
		OptdrMinEntry.grid(row=iFrow,column=5,  pady=nPady, padx=nPadx)


#		Mesh frame
		
		MeshFrame = LabelFrame(self.win_frame, text="Mesh")
		MeshFrame.grid(row=3, column=1,  pady=nPady, padx=nPadx, sticky=E+W+N)
		
		iFrow = 0
		xWidth = 12
	
		MeshAllBut = Button(MeshFrame, text="Generate Mesh", command = self.cmd_mesh_all)
		MeshAllBut.grid(row=iFrow, column=0, columnspan=2, pady=nPady, padx=nPadx, sticky=E+W)
		
		RegMeshBut = Button(MeshFrame, text="Regolarize", command = self.cmd_reg_mesh, width = xWidth)
		RegMeshBut.grid(row=iFrow, column=2, pady=nPady, padx=nPadx)
		
		TweakMeshBut = Button(MeshFrame, text="Tweak all", command = self.cmd_tweak_mesh, width = xWidth)
		TweakMeshBut.grid(row=iFrow, column=3, pady=nPady, padx=nPadx)
		
		TweakMeshBut = Button(MeshFrame, text="Tweak ortho", command = self.cmd_tweak_mesh_ortho, width = xWidth)
		TweakMeshBut.grid(row=iFrow, column=4, pady=nPady, padx=nPadx)

#		Penalization area
		iFrow += 1
		
		SoledgeMeshBut = Button(MeshFrame, text="Generate Soledge Mesh", command = self.cmd_soledge_mesh)
		SoledgeMeshBut.grid(row=iFrow, column=0, columnspan=2, pady=nPady, padx=nPadx, sticky=E+W)
		
		self.PenalizationCkButValue = IntVar()
		self.PenalizationCkButValue.set(1)
		ShowNeighCkBut = Checkbutton(MeshFrame, text="Penalization", variable = self.PenalizationCkButValue)
		ShowNeighCkBut.grid(row=iFrow,column=3, padx=nPadx,  pady=nPady)

		self.ReversePenalizationButValue = StringVar()
		self.ReversePenalizationButValue.set("Rev. Penaliz.")
		ReversePenalizationBut = Button(MeshFrame, textvariable=self.ReversePenalizationButValue, command = self.cmd_reverse_penalization, width = xWidth)
		ReversePenalizationBut.grid(row=iFrow, column=4, pady=nPady, padx=nPadx)

#		Split MagZones area
		iFrow += 1

		self.NumZonesLabelValue = StringVar()
		self.NumZonesLabelValue.set("N. zones: 0")
		NumZonesLabel = Label(MeshFrame, textvariable = self.NumZonesLabelValue, width = xWidth)
		NumZonesLabel.grid(row=iFrow,column=0,  pady=nPady, padx=nPadx)
		
		SplitVertBut = Button(MeshFrame, text="Split Vert.", command = self.cmd_split_vertical, width = xWidth)
		SplitVertBut.grid(row=iFrow, column=1, pady=nPady, padx=nPadx)

		SplitHorBut = Button(MeshFrame, text="Split Hor.", command = self.cmd_split_horizontal, width = xWidth)
		SplitHorBut.grid(row=iFrow, column=2, pady=nPady, padx=nPadx)

		SplitResBut = Button(MeshFrame, text="Split Reset", command = self.cmd_split_reset, width = xWidth)
		SplitResBut.grid(row=iFrow, column=3, pady=nPady, padx=nPadx)

#		Starting default values
		self.resize_mode	= True
		self.SkipEvents		= False
		self.sel_mode		= False				#segment selection mode
		self.sel_to_mesh	= False				#if true segment was selected and it is ready to be meshed
		self.align_psimode	= False
		self.force_ortho	= False
		self.zone_sel		= -1
		self.side_sel		= -1
		self.sub_sel		= -1
		self.old_plot_type	= -1
		self.InButton			   = False
		self.InReversePenalization = False
		
		modeAdjust		= 0
		self.RefinementSizeLeftEntryValue.set(0.)
		self.RefinementSizeRightEntryValue.set(0.)
		
		self.ZoneNeighbours = 0

#		self.OptdMinEntryValue.set("{:.3f}".format(dmin_ad))
#		self.OptLenghtEntryValue.set("{:.1f}".format(length_ad))
#		self.OptTargetEntryValue.set("{:.3f}".format(target_ad))
#		self.OptRangeEntryValue.set("{:.3f}".format(range_ad))
#		self.OptAlphaEntryValue.set("{:.3f}".format(alpha_ad))
#		self.OptdrMinEntryValue.set("{:.3f}".format(drmin_ad))


		self.PsiMinEntryValue.trace("w", self.cmd_psi_min_entry)
		self.PsiMaxEntryValue.trace("w", self.cmd_psi_max_entry)

		self.OptimizationEntryValue.trace("w", self.cmd_optimization_entry)
		self.OptdMinEntryValue.trace("w", self.cmd_opt_dmin_entry)
		self.OptLenghtEntryValue.trace("w", self.cmd_opt_lenght_entry)
		self.OptTargetEntryValue.trace("w", self.cmd_opt_target_entry)
		self.OptRangeEntryValue.trace("w", self.cmd_opt_range_entry)
		self.OptAlphaEntryValue.trace("w", self.cmd_opt_alpha_entry)
		self.OptdrMinEntryValue.trace("w", self.cmd_opt_drmin_entry)

		self.update_data()	
			
#	Internal functions
#===============================================================================

	def set_optmizazion_values(self):
		self.SkipEvents 	= True
		Optimization = self.root.Config.Optimization
		self.OptimizationScaleValue.set(1.-Optimization.orthoptim)
		self.OptimizationEntryValue.set("{:.3f}".format(1.-Optimization.orthoptim))
		self.OptdMinEntryValue.set(Optimization.dmin_ad)
		self.OptLenghtEntryValue.set(Optimization.length_ad)
		self.OptTargetEntryValue.set(Optimization.target_ad)
		self.OptRangeEntryValue.set(Optimization.range_ad)
		self.OptAlphaEntryValue.set(Optimization.alpha_ad)
		self.OptdrMinEntryValue.set(Optimization.drmin_ad)
		self.SkipEvents 	= False

	def cmd_donothing(self, empty):
		return

#	clear segments and grid

	def cmd_clear_segments(self):
		self.cmd_clear_grid()

		self.align_psimode = False
		
		MagZones		= self.root.Config.MagZones
		MagMegazones	= self.root.Config.MagMegazones
		MagPMegazones	= self.root.Config.MagPMegazones
		
		self.root.Config.OMPseg.ismeshed = False
		for k in range(len(MagMegazones)):
			MagMegazones[k].ismeshed = False

		for k in range(len(MagPMegazones)):
			MagPMegazones[k].ismeshed 			= False
			MagPMegazones[k].isaligned			= False
			MagPMegazones[k].meshchanged		= True
			if(hasattr(MagPMegazones[k], 'refpoints')):
				MagPMegazones[k].refpoints.R		= np.array([])
				MagPMegazones[k].refpoints.Z		= np.array([])
			if(hasattr(MagPMegazones[k], 'subrefpoints')):
				MagPMegazones[k].subrefpoints[0].R	= np.array([])
				MagPMegazones[k].subrefpoints[0].Z	= np.array([])
				MagPMegazones[k].subrefpoints[1].R	= np.array([])
				MagPMegazones[k].subrefpoints[1].Z	= np.array([])

		for k in range(len(MagZones)):
			MagZones[k].northaligned		= False
			MagZones[k].southaligned		= False
			if(hasattr(MagZones[k], 'north')): MagZones[k].north.ismeshed	= False
			if(hasattr(MagZones[k], 'south')): MagZones[k].south.ismeshed	= False
			if(hasattr(MagZones[k], 'east')):  MagZones[k].east.ismeshed 	= False
			if(hasattr(MagZones[k], 'west')):  MagZones[k].west.ismeshed 	= False

		self.root.set_Segments_OK(self.root.Config, False)	
		self.update_figure()

#	clear grid
	
	def cmd_clear_grid(self):
		MagZones		= self.root.Config.MagZones
		MagMegazones	= self.root.Config.MagMegazones
		MagPMegazones	= self.root.Config.MagPMegazones
		for k in range(len(MagZones)):
			MagZones[k].orthomeshchanged	= True
			if(hasattr(MagZones[k], 'gridR')):
				nPsi			= len(MagMegazones[MagZones[k].mz].refpoints.psi)
				nTheta  		= len(MagPMegazones[MagZones[k].pmz].refpoints.R)
				MagZones[k].gridR	= np.empty((nPsi, nTheta), dtype='f8')
				MagZones[k].gridZ	= np.empty((nPsi, nTheta), dtype='f8')

		self.root.set_MagGrid_OK(self.root.Config, False)
		self.update_figure()
			
	def cmd_show_neighbors(self):
		if(self.root.Config.MagZones_OK) and (self.ShowNeighCkButValue.get() and (self.ZoneNeighbours == 0)):
			self.ZoneNeighbours = show_zones_neighbours(self, self.root.Config.MagZones)
		else:
			if(self.ZoneNeighbours != 0): self.ZoneNeighbours = 0


	def cmd_select_segment(self):
		self.cid_mouse_button_press = self.Fig.canvas.mpl_connect('button_press_event', self.cmd_mouse_button_press_cut_selected)
		self.cursor = Cursor(self.Ax, useblit=False, color='red', linewidth=1)

	def cmd_mouse_button_press_cut_selected(self, event):
		self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press)
		self.cid_mouse_button_press = 0
		self.cursor=0

#		find closest edge
		
		MagZones		= self.root.Config.MagZones
		MagMegazones	= self.root.Config.MagMegazones
		MagPMegazones	= self.root.Config.MagPMegazones
		OMPseg			= self.root.Config.OMPseg

		self.sel_mode	= False
		Rpt			 	= event.xdata
		Zpt		 		= event.ydata
		dmin     		= 1e10
		self.zone_sel	= -1
		self.side_sel	= -1
		nMagZones			= len(MagZones)
		for k in range(nMagZones):
#			north
			d, i = find_closest_segment(Rpt, Zpt, MagZones[k].north.R, MagZones[k].north.Z)
			if(d < dmin):
				dmin 			= d
				self.zone_sel	= k
				self.side_sel	= 0
				
#			south		
			d, i = find_closest_segment(Rpt, Zpt, MagZones[k].south.R, MagZones[k].south.Z)
			if(d < dmin):
				dmin			= d
				self.zone_sel	= k
				self.side_sel	= 1
				
#			east
			d, i = find_closest_segment(Rpt, Zpt, MagZones[k].east.R, MagZones[k].east.Z)
			if(d < dmin):
				dmin			= d
				self.zone_sel	= k
				self.side_sel	= 2
				
#			west
			d, i = find_closest_segment(Rpt, Zpt, MagZones[k].west.R, MagZones[k].west.Z)
			if(d < dmin):
				dmin			= d
				self.zone_sel	= k
				self.side_sel	= 3

#		Check for Special OMP segment
		d, i = find_closest_segment(Rpt, Zpt, OMPseg.R, OMPseg.Z)
		if(d < dmin):
			dmin 			= d
			self.zone_sel	= len(MagZones)								# Conventional values for OMP special segment
			self.side_sel	= 2

		self.sel_to_mesh	= True
		self.sel_mode		= True
		self.Rpt			= Rpt
		self.Zpt			= Zpt		
		isSet				= False
		self.SkipEvents 	= True

		if(self.zone_sel < len(MagZones)):
			pmz	= MagZones[self.zone_sel].pmz
			mz	= MagZones[self.zone_sel].mz
			if(self.side_sel < 2):
				if(((self.side_sel == 0) and hasattr(MagZones[self.zone_sel],"northaligned") and MagZones[self.zone_sel].northaligned) or
				   ((self.side_sel == 1) and hasattr(MagZones[self.zone_sel],"southaligned") and MagZones[self.zone_sel].southaligned)):

#					self.SkipEvents = False
					self.AlignementTypeOptMenuValue.set(self.AlignementTypeOptions[1])
#					self.SkipEvents = True
#					self.set_segment_alignment()

					if(hasattr(MagPMegazones[pmz],"subrefpoints") and hasattr(MagPMegazones[pmz].subrefpoints[self.sub_sel],"nPoints")):
						AlignementType	= 1
						pt_num			= MagPMegazones[pmz].subrefpoints[self.sub_sel].nPoints
						refineType		= MagPMegazones[pmz].subrefpoints[self.sub_sel].RefineType 
						refineSide		= MagPMegazones[pmz].subrefpoints[self.sub_sel].RefineSide 
						adjustMode		= MagPMegazones[pmz].subrefpoints[self.sub_sel].AdjustMode 
						paramL			= MagPMegazones[pmz].subrefpoints[self.sub_sel].ParamL 
						paramR			= MagPMegazones[pmz].subrefpoints[self.sub_sel].ParamR 
						ForceOrtho		= MagPMegazones[pmz].ForceOrtho
						isSet			= True
				else:
					if(MagPMegazones[pmz].isaligned):
						messagebox.showwarning("MeshGen", "In the same poloidal zone there is\na different aligned border.\nYou must first remove such alignement")
						return
					elif(hasattr(MagPMegazones[pmz],"refpoints") and hasattr(MagPMegazones[pmz].refpoints,"nPoints")):
						AlignementType	= 0
						pt_num			= MagPMegazones[pmz].refpoints.nPoints
						refineType		= MagPMegazones[pmz].refpoints.RefineType 
						refineSide		= MagPMegazones[pmz].refpoints.RefineSide 
						adjustMode		= MagPMegazones[pmz].refpoints.AdjustMode 
						paramL			= MagPMegazones[pmz].refpoints.ParamL 
						paramR			= MagPMegazones[pmz].refpoints.ParamR 
						ForceOrtho		= MagPMegazones[pmz].ForceOrtho
						isSet		= True
			elif((self.side_sel > 1) and hasattr(MagMegazones[mz],"refpoints") and hasattr(MagMegazones[mz].refpoints,"nPoints")):
				AlignementType	= 0
				pt_num			= MagMegazones[mz].refpoints.nPoints
				refineType		= MagMegazones[mz].refpoints.RefineType 
				refineSide		= MagMegazones[mz].refpoints.RefineSide 
				adjustMode		= MagMegazones[mz].refpoints.AdjustMode 
				paramL			= MagMegazones[mz].refpoints.ParamL 
				paramR			= MagMegazones[mz].refpoints.ParamR 
				ForceOrtho		= 0 
				isSet		= True
		elif(hasattr(OMPseg,"refpoints")):
			AlignementType	= 0
			pt_num			= OMPseg.refpoints.nPoints
			refineType		= OMPseg.refpoints.RefineType 
			refineSide		= OMPseg.refpoints.RefineSide 
			adjustMode		= OMPseg.refpoints.AdjustMode 
			paramL			= OMPseg.refpoints.ParamL 
			paramR			= OMPseg.refpoints.ParamR 
			ForceOrtho		= 0 
			isSet			= True

		if(isSet):
			self.NumPointsEntryValue.set(pt_num)
			self.RefinementTypeOptMenuValue.set(self.RefinementTypeOptions[refineType])
			self.RefinementSideOptMenuValue.set(self.RefinementSideOptions[refineSide])
			self.RefinementModeOptMenuValue.set(self.RefinementModeOptions[adjustMode])
			self.RefinementSizeLeftEntryValue.set(paramL)
			self.RefinementSizeRightEntryValue.set(paramR)
			self.ForceOrthoCkButValue.set(ForceOrtho)
			self.SkipEvents = False
			self.AlignementTypeOptMenuValue.set(self.AlignementTypeOptions[AlignementType])
		else:
			if(self.AlignementTypeOptMenuValue.get() == self.AlignementTypeOptions[1]): self.set_segment_alignment()
			else:																		self.sub_sel	= -1
			self.SkipEvents = False
				
		self.update_figure()
			
			
	def set_segment_alignment(self):
		MagZones = self.root.Config.MagZones
		MagPMegazones = self.root.Config.MagPMegazones
		dmin   = 1e10
	#	north
		pmz		= MagZones[self.zone_sel].pmz
		if(self.side_sel == 0):
			if((not hasattr(self.root.Config.MagZones[self.zone_sel],'subNorth')) or (len(MagZones[self.zone_sel].subNorth[0].R) == 0)):
				subdivide_selected_segment(self.root.Config, MagZones, self.zone_sel, self.side_sel) 
				
				if(not hasattr(self.root.Config.MagZones[self.zone_sel],'subNorth')):
					messagebox.showwarning("MeshGen", "Wall alignment is impossible!\nSegment does not intersect wall")
					self.sel_to_mesh	= False
					self.sel_mode		= False
					self.side_sel		= -1 
					return
				
			for n in range(2):
				d, i = find_closest_segment(self.Rpt, self.Zpt, MagZones[self.zone_sel].subNorth[n].R, MagZones[self.zone_sel].subNorth[n].Z)
				if(d < dmin):
					dmin			= d
					self.sub_sel	= n
						
			self.psi_wall_segment = MagZones[self.zone_sel].pB.coord[0]
			if(MagZones[self.zone_sel].subNorth[0].ismeshed or MagZones[self.zone_sel].subNorth[1].ismeshed):
				psimin	= MagPMegazones[pmz].align_psimin
				psimax	= MagPMegazones[pmz].align_psimax
			else:
				psimin = 0.99*self.psi_wall_segment
				psimax = 1.01*self.psi_wall_segment

			vScaleMin  = (psimin-self.psi_wall_segment)/(self.root.Config.MagZones[self.root.Config.MagPMegazones[pmz].list[0]].pA.coord[0]-self.psi_wall_segment)
			vScaleMax = (psimax-self.psi_wall_segment)/(self.root.Config.MagZones[self.root.Config.MagPMegazones[pmz].list[-1]].pB.coord[0]-self.psi_wall_segment)
					
	#	south
		elif(self.side_sel == 1):
			if((not hasattr(self.root.Config.MagZones[self.zone_sel],'subSouth')) or (len(MagZones[self.zone_sel].subSouth[0].R) == 0)):
				subdivide_selected_segment(self.root.Config, MagZones, self.zone_sel, self.side_sel) 
					
				if(not hasattr(self.root.Config.MagZones[self.zone_sel],'subSouth')):
					messagebox.showwarning("Mesh Generator", "Wall alignment is impossible!\nSegment does not intersect wall")
					self.sel_to_mesh	= False
					self.sel_mode		= False
					self.side_sel		= -1 
					return
				
			for n in range(2):
				d, i = find_closest_segment(self.Rpt, self.Zpt, MagZones[self.zone_sel].subSouth[n].R, MagZones[self.zone_sel].subSouth[n].Z)
				if(d < dmin):
					dmin			= d
					self.sub_sel	= n

			self.psi_wall_segment = MagZones[self.zone_sel].pD.coord[0]
			if(MagZones[self.zone_sel].subSouth[0].ismeshed or MagZones[self.zone_sel].subSouth[1].ismeshed):
				psimin	= MagPMegazones[pmz].align_psimin
				psimax	= MagPMegazones[pmz].align_psimax
			else:
				psimin	= 0.99*self.psi_wall_segment
				psimax  = 1.01*self.psi_wall_segment

			vScaleMin  = (psimin-self.psi_wall_segment)/(self.root.Config.MagZones[self.root.Config.MagPMegazones[pmz].list[0]].pD.coord[0]-self.psi_wall_segment)
			vScaleMax = (psimax-self.psi_wall_segment)/(self.root.Config.MagZones[self.root.Config.MagPMegazones[pmz].list[-1]].pC.coord[0]-self.psi_wall_segment)

		else:
			messagebox.showwarning("Mesh Generator", "Error in program, check output and report")
			print("\tset_segment_alignment: self.side_sel    =",self.side_sel)
			print("\t                       self.sub_sel     =",self.sub_sel)
			print("\t                       self.sel_to_mesh =",self.sel_to_mesh)
			print("\t                       self.sel_mode    =",self.sel_mode)
			return

		OldSkipEvents	= self.SkipEvents
		self.SkipEvents = True
		self.PsiMinEntryValue.set("{:.4f}".format(psimin))
		self.PsiMaxEntryValue.set("{:.4f}".format(psimax))
		self.PsiMinScale.set(vScaleMin)
		self.PsiMaxScale.set(vScaleMax)

		self.SkipEvents = OldSkipEvents


	def cmd_deselect_all(self):
		self.sel_mode = False
		self.side_sel = -1
		self.zone_sel = -1
		self.sub_sel  = -1
		self.update_figure()


	def cmd_refinement_type(self, a, b, c):
		if(self.RefinementTypeOptMenuValue.get() !=	self.RefinementTypeOptions[1]):
			self.RefinementSizeLeftEntry.state(["disabled"])
			self.RefinementSizeRightEntry.state(["disabled"])
		else:
			self.RefinementSizeLeftEntry.state(["!disabled"])
			self.RefinementSizeRightEntry.state(["!disabled"])

	def cmd_refinement_mode(self, a, b, c):
		if(self.RefinementModeOptMenuValue.get() ==	self.RefinementModeOptions[0]):
			self.RefinementSizeLeftLabelValue.set("Left")								#Relative
			self.RefinementSizeRightLabelValue.set("Right")
		else:
			self.RefinementSizeLeftLabelValue.set("Left (cm)")							#Absolute
			self.RefinementSizeRightLabelValue.set("Right (cm)")


	def cmd_move_flux_min_scale(self, slider):
		if(self.SkipEvents): return
		vScaleMin = self.PsiMinScale.get()
		if(self.side_sel == 0):
			psival	= self.root.Config.MagZones[self.zone_sel].pB.coord[0]
			pmz		= self.root.Config.MagZones[self.zone_sel].pmz
			nz		= self.root.Config.MagPMegazones[pmz].list[0]
			psimin	= self.root.Config.MagZones[nz].pA.coord[0]							#era pA 
			psi		= psival+vScaleMin*(psimin-psival)
			self.PsiMinEntryValue.set("{:.3f}".format(psi))

		if(self.side_sel == 1):
			psival		= self.root.Config.MagZones[self.zone_sel].pD.coord[0]
			pmz			= self.root.Config.MagZones[self.zone_sel].pmz
			nz			= self.root.Config.MagPMegazones[pmz].list[0]
			psimin		= self.root.Config.MagZones[nz].pD.coord[0]
			psi			= psival+vScaleMin*(psimin-psival)		
			self.PsiMinEntryValue.set("{:.3f}".format(psi))


	def cmd_move_flux_max_scale(self, slider):
		if(self.SkipEvents): return
		vScaleMax = self.PsiMaxScale.get()
		if(self.side_sel == 0):
			psival	= self.root.Config.MagZones[self.zone_sel].pB.coord[0]
			pmz		= self.root.Config.MagZones[self.zone_sel].pmz
			nz		= self.root.Config.MagPMegazones[pmz].list[-1]
			psimax	= self.root.Config.MagZones[nz].pB.coord[0]
			psi	= psival+vScaleMax*(psimax-psival)
			self.PsiMaxEntryValue.set("{:.4f}".format(psi))

		if(self.side_sel == 1):
			psival	= self.root.Config.MagZones[self.zone_sel].pD.coord[0]
			pmz		= self.root.Config.MagZones[self.zone_sel].pmz
			nz		= self.root.Config.MagPMegazones[pmz].list[-1]
			psimax	= self.root.Config.MagZones[nz].pC.coord[0]
			psi		= psival+vScaleMax*(psimax-psival)
			self.PsiMaxEntryValue.set("{:.4f}".format(psi))

	def cmd_psi_min_entry(self, a, b, c):
		if(self.SkipEvents): return
		
		if(safe_eval(self.PsiMinEntryValue.get()) > self.psi_wall_segment): self.PsiMinEntryValue.set("{:.4f}".format(self.psi_wall_segment))
		self.update_figure()

	def cmd_psi_max_entry(self, a, b, c):
		if(self.SkipEvents): return
		if(safe_eval(self.PsiMaxEntryValue.get()) < self.psi_wall_segment): self.PsiMaxEntryValue.set("{:.4f}".format(self.psi_wall_segment))
		self.update_figure()

	def cmd_mesh_segment(self):
		if(not self.sel_mode): return

		MagZones = self.root.Config.MagZones
		MagPMegazones = self.root.Config.MagPMegazones

		mesh_pt_num = self.NumPointsEntryValue.get()
		if(mesh_pt_num < 3):
			messagebox.showerror("MeshGen", "Minimum number of segment points is 3")
			return
			
		if(self.RefinementTypeOptMenuValue.get() ==	self.RefinementTypeOptions[0]):
			RefineType = 0															#None
		elif(self.RefinementTypeOptMenuValue.get() ==	self.RefinementTypeOptions[1]):	
			RefineType = 1															#Linear
		elif(self.RefinementTypeOptMenuValue.get() ==	self.RefinementTypeOptions[2]):	
			RefineType = 2															#Custom

		if(self.RefinementSideOptMenuValue.get() == self.RefinementSideOptions[0]):
			RefineSide = 0															#Left
		elif(self.RefinementSideOptMenuValue.get() == self.RefinementSideOptions[1]):
			RefineSide = 1															#Right
		elif(self.RefinementSideOptMenuValue.get() == self.RefinementSideOptions[2]):
			RefineSide = 2															#Both

		if(self.RefinementModeOptMenuValue.get() ==	self.RefinementModeOptions[0]): #relative
			modeAdjust = 0
		elif(self.RefinementModeOptMenuValue.get() == self.RefinementModeOptions[1]):	#absolute
			modeAdjust = 1
		
		if(self.AlignementTypeOptMenuValue.get() == self.AlignementTypeOptions[0]):
			Alignement = False
			psimin	   = 0.
			psimax	   = 0.
		elif(self.AlignementTypeOptMenuValue.get() == self.AlignementTypeOptions[1]):
			Alignement = True
			try:
				psimin	   = eval(self.PsiMinEntryValue.get())
				psimax	   = eval(self.PsiMaxEntryValue.get())
			except:
				pass

		paramL		= self.RefinementSizeLeftEntryValue.get()
		paramR		= self.RefinementSizeRightEntryValue.get()
		ForceOrtho	= self.ForceOrthoCkButValue.get()
		meshed_segment = mesh_segment(self.root, self.root.Config, self.zone_sel, self.side_sel, self.sub_sel,
									mesh_pt_num, RefineType, RefineSide, modeAdjust, paramL, paramR, ForceOrtho, Alignement, psimin, psimax)
		if(meshed_segment):
			self.sel_to_mesh	= False

			self.root.set_MagGrid_OK(self.root.Config, False)
			self.update_figure()
		

	def cmd_align_type(self, a, b, c):
		if(self.SkipEvents): return
		if(self.AlignementTypeOptMenuValue.get() == self.AlignementTypeOptions[0]):
			self.PsiMinScale.state(["disabled"])					#Allign with flux
			self.PsiMaxScale.state(["disabled"])
			self.PsiMinEntry.state(["disabled"])
			self.PsiMaxEntry.state(["disabled"])
			self.align_psimode = False
			self.sub_sel  = -1
		else:
			self.PsiMinEntry.state(["!disabled"])					#Allign with wall
			self.PsiMaxEntry.state(["!disabled"])	
			self.PsiMinScale.state(["!disabled"])
			self.PsiMaxScale.state(["!disabled"])	
			self.align_psimode = True
			if(self.sel_to_mesh):	self.set_segment_alignment()		

		self.update_figure()

#	Optimization commands
#	=====================

#	edit optimization scale

	def cmd_move_optimization_scale(self, slider):
		if(self.SkipEvents): return
		opti = self.OptimizationScaleValue.get()
		self.SkipEvents = True
		self.OptimizationEntryValue.set("{:.3f}".format(opti))
		self.SkipEvents = False
		for k in range(len(self.root.Config.MagZones)): self.root.Config.MagZones[k].orthomeshchanged = True

	def cmd_optimization_entry(self, a, b, c):
		if(self.SkipEvents): return
		opti = safe_eval(self.OptimizationEntryValue.get())
		self.SkipEvents = True
		if(opti < 0.):
			self.OptimizationEntryValue.set("{:.3f}".format(0.))	
			self.OptimizationScaleValue.set(0.)
		elif(opti > 1.):
			self.OptimizationEntryValue.set("{:.3f}".format(1.))	
			self.OptimizationScaleValue.set(1.)
		else:
			self.OptimizationScaleValue.set(opti)
		self.SkipEvents = False


#	edit optimal dmin

	def cmd_opt_dmin_entry(self, a, b, c):
		for k in range(len(self.root.Config.MagZones)): self.root.Config.MagZones[k].orthomeshchanged = True

#	edit optimal range

	def cmd_opt_range_entry(self, a, b, c):
		for k in range(len(self.root.Config.MagZones)): self.root.Config.MagZones[k].orthomeshchanged = True

#	edit optimal lenght

	def cmd_opt_lenght_entry(self, a, b, c):
		for k in range(len(self.root.Config.MagZones)): self.root.Config.MagZones[k].orthomeshchanged = True
	
#	edit optimal alpha

	def cmd_opt_alpha_entry(self, a, b, c):
		for k in range(len(self.root.Config.MagZones)): self.root.Config.MagZones[k].orthomeshchanged = True

#	edit optimal target

	def cmd_opt_target_entry(self, a, b, c):
		for k in range(len(self.root.Config.MagZones)): self.root.Config.MagZones[k].orthomeshchanged = True

#	edit optimal drmin

	def cmd_opt_drmin_entry(self, a, b, c):
		for k in range(len(self.root.Config.MagZones)): self.root.Config.MagZones[k].orthomeshchanged = True

#	Mesh commands
#	=============

	def cmd_mesh_all(self):
		self.root.Config.Optimization			= types.SimpleNamespace()
		self.root.Config.Optimization.orthoptim	= (1.-self.OptimizationScaleValue.get())
		self.root.Config.Optimization.dmin_ad	= self.OptdMinEntryValue.get()
		self.root.Config.Optimization.range_ad	= self.OptRangeEntryValue.get()
		self.root.Config.Optimization.length_ad	= self.OptLenghtEntryValue.get()
		self.root.Config.Optimization.alpha_ad	= self.OptAlphaEntryValue.get()
		self.root.Config.Optimization.target_ad	= self.OptTargetEntryValue.get()
		self.root.Config.Optimization.drmin_ad	= self.OptdrMinEntryValue.get()
		mesh_all_zones(self.root, self.root.Config, self.root.Config.Optimization)
		self.NumZonesLabelValue.set("N. zones: {:d}".format(len(self.root.Config.MagZones)))


	def cmd_reg_mesh(self):
		if(not self.root.Config.MagGrid_OK):return

		regularize_mesh_size(self.root.Config)
		conform_mesh(self.root.Config)

		self.root.set_Mesh_OK(self.root.Config, False)
		self.update_figure()
		

	def cmd_tweak_mesh(self):
		if(not self.root.Config.MagGrid_OK):return

		dmin_ad		= self.OptdMinEntryValue.get()
		drmin_ad	= self.OptdrMinEntryValue.get()
		tweak_mesh(self.root.Config, dmin_ad, drmin_ad)
		conform_mesh(self.root.Config)

		self.root.set_Mesh_OK(self.root.Config, False)
		self.update_figure()

	def cmd_tweak_mesh_ortho(self):
		if(not self.root.Config.MagGrid_OK):return

		dmin_ad		= self.OptdMinEntryValue.get()
		drmin_ad	= self.OptdrMinEntryValue.get()
		tweak_mesh_orho(self.root.Config, dmin_ad, drmin_ad)
		conform_mesh(self.root.Config)

		self.root.set_Mesh_OK(self.root.Config, False)
		self.update_figure()


	def cmd_soledge_mesh(self):
		if((self.root.Config.MagGrid_OK) and (not self.root.Config.Mesh_OK)):
			self.root.Config.use_penalization = self.PenalizationCkButValue.get() 
			generate_soledge_mesh(self.root, self.root.Config)

			if(not hasattr(self.root.Config, "Split")): self.root.Config.Split = []

			if(len(self.root.Config.Split) > 0):
				Split = self.root.Config.Split
				for k in range(Split.shape[1]):
					if(Split[0,k] == 0):
						ZonesToSplit = np.copy(self.root.Config.PMegazones[self.root.Config.Zones[Split[1,k]-1].pmz].list)
						vertical_split(self.root.Config.Zones, ZonesToSplit, Split[2,k])
					else:
						ZonesToSplit = np.copy(self.root.Config.Megazones[self.root.Config.Zones[Split[1,k]-1].mz].list)
						horizontal_split(self.root.Config.Zones, ZonesToSplit, Split[2,k]);

					recompute_megazones(self.root.Config)

			self.post_split()
			
#	Reverse penalization command
#	==================

	def cmd_reverse_penalization(self):
		if(not self.root.Config.Mesh_OK): return
		if((not self.InButton) and (not self.InReversePenalization)):
			if(self.ShowInFigureOptMenuValue.get() != self.ShowInFigureOptions[4]):			#not Cent. Mesh
				messagebox.showinfo("Mesh Gen", "Sorry you must select 'Cent. Mesh' view for this function")
				return
			self.ReversePenalizationButValue.set("Cancel")
			self.cid_mouse_button_press_reverse = self.Fig.canvas.mpl_connect('button_press_event', self.cmd_button_press_reverse_penalization)
			self.cursor = Cursor(self.Ax, horizOn=True, vertOn=True, useblit=False, color='red', linewidth=1)
			self.InButton = self.InReversePenalization = True
		elif(self.InReversePenalization):
			self.Fig.canvas.mpl_disconnect(self.cid_mouse_button_press_reverse)
			self.cid_mouse_button_press_limit = 0
			self.cursor=0
			self.ReversePenalizationButValue.set("Rev. Penaliz.")
			self.InButton = self.InReversePenalization = False


	def cmd_button_press_reverse_penalization(self, event):
		R = event.xdata
		Z = event.ydata
		Config = self.root.Config
		kZone  = get_zone_of_point(Config, R, Z)
		Zone   = Config.Zones[kZone]
		for i in range(Zone.Nx):
			for j in range(Zone.Nz):
				MeshPath = Path(np.array([[Zone.gridR[i,j], Zone.gridR[i,j+1], Zone.gridR[i+1,j+1], Zone.gridR[i+1,j], Zone.gridR[i,j]],
										  [Zone.gridZ[i,j], Zone.gridZ[i,j+1], Zone.gridZ[i+1,j+1], Zone.gridZ[i+1,j], Zone.gridZ[i,j]]]).T, closed=True)
				if(MeshPath.contains_point([R, Z])): 
					Zone.Chi[i,j] = 1 - Zone.Chi[i,j]
					self.update_figure()
					return


			
#	Split zones commands
#	====================

	def cmd_split_vertical(self):
		split_par = ask_to_split_zones(self.root)
		if(len(split_par) == 2):
			if(len(self.root.Config.Split)==0): self.root.Config.Split = np.array([[0],[split_par[0]],[split_par[1]]])
			else:								self.root.Config.Split = np.append(self.root.Config.Split,np.array([[0],[split_par[0]],[split_par[1]]]),axis=1)

			ZonesToSplit = np.copy(self.root.Config.PMegazones[self.root.Config.Zones[split_par[0]-1].pmz].list)
			vertical_split(self.root.Config.Zones, ZonesToSplit, split_par[1])
			recompute_megazones(self.root.Config)
			self.post_split()


	def cmd_split_horizontal(self):
		split_par = ask_to_split_zones(self.root)
		if(len(split_par) == 2):
			if(len(self.root.Config.Split)==0): self.root.Config.Split = np.array([[1],[split_par[0]],[split_par[1]]])
			else:								self.root.Config.Split = np.append(self.root.Config.Split,np.array([[1],[split_par[0]],[split_par[1]]]),axis=1)

			ZonesToSplit = np.copy(self.root.Config.Megazones[self.root.Config.Zones[split_par[0]-1].mz].list)
			horizontal_split(self.root.Config.Zones, ZonesToSplit, split_par[1]);
			recompute_megazones(self.root.Config)
			self.post_split()

	def cmd_split_reset(self):
		self.root.Config.Split = []
		self.root.set_Mesh_OK(self.root.Config, False)
		self.cmd_soledge_mesh()

	def post_split(self):
		self.root.set_Mesh_OK(self.root.Config, True)
		self.NumZonesLabelValue.set("N. zones: {:d}".format(len(self.root.Config.Zones)))
		if(self.ShowInFigureOptMenuValue.get() == self.ShowInFigureOptions[11]):
			self.update_figure()
		elif(self.ShowInFigureOptMenuValue.get() == self.ShowInFigureOptions[12]): 
			self.resize_mode = False
			self.update_figure()

#	show in figure command from option menu
#	=======================================

	def cmd_show_in_figure(self, a, b, c):

		if((self.ShowInFigureOptMenuValue.get() == self.ShowInFigureOptions[1]) or 
		   (self.ShowInFigureOptMenuValue.get() == self.ShowInFigureOptions[2]) or
		   (self.ShowInFigureOptMenuValue.get() == self.ShowInFigureOptions[3])):
			if(not self.root.Config.MagGrid_OK): return

		elif(self.ShowInFigureOptMenuValue.get() != self.ShowInFigureOptions[0]):
			if(self.root.Config.MagGrid_OK):
				if(not self.root.Config.Mesh_OK): self.cmd_soledge_mesh()
			else:
				return
			
		self.update_figure()


#===============
#	Update data 
#===============

	def update_data(self):
		self.SkipEvents = True
		if(self.root.Config.xPoints_OK and not self.sel_mode):
			self.PsiMinEntryValue.set("{:.4f}".format(self.root.Config.X_points[0].psi))
			self.PsiMaxEntryValue.set("{:.4f}".format(self.root.Config.X_points[0].psi))

		if(self.root.Config.Mesh_OK):
			self.NumZonesLabelValue.set("N. zones: {:d}".format(len(self.root.Config.Zones)))
		else:	
			self.NumZonesLabelValue.set("N. zones: 0")

		self.SkipEvents = False
		self.cmd_show_neighbors()
		self.update_figure()
		
		
#==============
#	Main figure 
#==============

#update main figure

	def update_figure(self):

		if(not  self.root.Config.MagZones_OK): return

		Config			= self.root.Config
		MagZones		= self.root.Config.MagZones
		MagMegazones	= self.root.Config.MagMegazones
		MagPMegazones	= self.root.Config.MagPMegazones
		Zones			= self.root.Config.Zones
		Megazones		= self.root.Config.Megazones
		Fig				= self.Fig
		Ax				= self.Ax

		nMagZones		= len(MagZones)
		nMagMegazones	= len(MagMegazones)
		nMagPMegazones	= len(MagPMegazones)
		nZones			= len(Zones)
		nMegazones		= len(Megazones)
		
		Xlims = Ax.get_xlim()
		Ylims = Ax.get_ylim()
		if(self.Cb != 0):
			Fig.clf()
			self.Ax = Ax = Fig.add_subplot(111)
			self.Cb = 0
		else:
			Ax.clear()
		
		Ax.tick_params(axis='both', which='major', labelsize=8)
		Ax.autoscale(axis='both', tight=True)
		Ax.axis('on')

		PlotWall = True
		Ax.set_aspect(1.)
		self.InfoLabelValue.set("")
		
		self.new_plot_type	= 0
		if((not Config.MagGrid_OK) or (self.ShowInFigureOptMenuValue.get() == self.ShowInFigureOptions[0])):				#Show zones borders and segements

			for k in range(nMagZones):
				sMid = floor(len(MagZones[k].south.R)/2)
				nMid = floor(len(MagZones[k].north.R)/2)
				R	 = (MagZones[k].south.R[sMid] + MagZones[k].north.R[nMid])*0.5
				Z	 = (MagZones[k].south.Z[sMid] + MagZones[k].north.Z[nMid])*0.5
				if(MagZones[k].south.ismeshed):
					if(MagZones[k].southaligned):
						Ax.plot(MagZones[k].south.R,MagZones[k].south.Z,'g-', linewidth=3)
					else:
						Ax.plot(MagZones[k].south.R,MagZones[k].south.Z,'g-', linewidth=1)
	
				else:
					if(MagZones[k].southaligned): 
						Ax.plot(MagZones[k].south.R,MagZones[k].south.Z,'r-', linewidth=3)
					else:
						Ax.plot(MagZones[k].south.R,MagZones[k].south.Z,'r-', linewidth=1)
	
				if(MagZones[k].north.ismeshed):
					if(MagZones[k].northaligned):
						Ax.plot(MagZones[k].north.R,MagZones[k].north.Z,'g-', linewidth=3)
					else:
						Ax.plot(MagZones[k].north.R,MagZones[k].north.Z,'g-', linewidth=1)
	
				else:
					if(MagZones[k].northaligned):
						Ax.plot(MagZones[k].north.R,MagZones[k].north.Z,'r-', linewidth=3)
					else:
						Ax.plot(MagZones[k].north.R,MagZones[k].north.Z,'r-', linewidth=1)
	
				if(MagZones[k].east.ismeshed):
					Ax.plot(MagZones[k].east.R,MagZones[k].east.Z,'g-', linewidth=1)
				else:
					Ax.plot(MagZones[k].east.R,MagZones[k].east.Z,'r-', linewidth=1)
	
				if(MagZones[k].west.ismeshed):
					Ax.plot(MagZones[k].west.R,MagZones[k].west.Z,'g-', linewidth=1)
				else:
					Ax.plot(MagZones[k].west.R,MagZones[k].west.Z,'r-', linewidth=1)

				Ax.text(R,Z,"{:d}".format(k+1))
	
			OMPseg	= self.root.Config.OMPseg
			if(OMPseg.ismeshed):													# Plot special OMP segment
				Ax.plot(OMPseg.R,OMPseg.Z,'g-', linewidth=1)
			else:
				Ax.plot(OMPseg.R,OMPseg.Z,'r-', linewidth=1)
	
			if(self.sel_mode):
				zone_sel = self.zone_sel
				R = []; Z = []
				if(zone_sel < len(MagZones)):										#It is a normal zone segment that is selected
					if(self.side_sel == 0):
						if(self.sub_sel < 0):
							R = MagZones[zone_sel].north.R
							Z = MagZones[zone_sel].north.Z
						else:
							try:
								psimin = eval(self.PsiMinEntryValue.get())
								psimax = eval(self.PsiMaxEntryValue.get())
								Ax.contour(Config.r2D, Config.z2D, Config.flux2D, levels = [psimin], colors ='m', linestyles='dashed')
								Ax.contour(Config.r2D, Config.z2D, Config.flux2D, levels = [psimax], colors ='y', linestyles='dashed')
								R = MagZones[zone_sel].subNorth[self.sub_sel].R
								Z = MagZones[zone_sel].subNorth[self.sub_sel].Z
								if(MagZones[zone_sel].north.ismeshed):
									Ax.plot(MagZones[zone_sel].north.R,MagZones[zone_sel].north.Z,'g-', linewidth=3)
								else:
									Ax.plot(MagZones[zone_sel].north.R,MagZones[zone_sel].north.Z,'r-', linewidth=3)
							except:
								pass
		
					elif(self.side_sel == 1):
						if(self.sub_sel < 0):
							R = MagZones[zone_sel].south.R
							Z = MagZones[zone_sel].south.Z						
						else:
							try:
								psimin = eval(self.PsiMinEntryValue.get())
								psimax = eval(self.PsiMaxEntryValue.get())
								Ax.contour(Config.r2D,Config.z2D,Config.flux2D, levels = [psimin], colors ='m', linestyles='dashed')
								Ax.contour(Config.r2D,Config.z2D,Config.flux2D, levels = [psimax], colors ='y', linestyles='dashed')
								R = MagZones[zone_sel].subSouth[self.sub_sel].R
								Z = MagZones[zone_sel].subSouth[self.sub_sel].Z
								if(MagZones[zone_sel].south.ismeshed):
									Ax.plot(MagZones[zone_sel].south.R,MagZones[zone_sel].south.Z,'g-', linewidth=3)
								else:
									Ax.plot(MagZones[zone_sel].south.R,MagZones[zone_sel].south.Z,'r-', linewidth=3)
							except:
								pass
								
					elif(self.side_sel == 2):
						R = MagZones[zone_sel].east.R
						Z = MagZones[zone_sel].east.Z
							
					elif(self.side_sel == 3):
						R = MagZones[zone_sel].west.R
						Z = MagZones[zone_sel].west.Z
				elif(zone_sel == len(MagZones)):
					if(OMPseg.ismeshed):											# OMP segment
						Ax.plot(OMPseg.R,OMPseg.Z,'g-', linewidth=3)
					else:
						Ax.plot(OMPseg.R,OMPseg.Z,'r-', linewidth=3)
				if(len(R) > 0):														# for sub-segments overplot dashed line and Left/Right
					iLeft  = 10
					iRight = len(R)-11
					if(iLeft > iRight):
						iLeft   = 0
						iRight  = len(R)-1				
	
					Ax.plot(R, Z, 'b--', linewidth=2)
					Ax.text(R[iLeft],  Z[iLeft],  "L")
					Ax.text(R[iRight], Z[iRight], "R")
	
			for k in range(nMagMegazones):											# plot ref points
				if(MagMegazones[k].ismeshed):
					Ax.plot(MagMegazones[k].refpoints.R,MagMegazones[k].refpoints.Z,'k.')
					
			for k in range(nMagPMegazones):
				if(hasattr(MagPMegazones[k],'refpoints')):
					if(len(MagPMegazones[k].refpoints.R) > 0):
						Ax.plot(MagPMegazones[k].refpoints.R,MagPMegazones[k].refpoints.Z,'k.')

		elif(self.ShowInFigureOptMenuValue.get() == self.ShowInFigureOptions[1]):				#Show grid

			EdgeLineType	='b-'
			EdgeLineWidth	= 1
			for k in range(nMagZones):
				if(MagZones[k].meshortho):
					Ax.plot(MagZones[k].gridR,   MagZones[k].gridZ,   'g-')
					Ax.plot(MagZones[k].gridR.T, MagZones[k].gridZ.T, 'g-')
				else:
					Ax.plot(MagZones[k].gridR,   MagZones[k].gridZ,   '-', color=[1.,0.8,0.2])
					Ax.plot(MagZones[k].gridR.T, MagZones[k].gridZ.T, '-', color=[1.,0.8,0.2])

				Ax.plot(MagZones[k].gridR[0,:],  MagZones[k].gridZ[0,:],  EdgeLineType, linewidth=EdgeLineWidth)		#south
				Ax.plot(MagZones[k].gridR[:,-1], MagZones[k].gridZ[:,-1], EdgeLineType, linewidth=EdgeLineWidth)		#East
				Ax.plot(MagZones[k].gridR[-1,:], MagZones[k].gridZ[-1,:], EdgeLineType, linewidth=EdgeLineWidth)		#North
				Ax.plot(MagZones[k].gridR[:,0],  MagZones[k].gridZ[:,0],  EdgeLineType, linewidth=EdgeLineWidth)		#West

			self.InfoLabelValue.set("")
			
		elif(self.ShowInFigureOptMenuValue.get() == self.ShowInFigureOptions[2]):			#Show quality

			QualMin = np.empty(nMagZones, dtype='f8')
			QualMax = np.empty(nMagZones, dtype='f8')
			
			for k in range(nMagZones):
				m = MagZones[k].gridR.shape[0]
				if(m != 0):
					vx		= MagZones[k].gridR[0,1]-MagZones[k].gridR[0,0]
					vy		= MagZones[k].gridZ[0,1]-MagZones[k].gridZ[0,0]
					vx1		= MagZones[k].gridR[1,0]-MagZones[k].gridR[0,0]
					vy1		= MagZones[k].gridZ[1,0]-MagZones[k].gridZ[0,0]
					qual	= (vx1*vy-vy1*vx)/(np.sqrt(vx**2+vy**2)*np.sqrt(vx1**2+vy1**2))
					sens	= np.sign(qual)
				else:
					sens	= 1

				qual	= np.zeros(MagZones[k].gridR.shape)
				vx		= MagZones[k].gridR[:-1,1:] - MagZones[k].gridR[:-1,:-1]
				vy		= MagZones[k].gridZ[:-1,1:] - MagZones[k].gridZ[:-1,:-1]
				vx1		= MagZones[k].gridR[1:,:-1] - MagZones[k].gridR[:-1,:-1]
				vy1		= MagZones[k].gridZ[1:,:-1] - MagZones[k].gridZ[:-1,:-1]
				qual	= (vx1*vy-vy1*vx)/(np.sqrt(vx**2+vy**2)*np.sqrt(vx1**2+vy1**2))*sens
				QualMin = np.min(qual)
				QualMax = np.max(qual)
				MagZones[k].qual = qual

			QualMin = np.min(QualMin)	
			QualMax = np.max(QualMax)	
			for k in range(nMagZones):					
				qual = np.append(MagZones[k].qual, MagZones[k].qual[:,-1].reshape(MagZones[k].qual.shape[0],1), axis=1)
				qual = np.append(qual, qual[-1,:].reshape(1,qual.shape[1]), axis=0)
				inQual = Ax.pcolormesh(MagZones[k].gridR, MagZones[k].gridZ, qual, shading='gouraud', vmin=QualMin, vmax=QualMax, cmap=plt.get_cmap('RdYlGn'))
				MagZones[k].qual = 0
				
			self.Cb = Fig.colorbar(inQual, ax=Ax)
			self.Cb.ax.tick_params(labelsize=8)
			self.InfoLabelValue.set("Info: Minimun quality factor = {:.3f}".format(QualMin))

		elif(self.ShowInFigureOptMenuValue.get() == self.ShowInFigureOptions[3]):			#Show mesh size

			QualMin = np.empty(nMagZones, dtype='f8')
			QualMax = np.empty(nMagZones, dtype='f8')
			for k in range(nMagZones):
				vx		= MagZones[k].gridR[:-1,1:] - MagZones[k].gridR[:-1,:-1]
				vy		= MagZones[k].gridZ[:-1,1:] - MagZones[k].gridZ[:-1,:-1]
				vx1		= MagZones[k].gridR[1:,:-1] - MagZones[k].gridR[:-1,:-1]
				vy1		= MagZones[k].gridZ[1:,:-1] - MagZones[k].gridZ[:-1,:-1]
				qualr1	= np.abs((vx1*vy-vy1*vx)/(np.sqrt(vx**2+vy**2)))*1000
				qualz1	= np.abs((vx1*vy-vy1*vx)/(np.sqrt(vx1**2+vy1**2)))*1000
				qual	= np.where(qualr1 < qualz1, qualr1,qualz1 )
				QualMin = np.min(qual)
				QualMax = np.max(qual)
				MagZones[k].qual = qual
				
			QualMin = np.min(QualMin)	
			QualMax = np.max(QualMax)	
			for k in range(nMagZones):					
				qual = np.append(MagZones[k].qual, MagZones[k].qual[:,-1].reshape(MagZones[k].qual.shape[0],1), axis=1)
				qual = np.append(qual, qual[-1,:].reshape(1,qual.shape[1]), axis=0)
				inQual = Ax.pcolormesh(MagZones[k].gridR, MagZones[k].gridZ, qual, shading='gouraud', vmin=QualMin, vmax=QualMax, cmap=plt.get_cmap('hot'))
				MagZones[k].qual = 0
				
			self.Cb = Fig.colorbar(inQual, ax=Ax)
			self.Cb.ax.tick_params(labelsize=8)
			self.InfoLabelValue.set("Info: Minimun distance = {:.3f} mm.".format(QualMin))
			
		elif(self.ShowInFigureOptMenuValue.get() == self.ShowInFigureOptions[4]):			#Cent. Mesh

			for k in range(nZones):
				Ax.plot(Zones[k].gridR,   Zones[k].gridZ,   'b-')
				Ax.plot(Zones[k].gridR.T, Zones[k].gridZ.T, 'b-')
				iIn, jIn = np.where(Zones[k].Chi == 0)
				Ax.plot(Zones[k].gridRc[iIn, jIn],  Zones[k].gridZc[iIn, jIn],  'r.')
				iOut, jOut = np.where(Zones[k].Chi == 1)
				Ax.plot(Zones[k].gridRc[iOut, jOut],  Zones[k].gridZc[iOut, jOut],  'y.')
				iIn = 0; jIn = 0; iOut = 0; jOut = 0

		elif(self.ShowInFigureOptMenuValue.get() == self.ShowInFigureOptions[5]):			#Chi

			ChiMin = np.empty(nZones, dtype='f8')
			ChiMax = np.empty(nZones, dtype='f8')
			for k in range(nZones):
				ChiMin[k] = Zones[k].Chi.min()
				ChiMax[k] = Zones[k].Chi.max()
				
			ChiMin = ChiMin.min()
			ChiMax = ChiMax.max()
				
			for k in range(nZones):
				Chi = np.append(Zones[k].Chi, Zones[k].Chi[:,-1].reshape(Zones[k].Chi.shape[0],1), axis=1)
				Chi = np.append(Chi, Chi[-1,:].reshape(1,Chi.shape[1]), axis=0)
				InChi = Ax.pcolormesh(Zones[k].gridR, Zones[k].gridZ, Chi, vmin=ChiMin, vmax=ChiMax)

			Chi  = 0
			self.Cb = Fig.colorbar(InChi, ax=Ax)
			self.Cb.ax.tick_params(labelsize=8)

		elif(self.ShowInFigureOptMenuValue.get() == self.ShowInFigureOptions[6]):			#Br

			BrMin = np.empty(nZones, dtype='f8')
			BrMax = np.empty(nZones, dtype='f8')
			for k in range(nZones):
				BrMin[k] = np.min(Zones[k].Br)
				BrMax[k] = np.max(Zones[k].Br)
				
			BrMin = np.min(BrMin)
			BrMax = np.max(BrMax)

			for k in range(nZones):
				Br = np.append(Zones[k].Br, Zones[k].Br[:,-1].reshape(Zones[k].Br.shape[0],1), axis=1)
				Br = np.append(Br, Br[-1,:].reshape(1,Br.shape[1]), axis=0)
				InBr = Ax.pcolormesh(Zones[k].gridR, Zones[k].gridZ, Br, shading='gouraud', vmin=BrMin, vmax=BrMax)

			Br = 0
			self.Cb = Fig.colorbar(InBr, ax=Ax)
			self.Cb.ax.tick_params(labelsize=8)

		elif(self.ShowInFigureOptMenuValue.get() == self.ShowInFigureOptions[7]):			#Bz

			BzMin = np.empty(nZones, dtype='f8')
			BzMax = np.empty(nZones, dtype='f8')
			for k in range(nZones):
				BzMin[k] = np.min(Zones[k].Bz)
				BzMax[k] = np.max(Zones[k].Bz)
				
			BzMin = np.min(BzMin)
			BzMax = np.max(BzMax)

			for k in range(nZones):
				Bz = np.append(Zones[k].Bz, Zones[k].Bz[:,-1].reshape(Zones[k].Bz.shape[0],1), axis=1)
				Bz = np.append(Bz, Bz[-1,:].reshape(1,Bz.shape[1]), axis=0)
				InBz = Ax.pcolormesh(Zones[k].gridR, Zones[k].gridZ, Bz, shading='gouraud', vmin=BzMin, vmax=BzMax)

			Bz = 0
			self.Cb = Fig.colorbar(InBz, ax=Ax)
			self.Cb.ax.tick_params(labelsize=8)

		elif(self.ShowInFigureOptMenuValue.get() == self.ShowInFigureOptions[8]):			#Bphi

			BphiMin = np.empty(nZones, dtype='f8')
			BphiMax = np.empty(nZones, dtype='f8')
			for k in range(nZones):
				BphiMin[k] = np.min(Zones[k].Bphi)
				BphiMax[k] = np.max(Zones[k].Bphi)
				
			BphiMin = np.min(BphiMin)
			BphiMax = np.max(BphiMax)

			for k in range(nZones):
				Bphi = np.append(Zones[k].Bphi, Zones[k].Bphi[:,-1].reshape(Zones[k].Bphi.shape[0],1), axis=1)
				Bphi = np.append(Bphi, Bphi[-1,:].reshape(1,Bphi.shape[1]), axis=0)
				InBphi = Ax.pcolormesh(Zones[k].gridR, Zones[k].gridZ, Bphi, shading='gouraud', vmin=BphiMin, vmax=BphiMax)

			Bphi = 0
			self.Cb = Fig.colorbar(InBphi, ax=Ax)
			self.Cb.ax.tick_params(labelsize=8)
				
		elif(self.ShowInFigureOptMenuValue.get() == self.ShowInFigureOptions[9]):			#x mesh

			xcMin = np.empty(nZones, dtype='f8')
			xcMax = np.empty(nZones, dtype='f8')
			for k in range(nZones):
				xcMin[k] = Zones[k].x.min()
				xcMax[k] = Zones[k].x.max()
				
			xcMin = xcMin.min()
			xcMax = xcMax.max()

			for k in range(nZones):
				xc = np.append(Zones[k].x, Zones[k].x[:,-1].reshape(Zones[k].x.shape[0],1), axis=1)
				xc = np.append(xc, xc[-1,:].reshape(1,xc.shape[1]), axis=0)
				Inxc = Ax.pcolormesh(Zones[k].gridR, Zones[k].gridZ, xc, shading='gouraud', vmin=xcMin, vmax=xcMax)

			xc = 0
			self.Cb = Fig.colorbar(Inxc, ax=Ax)
			self.Cb.ax.tick_params(labelsize=8)

		elif(self.ShowInFigureOptMenuValue.get() == self.ShowInFigureOptions[10]):			#z mesh

			zcMin = np.empty(nZones, dtype='f8')
			zcMax = np.empty(nZones, dtype='f8')
			for k in range(nZones):
				zcMin[k] = Zones[k].z.min()
				zcMax[k] = Zones[k].z.max()
				
			zcMin = zcMin.min()
			zcMax = zcMax.max()

			for k in range(nZones):
				zc = np.append(Zones[k].z, Zones[k].z[:,-1].reshape(Zones[k].z.shape[0],1), axis=1)
				zc = np.append(zc, zc[-1,:].reshape(1,zc.shape[1]), axis=0)
				Inzc = Ax.pcolormesh(Zones[k].gridR, Zones[k].gridZ, zc, shading='gouraud', vmin=zcMin, vmax=zcMax)

			zc = 0
			self.Cb = Fig.colorbar(Inzc, ax=Ax)
			self.Cb.ax.tick_params(labelsize=8)

		elif(self.ShowInFigureOptMenuValue.get() == self.ShowInFigureOptions[11]):			#Zones sizes
			self.new_plot_type	= 1
			for k in range(nZones):
				zmin	= Zones[k].zb[0,0]
				zmax	= Zones[k].zb[0,-1]
				xmin	= Zones[k].xb[0,0]
				xmax	= Zones[k].xb[-1,0]
				Ax.plot([zmin, zmax],[xmin, xmin],'b-')
				Ax.plot([zmin, zmax],[xmax, xmax],'b-')
				Ax.plot([zmin, zmin],[xmin, xmax],'b-')
				Ax.plot([zmax, zmax],[xmin, xmax],'b-')
				zc = (zmin+zmax)/2.
				xc = (xmin+xmax)/2.
#				print("k+1, zmin, zmax=",k+1,zmin, zmax)
#				print("k+1, xmin, xmax=",k+1,xmin, xmax)
#				print("k+1, zc,   xc  =",k+1,zc,xc)
				Ax.text(zc,xc,"{:d}={:d}x{:d}".format(k+1,Zones[k].Nz,Zones[k].Nx), horizontalalignment ="center", verticalalignment ="center", size = "xx-small")
			Ax.autoscale(axis='both', tight=False)
			Ax.set_aspect("auto")
			Ax.axis('off')
			PlotWall = False

		elif(self.ShowInFigureOptMenuValue.get() == self.ShowInFigureOptions[12]):			#Zones bar plot
			self.new_plot_type	= 2

			Sizes  = np.zeros(nZones, dtype='f8')
			kSizes = np.arange(nZones, dtype='f8') + 1.5
			for k in range(nZones): Sizes[k] = Zones[k].Nx*Zones[k].Nz
			
			Ax.bar(kSizes, Sizes, 1)
			for k in range(nZones):		
				Ax.text(k+1.,Sizes[k]+1.,"{:d}=({:d},{:d})".format(k+1,Zones[k].Nz,Zones[k].Nx), rotation = "vertical", horizontalalignment ="center", verticalalignment ="bottom", size = "xx-small")

			Ax.set_aspect("auto")
			Ax.autoscale(enable=True, axis='x', tight=True)
			Ax.autoscale(enable=True, axis='y', tight=False)

#			self.InfoLabelValue.set("Min={:d}, Max={:d}, Mean={:d}, STDV={:d}".format(Sizes.min(), Sizes.max(), floor(np.mean(Sizes)+0.5), floor(sqrt(np.mean(Sizes**2)-np.mean(Sizes)**2))))
			self.InfoLabelValue.set("Min={:d}, Max={:d}, Mean={:d},".format(floor(Sizes.min()), floor(Sizes.max()), floor(np.mean(Sizes)+0.5)))
			PlotWall = False
			self.resize_mode = True

		elif(self.ShowInFigureOptMenuValue.get() == self.ShowInFigureOptions[13]):			#Splitted sizes
			for k in range(len(Zones)):
				EdgeLineType	='b-'
				EdgeLineWidth	= 1
				Ax.plot(Zones[k].gridR[0,:],  Zones[k].gridZ[0,:],  EdgeLineType, linewidth=EdgeLineWidth)		#south
				Ax.plot(Zones[k].gridR[:,-1], Zones[k].gridZ[:,-1], EdgeLineType, linewidth=EdgeLineWidth)		#East
				Ax.plot(Zones[k].gridR[-1,:], Zones[k].gridZ[-1,:], EdgeLineType, linewidth=EdgeLineWidth)		#North
				Ax.plot(Zones[k].gridR[:,0],  Zones[k].gridZ[:,0],  EdgeLineType, linewidth=EdgeLineWidth)		#West
	
				iCenter = floor(Zones[k].gridR.shape[1]/2)
				R=(Zones[k].gridR[0,iCenter]+Zones[k].gridR[-1,iCenter])*0.5
				Z=(Zones[k].gridZ[0,iCenter]+Zones[k].gridZ[-1,iCenter])*0.5
				Ax.text(R,Z,"{:d}".format(k+1), color='green', size = "x-small")

		elif(self.ShowInFigureOptMenuValue.get() == self.ShowInFigureOptions[14]):			#Mesh indexes

			for k in range(nZones):
				Ax.plot(Zones[k].gridR,   Zones[k].gridZ,   'b-')
				Ax.plot(Zones[k].gridR.T, Zones[k].gridZ.T, 'b-')
				iIn, jIn = np.where(Zones[k].Chi == 0)
				if(len(iIn) > 0):
					for i in range(len(iIn)):
						Ax.text(Zones[k].gridRc[iIn[i], jIn[i]],  Zones[k].gridZc[iIn[i], jIn[i]],"({:d},{:d},{:d})".format(k+1,iIn[i]+1, jIn[i]+1), color='red', size = "xx-small", horizontalalignment = "center", verticalalignment = "center")

				iOut, jOut = np.where(Zones[k].Chi == 1)
				if(len(iOut) > 0):
					for i in range(len(iOut)):
						Ax.text(Zones[k].gridRc[iOut[i], jOut[i]],  Zones[k].gridZc[iOut[i], jOut[i]],"({:d},{:d},{:d})".format(k+1,iOut[i]+1, jOut[i]+1), color='yellow', size = "xx-small", horizontalalignment = "center", verticalalignment = "center")

				iIn = 0; jIn = 0; iOut = 0; jOut = 0
		else:
			print("ERROR invalid option in update_figure()in mesh_gen.py check code ")
			print("\tself.ShowInFigureOptMenuValue.get()=",self.ShowInFigureOptMenuValue.get())
			print("\tself.ShowInFigureOptions=",self.ShowInFigureOptions)

#			if(self.align_psimode):
#				psimin = eval(self.PsiMinEntryValue.get())
#				psimax = eval(self.PsiMaxEntryValue.get())

#				Ax.contour(Config.r2D, Config.z2D, Config.flux2D, levels = [psimin, psimax], colors ='m', linestyles='dashed')
#				Ax.contour(Config.r2D, Config.z2D, Config.flux2D, levels = [psimin, psimax], colors ='m', linestyles='dashed')
		if(Config.wall_OK and PlotWall): plot2d_walls(Ax, Config.Walls, self.root.ShowPlasmaWallCheckValue.get(), self.root.ShowEireneWallCheckValue.get(), self.root.ShowExtraWallCheckValue.get())

		if(self.old_plot_type != self.new_plot_type): self.resize_mode = True
		self.old_plot_type = self.new_plot_type
		
		if((not ((Xlims[0] == 0.) and (Xlims[1] == 1.))) and (not self.resize_mode)):
			Ax.set_xlim(Xlims[0], Xlims[1])
			Ax.set_ylim(Ylims[0], Ylims[1])
		else:
			self.resize_mode = False
	
		Fig.canvas.draw()
