import types

def reset_config():
	Config = types.SimpleNamespace()

	Config.MeshFile				= False
	Config.in_equ_OK			= False
	Config.extrapol_OK			= False
	Config.wall_OK				= False
	Config.Frontiers_OK			= False
	Config.new_equ_OK			= False
	Config.equ_OK				= False
	Config.xPoints_OK			= False
	Config.MagZones_OK			= False
	Config.Segments_OK			= False
	Config.MagGrid_OK			= False
	Config.Mesh_OK				= False
	Config.transp_prof_OK		= False
	Config.transp_values_OK 	= False
	Config.feedback_transp_OK 	= False	
	Config.feedback_puffing_OK 	= False	
	Config.Eirene_OK			= False
	Config.CustomPlots_OK		= False

	Config.show_ortho		= False
	
	Config.X_points			= []
	Config.MagZones			= []
	Config.MagMegazones		= []
	Config.MagPMegazones	= []
	Config.Splits			= []
	Config.Zones			= []
	Config.Megazones		= []
	Config.Frontiers		= []
	Config.Eirene			= []
	Config.CustomPlots		= []
	Config.FeedbackTransp	= 0
	Config.FeedbackPuffing	= 0
	Optimization = types.SimpleNamespace()
	Optimization.orthoptim	= 0.5
	Optimization.dmin_ad	= 2e-3
	Optimization.range_ad	= 0.2
	Optimization.length_ad	= 1000.
	Optimization.alpha_ad	= 0.1
	Optimization.target_ad	= 0.75
	Optimization.drmin_ad	= 0.1
	Config.Optimization = Optimization


	return Config