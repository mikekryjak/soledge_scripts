BC_UKNOWN   = -1
BC_TRIANGLE = 0
BC_WALL		= 3
BC_CORE		= 1

CORE_NEIGHBOUR		= -2
NS_BORDER_NEIGHBOUR = -3
WALL_NEIGHBOUR		= -5
EW_BORDER_NEIGHBOUR = -6

DEBUG		= 2				#0 NO Message, 1= Few message, 2=many messages, 3=plot also

CIRCLE_RADIUS = -0.01		#Circle radius in surrounding X points (if negative factor of R) 

TRANS_MESH_DMAX = 0.5		#Max distance in searching triangles (typ. = 0.1-1.0)	

BALLOONING_TAB_NAMES	= ["Chi_e ", "D    ", "Nu   ", "Chi_i ", "Vpinch"] 
BALLOONING_NAMES		= ["Chie", "D", "Nu", "Chi", "Vpinch"] 					#Balloning names in mesh file
#BALLOONING_TAB_NAMES	= ["Chi_e ", "D    ", "Nu   ", "Chi_i "] 
#BALLOONING_NAMES		= ["Chie", "D", "Nu", "Chi"] 					#Balloning names in mesh file
DEF_BALLOONING_MODE		= 1													#deafult balloning mode 1/Bt

ZONE_MODES				= ("Single", "Connected", "Non Valid", "Same x")
DEF_ZONE_MODE			= 1													#deafult Connected mode
ZONE_MODE_SINGLE		= 0
ZONE_MODE_CONNECTED		= 1
ZONE_MODE_SAME_FLUX		= 2
ZONE_MODE_SAME_X		= 3

FEEDBACK_TAB_NAMES		= ["Ne  ", "Ni  ", "Te  ", "Ti  ","D   ", "Chi "]
FEEDBACK_NAMES			= ["Dense",   "Densi",    "Tempe",   "Tempi",   "D",    "Chi"]		#Feedback names

LENGTH_TYPES			= ("DIST", "LENGTH", "XM", "PSI_POL", "RHO_POL", "RHO_V")

CUSTOM_PLOTS_TYPES		= ("Parallel", "Perpendicular", "Temporal")

MAX_N_XPOINTS			= 6													#Max number of X points to find

KB						= 1.3806e-23										#boltzamn constant
EV						= 1.6022e-19										#Enegy associated to 1 eV in Joules proton charge
HPLANCK					= 6.626070040e-34									#Planck constant J.s
CSPEED					= 299792458e0										#speed of light
PI_TO_DEG				= 57.29577951308232

STEP_FOLLOW_GRAD		= 1e-3												#step in following gradient (usually 1e-3)

ASK_XPOINTS				= 0													#ask x points 0= no ask, 1= ask

MIN_AREA_TRIANGLES		= 1e-8												#Minimum area triangles
MIN_ORDER_TRIANGLES		= 1e-8												#Minimum order triangles

MIN_AUTO_VALUE			= -1e40
MAX_AUTO_VALUE			=  1e40

DRAWING_WALL			= 0
EXTERNAL_PLASMA_WALL	= 1													#Wall external to plasma domain
INTERNAL_PLASMA_WALL	= -1												#Wall internal to plasma domain
EXTERNAL_EIRENE_WALL	= 2
INTERNAL_EIRENE_WALL	= -2

WALL_TYPE_NAMES = ("External plasma boundary", "Internal plasma boundary", "External eirene boundary", "Internal eirene boundary", "Extra")
WALL_TYPE_VALUES  = (1,-1,2,-2,0)
WALL_TYPE_INDEX   = (4,0,2,3,1) 

FILE_TYPE_SOLEDGE2D		= 0
FILE_TYPE_SOLEDGE3X		= 1

#		https://xkcd.com/color/rgb/
#		https://matplotlib.org/3.1.1/gallery/color/named_colors.html

LINE_COLORS = ("Black", "Grey", "Brown", "Orange", "Gold", "Yellow", "Olive", "Green","Cyan" , "Blue", "Red", "Purple","Pink")
LINE_STYLES = ("solid", "dashed", "dashdot" "dotted", "loosely dotted", "densely dotted", "loosely dashed", "densely dashed", "dashdotted", "loosely dashdotted", "densely dashdotted", "dashdotdotted", "loosely dashdotdotted", "densely dashdotdotted")
LINE_WIDTH = ("0.5", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5", "4.0", "4.5", "5.0")

