
# Function definition is here

import h5py
import numpy					as np
from matplotlib.path 			import Path
import scipy.interpolate		as interp
from mesh.get_rz_core_sep		import get_rz_core_sep
from routines.h5_routines		import h5_read
from routines.globals			import DEBUG

#==============================================================================
# T
#==============================================================================

def get_rho_in_out_core_sep(Config, pts_r, pts_z, RZcore=None, RZsep=None, rho_type = "PSI_POL"):

	if(DEBUG > 1): print("get_rho_in_out_core_sep")

	RZcore, RZsep, CoreMegazone, SepMegazone, jSep = get_rz_core_sep(Config)
	f_psi		= interp.RectBivariateSpline(Config.r2D[0,:], Config.z2D[:,0], Config.flux2D.T)				#RectBivariateSpline want [x,y] order

	psi_sep		= np.mean(f_psi.ev(RZsep[:,0], RZsep[:,1]))
	nSep = RZsep.shape[0]
	if((RZsep[0,0] == RZsep[-1,0]) and (RZsep[0,1] == RZsep[-1,1])):
		nSepPath = nSep
	else:
		nSepPath = nSep + 1

	if((RZsep[0,0] == RZsep[-1,0]) and (RZsep[0,1] == RZsep[-1,1])):
		SepVerts = RZsep
	else:
		SepVerts = np.empty((nSepPath,2), dtype='f8')
		SepVerts[:-1,0] = RZsep[:,0]
		SepVerts[:-1,1] = RZsep[:,1]
		SepVerts[-1,:]  = SepVerts[0,:]

	SepCodes		= np.empty(nSepPath, dtype='i')
	SepCodes[1:-1]	= Path.LINETO
	SepCodes[0]		= Path.MOVETO
	SepCodes[-1]	= Path.CLOSEPOLY
	SepPath			= Path(SepVerts, SepCodes)

	psi1		= Config.flux2D.reshape(Config.flux2D.size)

	psi_rz		= np.empty((Config.flux2D.size,2), dtype='f8')
	psi_rz[:,0] = Config.r2D.reshape(Config.r2D.size)
	psi_rz[:,1] = Config.z2D.reshape(Config.r2D.size)
	psi_in		= np.where(SepPath.contains_points(psi_rz)); psi_in = psi_in[0]
	psi_axis	= np.min(psi1[psi_in])
	psi_rz		= 0
	psi1		= 0
	psi_in		= 0

	pts_rz = np.empty((len(pts_r),2), dtype='f8')
	pts_rz[:,0] = pts_r
	pts_rz[:,1] = pts_z

	pts_inout = np.where(SepPath.contains_points(pts_rz),1,0)
	pts_in	  = np.where(pts_inout == 1); pts_in  = pts_in[0]
	pts_out   = np.where(pts_inout == 0); pts_out = pts_out[0]
	pts_rz	  = 0
	pts_inout = 0

	if(rho_type.upper() == "PSI_POL"):	pts_rho = np.abs((f_psi.ev(pts_r, pts_z)-psi_axis)/(psi_sep-psi_axis))				#Normalized poloidal flux does not need sqrt
	else:								pts_rho = np.sqrt(np.abs((f_psi.ev(pts_r, pts_z)-psi_axis)/(psi_sep-psi_axis)))		#rho of poloidal flux needs sqrt

	if(DEBUG > 1): print("get_rho_in_out_core_sep: completed")

	return pts_rho, pts_in, pts_out, RZcore, RZsep



