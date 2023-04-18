import numpy 					as np
from matplotlib.path 			import Path

def get_in_out(ContourR, ContourZ, PtsR, PtsZ):

	PtsRZ	  = np.array([PtsR, PtsZ]).T
	if((ContourR[0] == ContourR[-1]) and (ContourZ[0] == ContourZ[-1])):
		ContourRZ = np.array([ContourR, ContourZ]).T
	else:
		ContourRZ = np.array([np.append(ContourR,ContourR[0]), np.append(ContourZ,ContourZ[0])]).T

	SepPath			= Path(ContourRZ, closed=True)

	pts_inout = np.where(SepPath.contains_points(PtsRZ),1,0)
	pts_in	  = np.where(pts_inout == 1)[0]
	pts_out   = np.where(pts_inout == 0)[0]

	return pts_in, pts_out
