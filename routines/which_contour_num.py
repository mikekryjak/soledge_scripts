import numpy as np

def which_contour_num(cin, p1):

	nArc = len(cin.arc)
	d = np.zeros(nArc, dtype='f8')
	for k in range(nArc):
		dist = np.sqrt((cin.arc[k].x - p1.x)**2 + (cin.arc[k].y - p1.y)**2)
		d[k] = np.min(dist)
		
	return np.argmin(d)
