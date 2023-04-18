import types
import numpy					as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot 		as pyp
try:
	from matplotlib	import _cntr as cntr
except:
	from legacycontour  import _cntr as cntr


def contour_better(x, y, f, val):


	c	  = types.SimpleNamespace()
	c.num = 0
	c.arc = []

	cs = cntr.Cntr(x, y, f)
	for level in val:
		res = cs.trace(level)
		nsegs = len(res) // 2
		segments, codes = res[:nsegs], res[nsegs:]
		
		c.num += nsegs
		for k in range(nsegs-1,-1,-1):
#		for k in range(nsegs):
			c.arc.append(types.SimpleNamespace())
			c.arc[-1].x	= segments[k][:,0]
			c.arc[-1].y	= segments[k][:,1]
		
	return c


def contour_better_plot(x, y, f, val):

	cs = pyp.contour(x, y, f, levels=val)
	pyp.show(block=False)
	pyp.close()
	
	contours, npaths, nsegs = get_contour_verts(cs)
	
	c = types.SimpleNamespace()
	if(len(contours) == 0):
		c.num = 0
		return

	c.num = npaths[0]
	c.arc = []
	for k in range(npaths[0]):
		c.arc.append(types.SimpleNamespace())
		c.arc[-1].x	= contours[0][k][:,0]
		c.arc[-1].y	= contours[0][k][:,1]
		
	return c

#=========================================================
# This routine generate mesh nodes from flux contour file
#=========================================================

def get_contour_verts(cn):
	contours = []
	npaths   = []
	nsegs	 = []
#	for each contour line
	for cc in cn.collections:
		paths = []
		nsegs1= []
		npath = 0
#		for each separate section of the contour line
		for pp in cc.get_paths():
			npath += 1
			xy = []
#			for each segment of that section
			nseg = 0
#			for vv in pp.iter_segments():
			for vv in pp.vertices:
				xy.append(vv)
				nseg += 1
			nsegs1.append(nseg)
			paths.append(np.vstack(xy))
		nsegs.append(nsegs1)
		npaths.append(npath)
		contours.append(paths)
	return contours, npaths, nsegs
		
		