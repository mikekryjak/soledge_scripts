import types
import numpy					as np
import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt

## Legacycontour broke... let's try just not having this
# try:
# 	from matplotlib	import _cntr as cntr
# except:
# 	from legacycontour  import _cntr as cntr


# def contour_better(x, y, f, val):


# 	c	  = types.SimpleNamespace()
# 	c.num = 0
# 	c.arc = []

# 	cs = cntr.Cntr(x, y, f)
# 	for level in val:
# 		res = cs.trace(level)
# 		nsegs = len(res) // 2
# 		segments, codes = res[:nsegs], res[nsegs:]
		
# 		c.num += nsegs
# 		for k in range(nsegs-1,-1,-1):
# #		for k in range(nsegs):
# 			c.arc.append(types.SimpleNamespace())
# 			c.arc[-1].x	= segments[k][:,0]
# 			c.arc[-1].y	= segments[k][:,1]
		
# 	return c


# def contour_better_plot(x, y, f, val):

# 	cs = pyp.contour(x, y, f, levels=val)
# 	pyp.show(block=False)
# 	pyp.close()
	
# 	contours, npaths, nsegs = get_contour_verts(cs)
	
# 	c = types.SimpleNamespace()
# 	if(len(contours) == 0):
# 		c.num = 0
# 		return

# 	c.num = npaths[0]
# 	c.arc = []
# 	for k in range(npaths[0]):
# 		c.arc.append(types.SimpleNamespace())
# 		c.arc[-1].x	= contours[0][k][:,0]
# 		c.arc[-1].y	= contours[0][k][:,1]
		
# 	return c

def contour_better(x, y, f, levels):
    """
    Returns a namespace with:
      - num : total number of contour segments
      - arc : list of objects, each having .x and .y arrays
    """
    # Matplotlib demands ascending levels
    levels = sorted(levels)

    # Create figure+axes so we can close cleanly later
    fig, ax = plt.subplots()
    cs = ax.contour(x, y, f, levels=levels)
    plt.close(fig)

    arcs = []
    # cs.allsegs is a list (over levels) of lists of Nx2 arrays
    for level_segs in cs.allsegs:
        for seg in level_segs:
            arcs.append(types.SimpleNamespace(x=seg[:, 0], y=seg[:, 1]))

    return types.SimpleNamespace(num=len(arcs), arc=arcs)


def contour_better_plot(x, y, f, levels):
    """
    Like contour_better, but also briefly shows the plot window (non-blocking).
    """
    levels = sorted(levels)
    fig, ax = plt.subplots()
    cs = ax.contour(x, y, f, levels=levels)
    plt.show(block=False)
    plt.close(fig)

    arcs = []
    for level_segs in cs.allsegs:
        for seg in level_segs:
            arcs.append(types.SimpleNamespace(x=seg[:, 0], y=seg[:, 1]))

    return types.SimpleNamespace(num=len(arcs), arc=arcs)

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
		
		