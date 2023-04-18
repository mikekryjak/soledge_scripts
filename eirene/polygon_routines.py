import numpy as np
from numpy							import asarray
#from shapely.geometry 				import asLineString, LineString
from routines.intersect_contour		import intersect_2contours
from routines.find_closest_segment	import find_closest_segment
from interfaces.plot_and_ask		import plot_and_ask

#=======================================================================================================================================
# Function definition is here

#========================== 
# check many points in polygon 
#========================== 

def polygon_points_in(x, y, PolyPath):

	nPoints = len(x)
	Points  = np.empty((nPoints,2), dtype='f8')
	Points[:,0] = x
	Points[:,1] = y
		
	return PolyPath.contains_points(Points)

#========================== 
# check a single point in polygon 
#========================== 

def polygon_point_in(x, y, PolyPath):
	return PolyPath.contains_point(np.array([x,y]))
	
#=======================================================================================================================================

#========================== 
# Count crossing with edges
#========================== 

def polygon_count_in_out(x1, y1, x2, y2, nsegs, PolyPath):
	
	Points  = np.empty((nsegs,2), dtype='f8')
	Points[:,0] = x1 + np.arange(nsegs)*(x2-x1)/(nsegs-1)
	Points[:,1] = y1 + np.arange(nsegs)*(y2-y1)/(nsegs-1)

	IsIn = np.where(PolyPath.contains_points(Points), 1, 0)

	return np.sum(np.abs(IsIn[1:nsegs]-IsIn[0:nsegs-1]))

#=======================================================================================================================================

#========================== 
# check segment crossing with edges 
#========================== 

def polygon_in_out(x1, y1, x2, y2, PolyPath):

	IsInOut = np.sum(np.where(PolyPath.contains_points(np.array([[x1,y1],[x2,y2]])), 1, 0))

	if IsInOut == 1: return True			#one inside the other outside
	else: return False						#both poins inside or aoutside								


#=======================================================================================================================================

#================================
# find aproximate segment crossing with closed edges 
#================================

def polygon_inout_intersect(x1, y1, x2, y2, PolyPath):

	tol = 1.e-5
	err = 1.e0
	while err > tol:

		x3 = (x1 + x2)/2.
		y3 = (y1 + y2)/2.
		if polygon_in_out(x1, y1, x3, y3, PolyPath):
			x2 = x3
			y2 = y3
		else:
			x1 = x3
			y1 = y3

		err = np.sqrt((x1-x2)**2+(y1-y2)**2)

	x = (x1 + x2)/2.
	y = (y1 + y2)/2.

	return x,y, 0

#=======================================================================================================================================

#========================== 
# find true crossing crossing with edges 
#========================== 
"""
def polygon_true_intersect(Root, x1, y1, x2, y2, Xc, Yc):

	Seg     = LineString([(x1, y1), (x2, y2)])
	Countur = asLineString(np.array([Xc, Yc]).T)
	try:
		xy  = asarray(Countur.intersection(Seg))
		ind = 0
	except:																						#Fail when segment or part of it is on contour
		ind = 3

	if((ind != 0) or (len(xy.shape) > 1) or (xy.shape[0] != 2)): 
		d1, i1 = find_closest_segment(x1, y1, Xc, Yc)
		d2, i2 = find_closest_segment(x2, y2, Xc, Yc)
		if((d1 < 1e-6) and (d2 < 1e-6)):														#All segment is in contour
#			print("\t All segment in contour")
			return 0., 0., 3
		elif(d1 < 1e-6):																		#Sub-segment (1) in contour
			iP = np.argmin(np.sqrt((Xc[i1:i1+2] - x2)**2 + (Yc[i1:i1+2] - y2)**2))+i1
			return Xc[iP], Yc[iP], 1
		elif(d2 < 1e-6):																		#Sub-segment (2) in contour
			iP = np.argmin(np.sqrt((Xc[i2:i2+2] - x1)**2 + (Yc[i2:i2+2] - y1)**2))+i2
			return Xc[iP], Yc[iP], 2
		else:
			print("\tError in polygon_intersect no crossing")
			print("\t\tfor segment:({:.4f},{:.4f})-({:.4f},{:.4f})".format(x1,y1,x2,y2))
			Title= "Error no crossing (Contour=black, line=red)!"

			Xarrs  = [Xc, np.array([x1,x2])]
			Yarrs  = [Yc, np.array([y1,y2])]
			lines  = ["k-", "r.-"]
			LinesData = [Xarrs, Yarrs, lines]
			choice = plot_and_ask(Root, LinesData=LinesData, title="Error multiple crossing: Contour=black, line=red")
			return 0., 0., 3
	else: 
		return xy[0], xy[1], 0
"""