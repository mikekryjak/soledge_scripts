
#====================================================================================================
# Find intersection between two segments
#====================================================================================================

import numpy as np 

def intersection(xa, ya, xb, yb):

#	print("intersection:xa, ya, xb, yb",xa, ya, xb, yb)
	xa_min = min(xa)
	ya_min = min(ya)
	xa_max = max(xa)
	ya_max = max(ya)


	xb_min = min(xb)
	yb_min = min(yb)
	xb_max = max(xb)
	yb_max = max(yb)

	x_max_min = max(xa_min, xb_min)
	x_min_max = min(xa_max, xb_max)

	y_max_min = max(ya_min, yb_min)
	y_min_max = min(ya_max, yb_max)

#	print("intersection:x_max_min,x_min_max,y_max_min,y_min_max",x_max_min,x_min_max,y_max_min,y_min_max)
	if (x_max_min > x_min_max) or (y_max_min > y_min_max) : return 0., 0., False		#non intersection between two containing squares	


	dx_a = xa[1] - xa[0]
	dx_b = xb[1] - xb[0]

	dy_a = ya[1] - ya[0]
	dy_b = yb[1] - yb[0]

	if (x_min_max - x_max_min) < (y_min_max - y_max_min) :
		m_a = np.nan
		m_b = np.nan
		if dx_a != 0. : m_a = dy_a/dx_a
		if dx_b != 0. : m_b = dy_b/dx_b
		if ((dx_a == 0.) and (dx_b == 0.)) or (m_a == m_b) : return 0., 0., False		#non intersection they are parallel
	else :
		InVm_a = np.nan
		InVm_b = np.nan
		if dy_a != 0. : InVm_a = dx_a/dy_a
		if dy_b != 0. : InVm_b = dx_b/dy_b
		if ((dy_a == 0.) and (dy_b == 0.)) or (InVm_a == InVm_b) : return 0., 0., False		#non intersection they are parallel

#	print("intersection: dx_a,dx_b,dy_a,dy_b",dx_a,dx_b,dy_a,dy_b)
	if (x_min_max - x_max_min) < (y_min_max - y_max_min) :
		if dx_a == 0. :
			xi = xa[0]
			if dy_b != 0. :
				yi = m_b*xi - m_b*xb[0] + yb[0]
			else :
				yi = yb[0]
		elif dx_b == 0. :
			xi = xb[0]
			if dy_a != 0. :
				yi = m_a*xi - m_a*xa[0] + ya[0]
			else :
				yi = ya[0]
		elif dy_a == 0. :
			yi = ya[0]
			xi = (yi - yb[0] + m_b*xb[0])/m_b
		elif dy_b == 0. :
			yi = yb[0]
			xi = (yi - ya[0] + m_a*xa[0])/m_a
		else :
			xi = (m_b*xb[0] - m_a*xa[0] - (yb[0] - ya[0]))/(m_b - m_a)
			yi = m_a*xi - m_a*xa[0] + ya[0]
		if (np.isnan(xi) or np.isnan(yi)) : print("1: xa,ya,xb,yb,dx_a,dy_a,dx_b,dy_b=",xa,ya,xb,yb,dx_a,dy_a,dx_b,dy_b)
	else :
		if dy_a == 0. :
			yi = ya[0]
			if dx_b != 0. :
				xi = m_b*yi - InVm_b*yb[0] + xb[0]
			else :
				xi = xb[0]
		elif dy_b == 0. :
			yi = yb[0]
			if dx_a != 0. :
				xi = InVm_a*yi - InVm_a*ya[0] + xa[0]
			else :
				xi = xa[0]
		elif dx_a == 0. :
			xi = xa[0]
			yi = (xi - xb[0] + InVm_b*yb[0])/InVm_b
		elif dx_b == 0. :
			xi = xb[0]
			yi = (xi - xa[0] + InVm_a*ya[0])/m_a
		else :
			yi = (InVm_b*yb[0] - InVm_a*ya[0] - (xb[0] - xa[0]))/(InVm_b - InVm_a)
			xi = InVm_a*yi - InVm_a*ya[0] + xa[0]
		if (np.isnan(xi) or np.isnan(yi)) : print("2: xa,ya,xb,yb,dx_a,dy_a,dx_b,dy_b=",xa,ya,xb,yb,dx_a,dy_a,dx_b,dy_b)

#	print("intersection:xi,x_max_min,x_min_max,yi,y_max_min,y_min_max",xi,x_max_min,x_min_max,yi,y_max_min,y_min_max)
	if ((xi < x_max_min) or (xi > x_min_max) or (yi < y_max_min) or (yi > y_min_max)): return 0., 0., False									#intersection is outside segments	
#	print("intersection: xi, yi=",xi, yi)

	if (np.isnan(xi) or np.isnan(yi)) : print("3:xa,ya,xb,yb,dx_a,dy_a,dx_b,dy_b=",xa,ya,xb,yb,dx_a,dy_a,dx_b,dy_b)
	return xi, yi, True


