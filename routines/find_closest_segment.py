import numpy as np

def find_closest_segment(xPt, yPt, xArr, yArr):
	vx	 = xArr[1:] - xArr[:-1]
	vy	 = yArr[1:] - yArr[:-1]
	vx1	= xPt - xArr[:-1]
	vy1	= yPt - yArr[:-1]
	vx2	= xPt - xArr[1:]
	vy2	= yPt - yArr[1:]
	p	= (vx*vx1+vy*vy1)/(vx**2+vy**2)
	d	= np.where((p > 0) & (p<1), \
								np.abs(vx*vy1-vy*vx1)/np.sqrt(vx**2+vy**2), \
								np.where(np.sqrt(vx1**2+vy1**2) < np.sqrt(vx2**2+vy2**2), \
										 np.sqrt(vx1**2+vy1**2),  np.sqrt(vx2**2+vy2**2)))
	i	= np.argmin(d)
	d	= d[i]
	return d, i
