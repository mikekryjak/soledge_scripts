
# ################################################
# Interpolate data on rectangular mesh
# #################################################

#from scipy import interpolate
import numpy						as np
import numpy.matlib					as mat
import itertools					as it
import matplotlib.mlab				as ml
import scipy.interpolate			as interpolate


def interp2d(x, y, z, xnew, ynew, kind='spline'):

# routines uses the slow scipy.interpolate.interp2d when:
#	1) kind is equal to 'linear', 'cubic', 'quintic'
#	2) numbero of points is lower than 1000

	if(((kind == 'linear') or (kind == 'cubic') or (kind == 'quintic')) and (z.shape[0] < 1000)):
		print(" doing very slow scipy.interpolate.interp2d")
		f = interpolate.interp2d(y, x, z, kind=kind)     		 	#exchange x&y because interp2d use C reversed order ==>z= z[y,x]
	elif ((kind == 'spline') and (z.shape[0] < 1000)):
		print(" doing interpolate.SmoothBivariateSpline")
		BivSp = interpolate.SmoothBivariateSpline(y, x, z)

	if(isinstance(xnew, np.ndarray)):
		if(len(xnew.shape) == 2):
			nx = xnew.shape[0]
			ny = xnew.shape[1]
		elif (xnew.shape[0] > 1):
			nx = xnew.shape[0]
			ny = xnew.shape[0]

			xnew = np.transpose(mat.repmat(xnew, ny, 1))
			ynew = mat.repmat(ynew, nx, 1)
	else:
		nx = int(xnew)
		ny = int(ynew)

		min_x = np.min(x)
		max_x = np.max(x)
		xnew  = np.arange(nx, dtype='f8')*(max_x-min_x)/(nx-1)+min_x
		xnew  = np.transpose(mat.repmat(xnew, ny, 1))

		min_y = np.min(y)
		max_y = np.max(y)
		ynew  = np.arange(ny, dtype='f8')*(max_y-min_y)/(ny-1)+min_y
		ynew  = mat.repmat(ynew, nx, 1)


	if(((kind == 'linear') or (kind == 'cubic') or (kind == 'quintic')) and (z.shape[0] < 1000)):
		print("doing very slow scipy.interpolate.interp2d f computation",yx.shape)
		x1 = np.reshape(xnew, nx*ny)
		y1 = np.reshape(ynew, nx*ny)
		yx = [y1,x1]    		 							#exchange x&y because interp2d use C reversed order ==>z= z[y,x]
		znew = np.reshape(np.array(list(it.starmap(f, yx))), (ny,nx))
		znew = znew.T
	elif ((kind == 'spline') and (z.shape[0] < 1000)):
		print(" doing interpolate.SmoothBivariateSpline evaluation ")
		x1 = np.reshape(xnew, nx*ny)
		y1 = np.reshape(ynew, nx*ny)
		znew = np.reshape(BivSp.ev(y1,x1), (ny,nx))
		znew = znew.T
	else:
		znew = interpolate.griddata(np.array([x,y]).T, z, (xnew, ynew), method='linear')
		ii, jj  = np.where(np.isnan(znew))
		znew[ii,jj] = 0.
#		znew = ml.griddata(x, y, z, xnew, ynew, interp='linear')


	return xnew, ynew, znew

