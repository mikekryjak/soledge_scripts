import numpy 			as np
import numpy.matlib		as mat
from matplotlib.path 	import Path

def extend_mat1(a):
	b = np.zeros((a.shape[0]+1, a.shape[1]+1), dtype=a.dtype)
	b[:-1,:-1] = a
	b[-1, :-1] = a[-1,:]
	b[:-1,-1]  = a[:,-1]
	b[-1,-1]   = b[-1,-2]
	return b

def extend_mat(a):
	b = np.zeros((a.shape[0]+2, a.shape[1]+2), dtype=a.dtype)
	b[1:-1,1:-1] = a
	b[0,1:-1]    = a[0,:]
	b[-1,1:-1]   = a[-1,:]
	b[:,0]		 = b[:,1]
	b[:,-1]		 = b[:,-2]
	return b


def append_RZ(ROld, ZOld, RAdd, ZAdd):
	if(len(RAdd) > 1):
		Nx = len(RAdd)
		if(len(ROld) == 0):
			RNew = np.copy(RAdd.reshape(Nx,1))
			ZNew = np.copy(ZAdd.reshape(Nx,1))
		else:
			RNew = np.append(ROld, RAdd.reshape(Nx,1), axis=1)
			ZNew = np.append(ZOld, ZAdd.reshape(Nx,1), axis=1)
	else:
		if(len(ROld) == 0):
			RNew = np.copy(RAdd)
			ZNew = np.copy(ZAdd)
		else:
			RNew = np.append(ROld, RAdd)
			ZNew = np.append(ZOld, ZAdd)
		
	return RNew, ZNew

def argsort_mat(a):
	b  = np.argsort(a.reshape(-1))
	ii = mat.repmat(np.arange(a.shape[0]),a.shape[1],1).T.reshape(-1)[b]
	jj = mat.repmat(np.arange(a.shape[1]),1,a.shape[0]).reshape(-1)[b]
	return ii,jj


#	define path

def def_contour_path(R, Z):

	ContourPath = Path(np.array([R,Z]).T, closed=True)

	return ContourPath


def safe_eval(a):
	try:
		v = eval(a)
	except:
		v = 0.
	return v

def min_and_arg(a, axis=None):

	iMin = np.argmin(a, axis=axis)
	return a[iMin], iMin



def max_and_arg(a, axis=None):

	iMax = np.argmax(a, axis=axis)
	return a[iMax], iMax

def get_clockwise(x,y):
	Area	= np.sum((y[1:]+y[:-1])*(x[1:]-x[:-1]))+(y[0]+y[-1])*(x[0]-x[-1])
	if(Area > 0):	return 1
	else:			return -1
	


