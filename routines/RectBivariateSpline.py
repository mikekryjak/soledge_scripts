import numpy		as np
from scipy			import interpolate

def RectBivariateSpline(x, y, z, x_out, y_out, s=0.):
	z_min  = z.min()
	z_max  = z.max()
	z_norm = (z-z_min)/(z_max-z_min)
	SpFun  = interpolate.RectBivariateSpline(x, y, z_norm, s=s*0.05)

	if(len(x_out.shape) == 1):
		out	= SpFun.ev(x_out, y_out)*(z_max-z_min) + z_min
	else:
		out = np.reshape(SpFun.ev(x_out.reshape(x_out.size), y_out.reshape(y_out.size))*(z_max-z_min) + z_min, reshape(x_out.shape))

	return out
