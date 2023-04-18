
# Function definition is here

import types
import os
import numpy 						as np
import scipy.io
import h5py
from routines.h5_routines		import h5_read

#=========================================================
# This routine read experimental data from HDF5 file
#=========================================================

def load_exp_data(exp, shot, time):

	data			= types.SimpleNamespace()
	data.shot		= shot

	data.ts			= types.SimpleNamespace()
	data.lp			= types.SimpleNamespace()
	data.sp			= types.SimpleNamespace()
	data.ts.ok		= False
	data.lp.ok		= False
	data.sp.ok		= False

	path = os.getenv(exp+"_path",'.')
	try:
		if_data = h5py.File(path+"/"+exp+"_{:d}.h5".format(shot), "r")
	except:
		print("\tNot found h45 file: "+path+"/"+exp+"_{:d}.h5".format(shot),)
		return data

	try:
		group	= "thomson"
		times	= h5_read(if_data, group+ '/time', order='C')
		rpos	= h5_read(if_data, group+ '/rpos', order='C')
		zpos	= h5_read(if_data, group+ '/zpos', order='C')
		te		= h5_read(if_data, group+ '/te',   order='C')
		dte		= h5_read(if_data, group+ '/dte',  order='C')
		ne		= h5_read(if_data, group+ '/ne',   order='C')
		dne		= h5_read(if_data, group+ '/dne',  order='C')
		rho		= h5_read(if_data, group+ '/rho',  order='C')
		psi		= h5_read(if_data, group+ '/psi',  order='C')

		i			= np.argmin(np.abs(times-time))
		data.ts.time	= times[i]
		data.ts.rpos	= np.copy(rpos)
		data.ts.zpos	= np.copy(zpos)
		data.ts.te		= np.copy(te[i,:])
		data.ts.dte		= np.copy(dte[i,:])
		data.ts.ne		= np.copy(ne[i,:])
		data.ts.dne		= np.copy(dne[i,:])
		data.ts.rho		= np.copy(rho[i,:])
		data.ts.psi		= np.copy(psi[i,:])
		data.ts.ok		= True
	except:
		pass


	try:
		group	= "lp"
		times		= h5_read(if_data, group+ "/time",			order='C')
		fit_params	= h5_read(if_data, group+ "/fit_params",	order='C')
		probe		= h5_read(if_data, group+ "/probe",		order='C')
		pos_rz		= h5_read(if_data, group+ "/pos_rz",		order='C')
		pos_l		= h5_read(if_data, group+ "/pos_l",		order='C')
		dens		= h5_read(if_data, group+ "/dens",			order='C')
		densg		= h5_read(if_data, group+ "/densg",		order='C')
		te			= h5_read(if_data, group+ "/te",			order='C')
		jsat		= h5_read(if_data, group+ "/jsat",			order='C')
		jsatp		= h5_read(if_data, group+ "/jsatp",		order='C')
		vfl			= h5_read(if_data, group+ "/vfl",			order='C')
		r2			= h5_read(if_data, group+ "/r2",			order='C')
		alpha		= h5_read(if_data, group+ "/alpha",		order='C')
		pperp		= h5_read(if_data, group+ "/pperp",		order='C')
		ang			= h5_read(if_data, group+ "/ang",			order='C')
		area		= h5_read(if_data, group+ "/area",			order='C')
#		jsat2		= h5_read(if_data, group+ "/jsat2",			order='F')

		i			= np.argmin(np.abs(times-time))

		data.lp.time		= times[i]
		data.lp.fit_params	= np.copy(fit_params)
		data.lp.probe		= np.copy(probe)
		data.lp.pos_rz		= np.copy(pos_rz)
		data.lp.pos_l		= np.copy(pos_l)

		data.lp.dens		= np.copy(dens[i,:])
		data.lp.densg		= np.copy(densg[i,:])
		data.lp.te			= np.copy(te[i,:])
		data.lp.jsat		= np.copy(jsat[i,:])
		data.lp.jsatp		= np.copy(jsatp[i,:])
		data.lp.vfl			= np.copy(vfl[i,:])
		data.lp.r2			= np.copy(r2[i,:])
		data.lp.alpha		= np.copy(alpha[i,:])
		data.lp.pperp		= np.copy(pperp[i,:])
		data.lp.ang			= np.copy(ang[i,:])
		data.lp.area		= np.copy(area[i,:])
#		data.lp.jsat2		= np.copy(jsat2[i,:])
		data.lp.ok			= True
	except:
		pass

	try:
		group	= "strike_pts"
		n_sp	= h5_read(if_data, group+"/n_sp",	order='F')
		times	= h5_read(if_data, group+ "/times_1",	order='F')
		sp		= np.zeros((n_sp.shape[0]), dtype=[("t", 'f8'), ("r",'f8',(n_sp.max())), ("z",'f8',(n_sp.max())), ("pos_l",'f8',(n_sp.max()))])
		sp		= sp.view(np.recarray)
		
		for k in range(n_sp.max()):
			rz 		= h5_read(if_data, group+ "/rz_{:d}".format(k+1),	order='F')
			pos_l	= h5_read(if_data, group+ "/pos_l_{:d}".format(k+1),	order='F')
			ii = np.where(n_sp > k); ii = ii[0]
			sp.r[ii,k]		= rz[:,0]
			sp.z[ii,k]	 	= rz[:,1]
			sp.pos_l[ii,k] 	= pos_l

		i	= np.argmin(np.abs(times-time))
		data.sp.time	= times[i]
		data.sp.n_sp	= n_sp[i]
		data.sp.r		= np.copy(sp.r[i,:n_sp[i]])
		data.sp.z		= np.copy(sp.z[i,:n_sp[i]])
		data.sp.pos_l	= np.copy(np.sort(sp.pos_l[i,:n_sp[i]]))
		data.sp.ok		= True
	except:
		pass
	

	if_data.close()


	return data

