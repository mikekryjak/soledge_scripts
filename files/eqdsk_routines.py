from tkinter 						import messagebox
from tkinter.messagebox 			import Message
import types
import numpy						as np
from scipy							import interpolate
from matplotlib.path				import Path
from routines.RectBivariateSpline	import RectBivariateSpline
from routines.globals				import DEBUG

#====================================================================================================
# Compute fields and flux from an EQDSK file
#====================================================================================================

def eqdsk_compute_fields(eqdsk, nr=0, nz=0, smooth=1., out_length=-1., IgnoreErrors=0):

	if(DEBUG > 0): print("eqdsk_compute_fields")

	r_arr_eqdsk = np.arange(eqdsk.nw, dtype='f8')*eqdsk.rdim/(eqdsk.nw-1) + eqdsk.rleft
	z_arr_eqdsk = np.arange(eqdsk.nh, dtype='f8')*eqdsk.zdim/(eqdsk.nh-1) + (eqdsk.zmid - 0.5*eqdsk.zdim)
#	psirz_eqdsk = eqdsk.psirz + eqdsk.sibry										#add flux at boundary

	psirz_eqdsk = np.copy(eqdsk.psirz) 											#add flux at boundary

	if(out_length > -1.):
		r_min = np.min(eqdsk.rlim) - out_length
		r_max = np.max(eqdsk.rlim) + out_length
		z_min = np.min(eqdsk.zlim) - out_length
		z_max = np.max(eqdsk.zlim) + out_length

		i_min = np.where(r_arr_eqdsk >= r_min)[0][0]
		i_max = np.where(r_arr_eqdsk <= r_max)[0][-1]
		j_min = np.where(z_arr_eqdsk >= z_min)[0][0]
		j_max = np.where(z_arr_eqdsk <= z_max)[0][-1]

		r_min = max(r_arr_eqdsk[i_min], r_min)
		r_max = min(r_arr_eqdsk[i_max], r_max)
		z_min = max(z_arr_eqdsk[j_min], z_min)
		z_max = min(z_arr_eqdsk[j_max], z_max)

		if(DEBUG > 1): print("\tMesh reduced from {:d}x{:d} to {:d}x{:d}".format(int(eqdsk.nw),int(eqdsk.nh),i_max-i_min+1,j_max-j_min+1))
		r_arr_eqdsk = r_arr_eqdsk[i_min:i_max+1]
		z_arr_eqdsk = z_arr_eqdsk[j_min:j_max+1]
		psirz_eqdsk = psirz_eqdsk[j_min:j_max+1, i_min:i_max+1]
	else:
		r_min = r_arr_eqdsk[0]
		r_max = r_arr_eqdsk[-1]
		z_min = z_arr_eqdsk[0]
		z_max = z_arr_eqdsk[-1]
		

	Reshape = False
	if((nr > r_arr_eqdsk.size) or (r_min != r_arr_eqdsk[0]) or (r_max != r_arr_eqdsk[-1])):
		r_arr	= np.linspace(r_min, r_max, nr)
		Reshape = True	
	else:
		r_arr 	= r_arr_eqdsk

	if((nz > z_arr_eqdsk.size)or (z_min != z_arr_eqdsk[0]) or (z_max != z_arr_eqdsk[-1])):
		z_arr = np.linspace(z_min, z_max, nz)
		Reshape = True
	else:
		z_arr = z_arr_eqdsk

	r_mat, z_mat = np.meshgrid(r_arr, z_arr)
	r_mat1	 	 = r_mat.reshape(r_mat.size)
	z_mat1 		 = z_mat.reshape(r_mat.size)
	if(Reshape or (smooth != 0.)):
		if(DEBUG > 1): print("\tsmoothing flux: smooth=",smooth)

		psirz = np.reshape(RectBivariateSpline(r_arr_eqdsk, z_arr_eqdsk, psirz_eqdsk.T,  r_mat1, z_mat1, s=smooth), r_mat.shape)		#RectBivariateSpline wants MAT[x.size, y.size]
	else:
		psirz = psirz_eqdsk

	if(r_arr[0] == 0.):				#skip first point
		psirz	= psirz[:,1:]
		r_mat	= r_mat[:,1:]
		z_mat	= z_mat[:,1:]
		r_arr	= r_arr[1:]
		nr	    = r_arr.size
		r_mat1	= r_mat.reshape(r_mat.size)
		z_mat1 	= z_mat.reshape(r_mat.size)

	nTot = r_mat.size
	Br	 = np.empty(psirz.shape, dtype='f8')
	Bz	 = np.empty(psirz.shape, dtype='f8')

#	Python indexes matrix Bz[z,r]

	Dr	 = r_arr[2] - r_arr[0]
	Dz	 = z_arr[2] - z_arr[0]
	Bz[:,1:-1]	 = -(psirz[:,2:] - psirz[:,:-2])/Dr/r_mat[:,1:-1]
	Br[1:-1,:]	 =  (psirz[2:,:] - psirz[:-2,:])/Dz/r_mat[1:-1,:]

	Bz[:,0]	 	 = -(psirz[:,1]  - psirz[:,0])/(0.5*Dr)/r_arr[0]
	Bz[:,-1]	 = -(psirz[:,-1] - psirz[:,-2])/(0.5*Dr)/r_arr[-1]
	Br[0,:]		 =  (psirz[1,:]  - psirz[0,:])/(0.5*Dz)/r_arr
	Br[-1,:]	 =  (psirz[-1,:] - psirz[-2,:])/(0.5*Dz)/r_arr

#	to compute toroidal field first we find points inside seapratrix

#	define boundary path

	nbbbs = eqdsk.rbbbs.shape[0]
	if((eqdsk.rbbbs[0] == eqdsk.rbbbs[nbbbs-1]) and (eqdsk.zbbbs[0] == eqdsk.zbbbs[nbbbs-1])):
		nBoundPath = nbbbs
	else:
		nBoundPath = nbbbs + 1

	BoundVerts = np.empty((nBoundPath,2),dtype='f8')
	BoundCodes = np.empty(nBoundPath,dtype='i')

	if((eqdsk.rbbbs[0] == eqdsk.rbbbs[nbbbs-1]) and (eqdsk.zbbbs[0] == eqdsk.zbbbs[nbbbs-1])):
		BoundVerts[:,0] = eqdsk.rbbbs
		BoundVerts[:,1] = eqdsk.zbbbs
	else:
		BoundVerts[0:nBoundPath-1,0] = eqdsk.rbbbs
		BoundVerts[0:nBoundPath-1,1] = eqdsk.zbbbs
		BoundVerts[nBoundPath-1,0]   = eqdsk.rbbbs[0]
		BoundVerts[nBoundPath-1,1]   = eqdsk.zbbbs[0]

	BoundCodes[0:nBoundPath] = Path.LINETO
	BoundCodes[0] = Path.MOVETO
	BoundCodes[nBoundPath-1] = Path.CLOSEPOLY
	BoundPath = Path(BoundVerts, BoundCodes); BoundCodes= 0; BoundVerts = 0

	rz_mat1  = np.empty((nTot,2), dtype='f8')
	rz_mat1[:,0] = r_mat1
	rz_mat1[:,1] = z_mat1
	r_mat1		 = 0; z_mat1 = 0
	
	IsIn	= BoundPath.contains_points(rz_mat1)									#find point iside wall contour
	IsIn	= IsIn.reshape(r_mat.shape)
	ii, jj	= np.where(IsIn)														#find indexes of points inside separatrix
	rz_mat1 = 0; IsIn = 0
	Bphi = eqdsk.bcentr*eqdsk.rcentr/r_mat											#define set vacuum toroidal field everywhere

#	Attention! x coordinate must be sorted ascending for InterpolatedUnivariateSpline or interp

	psirz_min  = np.min(psirz[ii,jj])
	psirz_max  = np.max(psirz[ii,jj])
	
	Continue_loop=True
	Test=-1
	while Continue_loop:
		Test +=1
		Continue_loop=False
		if(eqdsk.simag > eqdsk.sibry):
			psi_fpol = np.arange(eqdsk.nw, dtype='f8')*(eqdsk.simag - eqdsk.sibry)/(eqdsk.nw-1) + eqdsk.sibry
			fpol	 = np.copy(eqdsk.fpol[::-1])											#reverse order for Fpol
			if(abs(psirz_max-eqdsk.simag) > 0.05*(eqdsk.simag - eqdsk.sibry)):
				if(abs(psirz_max + eqdsk.sibry -eqdsk.simag) < 0.05*(eqdsk.simag - eqdsk.sibry)):
					if(DEBUG > 1): print("\tadd separatrix flux value to flux")
					psirz = psirz + eqdsk.sibry
				elif(abs(-psirz_min-eqdsk.simag) < 0.05*(eqdsk.sibry - eqdsk.simag)):
					if(DEBUG > 1): print("\treverse  flux")
					psirz = -psirz
					Br	  = -Br
					Bz	  = -Bz
				elif(abs(-psirz_min+eqdsk.sibry-eqdsk.simag) < 0.05*(eqdsk.sibry - eqdsk.simag)):
					if(DEBUG > 1): print("\treverse  flux and add separatrix flux value to flux")
					psirz = -psirz+eqdsk.sibry
					Br	  = -Br
					Bz	  = -Bz
				elif(Test == 0):
					if(DEBUG > 1): print("\tchanged sibry snd simag with sibry1 and simag1")
					eqdsk.sibry = eqdsk.sibry1
					eqdsk.simag = eqdsk.simag1
					Continue_loop=True
				else:
					print("\tError in EQDSK data flux values do not agree !!")
					print("\t\teqdsk.sibry                     =",eqdsk.sibry)
					print("\t\teqdsk.simag                     =",eqdsk.simag)
					print("\t\tepsirz_min                      =",psirz_min)
					print("\t\tepsirz_max                      =",psirz_max-eqdsk.sibry)
					if(IgnoreErrors != 1):
						exit()
					else:
						MessageString = "EQDSK fluxes do not agree!!\n" + \
										"  eqdsk.sibry = {:.3f}\n" + \
										"  eqdsk.simag = {:.3f}\n" + \
										"  epsirz_min  = {:.3f}\n" + \
										"  epsirz_max  = {:.3f}"
						messagebox.showwarning("Eqdsk fields computation fields",  MessageString.format(eqdsk.sibry, eqdsk.simag, psirz_min, psirz_max-eqdsk.sibry))
		else:
			psi_fpol = np.arange(eqdsk.nw, dtype='f8')*(eqdsk.sibry - eqdsk.simag)/(eqdsk.nw-1) + eqdsk.simag
			fpol	 = np.copy(eqdsk.fpol)
			if(abs(psirz_min-eqdsk.simag) > 0.05*(eqdsk.sibry - eqdsk.simag)):
				if(abs(psirz_min+eqdsk.sibry-eqdsk.simag) < 0.05*(eqdsk.sibry - eqdsk.simag)):
					if(DEBUG > 1): print("\tadd separatrix flux value to flux")
					psirz = psirz + eqdsk.sibry
					psirz_min = np.min(psirz[ii,jj])
				elif(abs(-psirz_max-eqdsk.simag) < 0.05*(eqdsk.sibry - eqdsk.simag)):
					if(DEBUG > 1): print("\treverse  flux")
					psirz = -psirz
					Br	  = -Br
					Bz	  = -Bz
				elif(abs(-psirz_max+eqdsk.sibry-eqdsk.simag) < 0.05*(eqdsk.sibry - eqdsk.simag)):
					if(DEBUG > 1): print("\treverse  flux and add separatrix flux value to flux")
					psirz = -psirz+eqdsk.sibry
					Br	  = -Br
					Bz	  = -Bz
				elif(Test == 0):
					if(DEBUG > 1): print("\tchanged sibry snd simag with sibry1 and simag1")
					eqdsk.sibry = eqdsk.sibry1
					eqdsk.simag = eqdsk.simag1
					Continue_loop=True
				else:
					print("\tError in EQDSK data flux values do not agree !!")
					print("\t\teqdsk.sibry                     =",eqdsk.sibry)
					print("\t\teqdsk.simag                     =",eqdsk.simag)
					print("\t\tepsirz_min                      =",psirz_min-eqdsk.sibry)
					print("\t\tepsirz_max                      =",psirz_max)
					if(IgnoreErrors != 1):
						exit()
					else:
						MessageString = "EQDSK fluxes do not agree!!\n" + \
										"  eqdsk.sibry = {:.3f}\n" + \
										"  eqdsk.simag = {:.3f}\n" + \
										"  epsirz_min  = {:.3f}\n" + \
										"  epsirz_max  = {:.3f}"
						messagebox.showwarning("Eqdsk fields computation fields",  MessageString.format(eqdsk.sibry, eqdsk.simag, psirz_min, psirz_max-eqdsk.sibry))
		
#	If necessary fix 2pi flux  computing plasma current 

	bsBz	= RectBivariateSpline(r_mat[0,:], z_mat[:,0], Bz.T, eqdsk.rbbbs, eqdsk.zbbbs, s=0.)		#RectBivariateSpline wants MAT[x.size, y.size]

	bsBr	= RectBivariateSpline(r_mat[0,:], z_mat[:,0], Br.T, eqdsk.rbbbs, eqdsk.zbbbs, s=0.)
	SpB		= 0

	bsBp	= np.sqrt(bsBz**2 + bsBr**2)
	bsBp	= 0.5*(bsBp[1:]+bsBp[:-1])
	bsDs	= np.sqrt((eqdsk.rbbbs[1:]-eqdsk.rbbbs[:-1])**2 + (eqdsk.zbbbs[1:]-eqdsk.zbbbs[:-1])**2)
	Icom	= np.sum(bsBp*bsDs)/(4.e-7*np.pi)
	bsBp	= 0; bsDs = 0.

	if(1-abs(eqdsk.current/Icom) > 0.5):
		if(DEBUG > 1): print("\tRescaling Br and Bz because computed I is too big ")
		if(DEBUG > 1): print("\t\tNew I/Icomp = ",abs(eqdsk.current*(2.*np.pi)/Icom))
		Br			/= 2.*np.pi
		Bz			/= 2.*np.pi
		psirz		/= 2.*np.pi
		psi_fpol	/= 2.*np.pi
	elif(1-abs(eqdsk.current/Icom) < -0.5):
		print("\eqdsk_compute_fields: Attention!! Computed Id is too small, there is an error in EQDSK file?")


#	print("eqdsk.bcentr*eqdsk.rcentr= ",eqdsk.bcentr*eqdsk.rcentr)
#	print("eqdsk.simag= ",eqdsk.simag," eqdsk.sibry= ",eqdsk.sibry)
#	print("psi_fpol.min()= ",psi_fpol.min()," psi_fpol.max()= ",psi_fpol.max())
#	print("fpol.min()= ",fpol.min()," fpol.max()= ",fpol.max())
#	print("np.min(psirz[ii,jj])= ",np.min(psirz[ii,jj])," np.max(psirz[ii,jj])= ",np.max(psirz[ii,jj]))

	SpFpol   = interpolate.InterpolatedUnivariateSpline(psi_fpol, fpol)				#spline of fpol
	Bphi[ii,jj] = SpFpol(psirz[ii,jj])/r_mat[ii,jj]									#Bphi inside separatrix from Fpol of psi

	Bphi_in_ave = np.average(Bphi[ii,jj])
	if(Bphi_in_ave*eqdsk.bcentr < 0.): 												#Check in sep toroidal field
		Bphi[ii,jj] = -Bphi[ii,jj]
		if(DEBUG > 1): print("\tChanged sign to toroidal field inside separatrix")

	rz_b_psi = np.empty((r_mat.shape[0], r_mat.shape[1], 6), dtype='f8')

	rz_b_psi[:,:,0] = r_mat
	rz_b_psi[:,:,1] = z_mat
	rz_b_psi[:,:,2] = Br
	rz_b_psi[:,:,3] = Bz
	rz_b_psi[:,:,4] = Bphi
	rz_b_psi[:,:,5] = psirz

	if(DEBUG > 0): print("eqdsk_compute_fields: Completed")

	return rz_b_psi

#====================================================================================================
# Read a EQDSK file
#====================================================================================================

def load_eqdsk_file(file, f_data_type='f8', i_data_type='i4'):

	if(DEBUG > 0): print("load_eqdsk_file")

	try:
		fid = open(file,'r')
	except OSError:
		print("\tNot found file:",file)
		return

	eqdsk = types.SimpleNamespace()
#	line n.1

	tmp  = fid.readline()
	"""
	len_val  = 45
	case_dtype = "<U{:d}".format(len_val)

	s_start  = 0
	case  = tmp[s_start : s_start + len_val]; s_start += len_val

	values  = eval(",".join(tmp[s_start :].split()))

	idum    = int(values[0])
	nw	    = int(values[1])
	nh      = int(values[2])
	"""
	values  = tmp.split()

	eqdsk.case		= values[0]
	eqdsk.idum    	= int(eval(values[-3]))
	eqdsk.nw	    = int(eval(values[-2]))
	eqdsk.nh      	= int(eval(values[-1]))
	if(DEBUG > 1): print("\tnw,nh",eqdsk.nw,eqdsk.nh)
#	line n.2

#====================================================================================
	OffsetR	= 0.0											#ATTENTION HAND TEMPORARY SETTING TOOFFSET PLASMA POSITION
	OffsetZ	= 0.0
#====================================================================================

	values = eqdsk_read_line(fid, 5, 16)
	eqdsk.rdim		= float(values[0])
	eqdsk.zdim		= float(values[1])
	eqdsk.rcentr	= float(values[2]) + OffsetR
	eqdsk.rleft		= float(values[3]) + OffsetR
	eqdsk.zmid		= float(values[4]) + OffsetZ

	if(DEBUG > 1): print("\trdim,zdim,rcentr,rleft,zmid=",eqdsk.rdim,eqdsk.zdim,eqdsk.rcentr,eqdsk.rleft,eqdsk.zmid)

#	line n.3

	values = eqdsk_read_line(fid, 5, 16)
	eqdsk.rmaxis	= float(values[0]) + OffsetR			#R of magnetic axis in meter
	eqdsk.zmaxis	= float(values[1]) + OffsetZ			#Z of magnetic axis in meter
	eqdsk.simag1	= float(values[2])						#poloidal flux at magnetic axis in Weber /rad
	eqdsk.sibry1	= float(values[3])						#poloidal flux at the plasma boundary in Weber /rad
	eqdsk.bcentr	= float(values[4])						#Vacuum toroidal magnetic field in Tesla at RCENTR

#	line n.4

	values = eqdsk_read_line(fid, 5, 16)
	eqdsk.current	= float(values[0])						#Plasma current in Ampere
	eqdsk.simag		= float(values[1])						#poloidal flux at magnetic axis in Weber /rad
	eqdsk.xdum		= float(values[2])
	eqdsk.rmaxis	= float(values[3]) + OffsetR			#R of magnetic axis in meter
	eqdsk.xdum		= float(values[4])

#	line n.5

	values = eqdsk_read_line(fid, 5, 16)
	eqdsk.zmaxis	= float(values[0]) + OffsetZ			#Z of magnetic axis in meter
	eqdsk.xdum		= float(values[1])
	eqdsk.sibry		= float(values[2])						#poloidal flux at the plasma boundary in Weber /rad
	eqdsk.xdum		= float(values[3])
	eqdsk.xdum		= float(values[4])


	eqdsk.fpol		= eqdsk_read_array(fid, eqdsk.nw, f_data_type,	5)
	eqdsk.pres		= eqdsk_read_array(fid, eqdsk.nw, f_data_type,	5)
	eqdsk.ffprim	= eqdsk_read_array(fid, eqdsk.nw, f_data_type,	5)
	eqdsk.pprime	= eqdsk_read_array(fid, eqdsk.nw, f_data_type,	5)
	eqdsk.psirz		= eqdsk_read_mat(fid, eqdsk.nw, eqdsk.nh, f_data_type,	5)
	eqdsk.qpsi		= eqdsk_read_array(fid, eqdsk.nw, f_data_type,	5)

	values = eqdsk_read_line(fid, 2, 5, can_skip=True)
	eqdsk.nbbbs		= int(values[0])
	eqdsk.limitr	= int(values[1])
	if(DEBUG > 1): print("\tnbbbs,limitr=",eqdsk.nbbbs,eqdsk.limitr)

#	read separatrix contour

	try:
		rzbbbs	= eqdsk_read_array(fid, 2*eqdsk.nbbbs,  f_data_type, 5, can_rise=True)
	except:
		rzbbbs	= eqdsk_read_array(fid, 2*eqdsk.nbbbs,  f_data_type, 2)

	rzbbbs	= rzbbbs.reshape((eqdsk.nbbbs,2))
	eqdsk.rbbbs	= rzbbbs[:,0]
	eqdsk.zbbbs	= rzbbbs[:,1]

#	read wall profile

	try:
		rzlim	= eqdsk_read_array(fid, 2*eqdsk.limitr, f_data_type, 5, can_rise=True)
	except:
		rzlim	= eqdsk_read_array(fid, 2*eqdsk.limitr, f_data_type, 2)

	rzlim	= rzlim.reshape((eqdsk.limitr,2))
	eqdsk.rlim	= rzlim[:,0]
	eqdsk.zlim	= rzlim[:,1]

	fid.close()

	if((eqdsk.rlim[0] != eqdsk.rlim[-1]) or (eqdsk.zlim[0] != eqdsk.zlim[-1])):				#close wall if not closed
		eqdsk.rlim	= np.append(eqdsk.rlim, eqdsk.rlim[0])
		eqdsk.zlim	= np.append(eqdsk.zlim, eqdsk.zlim[0])



	if(DEBUG > 1): print("\t\tPlama current(MA) = ",1e-6*abs(eqdsk.current))
	if(DEBUG > 1): print("\t\tBt on axis    (T) = ",abs(eqdsk.bcentr))
	if(DEBUG > 1): print("\t\tq95               = ",eqdsk.qpsi[int(eqdsk.qpsi.size*0.95-1)])

	if(DEBUG > 0): print("load_eqdsk_file: Completed")


	return eqdsk
			
# read a line and get values
#====================================================================================================

def eqdsk_read_line(fid, n_line, len_values, can_skip=False, can_rise=False):

	if(can_rise): position = fid.tell()

	tmp		= fid.readline()
	if(can_skip and (len(tmp) < 2)): tmp = fid.readline()

	values  = tmp.split()
	if(len(values) < n_line):
		if(len(tmp) >= n_line*len_values):
			values = []
			for k in range(n_line): values.append(eval(tmp[len_values*k:len_values*(k+1)]))
		else:
			if(can_rise):
				fid.seek(position, 0)
				raise
			else:
				print("\teqdsk_read_line: ERROR")
				print("\n\t\tUnable to read {:d} values".format(n_line))
				print("\t\tFrom line: ==>",tmp,"<==")
				print("\t\tlen(line)={:d}, needed={:d}".format(len(tmp)-1,n_line*len_values))
				exit()
	else:
		values  = eval(",".join(values))

	return values

# read array
#====================================================================================================

def eqdsk_read_array(file_id, nx, data_type, n_line, can_rise=False):

	Array = np.empty(nx, dtype=data_type)

	n_full_line = int(nx/n_line)
	n_last_line = nx - n_full_line*n_line
	
	k = 0
	for l in range(n_full_line):
		if(l == 0): values = eqdsk_read_line(file_id, n_line, 16, can_skip=True, can_rise=can_rise)
		else:		values = eqdsk_read_line(file_id, n_line, 16, can_skip=False, can_rise=can_rise)

		for v in range(n_line):
			Array[k] = float(values[v]); k += 1

	if(n_last_line > 0):
		values = eqdsk_read_line(file_id, n_last_line, 16, can_skip=False, can_rise=can_rise)
		if(n_last_line == 1): 
			Array[k] = float(values)
		else:
			for v in range(n_last_line):
				Array[k] = float(values[v]); k += 1

	return Array



# read matrix
#====================================================================================================

def eqdsk_read_mat(file_id, nx, ny, data_type, n_line, can_rise=False):

	Mat = np.empty((ny, nx), dtype=data_type)

	n_full_line = int(nx*ny/n_line)
	n_last_line = nx*ny - n_full_line*n_line
	
	i = 0
	j = 0
	for l in range(n_full_line):
		if(l == 0): values = eqdsk_read_line(file_id, n_line, 16, can_skip=True, can_rise=can_rise)
		else:		values = eqdsk_read_line(file_id, n_line, 16, can_skip=False, can_rise=can_rise)

		for v in range(n_line):
			Mat[j,i] = float(values[v]); i += 1
			if(i == nx): i=0; j += 1


	if(n_last_line > 0):
		values = eqdsk_read_line(file_id, n_last_line, 16, can_skip=False, can_rise=can_rise)
		s_start = 0
		if(n_last_line == 1): 
			Mat[j,i] = float(values)
		else:
			for v in range(n_last_line):
				Mat[j,i] = float(values[v]); i += 1
				if(i == nx): i=0; j += 1


	return Mat


#====================================================================================================
# Read a EQDSK file
#====================================================================================================

def save_eqdsk_file(file, eqdsk):

	if(DEBUG > 0): print("save_eqdsk_file")

	fid = open(file,'w')

#	line n.1
#	fid.write("{:s} {:4d} {:4d} {:4d}\n".format(eqdsk.case,0,eqdsk.nw,eqdsk.nh))
	fid.write("{:s} {:4d} {:4d} {:4d}\n".format(eqdsk.case, eqdsk.idum,eqdsk.nw,eqdsk.nh))

#	line n.2
	fid.write("{:16.9e}{:16.9e}{:16.9e}{:16.9e}{:16.9e}\n".format(eqdsk.rdim, eqdsk.zdim, eqdsk.rcentr, eqdsk.rleft, eqdsk.zmid))

#	line n.3
	fid.write("{:16.9e}{:16.9e}{:16.9e}{:16.9e}{:16.9e}\n".format(eqdsk.rmaxis, eqdsk.zmaxis, eqdsk.simag, eqdsk.sibry, eqdsk.bcentr))

#	line n.4
	fid.write("{:16.9e}{:16.9e}{:16.9e}{:16.9e}{:16.9e}\n".format(eqdsk.current, eqdsk.simag, 0., eqdsk.rmaxis, 0.))

#	line n.5
	fid.write("{:16.9e}{:16.9e}{:16.9e}{:16.9e}{:16.9e}\n".format(eqdsk.zmaxis, 0., eqdsk.sibry, 0., 0.))

	eqdsk_write_array(eqdsk.fpol, fid, 5)
	eqdsk_write_array(eqdsk.pres, fid, 5)
	eqdsk_write_array(eqdsk.ffprim, fid, 5)
	eqdsk_write_array(eqdsk.pprime, fid, 5)
	eqdsk_write_mat(eqdsk.psirz, fid, 5)
	eqdsk_write_array(eqdsk.qpsi, fid, 5)

	fid.write("{:4d} {:4d}\n".format(eqdsk.nbbbs, eqdsk.limitr))

	rzbbbs	= np.array([eqdsk.rbbbs, eqdsk.zbbbs]).T.reshape(-1)
	eqdsk_write_array(rzbbbs, fid, 5)

	rzlim	= np.array([eqdsk.rlim, eqdsk.zlim]).T.reshape(-1)
	eqdsk_write_array(rzlim, fid, 5)

	fid.close()

	if(DEBUG > 0): print("write_eqdsk_file: Completed")


# write array
#====================================================================================================

def eqdsk_write_array(Array, file_id, n_line):

	nx = len(Array)
	n_full_line = int(nx/n_line)
	n_last_line = nx - n_full_line*n_line
	
	k = 0
	for l in range(n_full_line):
		for i in range(n_line): 
			file_id.write("{:16.9e}".format(float(Array[k])))
			k += 1
		file_id.write("\n")

	if(n_last_line > 0):
		for l in range(n_last_line):
			file_id.write("{:16.9e}".format(float(Array[k])))
			k += 1
		file_id.write("\n")

	file_id.write("\n")

# write matrix
#====================================================================================================

def eqdsk_write_mat(Mat, file_id, n_line):

	ny =Mat.shape[0]
	nx =Mat.shape[1]

	n_full_line = int(nx*ny/n_line)
	n_last_line = nx*ny - n_full_line*n_line
	
	i = 0
	j = 0
	for l in range(n_full_line):
		for v in range(n_line): 
			file_id.write("{:16.9e}".format(float(Mat[j,i])))
			i += 1
			if(i == nx): i=0; j += 1
		file_id.write("\n")

	if(n_last_line > 0):
		for v in range(n_last_line):
			file_id.write("{:16.9e}".format(float(Mat[j,i])))
			i += 1
			if(i == nx): i=0; j += 1
		file_id.write("\n")

	file_id.write("\n")
