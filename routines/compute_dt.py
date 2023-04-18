import types
import numpy 	as np
from mesh.broadcast_mesh import broadcast_mesh

def compute_dt(Config, Plasma, Metric, Refpar, Transp, Dn_p, CFL=0.3):


	Zones	= Config.Zones
	nZones	= len(Zones)
	
	if(len(Plasma) == 0):
		Plasma		= set_default_plasma(Zones)
		RefPar		= types.SimpleNamespace()
		RefPar.T0eV	= 500.
		RefPar.n0	= 5.e19

	reverse_B = False

	eps			= 0.
	min_density = 0.							#global_variables%min_density

#	from physics.f90
#	################

	kb			= 1.3806e-23
	eV			= 1.6022e-19
	m_u			= 1.6605e-27
	m_e			= 9.109e-31
	epsilon_0	= 8.854187818e-12
	
#	from compute_reference_geometry
#	###############################

	R0		= 2.								#reference_parameters%geometry%R0
	rs0	 	= 0.1								#reference_parameters%geometry%rs0
	
	n0		= Refpar.n0							#reference_parameters%fields%n0
	T0eV	= Refpar.T0eV
	T0		= Refpar.T0eV*eV/kb					#reference_parameters%fields%T0
	c0		= np.sqrt(kb*T0/m_u)				#reference_parameters%fields%c0
	tau0	= 2.*np.pi*R0/c0					#reference_parameters%fields%tau0

#	phi0	= ...									#reference_parameters%fields%phi0

	broadcast_mesh(Zones)
	
	Dt =[]																#prepare structure
	for iIon in range(len(Plasma)-1):
		Dt.append([])
		for k in range(nZones):
			Dt[-1].append(types.SimpleNamespace())
			Dt[-1][-1].Adv  = np.ones_like(Zones[k].gridRc)
			Dt[-1][-1].Diff = np.ones_like(Zones[k].gridRc)
		
	for k in range(nZones):

		bc_B = np.sqrt(Zones[k].bc_Br**2. + Zones[k].bc_Bz**2. + Zones[k].bc_Bphi**2.)
		B	 = bc_B[1:-1,1:-1]

		if(len(Metric) == 0):								#Compute Metric

#			from compute_jacobian.f90
#			==========================================
		
			T = Zones[k].bc_z*2.*np.pi
			P = Zones[k].bc_x
			R = Zones[k].bc_gridRc
			Z = Zones[k].bc_gridZc

#			Jacobian_dRdP = np.empty_like(Z)
#			Jacobian_dRdT = np.empty_like(Z)
#			Jacobian_dzdP = np.empty_like(Z)
#			Jacobian_dzdT = np.empty_like(Z)
		
#			computation of psi and theta derivatives
#			Jacobian_dRdP[1:-1,1:-1], Jacobian_dRdT[1:-1,1:-1] = dpsi_dtheta_centre(R,P,T)
#			Jacobian_dzdP[1:-1,1:-1], Jacobian_dzdT[1:-1,1:-1] = dpsi_dtheta_centre(Z,P,T)
		
			Jacobian_dRdP, Jacobian_dRdT = dpsi_dtheta_centre(R,P,T)
			Jacobian_dzdP, Jacobian_dzdT = dpsi_dtheta_centre(Z,P,T)

#			from compute_jacobian.f90(compute_inverse_derivative)
#			=====================================================

			det = Jacobian_dRdT*Jacobian_dzdP - Jacobian_dzdT*Jacobian_dRdP
			ii, jj = np.where(np.abs(det) > 1e-30)
		
			Jacobian_dTdR = np.zeros_like(Jacobian_dzdP)
			Jacobian_dTdz = np.zeros_like(Jacobian_dzdP)
			Jacobian_dPdR = np.zeros_like(Jacobian_dzdP)
			Jacobian_dPdz = np.zeros_like(Jacobian_dzdP)

			Jacobian_dTdR[ii,jj] = Jacobian_dzdP[ii,jj]/det[ii,jj]
			Jacobian_dTdz[ii,jj] = Jacobian_dRdP[ii,jj]/det[ii,jj]
			Jacobian_dPdR[ii,jj] = Jacobian_dzdT[ii,jj]/det[ii,jj]
			Jacobian_dPdz[ii,jj] = Jacobian_dRdT[ii,jj]/det[ii,jj]
		
#			from compute_metric_coefficients.f90
#			====================================
		
#			metric_Jacobian = np.where(Jacobian_dPdR*Jacobian_dTdz -
#											 Jacobian_dPdz*Jacobian_dTdR == 0., 0., np.abs(Zones[k].R/
#											(Jacobian_dPdR*Jacobian_dTdz -
#											 Jacobian_dPdz*Jacobian_dTdR)))
		
			metric_ctt = (Jacobian_dTdz**2. + Jacobian_dTdR**2.)*(rs0/(2.*np.pi))**2.

			metric_cpp = (Jacobian_dPdz**2. + Jacobian_dPdR**2.)*rs0**2.

			metric_cpt = (Jacobian_dPdz*Jacobian_dTdz + 
									Jacobian_dPdR*Jacobian_dTdR)*rs0**2./(2.*np.pi)

			"""
			diff_c   = metric_ctt*metric_cpp - metric_cpt**2.
			ii, jj = np.where(diff_c != 0.)

			metric_c_tt = np.zeros_like(Jacobian_dzdP)
			metric_c_pt = np.zeros_like(Jacobian_dzdP)
			metric_c_pp  = np.zeros_like(Jacobian_dzdP)
			
			metric_c_tt[ii,jj] = 1./diff_c[ii,jj]*metric_cpp[ii,jj]
		
			metric_c_pt[ii,jj] = -1./diff_c[ii,jj]*metric_cpt[ii,jj]

			metric_c_pp[ii,jj]   = 1./diff_c[ii,jj]*metric_ctt[ii,jj]
			"""

			metric_G = np.where(B != 0., (Zones[k].Br*Jacobian_dTdR + Zones[k].Bz*Jacobian_dTdz)/B*R0, 1e-15)

#			from compute_perp_velocity.f90
#			==============================

			Jacobian	= Zones[k].gridRc/(Jacobian_dPdR*Jacobian_dTdz - Jacobian_dPdz*Jacobian_dTdR)

		else:
			metric_cpp			= Metric[k].cpp[1:-1,1:-1]
			metric_ctt			= Metric[k].ctt[1:-1,1:-1]
			metric_G				= Metric[k].G[1:-1,1:-1]
			Jacobian				= Metric[k].Jac[1:-1,1:-1]

#		Gradient/curvature drift

		dBdP,dBdT	  = gradient_centre(bc_B, Zones[k].bc_x,2.*np.pi*Zones[k].bc_z)

		drifts_uBp = -2.*T0eV/Plasma[0][0].charge*Plasma[0][k].Temp[1:-1,1:-1]*dBdT*Zones[k].gridRc*Zones[k].Bphi/ \
						((Jacobian + eps)*(B**3. + eps))/(np.sqrt(metric_cpp)/rs0 + eps)
		drifts_uBt = 2.*T0eV/Plasma[0][0].charge*Plasma[0][k].Temp[1:-1,1:-1]*dBdP*Zones[k].gridRc*Zones[k].Bphi/ \
						((Jacobian + eps)*(B**3. + eps))/(np.sqrt(metric_cpp)/rs0 + eps)

		drifts_uBp /= c0
		drifts_uBt /= c0
		
		if(reverse_B):
			drifts_uBt = -drifts_uBt
			drifts_uBp = -drifts_uBp

		drifts_uEt = np.zeros_like(drifts_uBt)
		
		"""
#		Electric drift


		dPhidP,dPhidT = gradient_centre(zone%electric_fields(1)%phi_smooth, Zones[k].x,2.*np.pi*Zones[k].z)

		drifts_uEp = -phi0*dPhidT*Zones[k].R[1:-1,1:-1]*Zones[k].Bphi/((Jacobian + eps)*(B**2. + eps))/(np.sqrt(metric_cpp)/rs0 + eps)
		drifts_uEt =  phi0*dPhidP*Zones[k].R[1:-1,1:-1]*Zones[k].Bphi/((Jacobian + eps)*(B**2. + eps))/(np.sqrt(metric_ctt)/rs0 + eps)

#		adim
		drifts_uEp /= c0
		drifts_uEt /= c0
		
		if(reverse_B):
			drifts_uBt = -drifts_uBt
			drifts_uEt = -drifts_uEt
		"""

#		From compute_dt.f90
#		===================
		
		ii, jj  = np.where(Zones[k].Chi == 0)
		if(len(ii) > 0):
			bc_ii   = ii + 1
			bc_jj   = jj + 1

			for iIon in range(len(Plasma)-1):						#Now for each ion
				iPlasma = iIon+1
			
#				from compute_velocity_and_Mach.f90
#				===================================

				velocity = np.where(Plasma[iPlasma][k].Dens > min_density*0.1, Plasma[iPlasma][k].Gamma/Plasma[iPlasma][k].Dens, 0.)

#				from compute_transport_coefficients.f90(compute_perp_DD_diffusivities)
#				======================================================================

				perpD_p	= Dn_p[iIon][k]/(rs0**2/tau0)
			
#			advection terms
			
				vel		= velocity[ii,jj] + (drifts_uEt[ii,jj] + drifts_uBt[ii,jj])*np.sqrt(metric_ctt[ii,jj])*(2.*np.pi*R0/rs0)/metric_G[ii,jj]
							
				csad				= np.sqrt(Plasma[iPlasma][0].charge*Plasma[0][k].Temp[ii+1,jj+1] + Plasma[iPlasma][k].Temp[ii+1,jj+1])/Plasma[iPlasma][0].mass
#				print("iPlasma,ii,jj=",iPlasma,ii,jj)
	#			print("Dt[iIon][k].Adv.shape,csad.shape,metric_G.shape=",Dt[iIon][k].Adv.shape,csad.shape,metric_G.shape)
	#			print("Zones[k].bc_zp.shape, Zones[k].bc_zm.shape,vel.shape=",Dt[iIon][k].Adv.shape,csad.shape,metric_G.shape)
				Dt[iIon][k].Adv[ii,jj]	= np.abs((Zones[k].bc_zp[bc_ii,bc_jj] - Zones[k].bc_zm[bc_ii,bc_jj])/(metric_G[ii,jj]*(csad + np.abs(vel))))
				Dt[iIon][k].Adv[ii,jj]   *= tau0*CFL
			
#		diffusion term

				Dt[iIon][k].Diff[ii,jj]	= (0.5*(Zones[k].bc_x[bc_ii+1,bc_jj]-Zones[k].bc_x[bc_ii-1,bc_jj]))**2./(metric_cpp[ii,jj]*perpD_p[ii,jj])
				Dt[iIon][k].Diff[ii,jj]   *= tau0*CFL

	return Dt


#	from compute_jacobian.f90
#	=========================

def dpsi_dtheta_centre(X,P,T):
	dXdP	= (X[2:, 1:-1]*   (P[1:-1, 1:-1] - P[:-2, 1:-1])**2. - \
			   X[:-2, 1:-1]*  (P[2:, 1:-1]   - P[1:-1, 1:-1])**2. + \
			   X[1:-1, 1:-1]*((P[2:, 1:-1]  - P[1:-1, 1:-1])**2. - (P[1:-1, 1:-1]-P[:-2, 1:-1])**2.)) / \
			  ((P[2:, 1:-1]-P[1:-1, 1:-1])*(P[1:-1, 1:-1]-P[:-2, 1:-1])*(P[2:, 1:-1]-P[:-2, 1:-1]))
	dXdT	= (X[1:-1, 2:]*   (T[1:-1, 1:-1] - T[1:-1, :-2])**2. - \
			   X[1:-1, :-2]*  (T[1:-1, 2:]   - T[1:-1, 1:-1])**2. + \
			   X[1:-1, 1:-1]*((T[1:-1, 2:]   - T[1:-1, 1:-1])**2. - (T[1:-1, 1:-1]-T[1:-1, :-2])**2.))/ \
			  ((T[1:-1, 2:]-T[1:-1, 1:-1])*(T[1:-1, 1:-1]-T[1:-1, :-2])*(T[1:-1, 2:]-T[1:-1, :-2]))

	return dXdP, dXdT

#	from compute_perp_velocity.f90
#	==============================

def gradient_centre(X,P,T):
	dXdP	= (X[2:, 1:-1]*   (P[1:-1, 1:-1] - P[:-2, 1:-1])**2. - \
			   X[:-2, 1:-1]*  (P[2:, 1:-1]   - P[1:-1, 1:-1])**2. + \
			   X[1:-1, 1:-1]*((P[2:, 1:-1]  - P[1:-1, 1:-1])**2. - (P[1:-1, 1:-1]-P[:-2, 1:-1])**2.)) / \
			  ((P[2:, 1:-1]-P[1:-1, 1:-1])*(P[1:-1, 1:-1]-P[:-2, 1:-1])*(P[2:, 1:-1]-P[:-2, 1:-1]))
	dXdT	= (X[1:-1, 2:]*   (T[1:-1, 1:-1] - T[1:-1, :-2])**2. - \
			   X[1:-1, :-2]*  (T[1:-1, 2:]   - T[1:-1, 1:-1])**2. + \
			   X[1:-1, 1:-1]*((T[1:-1, 2:]   - T[1:-1, 1:-1])**2. - (T[1:-1, 1:-1]-T[1:-1, :-2])**2.))/ \
			  ((T[1:-1, 2:]-T[1:-1, 1:-1])*(T[1:-1, 1:-1]-T[1:-1, :-2])*(T[1:-1, 2:]-T[1:-1, :-2]))

	return dXdP, dXdT


def set_default_plasma(Zones):
	Plasma = []
	
	masses  = [0.,20.014]
	charges = [-1.,1.]
	for i in range(2):
		Plasma.append([])
	
		for k in range(len(Zones)):
				
			Plasma[i].append(types.SimpleNamespace())
			Plasma[i][k].Dens	= np.ones((Zones[k].gridRc.shape[0]+2,Zones[k].gridRc.shape[1]+2), dtype='f8')
			Plasma[i][k].Temp	= np.ones_like(Plasma[i][k].Dens)
					
			Plasma[i][k].velocity  = np.zeros_like(Plasma[i][k].Dens)
		
			Plasma[i][k].Alpham = np.zeros_like(Zones[k].gridRc)
			Plasma[i][k].Alphap = np.zeros_like(Zones[k].gridRc)
			
		Plasma[i][0].charge	= charges[i]
		Plasma[i][0].mass	= masses[i]
		
	return Plasma
