import types
import numpy						as np
from routines.contour_better		import contour_better
from routines.intersect_contour		import intersect_contour
from routines.which_contour_num		import which_contour_num
from routines.part_contour2			import part_contour2
from routines.show_no_intersection	import show_no_intersection
from routines.globals				import DEBUG
from interfaces.progressbar 		import ProgressBar

def find_zone_corners(Root, Config):

	if(DEBUG > 0): print("\nfind_zone_corners")

	prog_bar = ProgressBar(Root, Title="Finding zones borders", Label="Please wait", Value=0.) 

	X_points  = Config.X_points
	Frontiers = Config.Frontiers
	
	NptX	= len(X_points)
	Nsteps	= NptX*4+len(Frontiers)
	Nseg	= 0
	Npt		= NptX
	p1 	   = types.SimpleNamespace()
	thresh = 1.e-10
	for k in range(NptX):									# Loop on arcs from X-points
		for n in range(4):									# for each arc, find intersections of interest
			c1			= types.SimpleNamespace()
			c1.arc		= [types.SimpleNamespace()]
			c1.num		= 1
			c1.arc[0].x = X_points[k].cut[n].R
			c1.arc[0].y = X_points[k].cut[n].Z

#			Initializes the list of psi values on the frontier with the value
#			at the X-point and the value at the end of the domain

			psi_list 	= np.array([X_points[k].cut[n].psilim, X_points[k].psi])
			psi_type	= -np.ones(len(X_points)+1, dtype='i4')
			ind		 	= 2

#			Complete the list of psi values with all the values of psi of
#			other X-points that are between the 2 previously defined extreme values
			
#			if(X_points[k].cut[n].psilim != Config.psicore):									#Skip others separatrix values for core region
			for k1 in range(NptX):
				if(k1 != k):
					if(np.min(np.abs(psi_list-X_points[k1].psi)) > thresh):		# the psi value for this X-point has not been listed yet within the threshold
						if((X_points[k1].psi-X_points[k].psi)*(X_points[k1].psi-X_points[k].cut[n].psilim)<0):
							psi_list	  = np.append(psi_list, X_points[k1].psi)
							psi_type[ind] = k1
							ind			 += 1

#			Sort the list and keep only unique values (in case of connected X-points)

			psi_list_red, b = np.unique(psi_list, return_index=True) 
			psi_type	 	= psi_type[b]
			X_points[k].cut[n].Pt_lim = []

#			For each of the psi values, finds the intersection of the contour with the frontier

			k1 = 0
			for k2 in range(len(psi_list_red)):
				X_points[k].cut[n].Pt_lim.append(types.SimpleNamespace())
				c2 = contour_better(Config.r2D, Config.z2D, Config.flux2D, [psi_list_red[k2]])
				
				if(abs(psi_list_red[k2] - X_points[k].psi) > thresh):									#If it is not the psi of the starting X-point, there should be an intersection
					X  = intersect_contour(c1,c2)
					if(len(X) == 0):																	#If no intersection try to extend cut
						c1.arc[0].x = np.append(2*c1.arc[0].x[0]-c1.arc[0].x[1], np.append(c1.arc[0].x, 2*c1.arc[0].x[-1]-c1.arc[0].x[-2]))
						c1.arc[0].y = np.append(2*c1.arc[0].y[0]-c1.arc[0].y[1], np.append(c1.arc[0].y, 2*c1.arc[0].y[-1]-c1.arc[0].y[-2]))
						X  = intersect_contour(c1,c2)

					if(len(X) == 0):
						print("\tnX_points, n_cut=",k+1, n+1)
						print("\tpsi_list_red[k2]=",psi_list_red[k2])
						show_no_intersection(Root, Config, c1, c2)
						return False

					if(psi_type[k2] >= 0):
						FoundArc = False
						for iX in range(NptX):
							if(abs(psi_list_red[k2] - X_points[iX].psi) <= thresh):							# check al X points with same psi (connected double null)
								X_arcs = np.array([X_points[iX].branch[n].arc for n in range(4)])
								ii = np.where(X_arcs == X[0].arc2)[0]
								if(len(ii) > 0):
									FoundArc = True
									break
						if(not FoundArc): continue     				 										#By P. innocente to skip flux limit between x point and core and to prevent unnecessary zones


					X_points[k].cut[n].Pt_lim[k1].R			= X[0].x
					X_points[k].cut[n].Pt_lim[k1].Z			= X[0].y
					X_points[k].cut[n].Pt_lim[k1].psi_numc	= X[0].arc2
					X_points[k].cut[n].Pt_lim[k1].num		= Npt
					Npt += 1
				else:																					# If it is the psi of the starting X-point, we take the position of the X-point
					X_points[k].cut[n].Pt_lim[k1].R			= X_points[k].R
					X_points[k].cut[n].Pt_lim[k1].Z			= X_points[k].Z
					p1.x									= X_points[k].cut[n].R[0]
					p1.y									= X_points[k].cut[n].Z[0]
					X_points[k].cut[n].Pt_lim[k1].psi_numc	= which_contour_num(c2,p1)
					X_points[k].cut[n].Pt_lim[k1].num		= k

				X_points[k].cut[n].Pt_lim[k1].psi = psi_list_red[k2]		
				if(DEBUG > 1):	
					print("\tX_points[k].cut[n].Pt_lim[k1].num     =", X_points[k].cut[n].Pt_lim[k1].num+1)
					print("\tX_points[k].cut[n].Pt_lim[k1].psi_numc=", X_points[k].cut[n].Pt_lim[k1].psi_numc+1)
					print("\tX_points[k].cut[n].Pt_lim[k1].psi     =", X_points[k].cut[n].Pt_lim[k1].psi)
				k1	+= 1

			X_points[k].cut[n].seg_lim = []
			n_contours = k1
			for k1 in range(n_contours-1):
				p1		= types.SimpleNamespace()
				p2		= types.SimpleNamespace()
				cin		= types.SimpleNamespace()
				p1.x	= X_points[k].cut[n].Pt_lim[k1].R
				p1.y	= X_points[k].cut[n].Pt_lim[k1].Z
				p2.x	= X_points[k].cut[n].Pt_lim[k1+1].R
				p2.y	= X_points[k].cut[n].Pt_lim[k1+1].Z
				cin.x	= c1.arc[0].x
				cin.y	= c1.arc[0].y
				cout	= part_contour2(cin, p1, p2)
				
				X_points[k].cut[n].seg_lim.append(types.SimpleNamespace())
				X_points[k].cut[n].seg_lim[k1].R = cout.x
				X_points[k].cut[n].seg_lim[k1].Z = cout.y
				X_points[k].cut[n].seg_lim[k1].num = Nseg
				Nseg += 1

			X_points[k].cut[n].nPt_lim = n_contours
			Nstep		=	k*4+(n+1)												#python index
			prog_bar.Update(Nstep/Nsteps)

	Nseg2 = 0

#	Loop on poloidal limits

	for k in range(len(Frontiers)):
		c1			= types.SimpleNamespace()
		c1.arc		= [types.SimpleNamespace()]
		c1.arc[0].x = Frontiers[k].R
		c1.arc[0].y = Frontiers[k].Z
		psi_list	= np.array([Frontiers[k].psimin, Frontiers[k].psimax])
		psi_type	= -np.ones(len(X_points)+2, dtype='i4')
		ind		 	= 2
		for k1 in range(NptX):
			if(np.min(np.abs(psi_list-X_points[k1].psi)) > thresh):											# the psi value for this X-point has not been listed yet within the threshold
				if((X_points[k1].psi > Frontiers[k].psimin) and (X_points[k1].psi < Frontiers[k].psimax)):
					psi_list = np.append(psi_list, X_points[k1].psi)
					psi_type[ind] = k1
					ind += 1

		psi_list_red, b = np.unique(psi_list, return_index=True) 
		psi_type	 	= psi_type[b]
		Frontiers[k].Pt_lim = []
		k1 = 0
		for k2 in range(len(psi_list_red)):				
			c2 = contour_better(Config.r2D, Config.z2D, Config.flux2D, [psi_list_red[k2]])
			X  = intersect_contour(c1,c2)

			if(len(X) == 0):
				print("\tFrontiers=",k+1)
				print("\tpsi_list_red[k2]=",psi_list_red[k2])
				show_no_intersection(Root, Config, c1, c2)
				return False

			if(psi_type[k2] >= 0):
				FoundArc = False
				for iX in range(NptX):
					if(abs(psi_list_red[k2] - X_points[iX].psi) <= thresh):							# check al X points with same psi (connected double null)
						X_arcs = np.array([X_points[iX].branch[n].arc for n in range(4)])
						ii = np.where(X_arcs == X[0].arc2)[0]
						if(len(ii) > 0):
							FoundArc = True
							break
				if(not FoundArc): continue     				 										#By P. innocente to skip flux limit between x point and core and to prevent unnecessary zones
				
			Frontiers[k].Pt_lim.append(types.SimpleNamespace())
			Frontiers[k].Pt_lim[k1].R			= X[0].x
			Frontiers[k].Pt_lim[k1].Z			= X[0].y
			Frontiers[k].Pt_lim[k1].psi			= psi_list_red[k2]
			Frontiers[k].Pt_lim[k1].psi_numc	= X[0].arc2
			Frontiers[k].Pt_lim[k1].num			= Npt
			if(DEBUG > 1):	
				print("\tFrontiers[k].Pt_lim[k1].num      =", Frontiers[k].Pt_lim[k1].num+1)
				print("\tFrontiers[k].Pt_lim[k1].psi_numc =", Frontiers[k].Pt_lim[k1].psi_numc+1)
				print("\tFrontiers[k].Pt_lim[k1].psi      =", Frontiers[k].Pt_lim[k1].psi)
			Npt 							   += 1
			k1								   += 1

		Frontiers[k].seg_lim = []
		n_contours = k1
		for k1 in range(n_contours-1):
			p1		= types.SimpleNamespace()
			p2		= types.SimpleNamespace()
			cin		= types.SimpleNamespace()
			p1.x	= Frontiers[k].Pt_lim[k1].R
			p1.y	= Frontiers[k].Pt_lim[k1].Z
			p2.x	= Frontiers[k].Pt_lim[k1+1].R
			p2.y	= Frontiers[k].Pt_lim[k1+1].Z
			cin.x	= c1.arc[0].x
			cin.y	= c1.arc[0].y
			cout	= part_contour2(cin,p1,p2)
			
			Frontiers[k].seg_lim.append(types.SimpleNamespace())
			Frontiers[k].seg_lim[k1].R = cout.x
			Frontiers[k].seg_lim[k1].Z = cout.y
			Frontiers[k].seg_lim[k1].num = Nseg
			Nseg	+= 1
			Nseg2	+= 1

		Frontiers[k].nPt_lim = n_contours
		Nstep = (k+1)+NptX*4
		prog_bar.Update(Nstep/Nsteps)

	Config.ZoneElements			= types.SimpleNamespace()
	Config.ZoneElements.Nseg	= Nseg
	Config.ZoneElements.Nseg2	= Nseg2
	Config.ZoneElements.Npt		= Npt
	
	if(DEBUG > 1):	
		print("\tfind_zone_corners: Nseg =",Nseg)
		print("\tfind_zone_corners: Nseg2=",Nseg2)
		print("\tfind_zone_corners: Npt  =",Npt)

	prog_bar = 0

	if(DEBUG > 0): print("find_zone_corners: completed\n")

	return True

