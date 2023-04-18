import types
from tkinter						import messagebox
from tkinter.messagebox 			import Message
import numpy						as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot			as pyp
from math 							import floor, fmod, sqrt
from routines.contour_better		import contour_better
from routines.which_contour_num		import which_contour_num
from routines.intersect_contour		import intersect_contour
from routines.part_contour_per		import part_contour_per
from routines.show_no_intersection	import show_no_intersection
from routines.globals				import DEBUG
from routines.find_closest_segment	import find_closest_segment

# Detects radial sides
# Output: array of structure "cote" containing for each index:
#	- p1, p2: structures describing extremities of side:
#		-- R = radial coordinate
#		-- Z = vertical coordinate
#		-- coord = size 4 vector:
#			--- coord[0] = psi value of side
#			--- coord[1] = X-point number or frontier number to which corresponds the segment to which this extremety belongs
#			--- coord[2] = X-point cut number (to which corresponds the segment to which this extremety belongs
#						in case it corresponds to an X-point cut (semi integer between 2 cuts if branch starting or ending at X-point)
#						= -1 if poloidal frontier
#			--- coord[3] = type of segment segment to which this extremety belongs (1 = X-point cut, 2 = poloidal frontier)
#		- R, Z: coordinates of points of side arc

def find_cote_psi(self, Config):

	if(DEBUG > 0): print("\nfind_cote_psi")

	cote_colors = [['g-','r--','c--'],['r-.','c-','b--'],['c-.','b-.','y-']]

	X_points	 = Config.X_points
	ZoneElements = Config.ZoneElements
	Segments	 = Config.Segments
	
	Nseg  = ZoneElements.Nseg        					# total number of Segments
	NsegX = ZoneElements.Nseg-ZoneElements.Nseg2		#  Segments from X points (two sided)

	ncote	   = 0
	nX_points  = len(X_points)
	X_points_R = np.array([X_point.R for X_point in X_points])
	X_points_Z = np.array([X_point.Z for X_point in X_points])
	Cotes	   = []

	c1			= types.SimpleNamespace()
	c2			= types.SimpleNamespace()

#	First build a list of all psi contours that are needed for cotes
	thresh = 1.e-10
	psi_list_tmp = np.array([Segments[0].psimin,Segments[0].psimax])
	for iseg in range(1,Nseg):
		if(np.min(np.abs(psi_list_tmp-Segments[iseg].psimin)) > thresh): psi_list_tmp = np.append(psi_list_tmp, Segments[iseg].psimin)
		if(np.min(np.abs(psi_list_tmp-Segments[iseg].psimax)) > thresh): psi_list_tmp = np.append(psi_list_tmp, Segments[iseg].psimax)

	psi_list_all = np.unique(psi_list_tmp)

#	For each psi value, find all the intersection with extremeties of
#	Segments except at X-points where it is the extremeties of branches that are used

	for ipsi in range(len(psi_list_all)):
		
#		Establish a list of corner points
		ipt = -1
		corners_x		 = np.empty(0, dtype='f8')
		corners_y		 = np.empty(0, dtype='f8')
		corners_frontier = np.empty(0, dtype='f8') 		# frontier number if lies on frontier, -1 if not a frontier
		corners_nX		 = np.empty(0, dtype='i4')			# X-point number if lies on X-point cut, -1 if not an X-point cut
		corners_nXcut	 = np.empty(0, dtype='f8')			# number of X-point cut if lies on X-point cut (semi-integer if branch between 2 cuts at X-point), 0 if not an X-point cut

#		Intersections with Segments away from X-points

		for iseg in range(Nseg):
			if (abs(Segments[iseg].psimin-psi_list_all[ipsi]) < thresh): 
				distX = np.min((Segments[iseg].R[0] - X_points_R)**2 + (Segments[iseg].Z[0] - X_points_Z)**2)

				if(len(corners_x) > 0): distCorners = np.min((Segments[iseg].R[0] - corners_x)**2  + (Segments[iseg].Z[0] - corners_y)**2)
				else:					distCorners =10*thresh

				if((distX > thresh) and (distCorners > thresh)):									# It is not an X-point
#				if ((distCorners > thresh) or (ipt == -1)):				# It has not yet been listed
					ipt += 1
					corners_x = np.append(corners_x, Segments[iseg].R[0])
					corners_y = np.append(corners_y, Segments[iseg].Z[0])
					if (Segments[iseg].type == 2):								# frontier
						corners_frontier	= np.append(corners_frontier, Segments[iseg].nX)
						corners_nX			= np.append(corners_nX, -1)
						corners_nXcut		= np.append(corners_nXcut, -1)
					else:														# X-point cut
						corners_frontier	= np.append(corners_frontier, -1)
						corners_nX			= np.append(corners_nX, Segments[iseg].nX)
						corners_nXcut		= np.append(corners_nXcut, Segments[iseg].nB)
			if (abs(Segments[iseg].psimax - psi_list_all[ipsi]) < thresh):
				distX		= np.min((Segments[iseg].R[-1] - X_points_R)**2 + (Segments[iseg].Z[-1] - X_points_Z)**2)

				if(len(corners_x) > 0): distCorners = np.min((Segments[iseg].R[-1] - corners_x)**2  + (Segments[iseg].Z[-1] - corners_y)**2)
				else:					distCorners =10*thresh

				if((distX > thresh) and (distCorners > thresh)):	# It is not an X-point
#				if ((distCorners > thresh) or (ipt == -1)):					# It has not yet been listed
					corners_x = np.append(corners_x, Segments[iseg].R[-1])
					corners_y = np.append(corners_y, Segments[iseg].Z[-1])
					if(Segments[iseg].type == 2):								# frontier
						corners_frontier	= np.append(corners_frontier, Segments[iseg].nX)
						corners_nX			= np.append(corners_nX, -1)
						corners_nXcut		= np.append(corners_nXcut, -1)
					else:														# X-point cut
						corners_frontier	= np.append(corners_frontier, -1)
						corners_nX			= np.append(corners_nX, Segments[iseg].nX)
						corners_nXcut		= np.append(corners_nXcut, Segments[iseg].nB)

#		 Intersections with first points of X-point branches
		for iX in range(nX_points):
			if (abs(X_points[iX].psi - psi_list_all[ipsi]) < thresh):
				for ibranch in range(4):
					corners_x			= np.append(corners_x, X_points[iX].branch[ibranch].R)
					corners_y			= np.append(corners_y, X_points[iX].branch[ibranch].Z)
					corners_frontier	= np.append(corners_frontier, -1)
					corners_nX			= np.append(corners_nX, iX)
					corners_nXcut		= np.append(corners_nXcut, ibranch-0.5)
					print("iX                            =",iX)
					print("X_points[iX].branch[ibranch].R=",X_points[iX].branch[ibranch].R)
					print("X_points[iX].branch[ibranch].Z=",X_points[iX].branch[ibranch].Z)
#		 Now determines which segment of the contour these points lie on and
#		define cotes
		cpsi = contour_better(Config.r2D,Config.z2D,Config.flux2D, [psi_list_all[ipsi]])
		for iarc in range(len(cpsi.arc)):
			xstart				= cpsi.arc[iarc].x[:-1]
			ystart				= cpsi.arc[iarc].y[:-1]
			xend				= cpsi.arc[iarc].x[1:]
			yend				= cpsi.arc[iarc].y[1:]
			intersec_segs		= np.empty(0, dtype='i4')
			intersec_x			= np.empty(0, dtype='f8')
			intersec_y			= np.empty(0, dtype='f8')
			intersec_frontier	= np.empty(0, dtype='f8')
			intersec_nX			= np.empty(0, dtype='i4')
			intersec_nXcut		= np.empty(0, dtype='f8')
			for ipt in range(len(corners_x)):
#				d, new_seg = find_closest_segment(corners_x[ipt], corners_y[ipt] , cpsi.arc[iarc].x, cpsi.arc[iarc].y)
#				if(d < 1e-3):
				seglength2		= (xend - xstart)*(xend - xstart) + (yend - ystart)*(yend - ystart)
				crossproduct	= (corners_y[ipt] - ystart)*(xend - xstart) - (corners_x[ipt] - xstart)*(yend - ystart)
				dotproduct		= (corners_x[ipt] - xstart)*(xend - xstart) + (corners_y[ipt] - ystart)*(yend - ystart)
				new_seg			= np.where((np.abs(crossproduct) < thresh) & (dotproduct >= 0) & (dotproduct <= seglength2))[0]
				if(len(new_seg) > 0):
					intersec_segs		= np.append(intersec_segs, new_seg)
					intersec_x			= np.append(intersec_x, corners_x[ipt])
					intersec_y			= np.append(intersec_y, corners_y[ipt])
					intersec_frontier	= np.append(intersec_frontier, corners_frontier[ipt])
					intersec_nX			= np.append(intersec_nX, corners_nX[ipt])
					intersec_nXcut		= np.append(intersec_nXcut, corners_nXcut[ipt])

			sort_indices  		= np.argsort(intersec_segs)
			intersec_segs 		= intersec_segs[sort_indices]
			intersec_x			= intersec_x[sort_indices]
			intersec_y			= intersec_y[sort_indices]
			intersec_frontier	= intersec_frontier[sort_indices]
			intersec_nX			= intersec_nX[sort_indices]
			intersec_nXcut		= intersec_nXcut[sort_indices]

#			From the intersection points on the arc, define cotes
			n_intersec_segs = len(intersec_segs)
			for iseg in range(n_intersec_segs):
				isegp = iseg + 1
				if(isegp == n_intersec_segs): isegp = 0

				if ((intersec_nX[iseg] != -1) or (intersec_nX[isegp] != -1)):
#					The 2 extremeties of this cote are not on external
#					frontier so the cote is in the domain and should be defined
					if(iseg != n_intersec_segs-1):
						ilist = np.arange(intersec_segs[iseg]+1,intersec_segs[isegp]+1)
					else:
						Npts_arc   =  len(cpsi.arc[iarc].x)
						dFirstLast = sqrt((cpsi.arc[iarc].x[-1]-cpsi.arc[iarc].x[0])**2 + (cpsi.arc[iarc].y[-1]-cpsi.arc[iarc].y[0])**2)			#check for closed contours
						if(dFirstLast > 1e-9):
							ilist = np.append(np.arange(intersec_segs[iseg]+1,Npts_arc), np.arange(intersec_segs[isegp]+1))
						else:
							ilist = np.append(np.arange(intersec_segs[iseg]+1,Npts_arc-1), np.arange(intersec_segs[isegp]+1))

					if (len(ilist) > 1):
						Rcote		= cpsi.arc[iarc].x[ilist]
						Zcote		= cpsi.arc[iarc].y[ilist]
#						We keep this cote only if it does NOT intersect a cut (eliminates spurious cotes around X-points)
						good_cote = 1
						c1.arc = [types.SimpleNamespace()]
						c2.arc = [types.SimpleNamespace()]
						c1.arc[0].x		= Rcote
						c1.arc[0].y		= Zcote
						for iseg2 in range(Nseg):
							c2.arc[0].x = Segments[iseg2].R
							c2.arc[0].y = Segments[iseg2].Z
							X = intersect_contour(c1,c2)
							if (len(X) != 0):
								self.Ax.plot(Rcote, Zcote, "k-.", linewidth=3)
								self.Ax.text(Rcote[int(len(Rcote)/2)],Zcote[int(len(Zcote)/2)],"Bad cote (ipsi={:d},iarc={:d},seg={:d},seg2={:d})".format(ipsi,iarc,iseg+1,iseg2+1))
								good_cote = 0

						if (good_cote == 1):
							Cotes.append(types.SimpleNamespace())
							Cotes[-1].p1 = types.SimpleNamespace()
							Cotes[-1].p2 = types.SimpleNamespace()

#							Defining first extremety
							if(intersec_frontier[iseg] != -1):									# The first extremety of the cote lies on a frontier
								p1_case 			= 0
								Cotes[-1].p1.R		= intersec_x[iseg]
								Cotes[-1].p1.Z		= intersec_y[iseg]
								Cotes[-1].p1.coord	= np.array([psi_list_all[ipsi], intersec_frontier[iseg], -1, 2])
								RaddBeg				= np.array([intersec_x[iseg]])
								ZaddBeg				= np.array([intersec_y[iseg]])
							elif(intersec_nXcut[iseg] == floor(intersec_nXcut[iseg])):			# The first extremet of the cote lies on an X-point cut but is not at an X-point
								p1_case 			= 1
								Cotes[-1].p1.R		= intersec_x[iseg]
								Cotes[-1].p1.Z		= intersec_y[iseg]
								Cotes[-1].p1.coord	= np.array([psi_list_all[ipsi], intersec_nX[iseg], intersec_nXcut[iseg], 1])
								RaddBeg				= np.array([intersec_x[iseg]])
								ZaddBeg				= np.array([intersec_y[iseg]])
							else:																# We are at an X-point (ie between 2 cuts)
								p1_case 			= 2
								Cotes[-1].p1.R		= X_points_R[intersec_nX[iseg]]
								Cotes[-1].p1.Z		= X_points_Z[intersec_nX[iseg]]
								Cotes[-1].p1.coord	= np.array([psi_list_all[ipsi], intersec_nX[iseg], intersec_nXcut[iseg], 1])
								RaddBeg				= np.array([X_points_R[intersec_nX[iseg]], intersec_x[iseg]])
								ZaddBeg				= np.array([X_points_Z[intersec_nX[iseg]], intersec_y[iseg]])

#							Defining second extremety
							if (intersec_frontier[isegp] != -1):								# The second extremety of the cote lies on a frontier
								p2_case 			= 0
								Cotes[-1].p2.R		= intersec_x[isegp]
								Cotes[-1].p2.Z		= intersec_y[isegp]
								Cotes[-1].p2.coord	= np.array([psi_list_all[ipsi], intersec_frontier[isegp], -1, 2])
								RaddEnd				= np.array([intersec_x[isegp]])
								ZaddEnd				= np.array([intersec_y[isegp]])
							elif (intersec_nXcut[isegp] == floor(intersec_nXcut[isegp])):		# The second extremety of the cote lies on an X-point cut but is not at an X-point
								p2_case 			= 1
								Cotes[-1].p2.R		= intersec_x[isegp]
								Cotes[-1].p2.Z		= intersec_y[isegp]
								Cotes[-1].p2.coord	= np.array([psi_list_all[ipsi], intersec_nX[isegp], intersec_nXcut[isegp], 1])
								RaddEnd				= np.array([intersec_x[isegp]])
								ZaddEnd				= np.array([intersec_y[isegp]])
							else:																# We are at an X-point (ie between 2 cuts)
								p2_case 			= 2
								Cotes[-1].p2.R		= X_points_R[intersec_nX[isegp]]
								Cotes[-1].p2.Z		= X_points_Z[intersec_nX[isegp]]
								Cotes[-1].p2.coord	= np.array([psi_list_all[ipsi], intersec_nX[isegp], intersec_nXcut[isegp], 1])
								RaddEnd				= np.array([intersec_x[isegp], X_points_R[intersec_nX[isegp]]])
								ZaddEnd				= np.array([intersec_y[isegp], X_points_Z[intersec_nX[isegp]]])

							Cotes[-1].R		= np.concatenate((RaddBeg, Rcote, RaddEnd))
							Cotes[-1].Z		= np.concatenate((ZaddBeg, Zcote, ZaddEnd))
							Cotes[-1].iCote = len(Cotes)-1

							self.Ax.plot(Cotes[-1].R ,Cotes[-1].Z, cote_colors[p1_case][p2_case], linewidth=3)
							self.Ax.text(Cotes[-1].R[int(len(Cotes[-1].R)/2)],Cotes[-1].Z[int(len(Cotes[-1].R)/2)],"C_{:d}=({:d},{:d})".format(Cotes[-1].iCote+1,p1_case,p2_case))							
							if(DEBUG > 1):	
								print("\tncote             =",len(Cotes))
								print("\tCotes[-1].p1.coord=",np.array(Cotes[-1].p1.coord) + [0,1,1,0])
								print("\tCotes[-1].p2.coord=",np.array(Cotes[-1].p2.coord) + [0,1,1,0])

#						# cote has more than 1 point
#					cote in the domain (ie not connecting 2 frontiers)
#				iseg
#			iarc
	
#		ipsi

	self.Fig.canvas.draw()
	Config.Cotes = Cotes
	
	if(DEBUG > 0): print("find_cote_psi: computed\n")

	return True

def find_cote_psi_old(self, Config):

	if(DEBUG > 0): print("\nfind_cote_psi")
	
	X_points	 = Config.X_points
	ZoneElements = Config.ZoneElements
	Segments	 = Config.Segments
	
	Nseg  = ZoneElements.Nseg        					#total number of Segments
	NsegX = ZoneElements.Nseg-ZoneElements.Nseg2		#  Segments from X points (two sided)

	ncote	  = 0
	nX_points = len(X_points)
	Cotes	  = []
	
	cleft = types.SimpleNamespace()
	cleftp = types.SimpleNamespace()
	c1 	   = types.SimpleNamespace()
	cprog  = types.SimpleNamespace()
	p1 	   = types.SimpleNamespace()
	p2 	   = types.SimpleNamespace()
	c11	   = types.SimpleNamespace()
	c21	   = types.SimpleNamespace()
	for k in range(nX_points):

		for n in range(4):

#           psi_list = zeros(1,1+X_points)

			psi_list = np.array([X_points[k].cut[n].psilim, X_points[k].psi])
			psi_type = -np.ones(len(X_points)+1, dtype='i4')
			ind		 = 2
#			if(X_points[k].cut[n].psilim != Config.psicore):									#Skip others separatrix values for core region
			for k1 in range(nX_points):
				if(k1 != k):
#                   if psi_X  <  psi  <  psilim  ou bien  psilim  <  psi  <  psi_X
					if((X_points[k1].psi-X_points[k].psi)*(X_points[k1].psi-X_points[k].cut[n].psilim) < 0):
						psi_list	  = np.append(psi_list, X_points[k1].psi)
						psi_type[ind] = k1
						ind			 += 1

			b			 = np.argsort(psi_list)
			psi_list_red = psi_list[b]
			psi_type	 = psi_type[b]

#			first find intersections
			c1.arc = [types.SimpleNamespace()]
			c1.arc[0].x		= X_points[k].cut[n].R
			c1.arc[0].y		= X_points[k].cut[n].Z
			for k2 in range(len(psi_list_red)):
				if(psi_list_red[k2] == X_points[k].psi):											#special treatment here

					c2 = contour_better(Config.r2D, Config.z2D, Config.flux2D, [psi_list_red[k2]])
					for nb in range(4):																#for each leg
						X		= [types.SimpleNamespace()]
						X[0].x	= X_points[k].branch[nb].R
						X[0].y	= X_points[k].branch[nb].Z
						
						p1.x	= X[0].x
						p1.y	= X[0].y
						
						number	= which_contour_num(c2, p1)
						cprog.x = c2.arc[number].x
						cprog.y = c2.arc[number].y
						d = np.sqrt((cprog.x-X[0].x)**2+(cprog.y-X[0].y)**2)
						a = np.argmin(d)
						
						if(a + 1 < len(d)):
							v1x	= cprog.x[a+1]-cprog.x[a]
							v1y	= cprog.y[a+1]-cprog.y[a]
							vx	= X[0].x-cprog.x[a]
							vy	= X[0].y-cprog.y[a]
							
							if(v1x*vx+v1y*vy > 0):
								cleft.x = np.roll(cprog.x[:-1],-a-1)
								cleft.y = np.roll(cprog.y[:-1],-a-1)
							else:
								cleft.x = np.roll(cprog.x[:-1],-a)
								cleft.y = np.roll(cprog.y[:-1],-a)
								
							cleft.x = np.append(cleft.x,cleft.x[0])
							cleft.y = np.append(cleft.y,cleft.y[0])
						else:
							cleft.x = np.copy(cprog.x)
							cleft.y = np.copy(cprog.y)


#						should we turn left or right ?
						distX	= np.sqrt((cleft.x[0] - X_points[k].R)**2 + (cleft.y[0] - X_points[k].Z)**2)
						distX2	= np.sqrt((cleft.x[1] - X_points[k].R)**2 + (cleft.y[1] - X_points[k].Z)**2)
						if(distX2 < distX):						#cleft NOK : cleft -> cright
							sens = -1
							cleft.x = cleft.x[::-1]
							cleft.y = cleft.y[::-1]
						else:
							sens = 1


#						then same as below

						dist	= np.array([])
						type	= np.array([], dtype='i4')
						fine	= np.array([])
						num		= np.array([], dtype='i4')
						for ns in range(Nseg):
							if(Segments[ns].psimin == psi_list_red[k2]):
								if(Segments[ns].psimin_numc == number):
									X		= [types.SimpleNamespace()]
									X[0].x	= Segments[ns].R[0]
									X[0].y	= Segments[ns].Z[0]
								else:
									X		= []
								
							else:
								if(Segments[ns].psimax == psi_list_red[k2]):
									if(Segments[ns].psimax_numc == number):
										X	   = [types.SimpleNamespace()]
										X[0].x = Segments[ns].R[-1]
										X[0].y = Segments[ns].Z[-1]
									else:
										X		= []
								else:
									X	= []
							if(len(X) == 1): 											#intersect
								d	= np.sqrt((cleft.x  -X[0].x)**2 + (cleft.y - X[0].y)**2)
								a1	= np.where(d == d.min()); a1 = a1[0]
								cleftp.x = np.roll(cleft.x[:-1],-1)
								cleftp.y = np.roll(cleft.y[:-1],-1)
								cleftp.x  = np.append(cleftp.x, cleftp.x[0])
								cleftp.y  = np.append(cleftp.y, cleftp.y[0])
								dist	  = np.append(dist,a1)
								num		  = np.append(num, (ns+1)*np.ones(len(a1), dtype='i4')-1)
								if(ns+1 <= NsegX):
									type = np.append(type,np.ones(len(a1), dtype='i4'))				# Xpoint
								else:
									type = np.append(type, 2*np.ones(len(a1), dtype='i4'))			# lim

								vx	 = X[0].x - cleft.x[a1]
								vy	 = X[0].y - cleft.y[a1]
								vx1  = cleftp.x[a1]-cleft.x[a1]
								vy1  = cleftp.y[a1]-cleft.y[a1]
								fine = np.append(fine, d.min()*np.sign(vx*vx1 + vy*vy1))

#						end for nsegment

						if(len(fine) == 0):
							messagebox.showerror("Find Segments on flux", "Unable to find intersection\ntry to increase fields smoothing")
							return False

						pasmax = np.max(np.abs(fine)) + 1.
						a	   = dist*pasmax+fine
						b	   = np.argsort(a)
						a	   = a[b]
						stopit = False
						ind	   = 0
						segpa  = np.array([])
						while (not stopit):
							p1.psi = psi_list_red[k2]
							p2.psi = psi_list_red[k2]
							if(ind == 0):											#start from X point surrounding
								p1.x		= X_points[k].branch[nb].R
								p1.y		= X_points[k].branch[nb].Z
								p1.theta	= k
								p1.branch	= -2
								p1.type		= 1
							else:
								if(Segments[num[b[ind-1]]].psimin == psi_list_red[k2]):
									p1.x = Segments[num[b[ind-1]]].R[0]
									p1.y = Segments[num[b[ind-1]]].Z[0]
								else:
									p1.x = Segments[num[b[ind-1]]].R[-1]
									p1.y = Segments[num[b[ind-1]]].Z[-1]

								p1.theta	= Segments[num[b[ind-1]]].nX
								if(Segments[num[b[ind-1]]].nX != k):
									p1.branch	= Segments[num[b[ind-1]]].nB
								else:
									p1.branch =-2

								p1.type		= Segments[num[b[ind-1]]].type

							if(ind + 1 > len(b)):
								p2.x	  = X_points[k].R
								p2.y	  = X_points[k].Z
								p2.theta  = k
								p2.branch =-2
								p2.type	  = 1
							else:
								if(Segments[num[b[ind]]].psimin == psi_list_red[k2]):
									p2.x = Segments[num[b[ind]]].R[0]
									p2.y = Segments[num[b[ind]]].Z[0]
								else:
									p2.x = Segments[num[b[ind]]].R[-1]
									p2.y = Segments[num[b[ind]]].Z[-1]

								p2.theta  = Segments[num[b[ind]]].nX
								if(Segments[num[b[ind]]].nX != k):
									p2.branch = Segments[num[b[ind]]].nB
								else:											#Modif
									p2.branch =-2
								p2.type   = Segments[num[b[ind]]].type
								segpa	  = np.append(segpa,num[b[ind]])

							if(np.sqrt((p1.x-p2.x)**2 + (p1.y-p2.y)**2) == 0):
								ind += 1
							else:
								Cotes.append(types.SimpleNamespace())
								ncote			   += 1
								Cotes[-1].p1	 	= types.SimpleNamespace()
								Cotes[-1].p2	 	= types.SimpleNamespace()
								Cotes[-1].p1.R 		= p1.x
								Cotes[-1].p1.Z		= p1.y
								Cotes[-1].p2.R		= p2.x
								Cotes[-1].p2.Z		= p2.y
								cout				= part_contour_per(cprog, p1, p2, sens)
								Cotes[-1].p1.coord	= np.array([p1.psi,p1.theta,p1.branch,p1.type])
								Cotes[-1].p2.coord	= np.array([p2.psi,p2.theta,p2.branch,p2.type])
								if(DEBUG > 1):	
									print("\t1.ncote             =",ncote)
									print("\t1.Cotes[-1].p1.coord=",np.array(Cotes[-1].p1.coord) + [0,1,1,0])
									print("\t1.Cotes[-1].p2.coord=",np.array(Cotes[-1].p2.coord) + [0,1,1,0])
								if(ind + 1 > len(b)):
									
#									remove and close
									dbr = np.zeros(4)
									for nb2 in range(4):
										if(nb2 != nb):
											dbr[nb2] = np.min(np.sqrt((cout.x - X_points[k].branch[nb2].R)**2+(cout.y - X_points[k].branch[nb2].Z)**2))
										else:
											dbr[nb2] = 1e10


									idbr   = np.where(dbr == dbr.min())[0][0]
									p2.x   = X_points[k].branch[idbr].R
									p2.y   = X_points[k].branch[idbr].Z
									dist2  = np.sqrt((cout.x - p2.x)**2 + (cout.y - p2.y)**2)
									a	   = np.where(dist2 == dist2.min())[0][0]
									cout.x = np.append(cout.x[:a-1], p2.x)
									cout.y = np.append(cout.y[:a-1], p2.y)

								if(sens == 1):
									self.Ax.plot(cout.x ,cout.y, 'g.-', linewidth=3)
								else:
									self.Ax.plot(cout.x, cout.y, 'm.-', linewidth=3)

								Cotes[-1].R = np.copy(cout.x)
								Cotes[-1].Z = np.copy(cout.y)
								ind += 1
								if(ind +1 > len(b)):
									stopit = True
								else:
									segpa_num = np.where(segpa == num[b[ind]]); segpa_num = segpa_num[0]
									if((type[b[ind]] == 2) or (len(segpa_num) != 0)):						#lim
										stopit = True
									if((p2.type == 1) and (p2.theta==k)): 
										stopit = True

				else:																						#end of special treatment
					c2 = contour_better(Config.r2D, Config.z2D, Config.flux2D, [psi_list_red[k2]])
					X  = intersect_contour(c1,c2)															#should intersect once and just once		

					if(len(X) == 0):																	#If no intersection try to extend cut
						c1.arc[0].x = np.append(c1.arc[0].x, 2*c1.arc[0].x[-1]-c1.arc[0].x[-2])
						c1.arc[0].y = np.append(c1.arc[0].y, 2*c1.arc[0].y[-1]-c1.arc[0].y[-2])
						X  = intersect_contour(c1,c2)

					if(len(X) == 0):
						print("\tnX_points, n_cut=",k+1, n+1)
						print("\tpsi_list_red[k2]=",psi_list_red[k2])
						show_no_intersection(Root, Config, c1, c2)
						return False
					
					if(psi_type[k2] >= 0):
						X_arcs = np.array([X_points[psi_type[k2]].branch[n].arc for n in range(4)])
						ii = np.where(X_arcs == X[0].arc2); ii = ii[0]
						if(len(ii) == 0): continue     				 										#By P. innocente to skip flux limit between x point and core and to prevent unnecessary zones
	
					cprog.x = c2.arc[X[0].arc2].x
					cprog.y = c2.arc[X[0].arc2].y
					number = X[0].arc2
					
#					all contour close normally if extrapolation correct

#					look left
					d = np.sqrt((cprog.x - X[0].x)**2 + (cprog.y-X[0].y)**2)
					a = np.where(d == min(d)); a = a[0][0]
					if(a + 1 < len(d)):
						v1x = cprog.x[a+1]-cprog.x[a]
						v1y = cprog.y[a+1]-cprog.y[a]
						vx  = X[0].x-cprog.x[a]
						vy  = X[0].y-cprog.y[a]
						if(v1x*vx+v1y*vy > 0):
							cleft.x = np.roll(cprog.x[:-1],-a-1)
							cleft.y = np.roll(cprog.y[:-1],-a-1)
						else:
							cleft.x = np.roll(cprog.x[:-1],-a)
							cleft.y = np.roll(cprog.y[:-1],-a)

						cleft.x = np.append(cleft.x,cleft.x[0])
						cleft.y = np.append(cleft.y,cleft.y[0])
					else:
						cleft.x = np.copy(cprog.x)
						cleft.y = np.copy(cprog.y)

					dist	= np.array([])
					type	= np.array([], dtype='i4')
					fine	= np.array([])
					num		= np.array([], dtype='i4')
					for ns in range(Nseg):
						if(Segments[ns].psimin == psi_list_red[k2]):
							if(Segments[ns].psimin_numc == number):
								X		= [types.SimpleNamespace()]
								X[0].x 	= np.array([Segments[ns].R[0]])
								X[0].y 	= np.array([Segments[ns].Z[0]])
							else:
								X		= []

						else:
							if(Segments[ns].psimax == psi_list_red[k2]):
								if(Segments[ns].psimax_numc == number):
									X		= [types.SimpleNamespace()]
									X[0].x	= np.array([Segments[ns].R[-1]])
									X[0].y	= np.array([Segments[ns].Z[-1]])
								else:
									X		= []

							else:
								X	= []

						if(len(X) == 1):											#intersect
							d		 = np.sqrt((cleft.x - X[0].x)**2 + (cleft.y - X[0].y)**2)
							a1		 = np.where(d == d.min()); a1 = a1[0]
							cleftp.x = np.roll(cleft.x[:-1], -1)
							cleftp.y = np.roll(cleft.y[:-1], -1)
							cleftp.x = np.append(cleftp.x, cleftp.x[0])
							cleftp.y = np.append(cleftp.y, cleftp.y[0])
							dist	 = np.append(dist, a1)
							num		 = np.append(num, (ns+1)*np.ones(len(a1), dtype='i4')-1)
							if(ns+1 <= NsegX):
								type = np.append(type, np.ones(len(a1), dtype='i4'))			# Xpoint
							else:
								type = np.append(type, 2*np.ones(len(a1), dtype='i4'))			# lim

							vx   = X[0].x-cleft.x[a1]
							vy   = X[0].y-cleft.y[a1]
							vx1  = cleftp.x[a1]-cleft.x[a1]
							vy1  = cleftp.y[a1]-cleft.y[a1]
							fine = np.append(fine,d.min()*np.sign(vx*vx1+vy*vy1))

#					end for nsegment

					if(len(fine) == 0):
						messagebox.showerror("Find Segments on flux", "Unable to find intersection\ntry to increase fields smoothing")
						return False

					pasmax = np.max(np.abs(fine))+1
					a = dist*pasmax+fine
					b = np.argsort(a)
					a  = a[b]
					stopit	= False
					ind		= 0
					segpa	= np.array([])
					if(num[b[0]] != num[b[-1]]): b = np.append(b[-1],b)

					if(len(num) == 1):
						num  = np.append(num, num)
						b	 = np.array([0,1])
						type = np.append(type, type)

					if(ind+2 > len(b)): stopit = True

					while(not stopit):
						p1.psi = psi_list_red[k2]
						p2.psi = psi_list_red[k2]
						if(Segments[num[b[ind]]].psimin == psi_list_red[k2]):
							p1.x = Segments[num[b[ind]]].R[0]
							p1.y = Segments[num[b[ind]]].Z[0]
						else:
							p1.x = Segments[num[b[ind]]].R[-1]
							p1.y = Segments[num[b[ind]]].Z[-1]

						p1.theta	= Segments[num[b[ind]]].nX
						p1.branch	= Segments[num[b[ind]]].nB
						p1.type		= Segments[num[b[ind]]].type
						if(Segments[num[b[ind+1]]].psimin == psi_list_red[k2]):
							p2.x = Segments[num[b[ind+1]]].R[0]
							p2.y = Segments[num[b[ind+1]]].Z[0]
						else:
							p2.x = Segments[num[b[ind+1]]].R[-1]
							p2.y = Segments[num[b[ind+1]]].Z[-1]

						p2.theta	= Segments[num[b[ind+1]]].nX
						p2.branch	= Segments[num[b[ind+1]]].nB
						p2.type		= Segments[num[b[ind+1]]].type
						if((np.sqrt((p1.x-p2.x)**2+(p1.y-p2.y)**2)==0) and (psi_list_red[k2] != Config.psicore)):
							ind += 1
						else:
							segpa= np.append(segpa,num[b[ind]])
							if(np.sqrt((p1.x-p2.x)**2+(p1.y-p2.y)**2) < 1e-6):
								cout  = types.SimpleNamespace()
								cout.x= np.append(np.append(p1.x, cleft.x[:-1]), p2.x)
								cout.y= np.append(np.append(p1.y, cleft.y[:-1]), p2.y)
							else:
								cout = part_contour_per(cprog,p1,p2,1)

							OK = True
							for k3 in range(nX_points):
								for n3 in range(4):
									c11.arc		 = [types.SimpleNamespace()]
									c11.arc[0].x = np.append(X_points[k3].R, X_points[k3].branch[n3].R)
									c11.arc[0].y = np.append(X_points[k3].Z, X_points[k3].branch[n3].Z)
									
									c21.arc		 = [types.SimpleNamespace()]
									c21.arc[0].x = cout.x
									c21.arc[0].y = cout.y
									X = intersect_contour(c11, c21)
									if(len(X) > 0):	OK = False

							if(OK):
								Cotes.append(types.SimpleNamespace())
								ncote			+= 1
								Cotes[-1].p1	  = types.SimpleNamespace()
								Cotes[-1].p2	  = types.SimpleNamespace()
								Cotes[-1].p1.R = p1.x
								Cotes[-1].p1.Z = p1.y
								Cotes[-1].p2.R = p2.x
								Cotes[-1].p2.Z = p2.y
								Cotes[-1].p1.coord = np.array([p1.psi, p1.theta, p1.branch, p1.type])
								Cotes[-1].p2.coord = np.array([p2.psi, p2.theta, p2.branch, p2.type])
								self.Ax.plot(cout.x, cout.y, 'g-', linewidth=3)
								Cotes[-1].R = cout.x
								Cotes[-1].Z = cout.y
								if(DEBUG > 1):	
									print("\t2.ncote             =",ncote)
									print("\t2.Cotes[-1].p1.coord=",np.array(Cotes[-1].p1.coord) + [0,1,1,0])
									print("\t2.Cotes[-1].p2.coord=",np.array(Cotes[-1].p2.coord) + [0,1,1,0])

							ind += 1

						segpa_num = np.where(segpa == num[b[ind]]); segpa_num = segpa_num[0]
						if((type[b[ind]]==2) or (len(segpa_num) !=0) or (ind+2 > len(b))): stopit = True	#lim

#------------------------------------------------------------
#					look right

					cright	 = types.SimpleNamespace()
					crightp	 = types.SimpleNamespace()
					
					cright.x = cleft.x[::-1]
					cright.y = cleft.y[::-1]
					dist	 = np.array([])
					type	 = np.array([], dtype='i4')
					fine	 = np.array([])
					num		 = np.array([], dtype='i4')
					for ns in range(Nseg):
						if(Segments[ns].psimin == psi_list_red[k2]):
							if(Segments[ns].psimin_numc == number):
								X		= [types.SimpleNamespace()]
								X[0].x	= [Segments[ns].R[0]]
								X[0].y	= [Segments[ns].Z[0]]
							else:
								X 	= []
						else:
							if(Segments[ns].psimax == psi_list_red[k2]):
								if(Segments[ns].psimax_numc == number):
									X		= [types.SimpleNamespace()]
									X[0].x	= [Segments[ns].R[-1]]
									X[0].y	= [Segments[ns].Z[-1]]
								else:
									X	= []
							else:
								X = []

						if(len(X) == 1):														#intersect
							d		= np.sqrt((cright.x-X[0].x)**2+(cright.y-X[0].y)**2)
							a1 = np.where(d == d.min()); a1 = a1[0]
							crightp.x	= np.roll(cright.x[:-1],-1)
							crightp.y	= np.roll(cright.y[:-1],-1)
							crightp.x	= np.append(crightp.x,crightp.x[0])
							crightp.y	= np.append(crightp.y,crightp.y[0])
							dist		= np.append(dist,a1)
							num			= np.append(num, (ns+1)*np.ones(len(a1), dtype='i4')-1)
							if(ns+1 <= NsegX):
								type	= np.append(type, np.ones(len(a1), dtype='i4'))					# Xpoint
							else:
								type	= np.append(type, 2*np.ones(len(a1), dtype='i4'))				# lim

							vx	 = X[0].x-cright.x[a1]
							vy	 = X[0].y-cright.y[a1]
							vx1	 = crightp.x[a1]-cright.x[a1]
							vy1	 = crightp.y[a1]-cright.y[a1]
							fine = np.append(fine, d.min()*np.sign(vx*vx1+vy*vy1))

#					end for nsegment

					if(len(fine) == 0):
						messagebox.showerror("Find Segments on flux", "Unable to find intersection\ntry to increase fields smoothing")
						return False

					pasmax = np.max(np.abs(fine))+1
					a = dist*pasmax+fine
					b = np.argsort(a)
					a = a[b]
					if(num[b[0]] != num[b[-1]]): b= np.append(b[-1], b)

					if(len(num) == 1):
						num  = np.append(num, num)
						b	 = np.array([0, 1])
						type = np.append(type, type)

					stopit	= False
					ind		= 0
					segpa	= np.array([])
					if((ind+2 > len(b)) and (ind != 0)): stopit = True

					while(not stopit):
						p1.psi = psi_list_red[k2]
						p2.psi = psi_list_red[k2]
						if(Segments[num[b[ind]]].psimin == psi_list_red[k2]):
							p1.x = Segments[num[b[ind]]].R[0]
							p1.y = Segments[num[b[ind]]].Z[0]
						else:
							p1.x = Segments[num[b[ind]]].R[-1]
							p1.y = Segments[num[b[ind]]].Z[-1]

						p1.theta	= Segments[num[b[ind]]].nX
						p1.branch	= Segments[num[b[ind]]].nB
						p1.type		= Segments[num[b[ind]]].type
						if(Segments[num[b[ind+1]]].psimin == psi_list_red[k2]):
							p2.x = Segments[num[b[ind+1]]].R[0]
							p2.y = Segments[num[b[ind+1]]].Z[0]
						else:
							p2.x = Segments[num[b[ind+1]]].R[-1]
							p2.y = Segments[num[b[ind+1]]].Z[-1]

						p2.theta	= Segments[num[b[ind+1]]].nX
						p2.branch	= Segments[num[b[ind+1]]].nB
						p2.type		= Segments[num[b[ind+1]]].type
						if((np.sqrt((p1.x-p2.x)**2 + (p1.y-p2.y)**2)==0) and (psi_list_red[k2] != Config.psicore)): 
#	                        if(np.sqrt((p1.x-p2.x)^2+(p1.y-p2.y)^2)==0)
							ind += 1
						else:
							segpa = np.append(segpa, num[b[ind]])
							if(np.sqrt((p1.x-p2.x)**2+(p1.y-p2.y)**2) < 1e-6):
								cout   = types.SimpleNamespace()						
								cout.x = np.append(np.append(p1.x, cright.x[1:]), p2.x)
								cout.y = np.append(np.append(p1.y, cright.y[1:]), p2.y)
							else:
								cout = part_contour_per(cprog,p1,p2,-1)

							OK = True
							for k3 in range(nX_points):
								for n3 in range(4):
									c11.num = 1
									c11.arc[0].x = np.append(X_points[k3].R, X_points[k3].branch[n3].R)
									c11.arc[0].y = np.append(X_points[k3].Z, X_points[k3].branch[n3].Z)
									c21.num		 = 1
									c21.arc[0].x = cout.x
									c21.arc[0].y = cout.y
									X			 = intersect_contour(c11,c21)
									if(len(X) > 0): OK = False
									
							if(OK):
								ncote			   += 1
								Cotes.append(types.SimpleNamespace())
								Cotes[-1].p1	 	= types.SimpleNamespace()
								Cotes[-1].p2	 	= types.SimpleNamespace()
								Cotes[-1].p1.R 		= p1.x
								Cotes[-1].p1.Z 		= p1.y
								Cotes[-1].p2.R 		= p2.x
								Cotes[-1].p2.Z 		= p2.y
								Cotes[-1].p1.coord	= np.array([p1.psi,p1.theta,p1.branch,p1.type])
								Cotes[-1].p2.coord	= np.array([p2.psi,p2.theta,p2.branch,p2.type])
								if(DEBUG > 1):	
									print("\t3.ncote             =",ncote)
									print("\t3.Cotes[-1].p1.coord=",np.array(Cotes[-1].p1.coord) + [0,1,1,0])
									print("\t3.Cotes[-1].p2.coord=",np.array(Cotes[-1].p2.coord) + [0,1,1,0])
								Cotes[-1].R			= cout.x
								Cotes[-1].Z			= cout.y
								self.Ax.plot(cout.x,cout.y,'m-', linewidth=3)

							ind += 1
							
						segpa_num = np.where(segpa == num[b[ind]]); segpa_num = segpa_num[0]
						if((type[b[ind]]==2) or (len(segpa_num) != 0) or (ind+2 > len(b))): stopit = True				#lim

#					end if not separatrice

	self.Fig.canvas.draw()
	Config.Cotes = Cotes
	
	if(DEBUG > 0): print("find_cote_psi: computed\n")

	return True

