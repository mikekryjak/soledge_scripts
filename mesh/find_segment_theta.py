import types
import numpy 				as np
import scipy.interpolate	as interp
from routines.globals		import DEBUG

# Define segments along theta
# Output: array of structure "Segments" containing for each index:
#	- psimin_numc = index of psi
#	- psimax = psimax
#	- type = type of segment (1= inside mesh, 2 = poloidal border)
#	- Xtype = Type of  Xpoint
#	- nX = Number of Xpoint
#	- nB = Number of branch of Xpoint
#	- p1, p2 = coordinate of segment extreme (psi, nX, nB, Type)
#	- R, Z: coordinates of segment

def find_segment_theta(self, Config):

	if(DEBUG > 0): print("\nfind_segment_theta")

	X_points  = Config.X_points
	Frontiers = Config.Frontiers

	f_psi	  = interp.RectBivariateSpline(Config.r2D[0,:], Config.z2D[:,0], Config.flux2D.T)					#RectBivariateSpline want [x,y] order

	nX_points = len(X_points)
	
	Segments  = []
	for k in range(Config.ZoneElements.Nseg): Segments.append(types.SimpleNamespace())
	
	for k in range(nX_points):
		if(DEBUG > 1):	
			print("\tX_points        ",k+1)
			print("\tX_points[k].psi ",X_points[k].psi)
		for n in range(4):
			for ns in range(X_points[k].cut[n].nPt_lim-1):
				wh_min = np.where(np.append(X_points[k].cut[n].Pt_lim[ns].psi, X_points[k].cut[n].Pt_lim[ns+1].psi) == \
							np.min(np.append(X_points[k].cut[n].Pt_lim[ns].psi, X_points[k].cut[n].Pt_lim[ns+1].psi)))[0][0]
				wh_max = np.where(np.append(X_points[k].cut[n].Pt_lim[ns].psi, X_points[k].cut[n].Pt_lim[ns+1].psi) == \
							np.max(np.append(X_points[k].cut[n].Pt_lim[ns].psi, X_points[k].cut[n].Pt_lim[ns+1].psi)))[0][0]

				iSeg = X_points[k].cut[n].seg_lim[ns].num
				Segments[iSeg].iSeg			= iSeg
				Segments[iSeg].type			= 1								# inside mesh
				Segments[iSeg].Xtype		= X_points[k].cut[n].type
				Segments[iSeg].nX			= k
				Segments[iSeg].nB			= n

				Segments[iSeg].psimin		= np.min([X_points[k].cut[n].Pt_lim[ns].psi, X_points[k].cut[n].Pt_lim[ns+1].psi])
				Segments[iSeg].psimin_numc	= X_points[k].cut[n].Pt_lim[ns+wh_min].psi_numc

				Segments[iSeg].psimax		= np.max(np.append(X_points[k].cut[n].Pt_lim[ns].psi, X_points[k].cut[n].Pt_lim[ns+1].psi))
				Segments[iSeg].psimax_numc	= X_points[k].cut[n].Pt_lim[ns+wh_max].psi_numc

				psi1 = float(f_psi.ev(X_points[k].cut[n].seg_lim[ns].R[0],  X_points[k].cut[n].seg_lim[ns].Z[0]))								#float to convert from numpy.ndarray class
				psi2 = float(f_psi.ev(X_points[k].cut[n].seg_lim[ns].R[-1], X_points[k].cut[n].seg_lim[ns].Z[-1]))								#float to convert from numpy.ndarray class
				if(psi1 < psi2):
					Segments[iSeg].R = X_points[k].cut[n].seg_lim[ns].R
					Segments[iSeg].Z = X_points[k].cut[n].seg_lim[ns].Z
				else:
					Segments[iSeg].R = X_points[k].cut[n].seg_lim[ns].R[::-1]
					Segments[iSeg].Z = X_points[k].cut[n].seg_lim[ns].Z[::-1]

				Segments[iSeg].p1_coord = [psi1, k, Segments[iSeg].nB, 1]
				Segments[iSeg].p2_coord = [psi2, k, Segments[iSeg].nB, 1]
				if(DEBUG > 1):	
					print("\tiSeg                      ",iSeg+1)
					print("\tSegments[iSeg].psimin_numc",Segments[iSeg].psimin_numc+1)
					print("\tSegments[iSeg].psimax_numc",Segments[iSeg].psimax_numc+1)
					print("\tSegments[iSeg].psimin     ",Segments[iSeg].psimin)
					print("\tSegments[iSeg].psimax     ",Segments[iSeg].psimax)
					print("\tSegments[iSeg].type       ",Segments[iSeg].type)
					print("\tSegments[iSeg].Xtype      ",Segments[iSeg].Xtype)
					print("\tSegments[iSeg].nX         ",Segments[iSeg].nX+1)
					print("\tSegments[iSeg].nB         ",Segments[iSeg].nB+1)
					print("\tSegments[iSeg].p1_coord   ",np.array(Segments[iSeg].p1_coord) + [0,1,1,0])
					print("\tSegments[iSeg].p2_coord   ",np.array(Segments[iSeg].p2_coord) + [0,1,1,0])

				self.Ax.plot(Segments[iSeg].R ,Segments[iSeg].Z, "m--", linewidth=3)
				self.Ax.text(Segments[iSeg].R[int(len(Segments[iSeg].R)/2)],Segments[iSeg].Z[int(len(Segments[iSeg].Z)/2)],"S_{:d}".format(iSeg+1))							

	for k in range (len(Frontiers)):
		for ns in range(Frontiers[k].nPt_lim-1):
			wh_min = np.where(np.append(Frontiers[k].Pt_lim[ns].psi, Frontiers[k].Pt_lim[ns+1].psi) == \
					 np.min(np.append(Frontiers[k].Pt_lim[ns].psi, Frontiers[k].Pt_lim[ns+1].psi)))[0][0]
			wh_max = np.where(np.append(Frontiers[k].Pt_lim[ns].psi, Frontiers[k].Pt_lim[ns+1].psi) == \
				     np.max(np.append(Frontiers[k].Pt_lim[ns].psi, Frontiers[k].Pt_lim[ns+1].psi)))[0][0]

			iSeg = Frontiers[k].seg_lim[ns].num
			Segments[iSeg].iSeg			= iSeg
			Segments[iSeg].type			= 2				 # poloidal border
			Segments[iSeg].Xtype		= 0
			Segments[iSeg].nX			= k
			Segments[iSeg].nB			= -1

			Segments[iSeg].psimin		= np.min(np.append(Frontiers[k].Pt_lim[ns].psi, Frontiers[k].Pt_lim[ns+1].psi))
			Segments[iSeg].psimin_numc	= Frontiers[k].Pt_lim[ns+wh_min].psi_numc

			Segments[iSeg].psimax		= np.max(np.append(Frontiers[k].Pt_lim[ns].psi, Frontiers[k].Pt_lim[ns+1].psi))
			Segments[iSeg].psimax_numc	= Frontiers[k].Pt_lim[ns+wh_max].psi_numc

			Segments[iSeg].R		= Frontiers[k].seg_lim[ns].R
			Segments[iSeg].Z		= Frontiers[k].seg_lim[ns].Z
			Segments[iSeg].p1_coord = [Segments[iSeg].psimin, k, Segments[iSeg].nB,2]
			Segments[iSeg].p2_coord = [Segments[iSeg].psimax, k, Segments[iSeg].nB,2]
			
			self.Ax.plot(Segments[iSeg].R ,Segments[iSeg].Z, "m-.", linewidth=3)
			self.Ax.text(Segments[iSeg].R[int(len(Segments[iSeg].R)/2)],Segments[iSeg].Z[int(len(Segments[iSeg].Z)/2)],"S_{:d}".format(iSeg+1))							

			if(DEBUG > 1):	
				print("\tk,ns                      ",k,ns)
				print("\tiSeg                      ",iSeg+1)
				print("\tSegments[iSeg].psimin_numc",Segments[iSeg].psimin_numc+1)
				print("\tSegments[iSeg].psimax_numc",Segments[iSeg].psimax_numc+1)
				print("\tSegments[iSeg].psimin     ",Segments[iSeg].psimin)
				print("\tSegments[iSeg].psimax     ",Segments[iSeg].psimax)
				print("\tSegments[iSeg].type       ",Segments[iSeg].type)
				print("\tSegments[iSeg].Xtype      ",Segments[iSeg].Xtype)
				print("\tSegments[iSeg].nX         ",Segments[iSeg].nX+1)
				print("\tSegments[iSeg].nB         ",Segments[iSeg].nB+1)
				print("\tSegments[iSeg].p1_coord   ",np.array(Segments[iSeg].p1_coord) + [0,1,1,0])
				print("\tSegments[iSeg].p2_coord   ",np.array(Segments[iSeg].p2_coord) + [0,1,1,0])


	self.Fig.canvas.draw()
	Config.Segments = Segments

	if(DEBUG > 0): print("find_segment_theta: completed\n")

	return True
