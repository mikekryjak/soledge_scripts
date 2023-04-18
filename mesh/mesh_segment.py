from tkinter							import messagebox
from tkinter.messagebox 				import Message
import types
import numpy							as np
from scipy								import interpolate
from math 								import floor, sqrt
from interfaces.custom_mesh_interface	import custom_mesh_interface
from mesh.mesh_segment_linear			import mesh_segment_linear


def mesh_segment(root, Config, zone_sel, side_sel, sub_sel, 
					pt_num, refineType, refineSide, adjustMode, paramL, paramR, ForceOrtho, Alignement, psimin, psimax):

	MagZones		= Config.MagZones
	X_points		= Config.X_points
	MagMegazones	= Config.MagMegazones
	MagPMegazones	= Config.MagPMegazones
	OMPseg			= Config.OMPseg

	if(zone_sel < len(MagZones)):
		if(side_sel == 0):
			if(Alignement):
				R			= MagZones[zone_sel].subNorth[sub_sel].R
				Z			= MagZones[zone_sel].subNorth[sub_sel].Z
			else:
				R		= MagZones[zone_sel].north.R
				Z		= MagZones[zone_sel].north.Z

			MagZones[zone_sel].northaligned	= Alignement
			if(MagZones[zone_sel].Neighbour.north > -1): MagZones[MagZones[zone_sel].Neighbour.north].southaligned = Alignement

		elif(side_sel == 1):
			if(Alignement):
				R = MagZones[zone_sel].subSouth[sub_sel].R
				Z = MagZones[zone_sel].subSouth[sub_sel].Z
			else:		
				R = MagZones[zone_sel].south.R
				Z = MagZones[zone_sel].south.Z
				
			MagZones[zone_sel].southaligned	= Alignement
			if(MagZones[zone_sel].Neighbour.south > -1): MagZones[MagZones[zone_sel].Neighbour.south].northaligned = Alignement

		elif(side_sel == 2):
			R = MagZones[zone_sel].east.R
			Z = MagZones[zone_sel].east.Z
			
		elif(side_sel == 3):
			R = MagZones[zone_sel].west.R
			Z = MagZones[zone_sel].west.Z
	elif(zone_sel == len(MagZones)):						# Special segment for OMP
		R = OMPseg.R
		Z = OMPseg.Z
	else:
		messagebox.showerror("Mesh segment", "selected zone number is larger than Nzones+1\n and so does not corespond to anything.")
		return False

#	Define normalized coordinate along the segment to  mesh and eliminate double points
#	Compute coordinate
	dist = np.append(0.,np.cumsum(np.sqrt(np.diff(R)**2+np.diff(Z)**2)))

	ddist = np.diff(dist)											#	Eliminate double  points
	iNonDouble = np.append(0,np.where(ddist != 0.)[0]+1)
	R = R[iNonDouble]
	Z = Z[iNonDouble]
	dist = dist[iNonDouble]

	longueur = dist[-1]												#	Normalize
	dist = dist/longueur

#	Now mesh

	if(zone_sel < len(MagZones)):									# Normal segment

		if(side_sel >= 2):											#poloidal
			mz = MagZones[zone_sel].mz
			for k in range(len(MagMegazones[mz].list)):
				MagZones[MagMegazones[mz].list[k]].east.ismeshed = True
				MagZones[MagMegazones[mz].list[k]].west.ismeshed = True

#			Define pt, the normalized distribution of mesh points between 0 and 1

			if(refineType == 0):									#none
				pt = np.linspace(0,1,pt_num)
			elif(refineType == 1):									#linear
				pt = mesh_segment_linear(refineSide, adjustMode, pt_num, paramL, paramR, longueur)
			else:
				bezier_def	= False
				if(MagMegazones[mz].ismeshed):
					if((pt_num == MagMegazones[mz].refpoints.nPoints) and (MagMegazones[mz].refpoints.RefineType == 2)):
						xBezier = MagMegazones[mz].refpoints.xBezier
						yBezier = MagMegazones[mz].refpoints.yBezier
						bezier_def	= True

				if(not bezier_def):
					xBezier = np.array([], dtype='f8')
					yBezier = np.array([], dtype='f8')

				CustomMesh = custom_mesh_interface(root, Config, X_points, pt_num, xBezier, yBezier, R, Z, CurveType = refineType-2)							#Bezier or Spline
				if(len(CustomMesh.d) == 0): return False
				pt = CustomMesh.d

			#check
			if(len(pt) != pt_num):
				print(" len(pt)=",len(pt)," pt_num=",pt_num)
				messagebox.showerror("Mesh segment", "Unknow problem #1. Check code!")
				return False

			if(np.min(np.abs(pt[1:]-pt[:-1])) == 0):
				messagebox.showerror("Mesh segment", "Unknow problem #2. Check code!")
				return False

			f_psi = interpolate.RectBivariateSpline(Config.r2D[0,:], Config.z2D[:,0], Config.flux2D.T, s=0)

			MagMegazones[mz].refpoints			= types.SimpleNamespace()
			MagMegazones[mz].refpoints.R		= np.interp(pt, dist, R)
			MagMegazones[mz].refpoints.Z		= np.interp(pt, dist, Z)
			MagMegazones[mz].refpoints.psi		= f_psi.ev(MagMegazones[mz].refpoints.R, MagMegazones[mz].refpoints.Z)
			MagMegazones[mz].refpoints.psi[0]	= MagZones[zone_sel].pA.coord[0]
			MagMegazones[mz].refpoints.psi[-1]	= MagZones[zone_sel].pB.coord[0]
			MagMegazones[mz].ismeshed			= True
			for nz in range(len(MagMegazones[mz].list)):
				MagZones[MagMegazones[mz].list[nz]].orthomeshchanged = True
				MagZones[MagMegazones[mz].list[nz]].gridR			 = 0
				MagZones[MagMegazones[mz].list[nz]].gridZ			 = 0
			
			MagMegazones[mz].refpoints.nPoints		= pt_num							#Mesher parameters
			MagMegazones[mz].refpoints.RefineType	= refineType
			MagMegazones[mz].refpoints.RefineSide	= refineSide
			MagMegazones[mz].refpoints.AdjustMode	= adjustMode
			MagMegazones[mz].refpoints.ParamL		= paramL
			MagMegazones[mz].refpoints.ParamR		= paramR
			if(refineType == 2):
				MagMegazones[mz].refpoints.xBezier		= CustomMesh.xBezier
				MagMegazones[mz].refpoints.yBezier		= CustomMesh.yBezier


		if(side_sel <= 1): 					#radial
			mz = MagZones[zone_sel].pmz
			if(refineType == 0):									#none
				pt = np.linspace(0,1,pt_num)
			elif(refineType == 1):									#linear
				pt =  mesh_segment_linear(refineSide, adjustMode, pt_num, paramL, paramR, longueur)
			else:														# custom
				bezier_def	= False
				if(MagPMegazones[mz].ismeshed):
					if(((side_sel == 0) and (MagZones[zone_sel].northaligned)) or ((side_sel == 1) and (MagZones[zone_sel].southaligned))):
						if(((side_sel == 0) and MagZones[zone_sel].subNorth[sub_sel].ismeshed) or
						   ((side_sel == 1) and MagZones[zone_sel].subSouth[sub_sel].ismeshed)):
							if((pt_num == MagPMegazones[mz].subrefpoints[sub_sel].nPoints) and (MagPMegazones[mz].subrefpoints[sub_sel].RefineType == 2)):
								xBezier		= MagPMegazones[mz].refpoints.xBezier
								yBezier 	= MagPMegazones[mz].refpoints.yBezier
								bezier_def	= True
					else:
						if((pt_num == MagPMegazones[mz].refpoints.nPoints) and (MagPMegazones[mz].refpoints.RefineType == 2)):
							xBezier = MagPMegazones[mz].refpoints.xBezier
							yBezier = MagPMegazones[mz].refpoints.yBezier
							bezier_def	= True

				if(not bezier_def):
					xBezier = np.array([], dtype='f8')
					yBezier = np.array([], dtype='f8')

				CustomMesh = custom_mesh_interface(root, Config, X_points, pt_num, xBezier, yBezier, R, Z, CurveType = refineType-2)														#exponential
				pt = CustomMesh.d

																			#check
			if(len(pt) != pt_num):
				messagebox.showerror("Mesh segment", "Unknow problem #3. Check code!")
				return False

			if(np.min(np.abs(pt[1:]-pt[:-1])) == 0):
				messagebox.showerror("Mesh segment", "Unknow problem #4. Check code!")
				return False

			mz = MagZones[zone_sel].pmz
			if(Alignement):
				if(side_sel == 0):
					MagZones[zone_sel].subNorth[sub_sel].ismeshed = True
					doit = True
					for n in range(2): doit = doit and MagZones[zone_sel].subNorth[n].ismeshed

				if(side_sel == 1):
					MagZones[zone_sel].subSouth[sub_sel].ismeshed = True
					doit = True
					for n in range(2): doit = doit and MagZones[zone_sel].subSouth[n].ismeshed

				for k in range(len(MagPMegazones[mz].list)):
					MagZones[MagPMegazones[mz].list[k]].north.ismeshed = doit
					MagZones[MagPMegazones[mz].list[k]].south.ismeshed = doit

				MagPMegazones[mz].ismeshed = doit
			else:
				for k in range(len(MagPMegazones[mz].list)):
					MagZones[MagPMegazones[mz].list[k]].north.ismeshed = True
					MagZones[MagPMegazones[mz].list[k]].south.ismeshed = True


				MagPMegazones[mz].ismeshed = True

#			for k in range(len(MagPMegazones[mz].list)):
#				MagZones[MagPMegazones[mz].list[k]].northaligned = Alignement
#				MagZones[MagPMegazones[mz].list[k]].southaligned = Alignement

			MagPMegazones[mz].isaligned = Alignement
			if(Alignement):	MagPMegazones[mz].iAlignWall = MagZones[zone_sel].iAlignWall

			MagPMegazones[mz].refpoints = types.SimpleNamespace()		
			if(Alignement):
				MagPMegazones[mz].isaligned = Alignement
				if(not hasattr(MagPMegazones[mz],'subrefpoints')):
					MagPMegazones[mz].subrefpoints	= [types.SimpleNamespace(), types.SimpleNamespace()]
					for n in range(2):
						MagPMegazones[mz].subrefpoints[n].R = np.array([])
						MagPMegazones[mz].subrefpoints[n].Z = np.array([])

				MagPMegazones[mz].subrefpoints[sub_sel].R 	= np.interp(pt, dist, R)
				MagPMegazones[mz].subrefpoints[sub_sel].Z	= np.interp(pt, dist, Z)
				MagPMegazones[mz].refpoints.R				= np.append(MagPMegazones[mz].subrefpoints[0].R[:-1], MagPMegazones[mz].subrefpoints[1].R)
				MagPMegazones[mz].refpoints.Z				= np.append(MagPMegazones[mz].subrefpoints[0].Z[:-1], MagPMegazones[mz].subrefpoints[1].Z)
				MagPMegazones[mz].align_psimin				= psimin
				MagPMegazones[mz].align_psimax				= psimax

			else:
				MagPMegazones[mz].refpoints.R				= np.interp(pt, dist, R)
				MagPMegazones[mz].refpoints.Z				= np.interp(pt, dist, Z)

			MagPMegazones[mz].refpoints.nz		= zone_sel
			MagPMegazones[mz].refpoints.nzB		= side_sel
			MagPMegazones[mz].meshchanged		= True
			for nz in range(len(MagPMegazones[mz].list)):
				MagZones[MagPMegazones[mz].list[nz]].orthomeshchanged = True
				MagZones[MagPMegazones[mz].list[nz]].gridR			= 0
				MagZones[MagPMegazones[mz].list[nz]].gridZ			= 0
			
			if(Alignement):
				MagPMegazones[mz].subrefpoints[sub_sel].nPoints			= pt_num							#Mesher parameters
				MagPMegazones[mz].subrefpoints[sub_sel].RefineType		= refineType
				MagPMegazones[mz].subrefpoints[sub_sel].RefineSide		= refineSide
				MagPMegazones[mz].subrefpoints[sub_sel].AdjustMode		= adjustMode
				MagPMegazones[mz].subrefpoints[sub_sel].ParamL			= paramL
				MagPMegazones[mz].subrefpoints[sub_sel].ParamR			= paramR
				MagPMegazones[mz].ForceOrtho							= ForceOrtho
				if(refineType == 2):
					MagPMegazones[mz].subrefpoints[sub_sel].xBezier	= CustomMesh.xBezier
					MagPMegazones[mz].subrefpoints[sub_sel].yBezier	= CustomMesh.yBezier
			else:
				MagPMegazones[mz].refpoints.nPoints		= pt_num											#Mesher parameters
				MagPMegazones[mz].refpoints.RefineType	= refineType
				MagPMegazones[mz].refpoints.RefineSide	= refineSide
				MagPMegazones[mz].refpoints.AdjustMode	= adjustMode
				MagPMegazones[mz].refpoints.ParamL		= paramL
				MagPMegazones[mz].refpoints.ParamR		= paramR
				MagPMegazones[mz].ForceOrtho			= ForceOrtho
				if(refineType == 2):
					MagPMegazones[mz].refpoints.xBezier		= CustomMesh.xBezier
					MagPMegazones[mz].refpoints.yBezier		= CustomMesh.yBezier

	elif(zone_sel == len(MagZones)):
#		Special OMP segment
#		There are possibly points that will have to be forced at some given positions
#		(the points at the crossing between the OMP segment and the psi contour of zones on the way)
#		First double check that there are enough points to mesh this segment
		Npts_min = len(OMPseg.Intersec_R)
		if(pt_num < Npts_min):									# There are not enough points -> complain and don't do anything
			messagebox.showerror("Mesh segment", "Not enough mesh points specified to mesh special OMP segment\n"+ \
													 "One needs to have at least NzonesOMP+1 points,\n"+ \
													 "where NzonesOMP is the number of zones in the OMP.\n"+ \
													 "In this specific case, one needs at least {:d} points.".format(Npts_min))
			return False
		else:													# there are enough points -> let's do
			Rforce = OMPseg.Intersec_R[1:-1]					# Position of points to force
			Zforce = OMPseg.Intersec_Z[1:-1]

#			 Locate them on the (R,Z) segment to define their normalized coordinate
#			To do so we look at the change of sign of some dot product
			d_force = np.empty(len(Rforce), dtype='f8')
			for k in range(len(Rforce)):
				v1R = Rforce[k]-R[:-1]
				v1Z = Zforce[k]-Z[:-1]
				v2R = Rforce[k]-R[1:]
				v2Z = Zforce[k]-Z[1:]
				dotprod = v1R*v2R + v1Z*v2Z
				iseg = np.where(dotprod <= 0)[0]
				if (len(iseg) != 1):
					messagebox.showerror("Mesh segment", "Problem with locating point to force on OMP special segment.")
					return False
				else:
					d_force[k] = dist[iseg] + sqrt((Rforce[k]-R[iseg])**2+(Zforce[k]-Z[iseg])**2)/longueur

#			Now switch according to the type of mesh requested
			if(refineType == 0) or (refineType == 1):			#none or linear
#				Define distribution of points a priori
				if(refineType == 0): 							#none
					pt=np.linspace(0, 1, pt_num)
				else:											#linear
					pt = mesh_segment_linear(refineSide, adjustMode, pt_num, paramL, paramR, longueur)

#				Correct mesh points to have one at the interface between zones
				for iforce in range(len(d_force)):
					iclosest = np.argmin(np.abs(pt-d_force[iforce]))			# Find closest mesh point
					iclosest = min(max(iclosest,1),pt_num-2)				# Correct if we found the first one or the last one (forbidden)

					pt[iclosest] = d_force[iforce]							# Move
			else:															# custom
				bezier_def	= False
				if(MagMegazones[mz].ismeshed):
					if((pt_num == MagMegazones[mz].refpoints.nPoints) and (MagMegazones[mz].refpoints.RefineType == 2)):
						xBezier = MagMegazones[mz].refpoints.xBezier
						yBezier = MagMegazones[mz].refpoints.yBezier
						bezier_def	= True
				
				if(not bezier_def):
					xBezier = np.array([], dtype='f8')
					yBezier = np.array([], dtype='f8')

				CustomMesh = custom_mesh_interface(root, Config, X_points, pt_num, xBezier, yBezier, R, Z, CurveType = refineType-2, d_force=d_force)														#exponential
				if(len(CustomMesh.d) == 0): return False
				pt = CustomMesh.d

#			No correction to make as this was done directly in the custom mesh interface

			if(len(pt) != pt_num):											# Check that pt is correct
				messagebox.showerror("Mesh segment", "Unknow problem #5. Check code!")
				return False

			if(np.min(np.abs(np.diff(pt))) == 0):
				messagebox.showerror("Mesh segment", "Unknow problem #6. Check code!")
				return False

			NzonesOMP = len(OMPseg.zonelist)								# Loop on OMP zones to distribute them the ref points
			for izoneOMP in range(NzonesOMP):
				izone = OMPseg.zonelist[izoneOMP]
				mz = MagZones[izone].mz

	#			Extract ptsub list of points from pt for this zone
				if(izoneOMP == 0):	dsubmin = 0
				else:				dsubmin = d_force[izoneOMP-1]

				if (izoneOMP == NzonesOMP-1):	dsubmax = 1
				else:							dsubmax = d_force[izoneOMP]

				ptsub = pt[np.where((pt >= dsubmin) & (pt <= dsubmax))[0]]

	#			 Define mesh points by interpolation along R and Z, and store
				f_psi = interpolate.RectBivariateSpline(Config.r2D[0,:], Config.z2D[:,0], Config.flux2D.T, s=0)

				MagMegazones[mz].refpoints			= types.SimpleNamespace()
				MagMegazones[mz].refpoints.R		= np.interp(ptsub, dist, R)
				MagMegazones[mz].refpoints.Z		= np.interp(ptsub, dist, Z)
				MagMegazones[mz].refpoints.psi		= f_psi.ev(MagMegazones[mz].refpoints.R, MagMegazones[mz].refpoints.Z)
				MagMegazones[mz].refpoints.psi[0]	= MagZones[izone].pA.coord[0]
				MagMegazones[mz].refpoints.psi[-1]	= MagZones[izone].pB.coord[0]		
				MagMegazones[mz].ismeshed			= True
				for iz in MagMegazones[mz].list:
					MagZones[iz].orthomeshchanged = True
					MagZones[iz].gridR			  = 0
					MagZones[iz].gridZ			  = 0
					MagZones[iz].east.ismeshed	  = True								# Flag relevant zones as meshed
					MagZones[iz].west.ismeshed	  = True

				MagMegazones[mz].refpoints.nPoints		= len(ptsub)
				MagMegazones[mz].refpoints.RefineType	= 0
				MagMegazones[mz].refpoints.RefineSide	= 0
				MagMegazones[mz].refpoints.AdjustMode	= 0
				MagMegazones[mz].refpoints.ParamL		= 0.
				MagMegazones[mz].refpoints.ParamR		= 0.

#			Check if at this stage the special OMP segment should be considered as meshed
#			(I have left here the generic test function for that, but if we get here it should return true)
			OMPseg.ismeshed = True
			for imz in range(len(OMPseg.mzlist)):
				mz = OMPseg.mzlist[imz]
				if(not MagMegazones[mz].ismeshed): OMPseg.ismeshed = False

			OMPseg.refpoints = types.SimpleNamespace()
			OMPseg.refpoints.nPoints	= pt_num											#Mesher parameters
			OMPseg.refpoints.RefineType	= refineType
			OMPseg.refpoints.RefineSide	= refineSide
			OMPseg.refpoints.AdjustMode	= adjustMode
			OMPseg.refpoints.ParamL		= paramL
			OMPseg.refpoints.ParamR		= paramR
			if(refineType >= 2):
				OMPseg.refpoints.xBezier		= CustomMesh.xBezier
				OMPseg.refpoints.yBezier		= CustomMesh.yBezier

#			end izoneOMP

	Segments_OK	= True
	for k in range(len(MagMegazones)):
		Segments_OK = Segments_OK and MagMegazones[k].ismeshed

	for k in range(len(MagPMegazones)):
		Segments_OK = Segments_OK and MagPMegazones[k].ismeshed
	
	root.set_Segments_OK(Config, Segments_OK)

	return True