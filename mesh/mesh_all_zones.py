import types
import numpy 									as np
from routines.globals							import *
from routines.utils_walls						import set_XYlines_walls
from routines.contour_better					import contour_better
from routines.intersect_contour 				import intersect_contour
from mesh.follow_grad							import follow_grad
from mesh.optimized_mesh_zone_not_ortho			import optimized_mesh_zone_not_ortho
from mesh.conform_mesh							import conform_mesh
from interfaces.plot_and_ask					import plot_and_ask
from interfaces.progressbar						import ProgressBar

def mesh_all_zones(Root, Config, Optimization):

	if(DEBUG > 0): print("mesh_all_zones")

	MagZones		= Config.MagZones
	MagMegazones	= Config.MagMegazones
	MagPMegazones	= Config.MagPMegazones
	 
	numwait		= 0
	nMagZones		= len(MagZones)
	nMagPMegazones = len(MagPMegazones)
	
	for k in range(nMagPMegazones):
		if(MagPMegazones[k].ismeshed):
			numwait = numwait + len(MagPMegazones[k].refpoints.R)

	if(numwait == 0): return
	
	for k in range(len(MagZones)):
		mz	 = MagZones[k].mz
		pmz	 = MagZones[k].pmz
		if(MagMegazones[mz].ismeshed and MagPMegazones[pmz].ismeshed):
			nPsi			= len(MagMegazones[mz].refpoints.psi)
			nTheta  		= len(MagPMegazones[pmz].refpoints.R)
			MagZones[k].gridR	= np.empty((nPsi, nTheta), dtype='f8')
			MagZones[k].gridZ	= np.empty((nPsi, nTheta), dtype='f8')

			t1 = MagZones[k].Xtype_east
			t2 = MagZones[k].Xtype_west
			if(((t1 == 1) and (t2 == 1)) or MagPMegazones[pmz].ForceOrtho):
				MagZones[k].meshortho = True
			else:
				MagZones[k].meshortho = False

			if(MagZones[k].northaligned or MagZones[k].southaligned):
				for nz in MagPMegazones[pmz].list: MagZones[nz].meshortho = False

	for k in range(nMagPMegazones):
		MagPMegazones[k].meshortho = False
		for nz in MagPMegazones[k].list:
			if(MagZones[nz].meshortho): MagPMegazones[k].meshortho = True

	prog_bar = ProgressBar(Root, Label="meshing... (coffee break!)", Value=0.) 
		
	step = 0
	for k in range(nMagPMegazones):

		if(MagPMegazones[k].ismeshed and MagPMegazones[k].meshortho):				#By P. Innocente
			orthomeshchanged = False
			for nz in MagPMegazones[k].list:
				if(MagZones[nz].orthomeshchanged):	orthomeshchanged = True

		if(MagPMegazones[k].ismeshed and MagPMegazones[k].meshortho and orthomeshchanged):

#			Reference points for meshing Pmegazone (ie zones radially connected)

			R = MagPMegazones[k].refpoints.R
			Z = MagPMegazones[k].refpoints.Z

#			Psimin and psimax of Pmegazone

			psimin = MagZones[MagPMegazones[k].list[0]].pA.coord[0]
			psimax = MagZones[MagPMegazones[k].list[-1]].pB.coord[0]

			if(DEBUG > 1):	
				print("\tmegazone= ",k+1)
				print("\tpsimin= ",psimin)
				print("\tpsimmax= ",psimax)		

#			First step = mesh west side of Pmegazone

			step = step+1                 #first one

#			Put together segments of west side of Pmegazone, eliminating double points
			if(MagPMegazones[k].meshchanged):
				MagPMegazones[k].mesh		= [types.SimpleNamespace()]
				MagPMegazones[k].mesh[0].R = MagZones[MagPMegazones[k].list[0]].west.R
				MagPMegazones[k].mesh[0].Z = MagZones[MagPMegazones[k].list[0]].west.Z
				for k2 in range(1,len(MagPMegazones[k].list)):
					MagPMegazones[k].mesh[0].R = np.append(MagPMegazones[k].mesh[0].R, MagZones[MagPMegazones[k].list[k2]].west.R)
					MagPMegazones[k].mesh[0].Z = np.append(MagPMegazones[k].mesh[0].Z, MagZones[MagPMegazones[k].list[k2]].west.Z)

			c1			= types.SimpleNamespace()
			c1.arc		= [types.SimpleNamespace()]
			c1.arc[0].x = MagPMegazones[k].mesh[0].R
			c1.arc[0].y = MagPMegazones[k].mesh[0].Z
			x1			= c1.arc[0].x
			y1			= c1.arc[0].y
			d			= np.sqrt((x1[1:]-x1[:-1])**2+(y1[1:]-y1[:-1])**2)
			ii			= np.where(d != 0); ii = ii[0] + 1
			c1.arc[0].x = np.append(x1[0],x1[ii])
			c1.arc[0].y = np.append(y1[0],y1[ii])

#			Find intersections between psi contours of mesh and west side of Pmegazone

			for k3 in range(len(MagPMegazones[k].list)):				# Loop on zones in Pmegazone
				nz = MagPMegazones[k].list[k3]
				if(MagMegazones[MagZones[nz].mz].ismeshed):
					psilist = MagMegazones[MagZones[nz].mz].refpoints.psi
					if(DEBUG > 1):	print("\tmesh_zone_ortho.1 (all psi, first theta points): mz, nz =",k+1,nz+1)

					for k4 in range(1,len(psilist)-1):
						c2 = contour_better(Config.r2D, Config.z2D, Config.flux2D, [psilist[k4]])
						X = intersect_contour(c1,c2)
						if(len(X) > 1):
							print("Warning double point in zone ",nz)
							
						MagZones[nz].gridR[k4,0]	= X[0].x
						MagZones[nz].gridZ[k4,0]	= X[0].y
						
					MagZones[nz].gridR[0,0]		= MagZones[nz].west.R[0]				#i=k4=0...len(psilist)-1
					MagZones[nz].gridR[-1,0]	= MagZones[nz].west.R[-1]				#j=k2=0...len(R)-1
					MagZones[nz].gridZ[0,0]		= MagZones[nz].west.Z[0]
					MagZones[nz].gridZ[-1,0]	= MagZones[nz].west.Z[-1]

			prog_bar.Update(step/numwait)

#			Middle step: follow the gradient until you cross relevant psi contours
						
			psim = np.min(Config.flux2D)
			psiM = np.max(Config.flux2D)
			dpsi = psimax-psimin

			for k2 in range(1,len(R)-1):
				step = step+1
				if(MagPMegazones[k].meshchanged):
					MagPMegazones[k].mesh.append(types.SimpleNamespace())
					
					R1, Z1	= follow_grad(Config, R[k2], Z[k2], max(psimin-abs(dpsi)*0.05,psim), 5e-3)
					R2, Z2	= follow_grad(Config, R[k2], Z[k2], min(psimax+abs(dpsi)*0.05,psiM), 5e-3)
					Rm		= np.append(R1[:0:-1],R2)
					Zm		= np.append(Z1[:0:-1],Z2)
					MagPMegazones[k].mesh[k2].R = Rm
					MagPMegazones[k].mesh[k2].Z = Zm
				else:
					Rm = MagPMegazones[k].mesh[k2].R
					Zm = MagPMegazones[k].mesh[k2].Z
					
				d			= np.sqrt((Rm[1:]-Rm[:-1])**2+(Zm[1:]-Zm[:-1])**2)
				ii			= np.where(d != 0); ii = ii[0] + 1
				
				c1			= types.SimpleNamespace()
				c1.arc		= [types.SimpleNamespace()]
				c1.arc[0].x = np.append(Rm[0],Rm[ii])
				c1.arc[0].y = np.append(Zm[0],Zm[ii])

#				if(False):							#!!!!!! Temporaneo !!!!!
				for k3 in range(len(MagPMegazones[k].list)):
					nz = MagPMegazones[k].list[k3]
					if(MagZones[nz].meshortho):
						if(MagMegazones[MagZones[nz].mz].ismeshed):
							psilist = MagMegazones[MagZones[nz].mz].refpoints.psi
							if((DEBUG > 1) and (k2 == 1)):  print("\tmesh_zone_ortho.2 (all psi, 2:end-1 theta points): mz, nz =",k+1,nz+1)
							for k4 in range(len(psilist)):
								c2	= contour_better(Config.r2D, Config.z2D, Config.flux2D, [psilist[k4]])
								X	= intersect_contour(c1,c2)
								if(len(X) == 0):
									print("\tERROR: No intersection between contour")

									Xarrs,  Yarrs, lines = set_XYlines_walls(Root.Config)
									Xarrs.extend([c1.arc[0].x, c2.arc[0].x])
									Yarrs.extend([c1.arc[0].y, c2.arc[0].y])
									lines.extend(["r-", "g-"])
									LinesData = [Xarrs, Yarrs, lines]
									choice = plot_and_ask(Root, LinesData=LinesData, title="No contour intersection")
									return
									

								MagZones[nz].gridR[k4,k2] = X[0].x
								MagZones[nz].gridZ[k4,k2] = X[0].y
								
							MagZones[nz].orthomeshchanged = False	  			#By P. Innocente				 
								
				prog_bar.Update(step/numwait)

#			Last step = mesh east side of Pmegazone

			step = step+1														#last one

#			Put together segments of east side of Pmegazone, eliminating double points

			if(MagPMegazones[k].meshchanged):
				MagPMegazones[k].mesh.append(types.SimpleNamespace())
				MagPMegazones[k].mesh[-1].R = MagZones[MagPMegazones[k].list[0]].east.R
				MagPMegazones[k].mesh[-1].Z = MagZones[MagPMegazones[k].list[0]].east.Z
				for k2 in range(1,len(MagPMegazones[k].list)):
					MagPMegazones[k].mesh[-1].R	= np.append(MagPMegazones[k].mesh[-1].R, MagZones[MagPMegazones[k].list[k2]].east.R)
					MagPMegazones[k].mesh[-1].Z	= np.append(MagPMegazones[k].mesh[-1].Z, MagZones[MagPMegazones[k].list[k2]].east.Z)

			c1			= types.SimpleNamespace()
			c1.arc		= [types.SimpleNamespace()]
			c1.arc[0].x = MagPMegazones[k].mesh[len(R)-1].R
			c1.arc[0].y = MagPMegazones[k].mesh[len(R)-1].Z
			x1			= c1.arc[0].x
			y1			= c1.arc[0].y
			d			= np.sqrt((x1[1:]-x1[:-1])**2+(y1[1:]-y1[:-1])**2)
			ii			= np.where(d != 0); ii = ii[0] + 1
			c1.arc[0].x = np.append(x1[0],x1[ii])
			c1.arc[0].y = np.append(y1[0],y1[ii])

#			Find intersections between psi contours of mesh and east side of Pmegazone

			for k3 in range(len(MagPMegazones[k].list)):				 # Loop on zones in Pmegazone
				nz = MagPMegazones[k].list[k3]
				if(MagMegazones[MagZones[nz].mz].ismeshed):
					if(DEBUG > 1): print("\tmesh_zone_ortho.3 (all psi, last theta point): mz, nz =",k+1,nz+1)

					psilist = MagMegazones[MagZones[nz].mz].refpoints.psi
					for k4 in range(1,len(psilist)-1):
						c2	= contour_better(Config.r2D, Config.z2D, Config.flux2D, [psilist[k4]])
						X	= intersect_contour(c1,c2)
						MagZones[nz].gridR[k4,-1]	= X[0].x
						MagZones[nz].gridZ[k4,-1]	= X[0].y
						
					MagZones[nz].gridR[0,-1]	= MagZones[nz].east.R[0]				#i=k4=0...len(psilist)-1
					MagZones[nz].gridR[-1,-1]	= MagZones[nz].east.R[-1]				#j=k2=0...len(R)-1
					MagZones[nz].gridZ[0,-1]	= MagZones[nz].east.Z[0]
					MagZones[nz].gridZ[-1,-1]	= MagZones[nz].east.Z[-1]
					
			prog_bar.Update(step/numwait)

		MagPMegazones[k].meshchanged = False

#	 non orthogonal mesh

	for k in range(nMagPMegazones):
		if(MagPMegazones[k].ismeshed):
			n_last	= -1
			n_first	= -1
			for k1 in range(len(MagPMegazones[k].list)):
				nz = MagPMegazones[k].list[k1]
				
				if(MagZones[nz].meshortho):
					if(n_last == -1):
						n_last = k1
					else:
						if(MagZones[MagPMegazones[k].list[k1-1]].meshortho): n_last = k1

					if(n_first == -1): n_first = k1

			if(n_last > -1):
				containsortho = True
			else:
				containsortho = False

			if(not containsortho):
				meshstart	= MagPMegazones[k].refpoints.nz
				meshstartB	= MagPMegazones[k].refpoints.nzB				#1 north 2south

				num = np.where(MagPMegazones[k].list == meshstart)[0][0]
				if(meshstartB == 0):									#north
					refR = MagPMegazones[k].refpoints.R
					refZ = MagPMegazones[k].refpoints.Z
					for k1 in range(num+1,len(MagPMegazones[k].list)):
						nz = MagPMegazones[k].list[k1]
						if(MagMegazones[MagZones[nz].mz].ismeshed):
							if(MagZones[nz].orthomeshchanged):
								if(DEBUG > 1): print("\tmesh_zone_not_ortho.1 (last psi, all theta points): mz, nz, k1 =",k+1,nz+1,k1+1)

								optimized_mesh_zone_not_ortho(Root, Config, MagZones, MagMegazones, MagPMegazones, nz, refR ,refZ, 1, Optimization)
								refR = MagZones[nz].gridR[-1,:]
								refZ = MagZones[nz].gridZ[-1,:]
								MagZones[nz].orthomeshchanged = False
						else:
							break

					refR = MagPMegazones[k].refpoints.R
					refZ = MagPMegazones[k].refpoints.Z
					for k1 in range(num,-1,-1):
						nz = MagPMegazones[k].list[k1]
						if(MagMegazones[MagZones[nz].mz].ismeshed):
							if(MagZones[nz].orthomeshchanged):
								if(DEBUG > 1): print("\tmesh_zone_not_ortho.2 (last psi, all theta points): mz, nz, k1 =",k+1,nz+1,k1+1)
								optimized_mesh_zone_not_ortho(Root, Config, MagZones, MagMegazones, MagPMegazones, nz, refR, refZ, -1, Optimization)
								refR = MagZones[nz].gridR[0,:]
								refZ = MagZones[nz].gridZ[0,:]
								MagZones[nz].orthomeshchanged = False

						else:
							break

				else:															#1south  2north
					refR = MagPMegazones[k].refpoints.R
					refZ = MagPMegazones[k].refpoints.Z
					for k1 in range(num,len(MagPMegazones[k].list)):
						nz = MagPMegazones[k].list[k1]
						if(MagMegazones[MagZones[nz].mz].ismeshed):
							if(MagZones[nz].orthomeshchanged):
								if(DEBUG > 1): print("\tmesh_zone_not_ortho.3 (last psi, all theta points): mz, nz, k1 =",k+1,nz+1,k1+1)
								optimized_mesh_zone_not_ortho(Root, Config, MagZones, MagMegazones, MagPMegazones, nz, refR, refZ, 1, Optimization)
								refR = MagZones[nz].gridR[-1,:]
								refZ = MagZones[nz].gridZ[-1,:]
								MagZones[nz].orthomeshchanged = False
						else:
							break

					refR = MagPMegazones[k].refpoints.R
					refZ = MagPMegazones[k].refpoints.Z
					for k1 in range(num-1,-1,-1):
						nz = MagPMegazones[k].list[k1]
						if(MagMegazones[MagZones[nz].mz].ismeshed):
							if(MagZones[nz].orthomeshchanged):
								if(DEBUG > 1): print("\tmesh_zone_not_ortho.4 (last psi, all theta points): mz, nz, k1 =",k+1,nz+1,k1+1)
								optimized_mesh_zone_not_ortho(Root, Config, MagZones, MagMegazones, MagPMegazones, nz, refR, refZ, -1, Optimization)
								refR = MagZones[nz].gridR[0,:]
								refZ = MagZones[nz].gridZ[0,:]
								MagZones[nz].orthomeshchanged = False
						else:
							break


			else:
				if(n_last + 1 < len(MagPMegazones[k].list)):
					num		= n_last+1
					nz		= MagPMegazones[k].list[n_last]
					refR	= MagZones[nz].gridR[-1,:]
					refZ	= MagZones[nz].gridZ[-1,:]
					for k1 in range(num, len(MagPMegazones[k].list)):
						nz = MagPMegazones[k].list[k1]
						if(MagMegazones[MagZones[nz].mz].ismeshed):
							if(MagZones[nz].orthomeshchanged):
								if(DEBUG > 1): print("\tmesh_zone_not_ortho.5.1 (last psi, all theta points): mz, nz, k1 =",k+1,nz+1,k1+1)
								optimized_mesh_zone_not_ortho(Root, Config, MagZones, MagMegazones, MagPMegazones, nz ,refR, refZ, 1, Optimization)
								refR = MagZones[nz].gridR[-1,:]
								refZ = MagZones[nz].gridZ[-1,:]
								MagZones[nz].orthomeshchanged = False	
						else:
							break
						
				if(n_first > 0):
					num		= n_first-1
					nz		= MagPMegazones[k].list[n_first]
					refR	= MagZones[nz].gridR[0,:]
					refZ	= MagZones[nz].gridZ[0,:]
					for k1 in range(num, -1,-1):
						nz = MagPMegazones[k].list[k1]
						if(MagMegazones[MagZones[nz].mz].ismeshed):
							if(MagZones[nz].orthomeshchanged):
								if(DEBUG > 1): print("\tmesh_zone_not_ortho.5.2 (last psi, all theta points): mz, nz, k1 =",k+1,nz+1,k1+1)
								optimized_mesh_zone_not_ortho(Root, Config, MagZones, MagMegazones, MagPMegazones, nz ,refR, refZ, -1, Optimization)
								refR = MagZones[nz].gridR[0,:]
								refZ = MagZones[nz].gridZ[0,:]
								MagZones[nz].orthomeshchanged = False
						else:
							break

			step = step+len(MagPMegazones[k].refpoints.R)
			prog_bar.Update(step/numwait)


		MagPMegazones[k].meshchanged = False

	for k in range(len(MagZones)):
		mz	 = MagZones[k].mz
		pmz	 = MagZones[k].pmz
		if(MagMegazones[mz].ismeshed and MagPMegazones[pmz].ismeshed):
			MagZones[k].gridRc = 0.25*(MagZones[k].gridR[:-1,:-1] + MagZones[k].gridR[:-1,1:] + MagZones[k].gridR[1:,1:] + MagZones[k].gridR[1:,:-1])
			MagZones[k].gridZc = 0.25*(MagZones[k].gridZ[:-1,:-1] + MagZones[k].gridZ[:-1,1:] + MagZones[k].gridZ[1:,1:] + MagZones[k].gridZ[1:,:-1])

	MagGrid_OK = True
	for k in range(len(MagMegazones)):  MagGrid_OK = MagGrid_OK and MagMegazones[k].ismeshed
	for k in range(len(MagPMegazones)): MagGrid_OK = MagGrid_OK and MagPMegazones[k].ismeshed
		
	if(MagGrid_OK):	conform_mesh(Config)
	
	Root.set_MagGrid_OK(Config, MagGrid_OK)
	
	prog_bar = 0

	if(DEBUG > 0): print("mesh_all_zones: completed")
