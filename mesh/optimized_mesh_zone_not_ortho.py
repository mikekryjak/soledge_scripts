import types
import numpy						as np
from routines.part_contour_per		import part_contour_per
from routines.part_contour			import part_contour
from routines.intersect_contour		import intersect_contour
from routines.contour_better		import contour_better
from mesh.optim_mesh				import optim_mesh
from interfaces.plot_and_ask		import plot_and_ask
from routines.globals				import DEBUG


def PolyArea(x,y):
	return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

def optimized_mesh_zone_not_ortho(Root, Config, MagZones, MagMegazones, MagPMegazones, nz, RefR, RefZ, direction, Optimization):

	if(DEBUG > 0):	print("optimized_mesh_zone_not_ortho")

#	Builds a normalized coordinate along the reference points (based on distance between points
	nRefR = len(RefR)
	dist  = np.append(0., np.cumsum(np.sqrt(np.diff(RefR)**2+np.diff(RefZ)**2)))
	dist = dist/dist[-1]

#	Determines the "directeness" of the quadrangle (S,W,N,E)
	vx	 = MagZones[nz].south.R[1] - MagZones[nz].south.R[0]
	vy	 = MagZones[nz].south.Z[1] - MagZones[nz].south.Z[0]
	vx1  = MagZones[nz].west.R[1]  - MagZones[nz].west.R[0]
	vy1  = MagZones[nz].west.Z[1]  - MagZones[nz].west.Z[0]
	qual = (vx1*vy-vy1*vx)/(np.sqrt(vx**2+vy**2)*np.sqrt(vx1**2+vy1**2))
	sens = np.sign(qual)

#	Reference values of psi for meshing radial points (ie on W or E side)
	psilist = MagMegazones[MagZones[nz].mz].refpoints.psi
	pmz	    = MagZones[nz].pmz

	if(direction == 1):						#south to north
#		Define normalized coordinate along the points of south side
		Rm = MagZones[nz].south.R
		Zm = MagZones[nz].south.Z
		dm = np.append(0., np.cumsum(np.sqrt((Rm[1:]-Rm[:-1])**2+(Zm[1:]-Zm[:-1])**2)))
		dm = dm/dm[-1]
		
#		Distribute grid point along south side so that they have the same normalized coordinates as the ref points 

		MagZones[nz].gridR[0,:] = np.interp(dist, dm, Rm)
		MagZones[nz].gridZ[0,:] = np.interp(dist, dm, Zm)

#		Make sure the extremeties of the grid on the south side correspond to the zone corners

		MagZones[nz].gridR[0,0] = MagZones[nz].south.R[0]
		MagZones[nz].gridZ[0,0] = MagZones[nz].south.Z[0]
		MagZones[nz].gridR[0, -1] = MagZones[nz].south.R[-1]
		MagZones[nz].gridZ[0, -1] = MagZones[nz].south.Z[-1]

#		Now loop on iso-ipsi grid lines
		for i in range(1,len(psilist)-1):

# 			Iso-psi contour based on reference points for this megazone
			c1			= contour_better(Config.r2D, Config.z2D, Config.flux2D, [psilist[i]])

#			Contour for east and west sides of zone	
			cE 	   		= types.SimpleNamespace()
			cE.arc 		= [types.SimpleNamespace()]
			cE.arc[0].x = MagZones[nz].east.R
			cE.arc[0].y = MagZones[nz].east.Z
			
			cW 	   		= types.SimpleNamespace()
			cW.arc 		= [types.SimpleNamespace()]
			cW.arc[0].x = MagZones[nz].west.R
			cW.arc[0].y = MagZones[nz].west.Z

#			Intersections between east and west side and reference iso-psi contour
			XE			= intersect_contour(c1,cE)
			XW			= intersect_contour(c1,cW)
			cnum		= XE[0].arc1
			
			cin 		= types.SimpleNamespace()
			cin.x		= c1.arc[cnum].x
			cin.y		= c1.arc[cnum].y

			p1	 		= types.SimpleNamespace()
			p2	 		= types.SimpleNamespace()
			p1.x		= XW[0].x
			p1.y		= XW[0].y
			p2.x		= XE[0].x
			p2.y		= XE[0].y

#			Extract piece of iso-psi  contour between east and west sides (2 possible directions)
			cout		= [part_contour_per(cin,p1,p2,1), part_contour_per(cin,p1,p2,-1)]
			
#			Choose arc among the 2 possible ones
			Areas = np.zeros(2, dtype='f8')
			for na in range(2):
				Areas[na]	= PolyArea(np.append(MagZones[nz].gridR[i-1,:], cout[na].x[::-1]), np.append(MagZones[nz].gridZ[i-1,:], cout[na].y[::-1]))

			na = np.where(Areas == np.min(Areas))[0][0]

			Rm = cout[na].x
			Zm = cout[na].y

			if(not MagPMegazones[pmz].isaligned):
#				In case of none aligned zone, we simply interpolate the points along the chosen contour based on their distribution in reference points
				Rin = MagZones[nz].gridR[i-1,:]
				Zin = MagZones[nz].gridZ[i-1,:]
				din = np.append(0., np.cumsum(np.sqrt((Rin[1:]-Rin[:-1])**2+(Zin[1:]-Zin[:-1])**2)))
				din = din/din[-1]
				
				dout,Rout,Zout = optim_mesh(Root, din, Rin, Zin, Rm, Zm, sens, Optimization)
				MagZones[nz].gridR[i,:] = Rout
				MagZones[nz].gridZ[i,:] = Zout
			else:
#				Case of aligned mesh: we need to stick to the wall between align_psimin and align_psimax
				if((psilist[i] <= MagPMegazones[pmz].align_psimax) and (psilist[i] >= MagPMegazones[pmz].align_psimin)):
#					We are in the psi range where alignment is requested
#					align --> two parts
#					First part

					iAlignWall = MagPMegazones[pmz].iAlignWall
					Rwall = Config.Walls[iAlignWall].Rwall 			# BEWARE, THIS WILL NOT WORK IF ALIGNEMNT WITH FLOATTING OBJECT IS REQUESTED
					Zwall = Config.Walls[iAlignWall].Zwall

#					Compute normalized curvilinear coordinate of grid points at previous ipsi for first sub segment

					Rin1 = MagZones[nz].gridR[i-1,:len(MagPMegazones[pmz].subrefpoints[0].R)]
					Zin1 = MagZones[nz].gridZ[i-1,:len(MagPMegazones[pmz].subrefpoints[0].R)]
					din1 = np.append(0., np.cumsum(np.sqrt((Rin1[1:]-Rin1[:-1])**2+(Zin1[1:]-Zin1[:-1])**2)))
					din1 = din1/din1[-1]

#					Get intersection between relevant wall and considered segment
					
					c1 			= types.SimpleNamespace()
					c1.arc		= [types.SimpleNamespace()]
					c1.arc[0].x = Rwall
					c1.arc[0].y = Zwall
					c2 			= types.SimpleNamespace()
					c2.arc		= [types.SimpleNamespace()]
					c2.arc[0].x = Rm
					c2.arc[0].y = Zm
					X			= intersect_contour(c1,c2)								# intersection between relevant wall and considered segment
					cin 		= types.SimpleNamespace()
					cin.x		= Rm
					cin.y		= Zm
					if(len(X) > 1):
#						If more than 1 intersection, retain only the one the closest to the start of the segment (or the end if the start is in the wall)						d = np.zeros(len(X), dtype='f8')
						d = np.zeros(len(X), dtype='f8')
						if(Config.Walls[iAlignWall].WallPath.contains_point([c2.arc[0].x[0], c2.arc[0].y[0]])):
							for k in range(len(X)):
								d[k]  = np.sqrt((X[k].x- c2.arc[0].x[0])**2 + (X[k].y-c2.arc[0].y[0])**2)
						else:
							for k in range(len(X)):
								d[k] = np.sqrt((X[k].x- c2.arc[0].x[-1])**2+(X[k].y-c2.arc[0].y[-1])**2)

						a 		= np.argmin(d)
						X1 		= [types.SimpleNamespace()]
						X1[0].x	= X[a].x
						X1[0].y	= X[a].y
						X		= X1

					if(len(X) == 1):
						p1.x = X[0].x
						p1.y = X[0].y
					else:
						print("\tERROR: No intersection between aligned psi segment and wall")

						Xarrs  = [Rwall, Rm, MagZones[nz].east.R, MagZones[nz].west.R]
						Yarrs  = [Zwall, Zm, MagZones[nz].east.Z, MagZones[nz].west.Z]
						lines  = ["k-", "r-", "g-", "g-"]
						LinesData = [Xarrs, Yarrs, lines]
						choice = plot_and_ask(Root, LinesData=LinesData, title="Allignement problem: psi=red, east&west=green")
						return

#					Cut initial segment in 2 pieces at the intersection
					cout	= part_contour(cin,p1)
					Rm1		= cout.arc[0].x
					Zm1		= cout.arc[0].y
					Rm2		= cout.arc[1].x
					Zm2		= cout.arc[1].y

#					Mesh the first part of the segment before the intersection
					dout, Rout1, Zout1 = optim_mesh(Root, din1, Rin1, Zin1, Rm1, Zm1, sens, Optimization)

#					Now deal with part on the other side of the wall
#					Compute normalized curvilinear coordinate of grid points at previous ipsi for second sub segment

#					Second part

					Rin2 = MagZones[nz].gridR[i-1, len(MagPMegazones[pmz].subrefpoints[0].R)-1:]
					Zin2 = MagZones[nz].gridZ[i-1, len(MagPMegazones[pmz].subrefpoints[0].R)-1:]
					din2 = np.append(0., np.cumsum(np.sqrt((Rin2[1:]-Rin2[:-1])**2+(Zin2[1:]-Zin2[:-1])**2)))
					din2 = din2/din2[-1]

#					Mesh the second part of the segment after the intersection

					dout, Rout2, Zout2 = optim_mesh(Root, din2, Rin2, Zin2, Rm2, Zm2, sens, Optimization)
					Rout	= np.append(Rout1[:-1], Rout2)
					Zout	= np.append(Zout1[:-1], Zout2)

#					Store
					MagZones[nz].gridR[i,:] = Rout
					MagZones[nz].gridZ[i,:] = Zout
				else:
#					We are in the psi range where alignment is not requested
#					%normal
					Rin = MagZones[nz].gridR[i-1,:]
					Zin = MagZones[nz].gridZ[i-1,:]
					din = np.append(0., np.cumsum(np.sqrt((Rin[1:]-Rin[:-1])**2+(Zin[1:]-Zin[:-1])**2)))
					din = din/din[-1]
					
					dout, Rout, Zout = optim_mesh(Root, din, Rin, Zin, Rm, Zm, sens, Optimization)
					MagZones[nz].gridR[i,:] = Rout
					MagZones[nz].gridZ[i,:] = Zout

#		last one
		i	= len(psilist)-1
		Rm	= MagZones[nz].north.R
		Zm	= MagZones[nz].north.Z
		if(not MagPMegazones[pmz].isaligned):
			Rin = MagZones[nz].gridR[i-1,:]
			Zin = MagZones[nz].gridZ[i-1,:]
			din = np.append(0., np.cumsum(np.sqrt((Rin[1:]-Rin[:-1])**2+(Zin[1:]-Zin[:-1])**2)))
			din = din/din[-1]
			
			dout, Rout, Zout = optim_mesh(Root, din, Rin, Zin, Rm, Zm, sens, Optimization)
			MagZones[nz].gridR[i,:] = Rout
			MagZones[nz].gridZ[i,:] = Zout
		else:
			if((psilist[i] <= MagPMegazones[pmz].align_psimax) and (psilist[i] >= MagPMegazones[pmz].align_psimin)):
#				align --> two parts
				print("\talign_psimin, align_psimax, psilist[i] =",MagPMegazones[pmz].align_psimin, MagPMegazones[pmz].align_psimax,psilist[i])
				iAlignWall = MagPMegazones[pmz].iAlignWall
				Rwall = Config.Walls[iAlignWall].Rwall 			# BEWARE, THIS WILL NOT WORK IF ALIGNEMNT WITH FLOATTING OBJECT IS REQUESTED
				Zwall = Config.Walls[iAlignWall].Zwall

				Rin1 = MagZones[nz].gridR[i-1,:len(MagPMegazones[pmz].subrefpoints[0].R)]
				Zin1 = MagZones[nz].gridZ[i-1,:len(MagPMegazones[pmz].subrefpoints[0].R)]
				din1 = np.append(0., np.cumsum(np.sqrt((Rin1[1:]-Rin1[:-1])**2+(Zin1[1:]-Zin1[:-1])**2)))
				din1 = din1/din1[-1]
				
				c1 			= types.SimpleNamespace()
				c1.arc 		= [types.SimpleNamespace()]
				c1.arc[0].x = Rwall
				c1.arc[0].y = Zwall
				
				c2 			= types.SimpleNamespace()
				c2.arc 		= [types.SimpleNamespace()]
				c2.arc[0].x = Rm
				c2.arc[0].y = Zm
				X			= intersect_contour(c1,c2)
				cin 		= types.SimpleNamespace()
				cin.x		= Rm
				cin.y		= Zm
				if(len(X) > 1):
					d = np.zeros(len(X), dtype='f8')
					if(Config.Walls[iAlignWall].WallPath.contains_point([c2.arc[0].x[0],c2.arc[0].y[0]])):
						for k in range(len(X)):
							d[k]=np.sqrt((X[k].x- c2.arc[0].x[0])**2+(X[k].y-c2.arc[0].y[0])**2)
					else:
						for k in range(len(X)):
							d[k] = np.sqrt((X[k].x- c2.arc[0].x[-1])**2+(X[k].y-c2.arc[0].y[-1])**2)

					a 		= np.argmin(d)
					X		= intersect_contour(c1,c2)
					X1		= [types.SimpleNamespace()]
					X1[0].x	= X[a].x
					X1[0].y	= X[a].y
					X		= X1

				if(len(X) == 1):
					p1.x = X[0].x
					p1.y = X[0].y
					
				cout	= part_contour(cin,p1)
				Rm1		= cout.arc[0].x
				Zm1		= cout.arc[0].y
				Rm2		= cout.arc[1].x
				Zm2		= cout.arc[1].y
				dout, Rout1, Zout1 = optim_mesh(Root, din1, Rin1, Zin1, Rm1 ,Zm1, sens, Optimization) 
				Rin2 = MagZones[nz].gridR[i-1,len(MagPMegazones[pmz].subrefpoints[0].R)-1:]
				Zin2 = MagZones[nz].gridZ[i-1,len(MagPMegazones[pmz].subrefpoints[0].R)-1:]
				din2 = np.append(0., np.cumsum(np.sqrt((Rin2[1:]-Rin2[:-1])**2+(Zin2[1:]-Zin2[:-1])**2)))
				din2 = din2/din2[-1]
				
				dout, Rout2, Zout2 = optim_mesh(Root, din2, Rin2, Zin2, Rm2, Zm2, sens, Optimization)
				Rout				 = np.append(Rout1[:-1],Rout2)
				Zout				 = np.append(Zout1[:-1],Zout2)
				MagZones[nz].gridR[i,:] = Rout
				MagZones[nz].gridZ[i,:] = Zout
			else:
#				normal
				Rin = MagZones[nz].gridR[i-1,:]
				Zin = MagZones[nz].gridZ[i-1,:]
				din = np.append(0., np.cumsum(np.sqrt((Rin[1:]-Rin[:-1])**2+(Zin[1:]-Zin[:-1])**2)))
				din = din/din[-1]
				
				dout, Rout, Zout = optim_mesh(Root, din, Rin, Zin, Rm, Zm, sens, Optimization)
				MagZones[nz].gridR[i,:] = Rout
				MagZones[nz].gridZ[i,:] = Zout


		MagZones[nz].gridR[-1, 0]	= MagZones[nz].north.R[0]
		MagZones[nz].gridZ[-1, 0]	= MagZones[nz].north.Z[0]
		MagZones[nz].gridR[-1, -1] = MagZones[nz].north.R[-1]
		MagZones[nz].gridZ[-1, -1] = MagZones[nz].north.Z[-1]

	else: 				#north to south
#		last one
		Rm	= MagZones[nz].north.R
		Zm	= MagZones[nz].north.Z
		dm = np.append(0., np.cumsum(np.sqrt((Rm[1:]-Rm[:-1])**2+(Zm[1:]-Zm[:-1])**2)))
		dm = dm/dm[-1]
		
		MagZones[nz].gridR[-1,:] = np.interp(dist, dm, Rm)
		MagZones[nz].gridZ[-1,:] = np.interp(dist, dm, Zm)

		for i in range(len(psilist)-2,0,-1):
			c1 = contour_better(Config.r2D,Config.z2D,Config.flux2D,[psilist[i]])
			
			cE 			= types.SimpleNamespace()
			cE.arc 		= [types.SimpleNamespace()]
			cE.arc[0].x = MagZones[nz].east.R
			cE.arc[0].y = MagZones[nz].east.Z

			cW 			= types.SimpleNamespace()
			cW.arc 		= [types.SimpleNamespace()]
			cW.arc[0].x = MagZones[nz].west.R
			cW.arc[0].y = MagZones[nz].west.Z
			XE 			= intersect_contour(c1,cE)
			XW			= intersect_contour(c1,cW)
			cnum		= XE[0].arc1
			
			cin			= types.SimpleNamespace()
			cin.x		= c1.arc[cnum].x
			cin.y		= c1.arc[cnum].y
			
			p1			= types.SimpleNamespace()
			p1.x		= XW[0].x
			p1.y		= XW[0].y

			p2			= types.SimpleNamespace()
			p2.x		= XE[0].x
			p2.y		= XE[0].y

			cout		= [part_contour_per(cin,p1,p2,1), part_contour_per(cin,p1,p2,-1)]
			
#			choose arc
			Areas = np.zeros(2, dtype='f8')
			for na in range(2):
				Areas[na] = PolyArea(np.append(MagZones[nz].gridR[i+1,:],cout[na].x[::-1]), np.append(MagZones[nz].gridZ[i+1,:],cout[na].y[::-1]))

			na = np.argmin(Areas)
			Rm = cout[na].x												#R & Z of flux contour zone
			Zm = cout[na].y

			if(not MagPMegazones[pmz].isaligned):
				Rin = MagZones[nz].gridR[i+1,:]
				Zin = MagZones[nz].gridZ[i+1,:]
				din = np.append(0., np.cumsum(np.sqrt((Rin[1:]-Rin[:-1])**2+(Zin[1:]-Zin[:-1])**2)))
				din = din/din[-1]
				
				dout, Rout, Zout = optim_mesh(Root, din, Rin, Zin, Rm, Zm, -sens, Optimization)
				MagZones[nz].gridR[i,:] = Rout
				MagZones[nz].gridZ[i,:] = Zout
			else:
				if((psilist[i] <= MagPMegazones[pmz].align_psimax) and (psilist[i] >= MagPMegazones[pmz].align_psimin)):
#					align --> two parts

					iAlignWall = MagPMegazones[pmz].iAlignWall
					Rwall = Config.Walls[iAlignWall].Rwall 			# BEWARE, THIS WILL NOT WORK IF ALIGNEMNT WITH FLOATTING OBJECT IS REQUESTED
					Zwall = Config.Walls[iAlignWall].Zwall

					Rin1 = MagZones[nz].gridR[i+1, :len(MagPMegazones[pmz].subrefpoints[0].R)]
					Zin1 = MagZones[nz].gridZ[i+1, :len(MagPMegazones[pmz].subrefpoints[0].R)]
					din1 = np.append(0., np.cumsum(np.sqrt((Rin1[1:]-Rin1[:-1])**2+(Zin1[1:]-Zin1[:-1])**2)))
					din1 	= din1/din1[-1]
					
					c1			= types.SimpleNamespace()
					c1.arc		= [types.SimpleNamespace()]
					c1.arc[0].x = Rwall
					c1.arc[0].y = Zwall
					
					c2			= types.SimpleNamespace()
					c2.arc		= [types.SimpleNamespace()]
					c2.arc[0].x = Rm
					c2.arc[0].y = Zm
					X = intersect_contour(c1,c2)							 #Intersection between flux contour and wall
					
					cin			= types.SimpleNamespace()
					cin.x		= Rm
					cin.y		= Zm
					if(len(X) > 1):
#						If more than 1 intersection, retain only the one the closest to the start of the segment (or the end if the start is in the wall)						d = np.zeros(len(X), dtype='f8')
						d = np.zeros(len(X), dtype='f8')
						if(Config.Walls[iAlignWall].WallPath.contains_point([c2.arc[0].x[0],c2.arc[0].y[0]])):
							for k in range(len(X)):
								d[k] = np.sqrt((X[k].x- c2.arc[0].x[0])**2+(X[k].y-c2.arc[0].y[0])**2)
						else:
							for k in range(len(X)):
								d[k] = np.sqrt((X[k].x- c2.arc[0].x[-1])**2+(X[k].y-c2.arc[0].y[-1])**2)
					
						a 		= np.argmin(d)
						X1	     = [types.SimpleNamespace()]
						X1[0].x  = X[a].x
						X1[0].y  = X[a].y
						X	     = X1

					if(len(X) == 1):
						p1.x = X[0].x
						p1.y = X[0].y

					cout = part_contour(cin,p1)
					Rm1 = cout.arc[0].x
					Zm1 = cout.arc[0].y
					Rm2 = cout.arc[1].x
					Zm2 = cout.arc[1].y
					dout, Rout1, Zout1 =optim_mesh(Root, din1, Rin1, Zin1, Rm1, Zm1, -sens, Optimization)
					Rin2 = MagZones[nz].gridR[i+1, len(MagPMegazones[pmz].subrefpoints[0].R)-1:]
					Zin2 = MagZones[nz].gridZ[i+1, len(MagPMegazones[pmz].subrefpoints[0].R)-1:]
					din2 = np.append(0., np.cumsum(np.sqrt((Rin2[1:]-Rin2[:-1])**2+(Zin2[1:]-Zin2[:-1])**2)))
					din2 = din2/din2[-1]
					
					dout, Rout2, Zout2 = optim_mesh(Root, din2,Rin2,Zin2,Rm2,Zm2,-sens, Optimization)
					Rout	= np.append(Rout1[:-1],Rout2)
					Zout	= np.append(Zout1[:-1],Zout2)
					MagZones[nz].gridR[i,:] = Rout
					MagZones[nz].gridZ[i,:] = Zout
				else:
#					normal
					Rin = MagZones[nz].gridR[i+1,:]
					Zin = MagZones[nz].gridZ[i+1,:]
					din = np.append(0., np.cumsum(np.sqrt((Rin[1:]-Rin[:-1])**2+(Zin[1:]-Zin[:-1])**2)))
					din = din/din[-1]
					
					dout, Rout, Zout = optim_mesh(Root, din,Rin,Zin,Rm,Zm,-sens, Optimization)
					MagZones[nz].gridR[i,:] = Rout
					MagZones[nz].gridZ[i,:] = Zout

#		first one
		i	= 0
		Rm	= MagZones[nz].south.R
		Zm	= MagZones[nz].south.Z
		if( not MagPMegazones[pmz].isaligned):
			Rin = MagZones[nz].gridR[i+1,:]
			Zin = MagZones[nz].gridZ[i+1,:]
			din = np.append(0., np.cumsum(np.sqrt((Rin[1:]-Rin[:-1])**2+(Zin[1:]-Zin[:-1])**2)))
			din = din/din[-1]
			
			dout, Rout, Zout = optim_mesh(Root, din,Rin,Zin,Rm,Zm,-sens, Optimization)
			MagZones[nz].gridR[i,:] = Rout
			MagZones[nz].gridZ[i,:] = Zout
		else:
			if((psilist[i] <= MagPMegazones[pmz].align_psimax) and (psilist[i] >= MagPMegazones[pmz].align_psimin)):
#				align --> two parts

				iAlignWall = MagPMegazones[pmz].iAlignWall
				Rwall = Config.Walls[iAlignWall].Rwall 			# BEWARE, THIS WILL NOT WORK IF ALIGNEMNT WITH FLOATTING OBJECT IS REQUESTED
				Zwall = Config.Walls[iAlignWall].Zwall

				Rin1 = MagZones[nz].gridR[i+1, :len(MagPMegazones[pmz].subrefpoints[0].R)]
				Zin1 = MagZones[nz].gridZ[i+1, :len(MagPMegazones[pmz].subrefpoints[0].R)]
				din1 = np.append(0., np.cumsum(np.sqrt((Rin1[1:]-Rin1[:-1])**2+(Zin1[1:]-Zin1[:-1])**2)))
				din1		= din1/din1[-1]
				
				c1			= types.SimpleNamespace()
				c1.arc		= [types.SimpleNamespace()]
				c1.arc[0].x = Rwall
				c1.arc[0].y = Zwall

				c2			= types.SimpleNamespace()
				c2.arc		= [types.SimpleNamespace()]
				c2.arc[0].x = Rm
				c2.arc[0].y = Zm
				X			= intersect_contour(c1,c2)
				
				cin			= types.SimpleNamespace()
				cin.x = Rm
				cin.y = Zm
				if(len(X) > 1):
					d = np.zeros(len(X), dtype='f8')
					if(Config.Walls[iAlignWall].WallPath.contains_point([c2.arc[0].x[0],c2.arc[0].y[0]])):
						for k in range(len(X)):
							d[k] = np.sqrt((X[k].x - c2.arc[0].x[0])**2+(X[k].y-c2.arc[0].y[0])**2)
					else:
						for k in range(len(X)):
							d = np.sqrt((X[k].x - c2.arc[0].x[-1])**2+(X[k].y-c2.arc[0].y[-1])**2)
				
					a 		= np.argmin(d)
					X1 	    = [types.SimpleNamespace()]
					X1[0].x = X[a].x
					X1[0].y = X[a].y
					X	    = X1
					
				if(len(X) == 1):
					p1.x = X[0].x
					p1.y = X[0].y

				cout = part_contour(cin,p1)
				Rm1 = cout.arc[0].x
				Zm1 = cout.arc[0].y
				Rm2 = cout.arc[1].x
				Zm2 = cout.arc[1].y
				dout, Rout1, Zout1 = optim_mesh(Root, din1, Rin1, Zin1, Rm1, Zm1, -sens, Optimization)
				Rin2 = MagZones[nz].gridR[i+1, len(MagPMegazones[pmz].subrefpoints[0].R)-1:]
				Zin2 = MagZones[nz].gridZ[i+1, len(MagPMegazones[pmz].subrefpoints[0].R)-1:]
				din2 = np.append(0., np.cumsum(np.sqrt((Rin2[1:]-Rin2[:-1])**2+(Zin2[1:]-Zin2[:-1])**2)))
				din2 = din2/din2[-1]
				
				dout,Rout2,Zout2 = optim_mesh(Root, din2, Rin2, Zin2, Rm2, Zm2, -sens, Optimization)
				Rout	= np.append(Rout1[:-1],Rout2)
				Zout	= np.append(Zout1[:-1],Zout2)
				MagZones[nz].gridR[i,:] = Rout
				MagZones[nz].gridZ[i,:] = Zout
			else:
#				normal
				Rin = MagZones[nz].gridR[i+1,:]
				Zin = MagZones[nz].gridZ[i+1,:]
				din = np.append(0., np.cumsum(np.sqrt((Rin[1:]-Rin[:-1])**2+(Zin[1:]-Zin[:-1])**2)))
				din = din/din[-1]
				
				dout, Rout, Zout = optim_mesh(Root, din, Rin, Zin, Rm, Zm, -sens, Optimization)
				MagZones[nz].gridR[i,:] = Rout
				MagZones[nz].gridZ[i,:] = Zout

		MagZones[nz].gridR[0,0]  = MagZones[nz].south.R[0]
		MagZones[nz].gridZ[0,0]  = MagZones[nz].south.Z[0]
		MagZones[nz].gridR[0,-1] = MagZones[nz].south.R[-1]
		MagZones[nz].gridZ[0,-1] = MagZones[nz].south.Z[-1]

	if(DEBUG > 0):	print("optimized_mesh_zone_not_ortho: completed")

	return
