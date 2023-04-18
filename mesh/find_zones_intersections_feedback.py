import types
import numpy as np
from routines.intersect_contour			import intersect_contour
from mesh.compute_mesh_intersections	import compute_mesh_intersections


def find_zones_intersections_feedback(line, Config):

	MagZones	= Config.MagZones

	c1		= types.SimpleNamespace()
	c1.arc	= [types.SimpleNamespace()]
	c2		= types.SimpleNamespace()
	c2.arc	= [types.SimpleNamespace()]
	nz		= np.empty((0), dtype='i4')
	r12	= np.empty((0,2), dtype='f8')
	z12	= np.empty((0,2), dtype='f8')
	
	c1.arc[0].x = line[:,0]
	c1.arc[0].y = line[:,1]
	for k in range(len(MagZones)):
		c2.arc[0].x = MagZones[k].gridR[0,:]				#south border
		c2.arc[0].y = MagZones[k].gridZ[0,:]	
		X1  = intersect_contour(c1,c2)
		if(len(X1)>0):
			c2.arc[0].x = MagZones[k].gridR[-1,:]			#North border
			c2.arc[0].y = MagZones[k].gridZ[-1,:]	
			X2  = intersect_contour(c1,c2)
			if(len(X2)>0):
				nz = np.append(nz,k)
				r12 = np.append(r12, [[X1[0].x, X2[0].x]],axis=0)
				z12 = np.append(z12, [[X1[0].y, X2[0].y]],axis=0)


	Inters		= types.SimpleNamespace()					#intersections
	Inters.Flux	= types.SimpleNamespace()					#Flux direction
	if(len(nz)== 0):
		Inters.Flux.nz = np.empty((0))
		return Inters
	

	d12 = np.sqrt((r12[:,:] -line[0,0])**2 + (z12[:,:] -line[0,1])**2)
	
	dir = np.where(d12[:,0] < d12[:,1], 0, 1)   			#South to north or north to south
	ii	= np.arange(len(dir))
	d	= d12[ii,dir[ii]]
	a   = np.argsort(d)
	d12	= d12 - d[a[0]]

	Inters.Flux.d12			= np.empty(d12.shape, dtype='f8')
	Inters.Flux.r12			= np.empty(r12.shape, dtype='f8')
	Inters.Flux.z12			= np.empty(z12.shape, dtype='f8')
	
	Inters.Flux.nz  		= nz[a]
	Inters.Flux.d12[ii,0]	= d12[a[ii],dir[a[ii]]]
	Inters.Flux.d12[ii,1]	= d12[a[ii],1-dir[a[ii]]]
	Inters.Flux.r12[ii,0]	= r12[a[ii],dir[a[ii]]] 
	Inters.Flux.r12[ii,1]	= r12[a[ii],1-dir[a[ii]]] 
	Inters.Flux.z12[ii,0]	= z12[a[ii],dir[a[ii]]] 
	Inters.Flux.z12[ii,1]	= z12[a[ii],1-dir[a[ii]]]
	Inters.Flux.dir			= dir[a[ii]]										#0=south to north, 1=north to south

	Inters.sel				= False
	Inters.Flux.Lengths		= compute_mesh_intersections(Config, Inters.Flux)	#Compute lengths for conversion between L, X, Rho_v, Rho_F

	return Inters