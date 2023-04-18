import numpy				as np
import numpy.matlib			as mat
import scipy.interpolate	as interp
from math					import floor
from routines.globals		import DEBUG

def generate_soledge_xz(Config):

	if(DEBUG > 0): print("generate_soledge_xz")

	Zones		= Config.Zones
	Megazones	= Config.Megazones
	PMegazones	= Config.PMegazones
	
	nZones		= len(Zones)
	nMegazones	= len(Megazones)
	for k in range(nZones):
		Zones[k].isx_defined = False
		Zones[k].isz_defined = False


	f_psi	= interp.RectBivariateSpline(Config.r2D[0,:], Config.z2D[:,0], Config.flux2D.T)									#RectBivariateSpline wants [x,y] order
	minx	=  1e10
	maxx	= -1e10
	for k in range(nZones):
		if(not Zones[k].isx_defined):
			R		= Zones[k].gridR[:,0]
			Z		= Zones[k].gridZ[:,0]
			psi		= f_psi.ev(R, Z)
			minx_	= psi.min()
			maxx_	= psi.max()
			if(minx > minx_):	minx = minx_
			if(maxx < maxx_):	maxx = maxx_

			mz = Zones[k].mz
			for k1 in range(len(Megazones[mz].list)):
				nz		= Megazones[mz].list[k1]
				Zones[nz].xb = mat.repmat(psi,Zones[nz].gridR.shape[1],1).T
				Zones[nz].x	 = (Zones[nz].xb[:-1,:-1] + Zones[nz].xb[1:,:-1])*0.5
				Zones[nz].isx_defined = True

	for k in range(nZones):
		Zones[k].x	= (Zones[k].x  - minx)/(maxx - minx)
		Zones[k].xb	= (Zones[k].xb - minx)/(maxx - minx)


	mzin = np.empty(0, dtype='i4')
	while(len(mzin) != nMegazones):

		l	= np.zeros(nMegazones - len(mzin), dtype='i4')
		l2	= np.zeros(nMegazones - len(mzin), dtype='i4')
		k1 = 0
		for k in range(nMegazones):
			if(len(np.where(k == mzin)[0]) == 0):
				l[k1]	= len(Megazones[k].list)
				l2[k1]	= k
				k1	   += 1

		d	 = np.where(l == l.max())[0]
		mz	 = l2[d[0]]
		mzin = np.append(mzin,mz)

		containsz = False
		for k in range(len(Megazones[mz].list)):
			containsz = containsz or Zones[Megazones[mz].list[k]].isz_defined

		if(not containsz):
			dist1 = np.empty(0, dtype='f8')
			for k in range(len(Megazones[mz].list)):
				nz	  = Megazones[mz].list[k]
				m2	  = floor(Zones[nz].gridR.shape[0]/2)
				dist1 = np.append(dist1, np.sqrt((Zones[nz].gridR[m2, 1:] - Zones[nz].gridR[m2 ,:-1])**2 + (Zones[nz].gridZ[m2, 1:] - Zones[nz].gridZ[m2, :-1])**2))

			dist1 = np.append(0.,np.cumsum(dist1))

			dist = dist1/dist1[-1]

#			broadcast perp
			for k in range(len(Megazones[mz].list)):
				nz	= Megazones[mz].list[k]
				pmz	= Zones[nz].pmz
				for k1 in range(len(PMegazones[pmz].list)):
					nz2 = PMegazones[pmz].list[k1]
					if(not Zones[nz2].isz_defined):

						Zones[nz2].zb = np.dot(np.ones((Zones[nz2].gridZ.shape[0],1)),
											   dist[:Zones[nz2].gridZ.shape[1]].reshape((1,Zones[nz2].gridZ.shape[1])))
						Zones[nz2].isz_defined = True

				dist = dist[Zones[nz2].gridZ.shape[1]-1:]

		else:

#			begin
			tomesh = np.empty(0, dtype='i4')
			for k1 in range(len(Megazones[mz].list)):
				nz2 = Megazones[mz].list[k1]
				if(not Zones[nz2].isz_defined):
					tomesh=np.append(tomesh,nz2)
				else:
					beg = nz2
					break
				
			if(len(tomesh) > 0):
				dist1 = np.empty(0, dtype='f8')
				for k2 in range(len(tomesh)):
					nz	  = tomesh[k2]
					m2	  = floor(Zones[nz].gridR.shape[0]/2)
					dist1 = np.append(dist1, np.sqrt((Zones[nz].gridR[m2, 1:] - Zones[nz].gridR[m2 ,:-1])**2 + (Zones[nz].gridZ[m2, 1:] - Zones[nz].gridZ[m2, :-1])**2))
				
				dist1 = np.append(0.,np.cumsum(dist1))
	
				dist = dist1/dist1[-1]
				if(len(dist) >= 2):
					dist = dist/(dist[-1]-dist[-2])*(Zones[beg].zb[0,1] -Zones[beg].zb[0,0])			#rescale
					dist = dist-dist[-1] + Zones[beg].zb[0,0]
				else:
					dist = dist*Zones[beg].zb[0,0]
					
#				broadcast perp
				for k2 in range(len(tomesh)):
					nz  = tomesh[k2]
					pmz = Zones[nz].pmz
					for k1 in range(len(PMegazones[pmz].list)):
						nz2	  = PMegazones[pmz].list[k1]
						if(not Zones[nz2].isz_defined):

							Zones[nz2].zb = np.dot(np.ones((Zones[nz2].gridZ.shape[0],1)),
												   dist[:Zones[nz2].gridZ.shape[1]].reshape((1,Zones[nz2].gridZ.shape[1])))
							Zones[nz2].isz_defined = True
	
					dist = dist[Zones[nz2].gridZ.shape[1]-1:]


			tomesh = np.empty(0, dtype='i4')
			for k1 in range(len(Megazones[mz].list)-1,-1,-1):
				nz2 = Megazones[mz].list[k1]
				if(not Zones[nz2].isz_defined):
					tomesh = np.append(tomesh,nz2)
				else:
					beg	 = nz2
					break

			if(len(tomesh) > 0):
				tomesh = tomesh[::-1]
				dist1 = np.empty(0, dtype='f8')
				for k2 in range(len(tomesh)):
					nz	  = tomesh[k2]
					m2	  = floor(Zones[nz].gridR.shape[0]/2)
					dist1 = np.append(dist1, np.sqrt((Zones[nz].gridR[m2, 1:] - Zones[nz].gridR[m2 ,:-1])**2 + (Zones[nz].gridZ[m2, 1:] - Zones[nz].gridZ[m2, :-1])**2))
				
				dist1 = np.append(0.,np.cumsum(dist1))
	
				dist = dist1/dist1[-1]
				if(len(dist) >= 2):
					dist = dist/(dist[1]-dist[0])*(Zones[beg].zb[0,-1] -Zones[beg].zb[0,-2])			#rescale
					dist = dist + Zones[beg].zb[0,-1]
				else:
					dist = Zones[beg].zb[0,-1]+dist*(1.-Zones[beg].zb[0,-1])

#				broadcast perp
				for k2 in range(len(tomesh)):
					nz	= tomesh[k2]
					pmz = Zones[nz].pmz
					for k1 in range(len(PMegazones[pmz].list)):
						nz2 = PMegazones[pmz].list[k1]
						if(not Zones[nz2].isz_defined):
							Zones[nz2].zb = np.dot(np.ones((Zones[nz2].gridZ.shape[0],1)),
												   dist[:Zones[nz2].gridZ.shape[1]].reshape((1,Zones[nz2].gridZ.shape[1])))
							Zones[nz2].isz_defined = True
	
					dist = dist[Zones[nz2].gridZ.shape[1]-1:]

	for k in range(nZones):
		Zones[k].z=(Zones[k].zb[:-1,:-1]+Zones[k].zb[:-1,1:])*0.5

	if(DEBUG > 0): print("generate_soledge_xz: Completed")

	return
