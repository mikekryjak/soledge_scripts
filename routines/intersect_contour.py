import types
import numpy 			as np
#from shapely.geometry 	import asLineString
from numpy				import asarray
TOL = 1e-6

def intersect_contour(c1, c2):
	x = []
	for i in range (len(c1.arc)):
		for j in range(len(c2.arc)):
			if(len(c1.arc[i].x) < len(c2.arc[j].x)):
				xi, yi, is1, is2 = intersect_2contours(c1.arc[i].x, c1.arc[i].y, c2.arc[j].x, c2.arc[j].y)
			else:
				xi, yi, is2, is1 = intersect_2contours(c2.arc[j].x, c2.arc[j].y, c1.arc[i].x, c1.arc[i].y)


			for k in range(len(xi)):
				x.append(types.SimpleNamespace())
				x[-1].x	   = xi[k]
				x[-1].y	   = yi[k]
				x[-1].arc1 = i
				x[-1].arc2 = j
				x[-1].is1  = is1[k]
				x[-1].is2  = is2[k]

	"""
#	check custom intersection routine
	XX = intersect_contour1(c1, c2)

	if(len(XX) != len(x)):
		print("error in intersect_contour")
		print("len(XX)=", len(XX)," len(x) =", len(x))
		for k in range(len(XX)):
			print("XX[k].x   =",XX[k].x)
			print("XX[k].y   =",XX[k].y)
			print("XX[k].arc1=",XX[k].arc1)
			print("XX[k].arc2=",XX[k].arc2)
			print("c1.arc[XX[k].arc1].x=",c1.arc[XX[k].arc1].x)
			print("c1.arc[XX[k].arc1].y=",c1.arc[XX[k].arc1].y)
			print("c2.arc[XX[k].arc2].x=",c2.arc[XX[k].arc2].x)
			print("c2.arc[XX[k].arc2].y=",c2.arc[XX[k].arc2].y)
		for k in range(len(x)):
			print("x[k].x   =",x[k].x)
			print("x[k].y   =",x[k].y)
			print("x[k].arc1=",x[k].arc1)
			print("x[k].arc2=",x[k].arc2)
			print("x[k].is1 =",x[k].is1)
			print("x[k].is2 =",x[k].is2)
			print("c1.arc[x[k].arc1].x[x[k].is1:x[k].is1+2]=",c1.arc[x[k].arc1].x[x[k].is1:x[k].is1+2])
			print("c1.arc[x[k].arc1].y[x[k].is1:x[k].is1+2]=",c1.arc[x[k].arc1].y[x[k].is1:x[k].is1+2])
			print("c2.arc[x[k].arc2].x[x[k].is2:x[k].is2+2]=",c2.arc[x[k].arc2].x[x[k].is2:x[k].is2+2])
			print("c2.arc[x[k].arc2].y[x[k].is2:x[k].is2+2]=",c2.arc[x[k].arc2].y[x[k].is2:x[k].is2+2])
		exit()
	else:
		for k in range(len(XX)):
			if((abs(XX[k].x - x[k].x) > TOL) or (abs(XX[k].y - x[k].y) > TOL)):
				print("XX[k].x=",XX[k].x, " x[k].x=",x[k].x)
				print("XX[k].y=",XX[k].y,  "x[k].y=",x[k].y)
				exit()
	"""

	return x


def intersect_contour_shapely(c1, c2):
	x = []
	for i in range (len(c1.arc)):
		for j in range(len(c2.arc)):
			xy1   		= np.empty((len(c1.arc[i].x),2), dtype='f8')
			xy1[:,0]	= c1.arc[i].x
			xy1[:,1]	= c1.arc[i].y
			ls1			= asLineString(xy1)
			
			xy2   		= np.empty((len(c2.arc[j].x),2), dtype='f8')
			xy2[:,0]	= c2.arc[j].x
			xy2[:,1]	= c2.arc[j].y
			ls2			= asLineString(xy2)
			xy			= asarray(ls2.intersection(ls1))
			if(len(xy.shape) == 2):
				for k in range(xy.shape[0]):
					x.append(types.SimpleNamespace())
					x[-1].x	   = xy[k,0]
					x[-1].y	   = xy[k,1]
					x[-1].arc1 = i
					x[-1].arc2 = j
			elif(xy.shape[0] > 0):
				x.append(types.SimpleNamespace())
				x[-1].x	   = xy[0]
				x[-1].y	   = xy[1]
				x[-1].arc1 = i
				x[-1].arc2 = j	
	return x

def intersect_2contours(xs1, ys1, xs2, ys2):

	xs1L = np.where(xs1[:-1] < xs1[1:], xs1[:-1], xs1[1:]) - TOL
	xs1H = np.where(xs1[:-1] > xs1[1:], xs1[:-1], xs1[1:]) + TOL

	ys1L = np.where(ys1[:-1] < ys1[1:], ys1[:-1], ys1[1:]) - TOL
	ys1H = np.where(ys1[:-1] > ys1[1:], ys1[:-1], ys1[1:]) + TOL

	xs2L = np.where(xs2[:-1] < xs2[1:], xs2[:-1], xs2[1:]) - TOL
	xs2H = np.where(xs2[:-1] > xs2[1:], xs2[:-1], xs2[1:]) + TOL

	ys2L = np.where(ys2[:-1] < ys2[1:], ys2[:-1], ys2[1:]) - TOL
	ys2H = np.where(ys2[:-1] > ys2[1:], ys2[:-1], ys2[1:]) + TOL

	xis  = np.empty(0, dtype='f8') 
	yis  = np.empty(0, dtype='f8')
	is1	 = np.empty(0, dtype='i4')
	is2	 = np.empty(0, dtype='i4')
	for k1 in range(len(xs1H)):
		k2 = np.where(np.logical_not((xs1H[k1] <  xs2L) | (xs1L[k1] >  xs2H) | (ys1H[k1] <  ys2L) | (ys1L[k1] >  ys2H)))[0]
		if(len(k2) > 0):
			a1 =   ys1[k1+1] - ys1[k1]
			b1 = -(xs1[k1+1] - xs1[k1])
			c1 = xs1[k1]*a1 +ys1[k1]*b1

			a2 =   ys2[k2+1] - ys2[k2]
			b2 = -(xs2[k2+1] - xs2[k2])
			c2 = xs2[k2]*a2 + ys2[k2]*b2

			Det = a1*b2 - b1*a2
			iNonPar = np.where(Det != 0.)[0]
			if(len(iNonPar) > 0):
				if(len(iNonPar) < len(k2)):
					a2  = a2[iNonPar]
					b2  = b2[iNonPar]
					c2  = c2[iNonPar]
					Det = Det[iNonPar]
					k2  = k2[iNonPar]

				xi = (c1*b2 - b1*c2)/Det
				yi = (a1*c2 - c1*a2)/Det
				
				InSeg = np.where(np.logical_not((xi < xs1L[k1]) | (xi > xs1H[k1]) | (yi < ys1L[k1]) | (yi > ys1H[k1]) | 
												(xi < xs2L[k2]) | (xi > xs2H[k2]) | (yi < ys2L[k2]) | (yi > ys2H[k2])))[0]
				"""
				print("k1         =",k1)
				print("len(InSeg)=",len(InSeg))
				print("xs1       =",xs1[k1:k1+2])
				print("ys1       =",ys1[k1:k1+2])
				print("xs2       =",xs2[k2], xs2[k2+1])
				print("ys2       =",ys2[k2], ys2[k2+1])
				print("xi        =",xi)
				print("yi        =",yi)
				"""

				if(len(InSeg) > 0):
					for l1 in InSeg:
						IsOut = True
						for l2 in range(len(xis)):
							if((abs(xis[l2] -xi[l1]) < TOL) and (abs(yis[l2] - yi[l1]) < TOL)):
								IsOut = False
								break
						if(IsOut):
							xis = np.append(xis, xi[l1])
							yis = np.append(yis, yi[l1])
							is1 = np.append(is1, k1)
							is2 = np.append(is2, k2[l1])
		

	return xis, yis, is1, is2