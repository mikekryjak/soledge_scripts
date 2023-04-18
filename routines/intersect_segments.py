import numpy 			as np
TOL = 1e-10

def intersect_segments(xs1, ys1, xe1, ye1, xs2, ys2, xe2, ye2, remove_same = True):

	xs1L = np.where(xs1 < xe1, xs1, xe1) - TOL
	xs1H = np.where(xs1 > xe1, xs1, xe1) + TOL

	ys1L = np.where(ys1 < ye1, ys1, ye1) - TOL
	ys1H = np.where(ys1 > ye1, ys1, ye1) + TOL

	xs2L = np.where(xs2 < xe2, xs2, xe2) - TOL
	xs2H = np.where(xs2 > xe2, xs2, xe2) + TOL

	ys2L = np.where(ys2 < ye2, ys2, ye2) - TOL
	ys2H = np.where(ys2 > ye2, ys2, ye2) + TOL

	xis  = np.empty(0, dtype='f8') 
	yis  = np.empty(0, dtype='f8')
	is1	 = np.empty(0, dtype='i4')
	is2	 = np.empty(0, dtype='i4')
	for k1 in range(len(xs1H)):
		k2 = np.where(np.logical_not((xs1H[k1] <  xs2L) | (xs1L[k1] >  xs2H) | (ys1H[k1] <  ys2L) | (ys1L[k1] >  ys2H)))[0]
		if(len(k2) > 0):
			a1 =   ye1[k1] - ys1[k1]
			b1 = -(xe1[k1] - xs1[k1])

			a2 =   ye2[k2] - ys2[k2]
			b2 = -(xe2[k2] - xs2[k2])

			Det = a1*b2 - b1*a2
			iNonPar = np.where(Det != 0.)[0]
			if(len(iNonPar) > 0):
				if(len(iNonPar) < len(k2)):
					a2  = a2[iNonPar]
					b2  = b2[iNonPar]
					c2  = c2[iNonPar]
					Det = Det[iNonPar]
					k2  = k2[iNonPar]

				c1 = xs1[k1]*a1 +ys1[k1]*b1
				c2 = xs2[k2]*a2 + ys2[k2]*b2

				xi = (c1*b2 - b1*c2)/Det
				yi = (a1*c2 - c1*a2)/Det
				
				InSeg = np.where(np.logical_not((xi < xs1L[k1]) | (xi > xs1H[k1]) | (yi < ys1L[k1]) | (yi > ys1H[k1]) | 
												(xi < xs2L[k2]) | (xi > xs2H[k2]) | (yi < ys2L[k2]) | (yi > ys2H[k2])))[0]

				if(len(InSeg) > 0):
					for l1 in InSeg:
						IsOut = True
						if(remove_same):
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

#	Intersection between a line (given at two points) and a segment (given at two extrema)

def intersect_lines_to_segments(xLs, yLs, xLe, yLe, xSs, ySs, xSe, ySe, remove_same = True):

	xSL = np.where(xSs < xSe, xSs, xSe) - TOL
	xSH = np.where(xSs > xSe, xSs, xSe) + TOL

	ySL = np.where(ySs < ySe, ySs, ySe) - TOL
	ySH = np.where(ySs > ySe, ySs, ySe) + TOL

	xis  = np.empty(0, dtype='f8') 
	yis  = np.empty(0, dtype='f8')
	iL	 = np.empty(0, dtype='i4')
	iS	 = np.empty(0, dtype='i4')
	for k1 in range(len(xLs)):
		a1 =   yLe[k1] - yLs[k1]
		b1 = -(xLe[k1] - xLs[k1])

		a2 =   ySe - ySs
		b2 = -(xSe - xSs)

		Det = a1*b2 - b1*a2
		k2 = np.where(Det != 0.)[0]
		if(len(k2) > 0):

			a2  = a2[k2]
			b2  = b2[k2]
			Det = Det[k2]

			c1 = xLs[k1]*a1 +yLs[k1]*b1
			c2 = xSs*a2 + ySs*b2

			xi = (c1*b2 - b1*c2)/Det
			yi = (a1*c2 - c1*a2)/Det
				
			InSeg = np.where(np.logical_not((xi < xSL[k2]) | (xi > xSH[k2]) | (yi < ySL[k2]) | (yi > ySH[k2])))[0]

			if(len(InSeg) > 0):
				for l1 in InSeg:
					IsOut = True
					if(remove_same):
						for l2 in range(len(xis)):
							if((abs(xis[l2] -xi[l1]) < TOL) and (abs(yis[l2] - yi[l1]) < TOL)):
								IsOut = False
								break
					if(IsOut):
						xis = np.append(xis, xi[l1])
						yis = np.append(yis, yi[l1])
						iL = np.append(iL, k1)
						iS = np.append(iS, k2[l1])
		

	return xis, yis, iL, iS


def intersect_line_to_segment(xs1, ys1, xe1, ye1, xs2, ys2, xe2, ye2):

	xs1L = min(xs1, xe1) - TOL
	xs1H = max(xs1, xe1) + TOL

	ys1L = min(ys1, ye1) - TOL
	ys1H = max(ys1, ye1) + TOL

	xs2L = min(xs2, xe2) - TOL
	xs2H = max(xs2, xe2) + TOL

	ys2L = min(ys2, ye2) - TOL
	ys2H = max(ys2, ye2) + TOL

	a1 =   ye1 - ys1
	b1 = -(xe1 - xs1)

	a2 =   ye2 - ys2
	b2 = -(xe2 - xs2)

	Det = a1*b2 - b1*a2
	if(Det !=0.):
		c1 = xs1*a1 +ys1*b1
		c2 = xs2*a2 + ys2*b2

		xi = (c1*b2 - b1*c2)/Det
		yi = (a1*c2 - c1*a2)/Det
		if(not ((xi < xs1L) or (xi > xs1H) or (yi < ys1L) or (yi > ys1H) or
				(xi < xs2L) or (xi > xs2H) or (yi < ys2L) or (yi > ys2H))): return xi,	yi, True
		
	return 0., 0., False