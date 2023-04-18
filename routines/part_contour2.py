import types
import numpy				as np
from routines.part_contour	import part_contour

def part_contour2(Cin, p1, p2):

	C1 = part_contour(Cin, p1)

	Cin2= types.SimpleNamespace()
	C	= types.SimpleNamespace()
	
	d	= np.sqrt((Cin.x-p2.x)**2+(Cin.y-p2.y)**2)
	ind = np.where(d == d.min()); ind = ind[0][0]
	if(ind > C1.ind):
		Cin2.x	= C1.arc[1].x
		Cin2.y	= C1.arc[1].y
		C2		= part_contour(Cin2, p2)
		C.x		= C2.arc[0].x
		C.y		= C2.arc[0].y
	else:
		if(ind == C1.ind):
			C.x	= np.append(p1.x, p2.x)
			C.y	= np.append(p1.y, p2.y)
		else:
			Cin2.x	= C1.arc[0].x
			Cin2.y	= C1.arc[0].y
			C2		= part_contour(Cin2, p2)
			C.x		= C2.arc[1].x
			C.y		= C2.arc[1].y

	return C

