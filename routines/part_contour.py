import types
import numpy as np

def part_contour(Cin, p1):

#	 input ordered contour
#	 output pieces of contour delimited by p1

	if((Cin.x[0] == Cin.x[-1]) and (Cin.y[0] == Cin.y[-1])):
		Cin.x = Cin.x[:-1]
		Cin.y = Cin.y[:-1]

	d		= np.sqrt((Cin.x-p1.x)**2+(Cin.y-p1.y)**2)
	ind		= np.where(d == d.min()); ind = ind[0][0]
	vx		= p1.x-Cin.x[ind]
	vy		= p1.y-Cin.y[ind]

	if(ind > 0):
		vmx = p1.x-Cin.x[ind - 1]
		vmy = p1.y-Cin.y[ind - 1]
		
	if(ind < len(d)-1):
		vpx = p1.x-Cin.x[ind + 1]
		vpy = p1.y-Cin.y[ind + 1]

	C		= types.SimpleNamespace()
	C.arc	= [types.SimpleNamespace(), types.SimpleNamespace()]
	if(ind == 0):
		if(vx*vpx+vy*vpy > 0):
			C.arc[0].x = p1.x
			C.arc[0].y = p1.y
			C.arc[1].x = Cin.x
			C.arc[1].y = Cin.y
		else:
			C.arc[0].x = np.append(Cin.x[0], p1.x)
			C.arc[0].y = np.append(Cin.y[0], p1.y)
			C.arc[1].x = np.append(p1.x,	 Cin.x[1:])
			C.arc[1].y = np.append(p1.y,	 Cin.y[1:])

	else:
		if(ind == len(d)-1):
			if(vx*vmx+vy*vmy > 0):
				C.arc[1].x = p1.x
				C.arc[1].y = p1.y
				C.arc[0].x = Cin.x
				C.arc[0].y = Cin.y
			else:
				C.arc[0].x = np.append(Cin.x[:-1], p1.x)
				C.arc[0].y = np.append(Cin.y[:-1], p1.y)
				C.arc[1].x = np.append(p1.x,	   Cin.x[-1])
				C.arc[1].y = np.append(p1.y,	   Cin.y[-1])

		else:																# general case
			if(vx*vpx+vy*vpy > 0):											# p1 between ind and ind-1 probably
				C.arc[0].x = np.append(Cin.x[:ind], p1.x)
				C.arc[0].y = np.append(Cin.y[:ind], p1.y)
				C.arc[1].x = np.append(p1.x,		Cin.x[ind:])
				C.arc[1].y = np.append(p1.y,		Cin.y[ind:])
			else:
				C.arc[0].x = np.append(Cin.x[:ind+1], p1.x)
				C.arc[0].y = np.append(Cin.y[:ind+1], p1.y)
				C.arc[1].x = np.append(p1.x,		  Cin.x[ind+1:])
				C.arc[1].y = np.append(p1.y,		  Cin.y[ind+1:])
	C.ind = ind

	return C
