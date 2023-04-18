import types
import numpy as np

# Extract sub-contour of contour between 2 points
# The extracted contour starts from the point of the initial contour the closest to p1
# One then moves in direction dir until either the end of the contour or p2 are met

def part_contour_per(cin, p1, p2, dir):

	if((cin.x[0] == cin.x[-1]) and (cin.y[0] == cin.y[-1])):
		cin.x = cin.x[:-1]
		cin.y = cin.y[:-1]
		
	d = np.sqrt((cin.x-p1.x)**2+(cin.y-p1.y)**2)
	a = np.where(d == d.min()); a = a[0][0]
	
	cleft = types.SimpleNamespace()
	cout  = types.SimpleNamespace()
	if(dir == 1):
		if(a < len(d)-1):
			v1x	= cin.x[a + 1]-cin.x[a]
			v1y	= cin.y[a + 1]-cin.y[a]
			vx	= p1.x-cin.x[a]
			vy 	= p1.y-cin.y[a]
			if(v1x*vx+v1y*vy > 0):
				cleft.x = np.roll(cin.x[:-1], -(a+1))
				cleft.y = np.roll(cin.y[:-1], -(a+1))
			else:
				cleft.x =  np.roll(cin.x[:-1], -a)
				cleft.y =  np.roll(cin.y[:-1], -a)

			cleft.x = np.append(cleft.x, cleft.x[0])
			cleft.y = np.append(cleft.y, cleft.y[0])
		else:
			cleft.x = cin.x
			cleft.y = cin.y

#		check if p2 in [p1, cleft[0]]

		d	= np.sqrt((cleft.x-p2.x)**2+(cleft.y-p2.y)**2)
		a2	= np.where(d == d.min()); a2 = a2[0][0]
		if((a2 == len(d)-1) or (a2 == len(d)-2) or (a2 == 0)):
			vx	= p2.x-p1.x
			vy	= p2.y-p2.x
			v1x	= cleft.x[-1]-cleft.x[-2]
			v1y	= cleft.y[-1]-cleft.y[-2]
			if(vx*v1x+vy*v1y > 0):
				cout.x	= np.array([])
				cout.y	= np.array([])
			else:
				cout.x = cleft.x[:-1]
				cout.y = cleft.y[:-1]

		else:
			v1x	= cleft.x[a2 + 1] - cleft.x[a2]
			v1y	= cleft.y[a2 + 1] - cleft.y[a2]
			vx	= p2.x - cleft.x[a2]
			vy	= p2.y - cleft.y[a2]
			if(v1x*vx+v1y*vy < 0):
				cout.x = cleft.x[:a2]
				cout.y = cleft.y[:a2]
			else:
				cout.x = cleft.x[:a2+1]
				cout.y = cleft.y[:a2+1]

		cout.x = np.append(np.append(p1.x, cout.x), p2.x)
		cout.y = np.append(np.append(p1.y, cout.y), p2.y)
	else:

		if(a < len(d)-1):
			v1x = cin.x[a + 1]-cin.x[a]
			v1y = cin.y[a + 1]-cin.y[a]
			vx	= p1.x - cin.x[a]
			vy	= p1.y - cin.y[a]
			if(v1x*vx+v1y*vy > 0):
				cleft.x = np.roll(cin.x[:-1], -(a+1))
				cleft.y = np.roll(cin.y[:-1], -(a+1))
			else:
				cleft.x = np.roll(cin.x[:-1], -a)
				cleft.y = np.roll(cin.y[:-1], -a)

			cleft.x = np.append(cleft.x, cleft.x[0])
			cleft.y = np.append(cleft.y, cleft.y[0])
		else:
			cleft.x = cin.x
			cleft.y = cin.y

		cright		= types.SimpleNamespace()
		cright.x	= np.append(cleft.x[-2::-1], cleft.x[-2])
		cright.y	= np.append(cleft.y[-2::-1], cleft.y[-2])

		d			= np.sqrt((cright.x-p2.x)**2+(cright.y-p2.y)**2)
		a2			= np.where(d == d.min()); a2 = a2[0][0]
		if((a2 == len(d)-1) or (a2 == len(d)-2) or (a2 == 0)):
			vx	= p2.x - p1.x
			vy	= p2.y - p1.y
			v1x = cright.x[-1] - cright.x[-2]
			v1y = cright.y[-1] - cright.y[-2]
			if(vx*v1x+vy*v1y > 0):
				cout.x = np.array([])
				cout.y = np.array([])
			else:
				cout.x = cright.x[:-1]
				cout.y = cright.y[:-1]

		else:
			v1x	= cright.x[a2 + 1] - cright.x[a2]
			v1y	= cright.y[a2 + 1]  -cright.y[a2]
			vx	= p2.x - cright.x[a2]
			vy	= p2.y - cright.y[a2]
			if(v1x*vx+v1y*vy < 0):
				cout.x = cright.x[:a2]
				cout.y = cright.y[:a2]
			else:
				cout.x = cright.x[:a2+1]
				cout.y = cright.y[:a2+1]


		cout.x = np.append(np.append(p1.x, cout.x), p2.x)
		cout.y = np.append(np.append(p1.y, cout.y), p2.y)

	return cout