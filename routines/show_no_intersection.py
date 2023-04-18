from interfaces.plot_and_ask		import plot_and_ask
from routines.utils_walls			import set_XYlines_walls
from routines.contour_better		import contour_better

def show_no_intersection(Root, Config, c1, c2):

	X_points  = Config.X_points
	Xarrs, Yarrs, lines = set_XYlines_walls(Root.Config)
	for kk in range(len(X_points)):
		cflux = contour_better(Config.r2D, Config.z2D, Config.flux2D, [X_points[kk].psi])
		for j in range(len(cflux.arc)): 
			Xarrs.append(cflux.arc[j].x)
			Yarrs.append(cflux.arc[j].y)
			lines.append("r-")

	Xarrs.append(c1.arc[0].x)
	Yarrs.append(c1.arc[0].y)
	lines.append("g.-")

	for j in range(len(c2.arc)): 
		Xarrs.append(c2.arc[j].x)
		Yarrs.append(c2.arc[j].y)
		lines.append("b.-")

	LinesData = [Xarrs, Yarrs, lines]
	choice = plot_and_ask(Root, LinesData=LinesData, title="No intersection!")

	return