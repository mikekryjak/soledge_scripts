import types
import numpy					as np
from interfaces.plot_and_ask	import plot_and_ask
from routines.utils_walls		import set_XYlines_walls
from routines.globals			import DEBUG

def optim_mesh(Root, din, Rin, Zin, cx, cy, sens, Optimization):

	if(DEBUG > 2): print("optim_mesh")

#	 return optimized mesh on piece of contour c
#	 generated from near points R1, Z1

	Npts = len(Rin)
	dout = din

	Dcxy		= np.sqrt((cx[1:]-cx[:-1])**2+(cy[1:]-cy[:-1])**2)
	dist	 	= np.append(0., np.cumsum(Dcxy))
	if(dist[-1] == 0.):
		print("ERROR!!! in optim_mesh, zero distance for:")
		print("\tRin =",Rin)
		print("\tZin =",Zin)
		print("\tcx  =",cx)
		print("\tcy  =",cy)
		print("\tdist=",dist)

		Xarrs, Yarrs, lines = set_XYlines_walls(Root.Config)
		Xarrs.append(Rin)
		Yarrs.append(Zin)
		lines.append("r.-")

		LinesData = [Xarrs, Yarrs, lines]
		choice = plot_and_ask(Root, LinesData=LinesData, title="Mesh optimization problem: segment=red")
		return din, Rin, Zin

#	First eliminate double points in (cx,cy) to avoid crash of interp

	thresh = 1.e-12
	iSDouble = np.where(Dcxy < thresh, 1, 0)[0]
	nDoubles = np.sum(iSDouble)
	if(nDoubles > 0):
		long 		= len(cx)-nDoubles
		cx_clean	= np.empty((long), dtype='f8')
		cy_clean	= np.empty((long), dtype='f8')
		cx_clean[0] = cx[0]
		cy_clean[0] = cy[0]
		k2 = 0
		for k in range(1,len(cx)):
			if(isdouble[k-1] == 0):
				k2 = k2+1
				cx_clean[k2] = cx[k]
				cy_clean[k2] = cy[k]
	else:
		cx_clean = cx
		cy_clean = cy

	delta_min	= Optimization.dmin_ad*1e-3;
	delta_min	= delta_min/dist[-1]

	dist		= dist/dist[-1]
	
#	first iteration
	Rout = np.interp(dout, dist, cx_clean)
	Zout = np.interp(dout, dist, cy_clean)

#	 figure[0]
#	 clf
#	 hold on
#	 plot([Rin;Rout],[Zin;Zout],'b-')
#	 plot([Rout],[Zout],'ro')
#	 plot([Rin;Rout]',[Zin;Zout]','b-')

	long	= len(Rin)
	vx		= Rin[2:]	 - Rin[1:-1]
	vy		= Zin[2:]	 - Zin[1:-1]
	vx1		= Rout[1:-1] - Rin[1:-1]
	vy1		= Zout[1:-1] - Zin[1:-1]
	qual	= (vx1*vy-vy1*vx)/(np.sqrt(vx**2+vy**2)*np.sqrt(vx1**2+vy1**2))		# +...
#			  (vx2*vy-vy2*vx)/(np.sqrt(vx**2+vy**2)*np.sqrt(vx2**2+vy2**2))

	l1		= dout[2:]	 - dout[1:-1]
	l2		= dout[1:-1] - dout[:-2]
	cost	= qual*sens*Optimization.orthoptim+((l1/l2+l2/l1)-1.)*(1.-Optimization.orthoptim)

	vx1 = 0; vy1 = 0; l1 = 0; l2 = 0
#	figure(2)
#	plot(cost,'b.-')
#	hold on
#	plot(mean(cost)*ones(size(cost)),'g-')
#	plot(min(cost)*ones(size(cost)),'r-')

	cost_er	= np.array([cost.min()])
	n		= 0
	while((cost_er[-1] < Optimization.target_ad) and (n < Optimization.length_ad)):

#		compute grad
		grad 		= np.zeros(long-2, dtype = 'f8')
		dd	 		= np.zeros(len(dout), dtype = 'f8')
		dout_tt		= np.empty((long-2, len(dout)), dtype = 'f8')
		for k in range(long-2):
			dd_			= 0.01*(dout[k+2]-dout[k+1])
			
			a			= min(Optimization.range_ad, (dout[k+1]-dout[0]))
			dd[:k+2]	= (dd_+dd_*(dout[:k+2]-dout[k+1])/a)
			dd[:k+2]	= np.where(dd[:k+2] > 0., dd[:k+2], 0.)
			
			a			= min([Optimization.range_ad, (dout[-1]-dout[k+1])])
			dd[k+2:]	= (dd_-dd_*(dout[k+2:]-dout[k+1])/a)
			dd[k+2:]	= np.where(dd[k+2:] > 0., dd[k+2:], 0. )
			
#	         figure(8)
#	         plot(dd,'r.-')
			dout_t				= dout + dd
			dout_t[0]			= dout[0]
			dout_t[-1]			= dout[-1]
			
#	         dout_t[:k+1]=dout[:k+1]*(dout[k+1]+0.1*(dout(k+2)-dout[k+1]))/dout[k+1]
#	         dout_t(k+2:length(dout))=dout[-1]-((dout[-1]-dout[k+2:])*((dout[-1]-(dout[k+1]+0.1*(dout(k+2)-dout[k+1])))...
#	             /(dout[-1]-dout[k+1])))

			dout_tt[k,:]		= dout_t
			Rout_t				= np.interp(dout_t, dist, cx_clean)
			Zout_t				= np.interp(dout_t, dist, cy_clean)
			cost_t				= np.zeros(long-1)
			
			vx1		= Rout_t[1:-1] - Rin[1:-1]
			vy1		= Zout_t[1:-1] - Zin[1:-1]
			qual	= (vx1*vy-vy1*vx)/(np.sqrt(vx**2+vy**2)*np.sqrt(vx1**2+vy1**2))     # +...
#					  (vx2*vy-vy2*vx)/(np.sqrt(vx**2+vy**2)*np.sqrt(vx2**2+vy2**2))

			l1		= dout_t[2:]   - dout_t[1:-1]
			l2		= dout_t[1:-1] - dout_t[:-2]
			cost_t  = qual*sens*Optimization.orthoptim+((l1/l2+l2/l1)-1.)*(1.-Optimization.orthoptim)

			vx1 = 0; vy1 = 0; i1 = 0; i2 = 0
			
			grad[k]	= (np.min(cost_t[1:])-np.min(cost[1:]))/0.01

#			plot(grad,'k.-')

		alpha		= Optimization.alpha_ad
		eps			= 1e-6
		dout_new 	= dout
		for k in range(long-2):
			dout_new = dout_new + grad[k]/np.sqrt(sum(grad**2)+eps)*(dout_tt[k,:]-dout)*alpha
#			s = s-grad[k]*alpha

#		for k in range(1,len(ddout)-1):
#	         ddout[k] = max(ddout[k],0.4*(dout[k+1]-dout[k]))
#	         ddout[k] = min(ddout[k],0.4*(dout[k-1]-dout[k]))
#		dout_new	 = dout+ddout

		delt		 = dout_new[1:] - dout_new[:-1]
		delt		 = np.where(delt < delta_min, delta_min, delt)
		dout_new[1:] = dout_new[:-1] + delt
		dout_new	 = dout_new/dout_new[-1]
#		dout_new	 = dout_new/s
		R_out_new	 = np.interp(dout_new, dist, cx_clean)
		Z_out_new	 = np.interp(dout_new, dist, cy_clean)

		vx1			= R_out_new[1:-1] - Rin[1:-1]
		vy1			= Z_out_new[1:-1] - Zin[1:-1]
		qual		= (vx1*vy-vy1*vx)/(np.sqrt(vx**2+vy**2)*np.sqrt(vx1**2+vy1**2))	#+...
#					  (vx2*vy-vy2*vx)/(np.sqrt(vx**2+vy**2)*np.sqrt(vx2**2+vy2**2))

		l1			= dout_new[2:]	 - dout_new[1:-1]
		l2			= dout_new[1:-1] - dout_new[0:-2]
		cost_new	= qual*sens*Optimization.orthoptim+((l1/l2+l2/l1)-1)*(1-Optimization.orthoptim)

		vx1 = 0; vy1 = 0; i1 = 0; i2 = 0


#		figure(4)
#		hold on
#		plot([Rin;R_out_new],[Zin;Z_out_new],'b-')
#		plot([R_out_new],[Z_out_new],'ro')
#		plot([Rin;R_out_new]',[Zin;Z_out_new]','b-')
#		figure(5)
#		plot(cost_new,'b.-')
#		hold on
#		plot(mean(cost_new)*ones(size(cost_new)),'g-')
#		plot(min(cost_new)*ones(size(cost_new)),'r-')
#

		dout	= dout_new
		Rout	= R_out_new
		Zout	= Z_out_new
		cost	= cost_new

		cost_er = np.append(cost_er, np.min(cost[1:]))
		n		= n+1
#	    figure(6)
#	    plot(cost_er,'k.-')

	vx = 0; vy = 0
	
#	figure(4)
#	clf
#	hold on
#	plot([Rin;Rout],[Zin;Zout],'b-')
#	plot([Rout],[Zout],'ro')
#	plot([Rin;Rout]',[Zin;Zout]','b-')
#	figure(5)
#	plot(cost,'b.-')
#	hold on
#	plot(mean(cost)*ones(size(cost)),'g-')
#	plot(min(cost)*ones(size(cost)),'r-')

	if(DEBUG > 2): print("optim_mesh: completed")

	return dout,Rout,Zout