from math import floor
import numpy as np

def mesh_segment_linear(refineSide, adjustMode, pt_num, paramL, paramR, longueur):

	if(refineSide == 0) :								#left
		if(adjustMode == 1):							#absolute
			M		= np.zeros((2,2), dtype='f8')
			RHS		= np.zeros((2), dtype='f8')
			M[0,0]	= (pt_num-1)*pt_num/2
			M[0,1]	= pt_num-1
			RHS[0]	= 1
			M[1,0]	= 1
			M[1,1]	= 1
			RHS[1]	= paramL/100/longueur
			A		= np.dot(np.linalg.inv(M),RHS)
			dpt		= A[0]*np.arange(1,pt_num)+A[1]
			dpt		= np.where(dpt > 1.e-4, dpt, 1.e-4)
		else:
			dpt		= 1+paramL*np.arange(1,pt_num)
							
		pt	= np.zeros((pt_num), dtype='f8')
		for k in range(1,pt_num): pt[k] = pt[k-1]+dpt[k-1]

		pt	= (pt-pt[0])/(pt[-1]-pt[0])

	elif(refineSide == 1):				#right
		if(adjustMode == 1):			#absolute
			M		= np.zeros((2,2), dtype='f8')
			RHS		= np.zeros((2), dtype='f8')
			M[0,0]	= (pt_num-1)*pt_num/2
			M[0,1]	= pt_num-1
			RHS[0]	= 1
			M[1,0]	= pt_num-1
			M[1,1]	= 1
			RHS[1]	= paramR/100/longueur
			A 		= np.dot(np.linalg.inv(M),RHS)
			dpt		= A[0]*np.arange(1,pt_num)+A[1]
			dpt		= np.where(dpt > 1.e-4, dpt, 1.e-4)
		else:
			dpt		= 1+paramR*np.arange(1,pt_num)
								
		pt		= np.zeros(pt_num, dtype='f8')
		for k in range(1,pt_num): pt[k]=pt[k-1]+dpt[k-1]

		pt		= (pt-pt[0])/(pt[-1]-pt[0])

		if(adjustMode == 0): pt	= 1 - pt[::-1]

	else:												# both
		if(adjustMode == 1):							# absolute
																			# part1
			M		= np.zeros((2,2), dtype='f8')
			RHS		= np.zeros((2), dtype='f8')
			M[0,0]	= (pt_num/2-1)*pt_num/2/2
			M[0,1]	= pt_num/2-1
			RHS[0]	= 0.5
			M[1,0]	= 1
			M[1,1]	= 1
			RHS[1]	= paramL/100/longueur
			A 		= np.dot(np.linalg.inv(M),RHS)
			dpt		= A[0]*np.arange(1,floor(pt_num/2))+A[1]
			dpt		= np.where(dpt > 1.e-4, dpt, 1.e-4)
			pt1		= np.zeros(1,floor(pt_num/2), dtype='f8')
			for k in range(1,floor(pt_num/2)): pt1[k]=pt1[k-1]+dpt[k-1]
																			# part2
			M		= np.zeros((2,2), dtype='f8')
			RHS		= np.zeros((2), dtype='f8')
			M[0,0]	= (pt_num/2-1)*pt_num/2/2
			M[0,1]	= pt_num/2-1
			RHS[0]	= 0.5
			M[1,0]	= pt_num/2-1
			M[1,1]	= 1
			RHS[1]	= paramR/100/longueur
			A 		= np.dot(np.linalg.inv(M),RHS)
			dpt		= A[0]*np.arange(1,pt_num-floor(pt_num/2)+1)+A[1]
			dpt		= np.where(dpt > 1.e-4, dpt, 1.e-4)
			pt2		= np.zeros(1,pt_num-floor(pt_num/2)+1, dtype='f8')
			for k in range(1,pt_num-floor(pt_num/2)+1): pt2[k]=pt2[k-1]+dpt[k-1]
								
			pt		= np.append(pt1,0.5+pt2[1:])
			pt		= (pt-pt[0])/(pt[-1]-pt[0])
		else:
			dpt		= 1+paramL*np.arange(1,floor(pt_num/2))
			pt1		= np.zeros(floor(pt_num/2), dtype='f8')
			for k in range(1,floor(pt_num/2)): pt1[k]=pt1[k-1]+dpt[k-1]

			pt1		= (pt1-pt1[0])/(pt1[-1]-pt1[0])*0.5
			dpt		= 1+paramR*np.arange(1,pt_num-floor(pt_num/2)+1)
			pt2		= np.zeros((pt_num-floor(pt_num/2)+1), dtype='f8')
			for k in range(1,pt_num-floor(pt_num/2)+1): pt2[k]=pt2[k-1]+dpt[k-1]

			pt2		= (pt2-pt2[0])/(pt2[-1]-pt2[0])
			pt2		= (1-pt2[::-1])*0.5
			pt		= np.append(pt1,0.5+pt2[1:])

	return pt