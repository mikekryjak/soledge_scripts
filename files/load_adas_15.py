import types
import os
import numpy 				as np
from routines.globals		import DEBUG, CSPEED, HPLANCK

def load_adas_15(file_name):
	
#	data

	cslash = "/"
	colon  = ":"
	cblnk8 = "        "
	ckey1  = "filmem"
	ckey2  = "type"
	ckey3  = "indm"
	ckey4  = "isel"
	ckey5  =  "pl"
	ckey6  = "ss"
	ckey7  = "pb"
	ckey8  = "pp"
	ckey9  = "sz"
	ckey10 = "ispb"
	ckey11 = "ispp"

	adas = types.SimpleNamespace()
	adas.has_data = False

#	start opening file

	if(DEBUG > 0): print("load_adas_15: reading from: ",file_name)

	AdasPath = os.path.dirname(os.path.abspath(__file__))[:-5]+"data/adas/adf15/"

	if(not os.path.isfile(AdasPath+file_name)):	return adas

	fid   = open(AdasPath+file_name,'r')

	head  = fid.readline()[:-1]

	heads = head.split()
	nbsel = eval(heads[0])
	if  (head[:14].find("+") > -1): lptn_temp= False
	elif(head[:14].find(":") > -1): lptn_temp= True
	else:
		print("\tIncorrect element/ion field: ",head[10:14])
		return adas


	adas.has_data = True

	ipos1	 = head.find("/")+1
	adas.ion_name = head[ipos1:ipos1+5]

	if((adas.ion_name[1:1] == "+") or (adas.ion_name[1:1] == ":")): e_sym = adas.ion_name[0]
	else:														 	  e_sym = adas.ion_name[:2]

	adas.iz0  = i4eiz0(e_sym)
	try:
		is0 = eval(adas.ion_name[3:5])
	except:
		is0 = eval(adas.ion_name[2:4])
	is1   = is0 + 1

#-----------------------------------------------------------------------
# input connection vector and partition block if appropriate
#-----------------------------------------------------------------------

	lres = False
	ncptn_stack = 0
	if(lptn_temp):
		lptn = True
		icnctv = np.empty(0, dtype = 'i4')
		cstrg  = fid.readline()[:-1]
		if(cstrg[1:10] == "---------"):
			cstrg  = fid.readline()[:-1]
			fid.seek(fid.tell() - len(cstrg) - 1)
			if(cstrg[:2] != "//"):
				lres = True
				while (True):
					cstrg  = fid.readline()[:-1]
					if(cstrg[1:10] == "---------"): break
					ifirsta , ilasta = xxword(cstrg , " ")
					for j in range(len(ifirsta)): icnctv = np.append(icnctv, eval(cstrg[ifirsta[j]:ilasta[j]]))
		else:
#			print("\tNo connection/partition block  present")
			fid.seek(fid.tell() - len(cstrg) - 1)
			lptn = False

		if(lptn):
			lresol = False
			nptn, nptnc, iptnla, iptna, iptnca, cptn_stack, lptn, cstrg = xxrptn(fid, lresol)
			iptnl = np.max(iptnla)

	else:
		cstrg  = fid.readline()[:-1].lower()
		fid.seek(fid.tell() - len(cstrg) - 1)
		try:
			i = eval(xxhkey(cstrg, ckey3, cslash))
			lres = True
		except:
			pass


#-----------------------------------------------------------------------
# input data for each of the data-blocks.  Evaluate largest indm to see
# if confirms unresolved.
#-----------------------------------------------------------------------
	max_indm = 0

	blocks = []
	for iblk in range(1,nbsel+1):

#-----------------------------------------------------------------------
# input title, wavelength and other information (check block exists)
#-----------------------------------------------------------------------

		blocks.append(types.SimpleNamespace())

		cstrg = fid.readline()[:-1].lower()
		if (cstrg[0] != "c"):
			ipos1 =  cstrg.find('a') + 1
			ipos2 =  cstrg.find('/') - 1

			if(ipos1 == 0):
				ipos1 = 10
				blocks[-1].cwavel = cstrg[0:10]
			elif (ipos1 > 9):
				blocks[-1].cwavel = cstrg[0:10]
			else:
				blocks[-1].cwavel = cstrg[0:ipos1]

			blocks[-1].cwavel = blocks[-1].cwavel.strip()
			if( blocks[-1].cwavel[-1] == 'a'): blocks[-1].wavel = eval(blocks[-1].cwavel[:-1])			#wavelenght  (nm)
			else:						  		blocks[-1].wavel = eval(blocks[-1].cwavel)

			blocks[-1].energy = CSPEED*HPLANCK/(blocks[-1].wavel*1e-9)

			heads = cstrg[ipos1:ipos2].split()
			blocks[-1].ida = eval(heads[0])
			blocks[-1].ita = eval(heads[1])

			blocks[-1].ctype  = xxhkey(cstrg, ckey2, cslash)
			lsup = False
			blocks[-1].cfile = xxhkey(cstrg, ckey1, cslash)
			if(blocks[-1].cfile == cblnk8):
				ipos1 = ipos2+2
				ipos2 = cstrg.find('/')+ipos1-1

				i = eval(xxhkey(cstrg[ipos1:ipos2+1], ckey5, colon))
				if(i != iptnl):
					print("\tERROR: partition level mismatch: isel = ",iblk)
					exit()

				i = eval(xxhkey(cstrg[ipos1:ipos2+1], ckey6, colon))
				if(i != is0):
					print("\tERROR: superstage mismatch: isel = ",iblk)
					exit()
				else:
					blocks[-1].isstgr = i+1

				try:
					blocks[-1].iszr = eval(xxhkey(cstrg[ipos1:ipos2+1], ckey9, colon))
				except:
					print("\tERROR: ionis. stage error: isel = ",iblk)
					exit()

				blocks[-1].iszr = i

				if(blocks[-1].ctype == 'excit'):
					try:
						blocks[-1].ispbr = eval(xxhkey(cstrg[ipos1:ipos2+1], ckey7, colon))
						blocks[-1].isppr = 0
					except:
						print("\tERROR: pb error: isel = ",iblk)
						exit()
				else:
					try:
						blocks[-1].ispbr = 0
						blocks[-1].isppr = eval(xxhkey(cstrg[ipos1:ipos2+1], ckey8, colon))
					except:
						print("\tERROR: pp error:' isel = ",iblk)
					exit()

				lsup = True

			if(not lptn_temp):
				blocks[-1].cindm = xxhkey(cstrg, ckey3 , cslash)
				try:
					max_indm = max(max_indm, eval(cindm[iblk]))
				except:
					max_indm = max(max_indm,1)

				if(blocks[-1].ctype == 'excit'):
					blocks[-1].ispbr = max_indm
					blocks[-1].isppr = 0
				elif(blocks[-1].ctype == 'recom'):
					blocks[-1].ispbr = 0
					blocks[-1].isppr = max_indm
				else:
					blocks[-1].ispbr = 0
					blocks[-1].isppr = 0

				blocks[-1].isstgr = is0 + 1
				blocks[-1].iszr   = is0

			elif(lptn_temp and (blocks[-1].ctype == 'excit')): blocks[-1].cindm = xxhkey(cstrg, ckey10 , cslash)
			elif(lptn_temp and (blocks[-1].ctype == 'recom')): blocks[-1].cindm = xxhkey(cstrg, ckey11 , cslash)
			elif(lptn_temp and (blocks[-1].ctype == 'chexc')): blocks[-1].cindm = xxhkey(cstrg, ckey11 , cslash)

			blocks[-1].isela = eval(xxhkey(cstrg, ckey4 , cslash))

			ndnum = blocks[-1].ida
			ntnum = blocks[-1].ita

#-----------------------------------------------------------------------
# input temperature, density and photon emissivities for block
#-----------------------------------------------------------------------
			blocks[-1].teda = adas_read_array(fid, ndnum, 'f8', 8)*1e6				#Convert density from 1/cm^3 to 1/m^3
			blocks[-1].teta = adas_read_array(fid, ntnum, 'f8', 8)

			blocks[-1].pec = np.empty((ndnum, ntnum), dtype = "f8")
			for itd in range(ndnum):
				blocks[-1].pec[itd, :] = adas_read_array(fid, ntnum, 'f8', 8)*1e-6	#Convert  photon emissivity coefficients from (cm^3/s) to (m^3/s)

			blocks[-1].pec_max = np.max(blocks[-1].pec[itd, 0])
		else:
			print("\tIinput ions./photon data set name: ",file_name)
			print("\tinconsistency in the number of wavelengths expected and read.")
			print("\t\tthe number expected = ",nbsel)
			print("\t\tthe number read in  = ",iblk-1)
			print("\t\ttherefore only ",blk-1," have been accepted")
			nbsel = iblk - 1
			break

	fid.close()

	adas.wavels   = np.array([blocks[i].wavel for i in range(len(blocks))])
	adas.energies = np.array([blocks[i].energy for i in range(len(blocks))])
	adas.blocks   = blocks

	return adas

#=======================================================00

def adas_read_array(file_id, n_points, data_type, n_on_lines):

	Array = np.empty(n_points, dtype=data_type)

	n_full_lines = int(n_points/n_on_lines)
	n_on_last_line = n_points - n_full_lines*n_on_lines

	k = 0
	for l in range(n_full_lines):
		tmp		= file_id.readline()
		if((l == 0) and (len(tmp) < 2)): tmp = file_id.readline()

		values  = eval(",".join(tmp.split()))
		for v in range(n_on_lines):
			Array[k] = float(values[v]); k += 1

	if(n_on_last_line > 0):
		tmp		= file_id.readline()
		values  = eval(",".join(tmp.split()))
		if(n_on_last_line == 1): 
			Array[k] = float(values)
		else:
			for v in range(n_on_last_line):
				Array[k] = float(values[v]); k += 1

	return Array

def xxhkey(ctext, ckey_in, cbreak):

	ckey   = ckey_in.rstrip()
	lentxt = len(ctext)
	lenkey = len(ckey)

	if(lenkey == 0): return ""

#-----------------------------------------------------------------------
# establish if a valid key can be found in 'ctext'
#-----------------------------------------------------------------------

	ipos1 = 0
	next_key = True
	while(next_key):
		ipos1 = ctext[ipos1:].find(ckey) + ipos1
		ipos2 = ipos1 + lenkey
		if(ipos2 > lentxt):	return ""
		elif(ipos1 == 0):	ipos1 = ipos2
		elif(ipos1 != -1):
			ipos1 = ipos1 - 1
			if((ctext[ipos1] == ' ') or (ctext[ipos1] == cbreak)):		#check if key preceeded by space or break
				ipos1 = ipos2
			else:
				ipos1=ipos2												#skip sub_key
				if(ctext[ipos1:].find(ckey) != -1): continue
				else:	return ""
		else:	return ""

		if((ctext[ipos2] != ' ') and (ctext[ipos2] != '=')):		#check if key delimited by space or break
			ipos1 = ipos2
			if(ctext[ipos1:].find(ckey) != -1): continue
			else:	return ""

		next_key = False
#-----------------------------------------------------------------------
# establish if equal sign exists after the key (seperated by blanks).
#-----------------------------------------------------------------------
		ipos2 = ctext[ipos1:].find('=')
		ipos3 = ipos2 + ipos1 + 1

		if (ipos3 >= lentxt):	return ""
		elif(ipos2 ==  0):		ipos2 = ipos3
		elif(ipos2 != -1):
			ipos2 = ipos3 - 2
			if(ctext[ipos1] == ' '): ipos2 = ipos3 + 1
			else: continue

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# see if valid key exists further down the 'ctext' string.
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		else:	return ""

#-----------------------------------------------------------------------
# find seperator character and identify response (if present).
#-----------------------------------------------------------------------

	ipos1 = -1
	ipos3 = ctext[ipos2:].find(cbreak)

	if(ipos3 == -1):	ipos3 = lentxt-1
	elif(ipos3 == 0):	return ""
	else:				ipos3 = ipos3 + ipos2 - 1

	for i in range(ipos2,ipos3+1):
		if((ipos1 == -1) and (ctext[i] != ' ' )): ipos1 = i
		if((ipos1 != -1) and (ctext[i] == ' ' )):
			ipos3 = i
			break

	if(ipos1 == -1): return ""

#-----------------------------------------------------------------------
# valid response found - set up 'cans'
#-----------------------------------------------------------------------

	cans = ctext[ipos1:ipos3+1].strip()

	return cans


#=======================================================00

def i4eiz0(e_sym):

	SYMBOLS = [	'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', \
				'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca', \
				'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', \
				'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr', \
				'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', \
				'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd', \
				'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', \
				'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', \
				'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', \
				'Pa','U ']

	ii = indexes_ignore_case(SYMBOLS, e_sym)
	if(len(ii) > 0): return ii[0]+1
	else:		     return -1

#=======================================================00

def indexes_ignore_case(list_strs, str):
	return [i for i, j in enumerate(list_strs) if(j.upper() == str.upper())]

#=======================================================00
#-----------------------------------------------------------------------
#
#  ****************** FORTRAN77 SUBROUTINE: XXWORD *********************
#
#  PURPOSE: TO EXTRACT THE Nfirst to (Nfirst+IWORDS-1) WORDS FROM AN
#           INPUT STRING. OUTPUTS THE FIRST AND LAST BYTE INDEXES OF
#           EACH WORD AS WELL AS THE TOTAL NUMBER OF WORDS FOUND.
#
#           A WORD = A STRING OF CHARACTERS SEPARATED BY ANY CHARACTER
#                    CONTAINED IN THE INPUT STRING CDELIM.
#
#  CALLING PROGRAM: GENERAL USE
#
#  SUBROUTINE:
#
#  INPUT : (C*(*)) CTEXT   = INPUT TEXT LINE CONTAINING STRING
#  INPUT : (C*(*)) CDELIM  = INPUT STRING CONTAINING DELIMITER CHARS.
#  INPUT : (I*4)   NFIRST  = THE INDEX NO. OF THE FIRST WORD TO EXTRACT.
#
#  I/O   : (I*4)   IWORDS  = INPUT : SIZE OF IFIRST, ILAST(ARRAYS)
#                                    (I.E. NUMBER OF WORDS TO EXTRACT)
#                          = OUTPUT: NUMBER OF REQUESTED WORDS FOUND
#
#  OUTPUT: (I*4)   IFIRST()= INDEX OF FIRST BYTE OF THE Nth WORD
#  OUTPUT: (I*4)   ILAST() = INDEX OF LAST  BYTE OF THE Nth WORD
#  OUTPUT: (I*4)   NWORDS  = THE TOTAL NUMBER OF WORDS FOUND IN CTEXT
#
#          (I*4)   LENTXT  = LENGTH IN BYTES OF 'CTEXT' STRING
#          (I*4)   IDELIM  = 0 => CTEXT CHARACTER IS NOT A DELIMITER
#                          > 0 => CTEXT CHARACTER IS A DELIMITER
#          (I*4)   ITOTAL  = NUMBER OF WORDS FOUND SO FAR
#          (I*4)   IINDEX  = IFIRST()/ILAST() INDEX OF CURRENT WORD
#          (I*4)   NLAST   = THE INDEX NO. OF THE LAST WORD TO EXTRACT
#          (I*4)   I       = GENERAL USE INDEX
#
#          (L*4)   LWORD   = .TRUE.  - PROCESSING AN IDENTIFIED WORD
#                            .FALSE. - PROCESSING SPACE BETWEEN WORDS
#
# ROUTINES: NONE
#
# NOTES:    IF THERE IS NO Nfirst WORD OR NO WORDS ARE FOUND
#           (I.E. INPUT STRING IS BLANK) THEN IWORDS=0
#
def xxword(ctext, cdelim, nfirst=0):

#-----------------------------------------------------------------------
# find the required words
#-----------------------------------------------------------------------

	ifirst = np.empty(0, dtype = 'i4')
	ilast = np.empty(0, dtype = 'i4')

	lentxt = len(ctext)
	itotal = -1
	lword  = False
	for i in range(lentxt):
		idelim = cdelim.find(ctext[i])
		if(lword):
			if (idelim > -1):
				lword  = False
				if(itotal >= nfirst): ilast[-1] = i				#in python last must mu be +1 to build a range
		else:
			if(idelim == -1):
				lword  = True
				itotal = itotal + 1
				if(itotal >= nfirst):
					ifirst = np.append(ifirst, i)
					ilast = np.append(ilast, lentxt)

	return ifirst, ilast

#-----------------------------------------------------------------------
#
#
#  Purpose:  To read and analyse a partition block in a datafile header
#
#  Calling program: adas416
#
#  Notes:  (1) Partition levels, partitions and partition components are 
#              labelled starting at 0 (but see (2)).
#          (2) Partition level 0 labels the resolved root partition level
#              partition level 1 labels the unresolved root partition
#              level.
#          (3) For an unresolved (standard) file, the partitions are each 
#              ionisation stage from the neutral to the bare nucleus and
#              they are labelled by the ion charge.  Each partition has 
#              just the one component.
#          (4) Distinguish the indexing (starting at 1) from the label
#              (starting at 0) .
#
#  Subroutine:
#
#  input : (i*4)  iunit     = unit to which input file is allocated
#  input : (i*4)  ndstack   = maximum no. of text lines in partition block 
#
#  input : (i*4)  ndptnl    = maximum level of partitions
#  input : (i*4)  ndptn     = maximum no. of partitions in one level
#  input : (i*4)  ndptnc    = maximum no. of components in a partition
#  input : (l*4)  lresol    = .true.  => resolved root partition
#                           = .false. => standard root partition
#
#  output: (i*4)  nptn()    = number of partitions in partition level
#                             1st dim: partition level 
#  output: (i*4)  nptnc(,)  = number of components in partition
#                             1st dim: partition level 
#                             2nd dim: member partition in partition level 
#  output: (i*4)  iptnla()  = partition level label (0=resolved root,1=
#                                                      unresolved root)
#                             1st dim: partition level index 
#  output: (i*4)  iptna(,)  = partition member label (labelling starts at 0)
#                             1st dim: partition level index 
#                             2nd dim: member partition index in partition 
#                             level 
#  output: (i*4)  iptnca(,,)= component label (labelling starts at 0)
#                             1st dim: partition level index 
#                             2nd dim: member partition index in partition level
#                             3rd dim: component index of member partition 
#  output: (l*4)  lptn      = .true.  => partition block present
#                           = .false. => partition block not present
#  output: (c*80) cstrg     = string marking end of partition block
#  output: (c*80) cptn_stack()=text lines of partition block       
#                              1st dim: text line pointer 
#

def xxrptn(fid, lresol):
 
#-----------------------------------------------------------------------
#  read partition data if present
#-----------------------------------------------------------------------

	lptn		= False
	cptn_stack	= []
	cptna		= []
	while (True):
		cstrg = fid.readline()[:-1]
		if(cstrg[:3] == "//#"):
			lptn	= True
			lptn_in = True
			ifirst,ilast = xxslen(cstrg)
			cptn_stack.append(cstrg[ifirst:ilast])
			cptna.append(cstrg[ifirst:ilast])
		elif((cstrg[:3] == "   ") and lptn_in):
			ifirst,ilast = xxslen(cstrg)

			cptn_stack.append("      " + cstrg[ifirst:ilast])
			cptna[-1] = cptna[-1]+ " "  + cstrg[ifirst:ilast]
		elif(cstrg[:3] == "---"):
			lptn_in = False
			break

#-----------------------------------------------------------------------
#  preliminary checks
#-----------------------------------------------------------------------

	iptnla = np.zeros(len(cptna), dtype='i4')
	if(lptn):
		for i in range(len(cptna)):
			if(cptna[i][3] == "0"): iptnla[i] = eval(cptna[i][4]) - 1
			else:				  	 iptnla[i] = eval(cptna[i][3:5]) - 1								#as index reduce to python index

		if((iptnla[i] > 1) and (len(cptna) == 1) and (not lresol)):
			print("\t*** xxrptn error *** partition level",i," in error")
			print("\t\tiptnla[i],len(cptna),lresol=",iptnla[i],len(cptna),lresol)
			exit()

		for i in range(len(cptna)):
			ifirst, ilast = xxslen(cptna[i])
			if(cptna[i][ilast-1] != "/"):
				print("\t*** xxrptn error *** partition level",i," not terminated")
				print("\t\tcptna[i][ilast-1]=",cptna[i][ilast-1])
				print("\t\tcptna[i][ilast-2]=",cptna[i][ilast-2])
				print("\t\tcptna[i]=",cptna[i])
				exit()

#-----------------------------------------------------------------------
#  analyse partitions
#-----------------------------------------------------------------------

	nptn  = np.empty(len(cptna), dtype = 'i4')
	ncpta = np.empty(len(cptna), dtype = 'i4')
	nptncM  = []
	iptnaM  = []
	iptncaM = []
	for i in range(len(cptna)):
		ifirst,ilast = xxslen(cptna[i])
		ifirsta, ilasta = xxword(cptna[i][ifirst+5:ilast] , '/')
		if((len(ifirsta)-2*int(len(ifirsta)/2)) != 0):
			print("\t*** xxrptn error *** partition level",i," partition count wrong")
			print("\t\tlen(ifirsta)=",len(ifirsta))
			exit()

		nptn[i]=int(len(ifirsta)/2)
		ic = 0
		iptnaM.append([])
		nptncM.append([])
		iptncaM.append([])
#		print("xxrptn: i,nwords,nptn(i)=",i,len(ifirsta),nptn[i])
		for j in range(nptn[i]):
			strg = cptna[i][ifirsta[2*j]+5:ilasta[2*j]+5]
#			print("\txxrptn-xxslen:strg="+strg)
			ifirst,ilast = xxslen(strg)
			if(strg[ifirst+1] != "0"): iptnaM[-1].append(eval(strg[ifirst+1:ilast]  ))
			else:					    iptnaM[-1].append(eval(strg[ifirst+2:ilast]))
			strg		= cptna[i][ifirsta[2*j+1]+5:ilasta[2*j+1]+5]

#			print("xxrptn-xxword:strg="+strg)
			ifirstac , ilastac = xxword(strg , ' ')
#			print("\txxrptn:nwordsc=",len(ifirstac))
			nptncM[-1].append(len(ifirstac))
#			print("\txxrptn: i,j,nptnc(i,j)=",i,j,nptncM[-1][-1])
#			print("\txxrptn: ifirstac(1),ilastac(1)=",ifirstac[0],ilastac[0])
#			print("\txxrptn:",strg[ifirstac[0]:ilastac[0]])

			iptncaM[-1].append([])
			for k in range(nptncM[-1][-1]):
				if(strg[ifirstac[k]] != "0"):	iptncaM[-1][-1].append(eval(strg[ifirstac[k]  :ilastac[k]]))
				else:							iptncaM[-1][-1].append(eval(strg[ifirstac[k]+1:ilastac[k]]))
#------------------------------------------------
#  components should be sequential and contiguous. 
#-------------------------------------------------
				if(iptncaM[-1][-1][-1] != ic):
					print("\t*** xxrptn error *** partition count",ic," out of sequence")
					print("\t\tiptnca[i,j,k],ic=",iptncaM[-1][-1][-1],ic)
					exit()
				else:
					ic += 1

			ncpta[i] = ic

#		print("\ti,ncpta(i)=",i,ncpta[i])

#----------------------------------------------------------------------------------------------------
#  final checks. level i+1 partition components should span partition i p-values
#----------------------------------------------------------------------------------------------------

	if(len(nptn) >= 2):
		for i in range(1,len(nptn)):
			if(ncpta[i-1] != nptn[i]):
				print("\t*** xxrptn error *** partition level",i," count incorrect")
				print("\t\tncpta[i-1],nptn[i]=",ncpta[i-1],nptn[i])
				exit()

#	Transfer python matrixes to numpy matrixes

	dim1  = len(cptna)
	dim2  = np.max(nptn)
	iptna = np.zeros((dim1,dim2), dtype = 'i4')
	nptnc = np.zeros((dim1,dim2), dtype = 'i4')
	for i in range(dim1):
		for j in range(dim2):
			iptna[i,j] = iptnaM[i][j] 
			nptnc[i,j] = nptncM[i][j]

	dim3 = np.max(nptnc)
	iptnca = np.zeros((dim1,dim2, dim3), dtype = 'i4')
	for i in range(dim1):
		for j in range(dim2):
			for k in range(dim3): iptnca[i,j,k] = iptncaM[i][j][k]
	
	return nptn, nptnc, iptnla, iptna, iptnca, cptn_stack, lptn, cstrg


#-----------------------------------------------------------------------

#  ****************** FORTRAN77 SUBROUTINE: XXSLEN *********************

#  PURPOSE: TO IDENTIFY THE FIRST AND LAST NON-BLANK CHARACTER IN A
#           STRING. (IF INPUT STRING IS BLANK IFIRST=ILAST=0)

#  CALLING PROGRAM: GENERAL USE

#  SUBROUTINE:

#  INPUT : (C*(*)) CSTRNG   = INPUT STRING FOR INTERROGATION

#  OUTPUT: (I*4)   IFIRST   = BYTE POSITION OF FIRST NON-BLANK CHARACTER
#                             IN INPUT STRING.
#  OUTPUT: (I*4)   ILAST    = BYTE POSITION OF LAST  NON-BLANK CHARACTER
#                             IN INPUT STRING.

#          (I*4)   I        = GENERAL USE
#          (I*4)   ILEN     = LENGTH OF 'CSTRNG' STRING IN BYTES

# ROUTINES: NONE

# NOTE:


# AUTHOR:  PAUL E. BRIDEN (TESSELLA SUPPORT SERVICES PLC)

def xxslen(cstrng):
	ifirst = -1
	ilast  = len(cstrng)
	for i in range(len(cstrng)):
		if(cstrng[i] != " "):
			if (ifirst == -1):ifirst = i
			ilast = i

	return ifirst, ilast+1							#move forward ilast for python range