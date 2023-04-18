import numpy as np

def reverse_megazone(MagZones, MagMegazones, nmz):
	
	for k in range(len(MagMegazones[nmz].list)):
#		return zone
		num						= MagMegazones[nmz].list[k]
		R						= MagZones[num].east.R
		Z						= MagZones[num].east.Z
		MagZones[num].east.R	= MagZones[num].west.R
		MagZones[num].east.Z	= MagZones[num].west.Z
		MagZones[num].west.R	= R
		MagZones[num].west.Z	= Z
		MagZones[num].north.R	= MagZones[num].north.R[::-1]
		MagZones[num].north.Z	= MagZones[num].north.Z[::-1]
		MagZones[num].south.R	= MagZones[num].south.R[::-1]
		MagZones[num].south.Z	= MagZones[num].south.Z[::-1]
		pA						= MagZones[num].pA
		pB						= MagZones[num].pB
		MagZones[num].pA		= MagZones[num].pD
		MagZones[num].pB		= MagZones[num].pC
		MagZones[num].pC		= pB
		MagZones[num].pD		= pA

	MagMegazones[nmz].list = MagMegazones[nmz].list[::-1]

