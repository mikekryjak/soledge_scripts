import types
import numpy	as np

def copyP(P):
	
	NewP		= types.SimpleNamespace()
	NewP.Br		= np.copy(P.Br)
	NewP.Bz		= np.copy(P.Bz)
	NewP.Bphi	= np.copy(P.Bphi)
	NewP.R		= np.copy(P.R)
	NewP.Z		= np.copy(P.Z)
	
	return NewP

def zerosP_like(P):
	
	NewP		= types.SimpleNamespace()
	NewP.Br		= np.zeros_like(P.Br)
	NewP.Bz		= np.zeros_like(P.Bz)
	NewP.Bphi	= np.zeros_like(P.Bphi)
	NewP.R		= np.zeros_like(P.R)
	NewP.Z		= np.zeros_like(P.Z)
	
	return NewP

	
