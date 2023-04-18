import numpy as np
import types

#	############
#	builds a 0d signal
#	############

def build_signal_0d(Values, Descr="", Units=""):

	signal0d = types.SimpleNamespace()
	signal0d.type	= 0
	signal0d.values = np.copy(Values)
	signal0d.descr  = Descr
	signal0d.units  = Units

	return signal0d

#	############
#	builds a 1d signal
#	############

def build_signal_1d(Values, Times, Descr="", Units=""):

	signal1d = types.SimpleNamespace()
	signal1d.type	= 1
	signal1d.values = np.copy(Values)
	signal1d.times  = np.copy(Times)
	signal1d.descr  = Descr
	signal1d.units  = Units

	return signal1d


#	############
#	builds a 2d signal
#	############

def build_signal_2d(Values, Times, xValues, Descr="", Units=""):

	signal2d = types.SimpleNamespace()
	signal2d.type	 = 2
	signal2d.values  = np.copy(Values)
	signal2d.times   = np.copy(Times)
	signal2d.xvalues = np.copy(xValues)
	signal2d.descr   = Descr
	signal2d.units   = Units

	return signal2d
