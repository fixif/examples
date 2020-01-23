# some tests

from fixif.LTI import Butter
from fixif.Structures import DFI,DFII,rhoDFII
from numpy import ceil, floor, log2, array, r_, vectorize
from fixif.FxP import Constant, FPF
from studySoP import SoP

def mulFPF(F1, F2):
	"""FPF of a multiplication"""
	return FPF(wl=F1.wl+F2.wl, msb=F1.msb+F2.msb)


# realization
F = Butter(9, 0.41)
R = DFII(F)

# ubar, wordlength
u_bar = 12
w_txyu = 16
w_z = 16




# SoP
SoP(R, u_bar, w_txyu, w_z)

#print(R.algorithmTxt('%.4e'))

#lsb, error_budget_y = R._compute_LSB(-12)
#print(lsb)
