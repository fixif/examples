# coding=utf-8

from fixif.LTI import Filter
from fixif.Structures import State_Space
from fixif.FXPF import FXPF_ABCD

from loadmat import loadmat
from numpy import matrix as mat


def testFxPF(S):
	# display
	print("<"+"-"*50+">")
	print(S)

	# naive MSB spec
	#TODO:

	# use FxPF
	wl = mat([w, ] * S.n)
	try:
		FXPF_ABCD(S.dSS.A, S.dSS.B, S.dSS.C, S.dSS.D, u_bound, wl)
	except:
		pass



# constant
u_bound = mat(10)
w = 8

# ---------------
# Xilinx example.
#This example is a 5\textsuperscript{th} order lowpass Elliptic filter. Its passband is between 0 and 1e-3 (normalized frequencies) with 0.1dB as allowed ripples, whereas its stopband is between 0.0018 and 1 (normalized frequencies) with 40dB as minimum attenuation. It has been designed with Filter Designer toolbox from Matlab.

F1 = Filter(den=[1, -4.989216395071318, 9.956976990745105, -9.935631971923312, 4.957198554483321, -0.9893271782278296],
			num=[0.0001878726842913545, -0.0005635670357698394, 0.0003756943544614141, 0.0003756943544614141, -0.0005635670357698394, 0.0001878726842913545], name="Xilinx")
#F1 = Elliptic(n=5, rp=0.1, rs=40, Wn=[50.0/50000, 90.0/50000], etype='lowpass', name="Xilinx") # do not work, since I cannot reproduce the same Elliptic filter design in scipy as in Matlab/fdatool
S1 = State_Space(F1)



# -------------------------
# Damien Lefebvre's example
# This large system ($n=10$) comes from control theory: the filter is used as a controller for an active control of vehicle longitudinal oscillation~\cite{Lefe03}
d = loadmat('exDL.mat')['DL_Cor']
A, B, C, D = [mat(d[x]) for x in ('a', 'b', 'c', 'd')]
S2 = State_Space(Filter(A=A, B=B.transpose(), C=C, D=D, name='longitudinal oscillation controller DL'))








