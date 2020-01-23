# coding=utf-8

from fixif.LTI import random_dSS, Filter, random_Elliptic
from fixif.Structures import State_Space, iterStructuresAndOptions
from fixif.LTI import Butter
from fixif.Structures import DFI,DFII, rhoDFII

from numpy import ones
from numpy.random import seed
from loadmat import loadmat
from numpy import matrix as mat

from studySoP import SoP


# ---------------
# Xilinx example.
#This example is a 5\textsuperscript{th} order lowpass Elliptic filter. Its passband is between 0 and 1e-3 (normalized frequencies) with 0.1dB as allowed ripples, whereas its stopband is between 0.0018 and 1 (normalized frequencies) with 40dB as minimum attenuation. It has been designed with Filter Designer toolbox from Matlab.

F1 = Filter(den=[1, -4.989216395071318, 9.956976990745105, -9.935631971923312, 4.957198554483321, -0.9893271782278296],
			num=[0.0001878726842913545, -0.0005635670357698394, 0.0003756943544614141, 0.0003756943544614141, -0.0005635670357698394, 0.0001878726842913545], name="Xilinx")
S1 = State_Space(F1)

# -------------------------
# Damien Lefebvre's example
# This large system ($n=10$) comes from control theory: the filter is used as a controller for an active control of vehicle longitudinal oscillation~\cite{Lefe03}
d = loadmat('exDL.mat')['DL_Cor']
A, B, C, D = [mat(d[x]) for x in ('a', 'b', 'c', 'd')]
S2 = State_Space(Filter(A=A, B=B.transpose(), C=C, D=D, name='longitudinal'))   # 'longitudinal oscillation controller DL'


# SDR example (from Fig. 15 "Software-Defined Radio FPGA Cores: Building towards a Domain-Specific Language")
# Fs=10kHz, Fstop1=2.190, Fpass1=2.1972, Fpass2=2.1974, Fstop2=2.210, Astop1=200dB, Apsass=0.1dB, Astop2=200dB
# This filter comes from a testbench in Software-Defined Radio system~\cite[Fig. 15]{}. It a 6th order Butterworth filter designed with the following parameters: sampling frequency = 10 kHz, lower cutoff frequency = 2.190 kHz, higher cutoff frequency = 2.210 kHz, passband ripple = 0.1 dB, stopband attenuation = 200 dB.
d = loadmat('SDR.mat')
#WARNING: this filter is designed with SOS-structure. When converted to state-space (by matlab), I am not sure if the spectral radius is lower than 1
#(when fdatool makes the single-structure conversion, the filter is not stable anymore)
A, B, C, D = [mat(d[x]) for x in ('A', 'B', 'C', 'D')]
S3 = State_Space(Filter(A=A, B=B.transpose(), C=C, D=D, name='SDR'))            # Software-Defined Radio


# random Elliptic filter
F4 = random_Elliptic(n=(2,3), seed=12345, quant=32)
print(F4.details())
S4 = State_Space(F4)

# Exemple Hinaomoto ISCAS2006
A = mat([[0,0,0.072, 0, 1.5], [1, 0, 0.300, 0, 0.20],[0,1,-0.1, 0,0.9], [0,0,0,0,0.05],[0,0,0,1,0.4]])
B = mat([[1, 0],[0,0], [0,0], [0,1],[0,0]])
C = mat([[1.1,2.7,0.9,0.4,1.5],[2.1,3.1,0.3,0.2,0.1],[5.4,1.6,-1.7, -6.6, 3.0]])
D = mat([[1,0.8],[0.3,0.6],[0.5,0.4]])
S4 = State_Space(Filter(A=A,B=B,C=C,D=D,name='ISCAS2006'))


# exemple PC (origine ?)
A= mat("""
[ 0.94961650496175  -0.00086079266447   0.05695647856580  -0.00206512941531;
  0.00243644867557   0.97324248160060  -0.05014627397097   0.00221130704225;
  0.00003538548664   0.00072499225469   0.99953199828010  -0.00053267821142;
  0.07016414190869  -0.00232763195634  -0.06768748276232   0.98860602697171 ]""")

B= mat("""
[ 0.08390894871690  -0.00001461862848  -0.00079960496498   0.00048214761580   0.00209408900920;
  0.07052658459684  -0.00000001872741   0.00123734566619  -0.00072744090143  -0.00223686575539;
 -0.00005903324913   0.00000001232224  -0.00072727303651   0.00043455865711   0.00153266700275;
 -0.11711256225159   0.00002446893748  -0.00223443371381   0.00152963752289   0.01136052996550 ]""")

C= 1024*mat("""
[                0   0.99876617000543   0.00072504533047   0.00222099952949;
                 0   0.00072504533047   0.99956696682531  -0.00152133486438;
                 0   0.00222099952949  -0.00152133486438   0.98869654855310;
  1.00000000000000   0.00079594908648  -0.00048022213387  -0.00210894664156;
                 0   0.99876617000543   0.00072504533047   0.00222099952949;
                 0   0.00072504533047   0.99956696682531  -0.00152133486438;
                 0   0.00222099952949  -0.00152133486438   0.98869654855310 ]""")
D= 1024*mat("""
[                0                  0   0.00123382999457  -0.00072504533047  -0.00222099952949;
                 0                  0  -0.00072504533047   0.00043303317469   0.00152133486438;
                 0                  0  -0.00222099952949   0.00152133486438   0.01130345144690;
                 0                  0  -0.00079594908648   0.00048022213387   0.00210894664156;
                 0                  0   0.00123382999457  -0.00072504533047  -0.00222099952949;
                 0                  0  -0.00072504533047   0.00043303317469   0.00152133486438;
                 0                  0  -0.00222099952949   0.00152133486438   0.01130345144690 ]""")
S5 = State_Space(Filter(A=A,B=B,C=C,D=D,name='exPC'))


lS = [State_Space(Filter(ss=S5.dSS[i, :])) for i in range(7)]
lS.append(S5)

#
# realization
F = Butter(7, 0.1)
S6 = rhoDFII(F)

# ubar, wordlength
u_bar = 12
w_txyu = 16
w_z = 16

for R in (S1, S2, S5, S4, S5, S6):
		print("----- " + R.name + " -----")
		#print(R.Z)

		SoP(R, u_bar, w_txyu, w_z)

		#print(R.algorithmTxt('%.4e'))





