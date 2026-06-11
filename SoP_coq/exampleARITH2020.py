# example for ARITH2020

from fixif.LTI import Butter
from fixif.Structures import DFII
from numpy import ceil, floor, log2, array, r_, vectorize
from fixif.FxP import Constant, FPF
from fixif.SoP import FxPSoP, SoP


def figure(sop, k, filename):
	print(sop)


	# # rouding using general odd-rounding
	# lsb = sorted(list(list(p.lsb for p in sop._productFPF)))
	# k = [0]*len(lsb)
	# k[-1] = sop._resFPF.lsb-2
	# for i in range(len(lsb)-1, 0, -1):
	# 	if lsb[i] == lsb[i-1]:
	# 		k[i-1] = k[i]
	# 	else:
	# 		k[i-1] = min(k[i], lsb[i])-1

	# # optimized algorithm
	# lsb = sorted(list(p.lsb for p in sop._productFPF))
	# lf = sop._resFPF.lsb
	#
	# for m in range(1,len(lsb)):
	# 	if lsb[m-1] <= lf-2 <= lsb[m]:
	# 		break
	# else:
	# 	m = len(lsb)
	#
	# k = [0] * len(lsb)
	# k[0] = lsb[1]-1
	# for i in range(0,m-3+1):
	# 	k[i+1] = lsb[i+2]-1
	# k[m-1] = lf-2
	# for i in range(len(lsb)-1,m+1-1,-1):
	# 	k[i-1] = lsb[i-1]



	with open(filename, 'w') as f:
		f.write(r"""
		\documentclass[tikz]{standalone}
		\begin{document}
		\begin{tikzpicture}""")
		# bits
		f.write(sop.sumLaTeX(colors=('blue!20', 'blue!20', 'blue!20'), axis=True, sort='lsb'))
		# rounding bars
		l = []
		for y, x in enumerate(k):
			l.extend([(-x, -y*1.3+1.15), (-x, -y*1.3-.15)])
		f.write(r"  \draw" + "--".join(["(%.2f,%.2f)" % p for p in l]) + "[color=orange,line width=2pt, line cap=round];")
		l = [l[-1],(-sop._resFPF.lsb, -len(sop._productFPF)*1.3+1.15), (-sop._resFPF.lsb, -len(sop._productFPF)*1.3-.45)]
		f.write(r"  \draw" + "--".join(["(%.2f,%.2f)" % p for p in l]) + "[color=purple,line width=2pt, line cap=round];")
		f.write(r"""
		\end{tikzpicture}
		\end{document}""")



# 1st SoP
F = Butter(9, 0.41)
R = DFII(F)

# ubar, wordlength
u_bar = 12
w_txyu = 16
w_z = 16

R.computeNaiveMSB(u_bar)
# FPF of t, x, u and y
zeta_bar = R.Hzeta.WCPG() * u_bar
FPF_txy = array([Constant(value=z[0, 0], wl=w_txyu).FPF for z in zeta_bar])
FPF_txu = r_[FPF_txy[:-1], Constant(value=u_bar, wl=w_txyu).FPF]

# Coefficients
Z = vectorize(lambda x: Constant(value=x, wl=w_z))(R.Zcomp)

# var Names
varT = [v.toStr() for v in R._varNameT]
varX = [v.toStr() for v in R._varNameX]
varU = [v.toStr() for v in R._varNameU]
varTXU = varT + varX + varU
varXp1 = [v.toStr(suffix='p' if R.isPnut() else '') for v in R._varNameX]
varY = [v.toStr() for v in R._varNameY]
varTXY = varT + varXp1 + varY

# list of SoP
lsop = []
for i in range(R.Z.shape[0]):
	lsop.append(FxPSoP(Z[i, :].tolist()[0], varTXU, FPF_txu, varTXY[i], FPF_txy[i]))
# consider lsop[0]
sop = lsop[0]
# force the MSB to 6 (arbitrary) and signed=False
sop._productFPF = [FPF(msb=6, lsb=v.lsb-1, signed=False) for v in sop._productFPF]      # HACK: lsb-1
sop._productFPF[-1]._lsb=-11    # HACK dernier terme
sop._resFPF = FPF(msb=6, lsb=sop._resFPF.lsb, signed=False)
figure(sop, [-31,-28,-26,-24,-23,-22,-22,-22,-12,-9],'sop1.tex')
for F in sop._productFPF:
	print(F)


# 2nd SoP (toy example)
size = 9
#FxPSoP should be used with constant, but here we just want to set the lsb...
sop = FxPSoP(Z[i, 0:size].tolist()[0], varTXU[0:size], FPF_txu[0:size], varTXY[0], FPF_txy[0])
lsb = [0,4,4,5,8,9,9,12,20]
msb = 30
# force the MSB to 6 (arbitrary) and signed=False
sop._productFPF = [FPF(msb=30, lsb=l, signed=False) for l in lsb]
sop._resFPF = FPF(msb=30, lsb=10, signed=False)
figure(sop, [3,4,4,7,8,8,9,9,12],'sop2.tex')

