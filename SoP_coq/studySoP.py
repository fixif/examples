
from numpy import ceil, floor, log2, array, r_, vectorize
from fixif.FxP import Constant, FPF


def mulFPF(F1, F2):
	"""FPF of a multiplication"""
	return FPF(wl=F1.wl+F2.wl, msb=F1.msb+F2.msb)


def strSOP( dest, FPF_dest, coefs, FPF_coef, var, FPF_var):

	sop = []
	dp = []
	for var, fv, co, fc in zip(var, FPF_var, coefs, FPF_coef):
		if co == 1:
			dp.append(' %s <%s>' % (var, repr(fv)))
			sop.append(fv)
		elif co == -1:
			dp.append('-%s <%s>' % (var, repr(fv)))
			sop.append(fv)
		elif co:
			dp.append( '%4g*%s <%s>' %  (co, var, repr(mulFPF(fv,fc))))
			sop.append(mulFPF(fv,fc))

	S = " + ".join(dp) if dp else '0'
	return '%s <%s> <- ' % (dest, repr(FPF_dest)) + S, sop


def SoP(R, u_bar, w_txyu, w_z, latex=None):

	# FPF of t, x, u and y
	zeta_bar = R.Hzeta.WCPG() * u_bar
	FPF_txy = array([Constant(value=z[0, 0], wl=w_txyu).FPF for z in zeta_bar])
	FPF_txu = r_[FPF_txy[:-1], Constant(value=u_bar, wl=w_txyu).FPF]

	# Coefficients
	FPF_Z = vectorize(lambda x: Constant(value=x, wl=w_z).FPF)(R.Zcomp)
	# print(FPF_Z)

	# var Names
	varT = [v.toStr() for v in R._varNameT]
	varX = [v.toStr() for v in R._varNameX]
	varU = [v.toStr() for v in R._varNameU]
	varTXU = varT + varX + varU
	varXp1 = [v.toStr(suffix='p' if R.isPnut() else '')	for v in R._varNameX]
	varY = [v.toStr() for v in R._varNameY]
	varTXY = varT + varXp1 + varY

	# SoP
	for i in range(R.Z.shape[0]):
		str, fsop = strSOP(varTXY[i], FPF_txy[i], R.Zcomp[i, :].tolist()[0], FPF_Z[i, :].tolist()[0], varTXU, FPF_txu)
		print(str)
		print('\n'.join(f.ParenthesisNotation() for f in fsop))



