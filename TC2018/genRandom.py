# coding=utf-8

from fixif.WCPG import WCPG_ABCD
from fixif.Structures import State_Space
from fixif.LTI import dSS
import numpy as np

quant = np.vectorize(lambda x,q: round(x*q)/q)


def random_dSS(seed):
	"""
	generate a 'random' ABCD 3x3 with good properties
	return a dSS, its WCPG, the result dictionary and a boolean (True if good properties)
	"""
	# set seed
	np.random.seed(seed)
	# build matrices
	A = np.diag([0.99, 0.8, 0.3])
	B = np.matrix('1; 0; 0')
	C = np.matrix('1 0 0')
	D = np.matrix(0)
	T = np.random.rand(3,3)
	Ap = np.linalg.inv(T)@A@T
	Bp = np.linalg.inv(T)@B
	Cp = C@T
	# simplified them
	App = quant(Ap,8)
	Bpp = quant(Bp,8)
	Cpp = quant(Cp,8)


	l = np.linalg.eig(App)[0]
	if all(np.vectorize(lambda x: abs(x)<(1-1e-8))(l)):
		res = {}
		try:
			W = WCPG_ABCD(App, Bpp, Cpp, D, res)
		except ValueError:
			return dSS(App, Bpp, Cpp, D), 0, {}, False
		else:
			return dSS(App, Bpp, Cpp, D), W, res, res['N']>5000
	else:
		return dSS(App, Bpp, Cpp, D), 0, {}, False

