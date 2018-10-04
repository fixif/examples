# coding=utf-8

from fixif.WCPG import WCPG_ABCD_res
import numpy as np

quant = np.vectorize(lambda x,q: round(x*q)/q)


def iter_randomABCD(N, initSeed=0):
	"""
	generate N 'random' ABCD 3x3 with good properties
	"""
	seed = initSeed
	for _ in range(N):
		# find a good random ABCD
		good = False
		while not good:
			A,B,C,D,good = randomABCD(seed)
			seed += 1
		# send it
		yield A,B,C,D, seed


def randomABCD(seed):
	"""
	generate a 'random' ABCD 3x3 with good properties
	return A, B, C, D and a boolean (True if good properties)
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
		try:
			W, res = WCPG_ABCD_res(App, Bpp, Cpp, D)
		except ValueError:
			return App, Bpp, Cpp, D, False
		else:
			return App, Bpp, Cpp, D, res['N']>5000
	else:
		return App, Bpp, Cpp, D, False
