# coding=utf-8


"""
Among some small simple matrices (3X3 matrices, coefficients with few bits)
find some with log2(WCPG) != log2(WCPG_naive)
"""


from genRandom import random_dSS
from fixif.LTI import random_Elliptic
import numpy as np





def naive_double_WCPG(S, nit):
	"""naive floating-point WCPG evaluation of the dSS S with only `nit iterations"""

	acc = np.zeros(S.D.shape)
	acc += np.absolute(S.D)
	powerA = np.matrix(np.eye(S.A.shape[1]))

	for i in range(0, nit):
		acc += np.absolute(S.C * powerA * S.B)
		powerA *= S.A
	return acc



def iter_somedSS(func, N, initSeed=0):
	"""
	iter some dSS with good properties
	"""
	seed = initSeed
	for _ in range(N):
		# find a good dSS
		good = False
		while not good:
			S, W, res, good = func(seed)
			seed += 1
		# send it
		yield S, W, res, seed


for f, t in [(random_dSS, 'random 3x3'), (random_Elliptic, 'elliptic')]:
	print("type=", t)
	for S, W, res, seed in iter_somedSS(f, 10, 12345):
		# exact WCPG
		#print("WCPG")
		l = np.ceil(np.log2(W[0, 0]))
		# dummy naive WCPG
		#print("very naive WCPG")
		W1 = naive_double_WCPG(S.dSS, 500)     # 500 iterations for the very naive implementation
		l1 = np.ceil(np.log2(W1[0, 0]))
		# double WCPG with good number of iterations
		#print("naive WCPG")
		W2 = naive_double_WCPG(S.dSS, res["N"])
		l2 = np.ceil(np.log2(W2[0, 0]))

		if len({l,l1,l2})>1:
			print("seed=", seed)
			print("S=",S)
			print("W=", W, "  log2=", l)
			print("W1=", W1, "  log2=", l1)
			print("W2=", W2, "  log2=", l2)

