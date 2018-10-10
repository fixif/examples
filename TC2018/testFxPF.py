# coding=utf-8


"""
Among some small simple matrices (3X3 matrices, coefficients with few bits)
find some with log2(WCPG) != log2(WCPG_naive)
"""


from genRandom import random_dSS
from fixif.LTI import random_Elliptic
from fixif.FXPF import FXPF_ABCD
import numpy as np
from numpy import matrix as mat



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




def testFxPF(S):
	# display
	print(S)

	# naive MSB spec
	#TODO:

	# use FxPF
	wl = mat([w, ] * S.n)
	try:
		FXPF_ABCD(S.dSS.A, S.dSS.B, S.dSS.C, S.dSS.D, u_bound, wl)
	except:
		pass




u_bound = mat(10)
w = 8


# test on several random dSS
# for f, t in [(random_dSS, 'random 3x3'), (random_Elliptic, 'elliptic')]:
# 	print("type=", t)
# 	for S, W, res, seed in iter_somedSS(f, 10, 12345):
# 		testFxPF(S)

# give an error
#  ERROR: Could not compute MSB: log2(1 - 2^(1-w)) is NaN
testFxPF(random_dSS(13289)[0])
