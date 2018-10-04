# coding=utf-8


"""
Among some small simple matrices (3X3 matrices, coefficients with few bits)
find some with log2(WCPG) != log2(WCPG_naive)
"""


from fixif.WCPG import WCPG_ABCD_res
from fixif.LTI import dSS
from genRandom import iter_randomABCD
import numpy as np





def naive_double_WCPG(A, B, C, D, nit):
	""""""

	acc = np.zeros(D.shape)
	acc += np.absolute(D)
	powerA = np.matrix(np.eye(A.shape[1]))

	for i in range(0, nit):
		acc += np.absolute(C * powerA * B)
		powerA *= A
	return acc



for A, B, C, D, seed in iter_randomABCD(1000,12345):
	# exact WCPG
	W, res = WCPG_ABCD_res(A, B, C, D)
	l = np.ceil(np.log2(W[0,0]))
	# dummy naive WCPG
	W1 = naive_double_WCPG(A, B, C, D, 500)
	l1 = np.ceil(np.log2(W1[0,0]))
	# double WCPG with good number of iterations
	W2 = naive_double_WCPG(A, B, C, D, res["N"])
	l2 = np.ceil(np.log2(W2[0,0]))

	if len({l,l1,l2})>1:
		print("seed=", seed)
		print("A=",A)
		print("B=",B)
		print("C=",C)
		print("D=",D)
		print("W=", W, "  log2=", l)
		print("W1=", W1, "  log2=", l1)
		print("W2=", W2, "  log2=", l2)

