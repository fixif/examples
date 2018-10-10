# coding=utf-8

import numpy as np
from scipy.signal import ellip
from fixif.LTI import Filter
from fixif.Structures import State_Space

quant = np.vectorize(lambda x, q: round(x*q)/q)


def randBetween(a, b):
	return np.random.rand()*(b-a)+a


def randomElliptic(seed, display=False):
	"""
	generate a 'random' ABCD 3x3 with good properties
	return a dSS, its WCPG, the result dictionary and a boolean (True if good properties)
	"""
	# set the seed
	np.random.seed(seed)
	N = np.random.randint(3, 8)
	rp = quant(randBetween(0, 10), 8)
	rs = np.random.randint(40, 80)
	Wn = [quant(np.random.rand(), 256), quant(np.random.rand(), 256)]
	typ = np.random.choice(['bandpass', 'bandstop'])
	# get the tf
	try:
		b, a = ellip(N, rp, rs, Wn, typ)
	except ValueError:
		return None, 0, {}, False
	# build a State-Space
	F = Filter(num=b, den=a)
	S = State_Space(F)
	# display if necessary
	if display:
		print("N=", N)
		print("rp=", rp)
		print("rs=", rs)
		print("Wn=", Wn)
		print("type=", typ)
		print("A=", S.dSS.A)
		print("B=", S.dSS.B)
		print("C=", S.dSS.C)
		print("D=", S.dSS.D)
	# try to compute the WCPG
	l = np.linalg.eig(S.A)[0]
	if all(np.vectorize(lambda x: abs(x)<(1-1e-8))(l)):
		res = {}
		try:
			W = S.dSS.WCPG(res)
		except ValueError:
			return S, 0, {}, False
		else:
			return S, W, res, res['N']>5000
	else:
		return S, 0, {}, False