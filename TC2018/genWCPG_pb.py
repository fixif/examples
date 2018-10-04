# coding=utf-8

""" Generate some A,B,C,D matrices
with spectral radius < 1 but where
 the WCPG cannot be computed (Gerhgorin's circle largen than 1 and norm(T)_2>1
 """


from fixif.WCPG import WCPG_ABCD_res
import numpy as np
from numpy import matrix as mat

quant = np.vectorize(lambda x: round(x*8)/8)

for seed in range(10000):
	np.random.seed(seed)

	A = np.diag([0.99, 0.8, 0.3])
	B = np.matrix('1; 0; 0')
	C = np.matrix('1 0 0')
	D = np.matrix(0)
	T = np.random.rand(3,3)
	Ap = np.linalg.inv(T)@A@T
	Bp = np.linalg.inv(T)@B
	Cp = C@T
	App = quant(Ap)
	Bpp = quant(Bp)
	Cpp = quant(Cp)


	l = np.linalg.eig(App)[0]
	if all(np.vectorize(lambda x: abs(x)<(1-1e-8))(l)):
		try:
			W, res = WCPG_ABCD_res(App, Bpp, Cpp, D)
		except ValueError:

			print("A=", App)
			print("B=", Bpp)
			print("C=", Cpp)
			print("D=", D)

		# if res['N']>5000:
		# 	print('seed=', seed)
		# 	print("N=", res['N'])
		# 	print(l)
		# 	print("App=", App)
