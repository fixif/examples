# coding=utf-8

from fixif.LTI import random_dSS, Filter
from fixif.Structures import State_Space

from string import Template
#from mpmath import floor, log, ceil
from numpy import zeros, ones, matrix, power, ndenumerate, kron, multiply, nditer
from numpy.random import seed
from matplotlib import pyplot as plt

from numpy import floor, log2, ceil


# functions to generate AMPL parameters
def generateAMPLParam(name, val, format="%f"):
	"""generate AMPL code for the parameter"""
	if isinstance(val, int):
		return "param %s := %d;" % (name, val)
	if val.shape[0] == 1 or val.shape[1] == 1:
		# 1d vector
		return "param %s := \n%s;" % (name, "\n".join(("%d "+format) % (i[0]+1, wm) for i, wm in ndenumerate(val)))
	else:
		# 2d array
		col = " ".join(str(x+1) for x in range(val.shape[1]))
		data = "\n".join(str(i+1)+" "+" ".join(format%x for x in nditer(val[i,:])) for i in range(val.shape[0]))
		return "param %s : %s := \n%s;" % (name, col, data)







# Define the system
seed(124)        # random I want to deal with
n, p, q = 5, 3, 2
S = random_dSS(n, p, q)
R = State_Space(Filter(ss=S))
ubar = 10*ones((q,1))

# constraints
wmax = 64*ones((R.l+R.n+R.p,1))
eps_final = 2**-6*ones((R.p,1))

# determining the MSB
# m = m_tilde + (1 if w<w_tilde else 0) component-wise
m_tilde = matrix(R.computeNaiveMSB(ubar)).transpose()
w_tilde = matrix(R.w_tilde(ubar)).transpose()

#print(m_tilde)
print(w_tilde)

# error
Weps = R.Hepsilon.WCPG()
E = multiply(Weps, kron(ones((R.p, 1)), power(2, m_tilde.transpose())))

# generate the .dat file
AMPLcode = {}
AMPLcode['def_p'] = generateAMPLParam('p', R.p)
AMPLcode['def_np'] = generateAMPLParam('np', R.p + R.n + R.l)
AMPLcode['def_u'] = generateAMPLParam('u', wmax, "%d")
AMPLcode['def_E'] = generateAMPLParam('E', E)
AMPLcode['def_eps'] = generateAMPLParam('eps', eps_final)
AMPLcode['def_wtilde'] = generateAMPLParam('wtilde', w_tilde, "%d")


with open("xmpl.dat.template") as f:
	AMPL = Template(f.read())
with open("xmpl.dat",'w') as f:
	f.write(AMPL.substitute(AMPLcode))



def delta(w):
	d = zeros(w.shape)
	for i in range(w.size):
		d[i] = 1 if w[i] < w_tilde[i] else 0
	return d


# determining the minimal word-length with uniform scheme
w_guess = int(max(ceil(log2(E*ones((R.l+R.n+R.p,1))) - log2(eps_final)))[0, 0])
wmin = w_guess if all( E*power(2, -w_guess + delta(w_guess*ones((R.l+R.n+R.p, 1)))) < eps_final ) else w_guess + 1
print(wmin)





# plot delta_y wrt to w
WL = range(3,32)
deltay = matrix(zeros((p,len(WL))))
dy = matrix(zeros((1,len(WL))))

for wl in WL:
	w = wl * ones((n + p, 1))
	# compute eps
	l = m_tilde + delta(w) - w + ones((n+p, 1))
	eps = power(2, l)
	# compute deltay
	deltay[:,wl-WL[0]] = Weps*eps
	dy[:,wl-WL[0]] = ones((1,p))*Weps*eps

#plt.semilogy(WL,deltay.transpose())
plt.semilogy(WL, dy.transpose())
plt.show()






