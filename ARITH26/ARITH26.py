# coding=utf-8

from fixif.LTI import random_dSS, Filter
from fixif.Structures import State_Space

from numpy import ones
from numpy.random import seed


# Define the system
seed(124)        # random I want to deal with
n, p, q = 5, 3, 2
S = random_dSS(n, p, q)
R = State_Space(Filter(ss=S))
ubar = 10*ones((q,1))

# constraints
eps_final = 2**-6*ones((R.p,1))



print(R.optimalUniformWL(ubar, eps_final))
R.optimalWL(ubar, eps_final, 'random.dat')

