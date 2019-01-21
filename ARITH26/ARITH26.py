# coding=utf-8

from fixif.LTI import random_dSS, Filter
from fixif.Structures import State_Space

from mpmath import floor, log, ceil
from numpy import ones


# system
n, p, q = 5, 3, 2
S = random_dSS(n, p, q)
R = State_Space(Filter(ss=S))
ubar = 10*ones((q,1))

# determining the MSB
w_tilde = R.w_tilde(ubar)

print(w_tilde)

