#TODO: give full description and move to examples/

from fixif.SIF import Realization

a1 = [1, 0.355407, 1.0]
a2 = [1, -1.091855, 1.0]
b1 = [1.66664, -0.75504]
b2 = [1.58007, -0.92288]

J = [[1, 0, 0], [-a1[0], 1, 0], [0, -1, 1]]
K = [[1, 0, 0], [0, 0, 0], [0, 0, 1], [0, 0, 0]]
L = [[0, 0, a2[0]]]
M = [[b1[0], b1[1], 0, 0], [a1[1], a1[2], 0, 0], [0, 0, b2[0], b2[1]]]
N = [[0.01], [0], [0]]
P = [[0, 0, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 0]]
Q = [[0], [0], [0], [0]]
R = [[0, 0, a2[1], a2[2]]]
S = [[0]]


r = Realization(None, (J, K, L, M, N, P, Q, R, S), surnameVarT=['t1', 'y1', 't2'], surnameVarX=['d1[0]', 'd1[1]', 'd2[0]', 'd2[1]'])
#print(r.algorithmTxt(coefFormat=None, withSurname=True, withTime=False))

#H = r.to_dTF()
#print(H)

u = [1,2,5,3,5,12,8,]

print(r.runCdouble(u))
