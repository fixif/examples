from fixif import *

import numpy

from fixif.FXPF import FXPF_ABCD, FXPF_readme

from fixif.SIF import Realization, SIF

from fixif.Structures import State_Space, DFI, DFII, LWDF, rhoDFIIt, iterStructuresAndOptions

from fixif.LTI import Gabarit, iter_random_Gabarit, random_Gabarit, Filter, random_Filter
from fixif.LTI import dTF, dSS, dTFmp, dSSmp

import sollya

import matplotlib.pyplot as plt

from pytest import mark

fakeFilter = random_Filter(2, 1, 1)

structlist = list()
resultlist = list()

@mark.parametrize("type", ('butter'))
@mark.parametrize("wl", (16, ))
@mark.parametrize("SandO", iterStructuresAndOptions(fakeFilter), ids=lambda x: x[0].name)
def test_FXPF(SandO, wl, type, method='python'):
    # create the initial transfer function (designMargin = 0)

    g = Gabarit(48000, [(0, 9600), (12000, None)], [(0.95, 1.05), -20])
    #g = Gabarit(48000, [(0, 9600), (12000, None)], [(0, -1), -20])
    H = g.to_dTF(method=method, ftype="butter", designMargin=0)
    filt = Filter(tf=H)
    # build the realizatoin from the structure and options
    st, options = SandO
    if options:
        R = st.makeRealization(filt, **options)
    else:
        R = st.makeRealization(filt)

    structlist.append(R.structureName)
    Rapprox = R.quantize(wl)



    # check if realization in Gabarit
    print('------> Implementing Realization: %s with coeff quantization %s' % (R.structureName, wl))

    SS = Rapprox.dSS
    order = SS.A.shape[0]
    outputs, inputs = SS.D.shape
    u_bound = numpy.array([1.0])

    wordlengths = list()
    errors = list()

    for w in range(16, 10, -1):
        try:
            wl = w * numpy.ones([order + outputs, ], dtype=numpy.int64)
            msb, lsb, error, additionalSteps = FXPF_ABCD(SS.A, SS.B, SS.C, SS.D, u_bound, wl)
            wordlengths.append(w)
            errors.append(error)
            print("\n FXPF for wl=%s are determined using %s additional steps.") % (wl, additionalSteps)
            #print ("Wordlengths: %s") % (wl)
            #print ("MSBs: %s") % (msb)
            print("LSBs: %s") % (lsb)
            print("Errors: %s") % (error)
        except:
            pass

    output_errors = [err[0,-1] for err in errors]
    resultlist.append(zip(wordlengths, output_errors))

    print("title= '' ")
    print ("structs=")%(structlist)
    print (" ")
    print("results=")%(resultlist)
        #for state in range(0, order):
        #plt.plot(wordlengths, [err[0, state] for err in errors], label="state_" + str(state + 1))

        #plt.plot(wordlengths, [err[0, -1] for err in errors], label="output", linewidth=2, marker='o')
        #plt.title("Realization " + R.structureName)
        #plt.legend()
        #plt.show()

