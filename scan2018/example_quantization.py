from fixif import *

import numpy



from fixif.FXPF import FXPF_ABCD, FXPF_readme

from fixif.SIF import Realization, SIF

from fixif.Structures import State_Space, DFI, DFII, LWDF, rhoDFII, LGS

from fixif.LTI import Gabarit, iter_random_Gabarit, random_Gabarit, Filter
from fixif.LTI import dTF, dSS, dTFmp, dSSmp
from scipy.signal import  freqz
import sollya

import matplotlib.pyplot as plt
from matplotlib2tikz import save as tikz_save
import numpy


def plotQuantizationErrorEffects(Realizations, g):
    for R in Realizations:

        for w in range(16, 7, -1):
            minG = -1
            Rq = R.quantize(w)

            Hq = Rq.to_dTF()
            # try:
            #     print ("Realization %s quantized to %d bits: spectral radius is %s \n") % (
            #     R.structureName, w, numpy.max(numpy.abs(numpy.linalg.eigvals(Rq.dSS.A))))
            #     print ("Realization %s quantized to %d bits: WCPG = %f") % (R.structureName, w, Rq.dSS.WCPG())
            # except:
            #     pass
            #

            fig = plt.figure(w)
            ax = fig.add_subplot(111)

            omega, h = freqz(Hq.num.transpose(), Hq.den.transpose())
            if type(R) == rhoDFII:
                ax.plot((g._Fs * 0.5 / numpy.pi) * omega, 20 * numpy.log10(abs(h)), label="rhoDFII (transposed:True)")
            else:
                ax.plot([(g._Fs * 0.5 / numpy.pi) * omega[i] for i in range(0, 206)],
                        [20 * numpy.log10(abs(h[i])) for i in range(0, 206)], label=R.structureName + " " + str(w))

        minG = min(minG, min(20 * numpy.log10(abs(h))))
        plt.ylim((0.9, 1.1))
        currentAxis = plt.gca()
        currentAxis.add_patch(g._bands[0].Rectangle(True, minG))
        # for b in g._bands:
        #    currentAxis.add_patch(b.Rectangle(True, minG))

        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1, 1))
        ax.grid('on')
        fig.savefig("quantization" + str(w), bbox_extra_artists=(lgd,), bbox_inches='tight')

        from matplotlib2tikz import save as tikz_save

        tikz_save("quantization" + str(w) + ".tex")


def plotRoundingErrorEffects(Realizations):
    u_bar = [1.0]
    N = 100

    fig = plt.figure()
    ax = fig.add_subplot(111)


    #u = numpy.random.random([1, N])

    for R in Realizations:


        u = R.generate_inputs(u_bar, N)
        exact = R.simulate(u)

        #ax.plot(range(0, N, 1), exact.transpose(), label=R.structureName + "exact")


        for w in range(15, 14, -1):
            w = 10
            #Rq = R.quantize(w)
            Rq=R
            u = Rq.generate_inputs(u_bar, N)

            print("Plotting Realization %s for %d bits") % (Rq.structureName, w)

            rounded = Rq.dSS.to_dSSmp().simulate_rounded(numpy.matrix(u), prec=w)
            rounded = numpy.array([float(r) for r in rounded])

            ax.plot(range(0, N, 1), [exact[0,i] - rounded[i] for i in range(0, rounded.shape[0])], label=R.structureName+str(w))
            # handles, labels = ax.get_legend_handles_labels()
            # lgd = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1, 1))
            # ax.grid('on')
            # fig.savefig("rounding" + str(w) + R.structureName, bbox_extra_artists=(lgd,), bbox_inches='tight')

    print("Doing the handles")
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1, 1))
    ax.grid('on')
    fig.savefig("rounding"+str(w), bbox_extra_artists=(lgd,), bbox_inches='tight')

    tikz_save("rounding"+str(w) + ".tex")



g = Gabarit(48000, [(0, 9600), (12000, None)], [(0.95, 1.05), -20])
print("Frequency specifications: \n %s") % g

H = g.to_dTF(ftype="butter", method="scipy")
print ("\nTransfer function order is %d \n %s") % (H.order, H)

#g.plot(H)

#testing a few realizations
R_SS = State_Space(Filter(tf=H))
R_SS_balanced = State_Space(Filter(tf=H), form="balanced")
R_DFII = DFII(Filter(tf=H), transposed=True)
R_rhoDFII = rhoDFII(Filter(tf=H), transposed=True, scaling="l2")
R_LGS = LGS(Filter(tf=H), transposed=True)

Realizations = [R_SS, R_SS_balanced, R_DFII,  R_LGS]


plotRoundingErrorEffects(Realizations)
print ("Yippi")










