"""
Author: Anastasia Volkova
Date: 22 August 2018

The purpose of this script is to demonstrate the workflow of the FiXiF project for the SCAN2018 conference.

Titile: "FiXiF toolbox : validated numerics for sound digital filter implementations"

The idea is to show that the combination of different arithmetics and their rigorous usage can be applied to
signal processing. That it is important to consider everything to have optimal or optimized implementations
in the contexts that require strict guarantees on the system safety.
Make a disclaimer for the worst-case analysis.
Talk about the IA and AA that are applied directly

Point: we are exploiting domain-specific knowledge to do the implementations.

We are not only providing the results that are numerically reliable by construction,
but also ensure that our tool is based on the proofs.

An important point is the verificaiton of digital filters. We provided a good method for that.


A possible workflow to show is:

I. TF to Code
    1. Transfer function
    2. Different realizations
    3. Fixed-Point algorithm (using the Asilomar code ?)
    4. C code generation
    5. Verification of the filter specs?

    Decisions:
    a) Do I show the overall error? With the coefficient quantization and the rounding erorrs?
    It would be good, yes.


    b) Do I show the FloPoCo? Do I show the FIRopt?
    Yes, i go back to the compiler slide and show that each component is used in two more projects.
    Then, I have one slide that shows what those projects do. Not more than one high-level slide per project.
    Try to condense it into one slide.

    c) Do I show Simulink?
    Maybe it is better to emphasize the design space exploration with the TF and just mention the filter conversion.

    d) Do I show LWDF?
    Will need to confirm that the lwdf toolbox works

    e) Do I add the SOS support? It should not be difficult but very useful for different works...
    Maybe, depends on the time. Should first start doing other structures, and only then perhaps add the SOS.

"""

from fixif import *

import numpy



from fixif.FXPF import FXPF_ABCD, FXPF_readme

from fixif.SIF import Realization, SIF

from fixif.Structures import State_Space, DFI, DFII, LWDF, rhoDFIIt, LGS

from fixif.LTI import Gabarit, iter_random_Gabarit, random_Gabarit, Filter
from fixif.LTI import dTF, dSS, dTFmp, dSSmp

import sollya

import matplotlib.pyplot as plt


#Creating frequency specifications
# gg = Gabarit(48000,[ (0,9600), (16400,None) ], [(0,-1), -20])
# bandstop
#g = Gabarit(48000, [(0, 9600), (12000, 14000), (16400, None)], [(0, -1), -20, (0, -1)])
g = Gabarit(48000, [(0, 9600), (12000, None)], [(0.95, 1.05), -20])
print(g)


#designing corresponding transfer function
H = g.to_dTF(ftype='butter', method='scipy')
print(H)
#we can instantly compute the WCPG that corresponds to this "ideal" transfer function
print ("Initial TF WCPG =  %s") % (H.WCPG())


structlist=list()
resultlist = list()
#creating different realizations of the transfer function H
R_SS = State_Space(Filter(tf=H))
R_SS_balanced = State_Space(Filter(tf=H), form="balanced")
R_DFII = DFII(Filter(tf=H), transposed=True)
R_LGS = LGS(Filter(tf=H), transposed=True)


Realizations = [R_SS_balanced, R_SS, R_DFII, R_LGS]

for R in Realizations:
    print('\n>>>>>>>>>>>>>>>>Implementing Realization: %s' % (R.structureName))
    structlist.append(R.structureName)
    #quantizing coefficients of the sutrcture
    Rapprox = R.quantize(16)
    #extracting corresponding State-space
    SS = Rapprox.dSS

    order = SS.A.shape[0]
    outputs, inputs = SS.D.shape

    #here we set the bound on the input signals to 1.0
    u_bound = numpy.array([1.0])

    wordlengths = list()
    errors = list()

    #testing how much one can decrease the wordlength while guaranteeing that no overflow occurs
    w = 16
    while True:
        w = w - 1
        try:
            #here we use uniform worldengths but this is not mandatory
            wl = w * numpy.ones([order + outputs, ], dtype=numpy.int64)

            msb, lsb, error, additionalSteps = FXPF_ABCD(SS.A, SS.B, SS.C, SS.D, u_bound, wl)
            wordlengths.append(w)
            errors.append(error)
            print("\n\nWCPG  %s") % (SS.WCPG())
            print("\nFXPF for wl=%s are determined using %s additional steps.") % (wl, additionalSteps)
            print ("Wordlengths: %s") % (wl)
            print ("MSBs: %s") % (msb)
            print("LSBs: %s") % (lsb)
            print("Errors: %s") % (error)
        except:
            break

    output_errors = [err[0, -1] for err in errors]
    #saving the error bounds that correspond to the output of the filter
    #we could have saved errors on the state variables as well
    resultlist.append(zip(wordlengths, output_errors))

    #print("title= '' ")
    #print ("structs= %s") % (structlist)
    #print (" ")
    #print("results= %s") % (resultlist)


#plotting the results
wl = range(16, 10, -1)
plt.xlim(10, 17)
plt.ylim(0, 1)

fig = plt.figure(1)
ax = fig.add_subplot(111)
for i in range(0, len(structlist)-1):
    ax.plot([r[0] for r in resultlist[i]], [r[1] for r in resultlist[i]], label=structlist[i], marker='o', linestyle='')

handles, labels = ax.get_legend_handles_labels()
lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5,-0.1))
ax.grid('on')
#fig.show()
fig.savefig('samplefigure', bbox_extra_artists=(lgd,), bbox_inches='tight')



#this code snippet plots the error bounds for states and outputs for each realization in Realizations
#for implementations with w = 16...10 bits
# for R in Realizations:
#     SS = R.dSS
#     order = SS.A.shape[0]
#     outputs, inputs = SS.D.shape
#     u_bound = numpy.array([1.0])
#
#     wordlengths = list()
#     errors = list()
#
#     for w in range(16, 10, -1):
#         try:
#             wl = w * numpy.ones([order + outputs, ], dtype=numpy.int64)
#             msb, lsb, error, additionalSteps = FXPF_ABCD(SS.A, SS.B, SS.C, SS.D, u_bound, wl)
#             wordlengths.append(w)
#             errors.append(error)
#             print("\n FXPF are determined using %s additional steps.") % ((additionalSteps))
#             print ("Wordlengths: %s") % (wl)
#             print ("MSBs: %s") % (msb)
#             print("LSBs: %s") % (lsb)
#             print("Errors: %s") % (error)
#         except:
#             pass
#
#
#     for state in range(0, order):
#         plt.plot(wordlengths, [err[0, state] for err in errors], label="state_" + str(state + 1))
#
#     plt.plot(wordlengths, [err[0, -1] for err in errors], label="output", linewidth=2, marker='o')
#     plt.title("Realization " + R.structureName)
#     plt.legend()
#     plt.show()




