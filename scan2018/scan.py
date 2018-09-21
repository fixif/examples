"""
Author: Anastasia Volkova
Date: 22 August 2018

The purpose of this script is to demonstrate the workflow of the FiXiF project for the SCAN2018 conference.

Titile: "FiXiF toolbox : validated numerics for sound digital filter implementations"

The idea is to show that the combination of different arithmetics and their rigorous usage can be applied to
signal processing. That it is important to consider everything to have optimal or optimized implementations
in the contexts that require strict guarantees on the system safety.
Make a disclaimer for the worst-case analysis.
Talk about the IA and AA that are applied directly (use the biblio from the LTICJR paper)

Point: we are exploiting domain-specific knowledge to do the implementations.

We are not only providing the results that are numerically reliable by construction,
but also ensure that our tool is based on the proofs.

An important point is the verificaiton of digital filters. We provided a good method for that.


IMPORTANT: I need to combine all sources of errors for the demonstration.
For rounding errors show the Asilomar results.
For quantization errors show the ARITH2017 paper.

DO NOT FORGET:
* references on the biblio starting from Balakrishnian.
* cite Silviu and his work for the TF for the FIR filters, and that we are working on the IIR cases now.

TODO:
* add the FXPF support
    => almost done, need to take care of the integration into the FiXiF toolbox.




The possible workflows to show are:

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
    Not sure because it did not work.
    Maybe it is better to emphasize the design space exploration with the TF and just mention the filter conversion.

    d) Do I show LWDF?
    Will need to confirm that the lwdf toolbox works

    e) Do I add the SOS support? It should not be difficult but very useful for different works...
    Maybe, depends on the time. Should first start doing other structures, and only then perhaps add the SOS.


II. Simulink to Code

    See answer c)



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



# gg = Gabarit(48000,[ (0,9600), (16400,None) ], [(0,-1), -20])
# bandstop
#g = Gabarit(48000, [(0, 9600), (12000, 14000), (16400, None)], [(0, -1), -20, (0, -1)])

g = Gabarit(48000, [(0, 9600), (12000, None)], [(0.95, 1.05), -20])

print(g)



H = g.to_dTF(ftype='butter', method='scipy')
print(H)
print ("Initial TF WCPG =  %s") % (H.WCPG())


#structures=list()

#structs = ['Li-Chu-Wu (transposed:False)', 'Li-Chu-Wu (transposed:False)', 'Li-Chu-Wu (transposed:True)', 'Li-Chu-Wu (transposed:True)', 'Direct Form II (transposed:False)', 'Direct Form II (transposed:False)', 'Direct Form II (transposed:True)', 'Direct Form II (transposed:True)', 'rho Direct Form II (equiv_dSS:False, transposed:True, scaling:None)', 'rho Direct Form II (equiv_dSS:False, transposed:True, scaling:None)', 'rho Direct Form II (equiv_dSS:False, transposed:True, scaling:l2)', 'rho Direct Form II (equiv_dSS:False, transposed:True, scaling:l2)', 'rho Direct Form II (equiv_dSS:False, transposed:True, scaling:l2-relaxed)', 'rho Direct Form II (equiv_dSS:False, transposed:True, scaling:l2-relaxed)', 'rho Direct Form II (equiv_dSS:False, transposed:False, scaling:None)', 'rho Direct Form II (equiv_dSS:False, transposed:False, scaling:None)', 'rho Direct Form II (equiv_dSS:False, transposed:False, scaling:l2)', 'rho Direct Form II (equiv_dSS:False, transposed:False, scaling:l2)', 'rho Direct Form II (equiv_dSS:False, transposed:False, scaling:l2-relaxed)', 'rho Direct Form II (equiv_dSS:False, transposed:False, scaling:l2-relaxed)', 'rho Direct Form II (equiv_dSS:True, transposed:True, scaling:None)', 'rho Direct Form II (equiv_dSS:True, transposed:True, scaling:None)', 'rho Direct Form II (equiv_dSS:True, transposed:True, scaling:l2)', 'rho Direct Form II (equiv_dSS:True, transposed:True, scaling:l2)', 'rho Direct Form II (equiv_dSS:True, transposed:True, scaling:l2-relaxed)', 'rho Direct Form II (equiv_dSS:True, transposed:True, scaling:l2-relaxed)', 'rho Direct Form II (equiv_dSS:True, transposed:False, scaling:None)', 'rho Direct Form II (equiv_dSS:True, transposed:False, scaling:None)', 'rho Direct Form II (equiv_dSS:True, transposed:False, scaling:l2)', 'rho Direct Form II (equiv_dSS:True, transposed:False, scaling:l2)', 'rho Direct Form II (equiv_dSS:True, transposed:False, scaling:l2-relaxed)', 'rho Direct Form II (equiv_dSS:True, transposed:False, scaling:l2-relaxed)', 'State-Space (form:None)', 'State-Space (form:None)', 'State-Space (form:ctrl)', 'State-Space (form:ctrl)', 'State-Space (form:obs)', 'State-Space (form:obs)', 'Direct Form I (transposed:False, nbSum:1)', 'Direct Form I (transposed:False, nbSum:1)', 'Direct Form I (transposed:False, nbSum:2)', 'Direct Form I (transposed:False, nbSum:2)', 'Direct Form I (transposed:True, nbSum:1)', 'Direct Form I (transposed:True, nbSum:1)', 'Direct Form I (transposed:True, nbSum:2)', 'Direct Form I (transposed:True, nbSum:2)', 'Li-Gevers-Sun (transposed:False)', 'Li-Gevers-Sun (transposed:False)', 'Li-Gevers-Sun (transposed:True)', 'Li-Gevers-Sun (transposed:True)']

#results = [[(16, 8.429448444553664e-05), (15, 0.00016858896889107328), (14, 0.00033717793778214656), (13, 0.00067435587556429312), (12, 0.0013487117511285862), (11, 0.0026974235022571725)], [(16, 0.00089613637153640566), (15, 0.0017922727430728113), (14, 0.0035845454861456226), (13, 0.0071690909722912453), (12, 0.014338181944582491), (11, 0.028676363889164981)], [(16, 8.429448444553664e-05), (15, 0.00016858896889107328), (14, 0.00033717793778214656), (13, 0.00067435587556429312), (12, 0.0013487117511285862), (11, 0.0026974235022571725)], [(16, 0.00089613637153640566), (15, 0.0017922727430728113), (14, 0.0035845454861456226), (13, 0.0071690909722912453), (12, 0.014338181944582491), (11, 0.028676363889164981)], [(16, 2.2610880833704232e-07), (15, 4.5221761667408464e-07), (14, 9.0443523334816928e-07), (13, 1.8088704666963386e-06), (12, 3.6177409333926771e-06), (11, 7.2354818667853542e-06)], [(16, 0.0026788384336079779), (15, 0.0053576768672159557), (14, 0.010715353734431911), (13, 0.021430707468863823), (12, 0.042861414937727646), (11, 0.085722829875455292)], [(16, 0.0019252460014695708), (15, 0.0038504920029391416), (14, 0.0077009840058782832), (13, 0.015401968011756566), (12, 0.030803936023513133), (11, 0.061607872047026266)], [(16, 0.00048743995025307262), (15, 0.00097487990050614524), (14, 0.0019497598010122905), (13, 0.003899519602024581), (12, 0.0077990392040491619), (11, 0.015598078408098324)], [(16, 0.019739012814791186), (15, 0.039478025629582372), (14, 0.078956051259164745), (13, 0.15791210251832949), (12, 0.31582420503665898), (11, 0.63164841007331796)], [(16, 0.011785992762994472), (15, 0.023571985525988944), (14, 0.047143971051977887), (13, 0.1018806373587692), (12, 0.20376127471753841), (11, 0.45130920104281164)], [(16, 0.0093284387694961956), (15, 0.018656877538992391), (14, 0.037313755077984782), (13, 0.074627510155969565), (12, 0.14925502031193913), (11, 0.29851004062387826)], [(16, 0.0017838113893133376), (15, 0.0035676227786266752), (14, 0.0071352455572533504), (13, 0.014270491114506701), (12, 0.028540982229013401), (11, 0.057081964458026803)], [], [(16, 0.01214209250392104), (15, 0.024284185007842081), (14, 0.048568370015684162), (13, 0.097136740031368324), (12, 0.21271566690329441), (11, 0.42543133380658882)], [], [(16, 0.0077037936319835493), (15, 0.015407587263967099), (14, 0.030815174527934197), (13, 0.061630349055868394), (12, 0.12326069811173679), (11, 0.24847452122347355)], [(16, 6.1999729941143338)], [(16, 0.0085638787449650311), (15, 0.017127757489930062), (14, 0.034255514979860124), (13, 0.068511029959720249), (12, 0.18467773852914851), (11, 0.50603460956328172)], [], [(16, 0.012036735374371672), (15, 0.024073470748743345), (14, 0.048146941497486689), (13, 0.096293882994973379), (12, 0.20835355374559936), (11, 0.41670710749119871)], [], [(16, 0.011785992762994472), (15, 0.023571985525988944), (14, 0.047143971051977887), (13, 0.1018806373587692), (12, 0.20376127471753841), (11, 0.45130920104281164)], [], [(16, 0.0017842707656194989), (15, 0.0035685415312389978), (14, 0.0071370830624779957), (13, 0.014274166124955991), (12, 0.028548332249911983), (11, 0.057096664499823965)], [], [(16, 0.01214209250392104), (15, 0.024284185007842081), (14, 0.048568370015684162), (13, 0.097136740031368324), (12, 0.21271566690329441), (11, 0.42543133380658882)], [], [(16, 0.0077050927212621561), (15, 0.015410185442524312), (14, 0.030820370885048624), (13, 0.061640741770097249), (12, 0.1232814835401945), (11, 0.24851609208038902)], [(16, 2.0532256664648791), (15, 4.1064513329297583), (14, 4.599926259393051), (13, 9.3653057758443499), (12, 18.855815599994443), (11, 35.66030433969604)], [(16, 0.0085634822052469521), (15, 0.017126964410493904), (14, 0.034253928820987808), (13, 0.068507857641975617), (12, 0.18466857621677343), (11, 0.50600902803290559)], [(16, 1.7499392150541522), (15, 3.5705179622343679), (14, 7.1420124869687367)], [(16, 0.012039819305284222), (15, 0.024079638610568443), (14, 0.048159277221136887), (13, 0.096318554442273774), (12, 0.20840620915618713), (11, 0.41681241831237426)], [(16, 2.2610373839394147e-07), (15, 4.5220747678788294e-07), (14, 9.0441495357576589e-07), (13, 1.8088299071515318e-06), (12, 3.6176598143030635e-06), (11, 7.2353196286061271e-06)], [(16, 0.0026786754398938524), (15, 0.0053573508797877048), (14, 0.01071470175957541), (13, 0.021429403519150819), (12, 0.042858807038301638), (11, 0.085717614076603277)], [(16, 2.2610373839394147e-07), (15, 4.5220747678788294e-07), (14, 9.0441495357576589e-07), (13, 1.8088299071515318e-06), (12, 3.6176598143030635e-06), (11, 7.2353196286061271e-06)], [(16, 0.0026786754398938524), (15, 0.0053573508797877048), (14, 0.01071470175957541), (13, 0.021429403519150819), (12, 0.042858807038301638), (11, 0.085717614076603277)], [(16, 0.0010311997176649108), (15, 0.0020623994353298216), (14, 0.0041247988706596432), (13, 0.0082495977413192863), (12, 0.016499195482638573), (11, 0.4611723523669734)], [(16, 0.0026786754398938654), (15, 0.0053573508797877308), (14, 0.010714701759575462), (13, 0.021429403519150923), (12, 0.042858807038301847), (11, 0.085717614076603693)], [], [], [], [], [], [(16, 0.000244140625), (15, 0.00048828125), (14, 0.0009765625), (13, 0.001953125), (12, 0.00390625), (11, 0.0078125)], [], [(16, 0.000244140625), (15, 0.00048828125), (14, 0.0009765625), (13, 0.001953125), (12, 0.00390625), (11, 0.0078125)], [(16, 0.0030667565331882474), (15, 0.0061335130663764947), (14, 0.012267026132752989), (13, 0.024534052265505979), (12, 0.052739065133003138), (11, 0.0056072248385281888)], [(16, 0.00048376125483948586), (15, 0.00096752250967897173), (14, 0.0056743211658752887), (13, 0.011348642331750577), (12, 0.022697284663501155), (11, 0.04539456932700231)], [(16, 0.0030667565331882474), (15, 0.0061335130663764947), (14, 0.012267026132752989), (13, 0.024534052265505979), (12, 0.052739065133003138), (11, 0.0056072248385281888)], [(16, 0.00048376125483948586), (15, 0.00096752250967897173), (14, 0.0056743211658752887), (13, 0.011348642331750577), (12, 0.022697284663501155), (11, 0.04539456932700231)]]
structlist=list()
resultlist = list()
R_SS = State_Space(Filter(tf=H))
R_SS_balanced = State_Space(Filter(tf=H), form="balanced")
R_DFII = DFII(Filter(tf=H), transposed=True)
R_LGS = LGS(Filter(tf=H), transposed=True)

Realizations = [R_SS_balanced]

for R in Realizations:
    print('\n>>>>>>>>>>>>>>>>Implementing Realization: %s' % (R.structureName))
    structlist.append(R.structureName)
    Rapprox = R.quantize(16)
    SS = Rapprox.dSS
    order = SS.A.shape[0]
    outputs, inputs = SS.D.shape
    u_bound = numpy.array([1.0])

    wordlengths = list()
    errors = list()

    w = 16
    while True:
        w = w - 1
        try:
            wl = w * numpy.ones([order + outputs, ], dtype=numpy.int64)
            msb, lsb, error, additionalSteps = FXPF_ABCD(SS.A, SS.B, SS.C, SS.D, u_bound, wl)
            wordlengths.append(w)
            errors.append(error)
            print("WCPG  %s") % (SS.WCPG())
            print("\nFXPF for wl=%s are determined using %s additional steps.") % (wl, additionalSteps)
            # print ("Wordlengths: %s") % (wl)
            # print ("MSBs: %s") % (msb)
            print("LSBs: %s") % (lsb)
            print("Errors: %s") % (error)
        except:
            break

    output_errors = [err[0, -1] for err in errors]
    resultlist.append(zip(wordlengths, output_errors))

    #print("title= '' ")
    #print ("structs= %s") % (structlist)
    #print (" ")
    #print("results= %s") % (resultlist)


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




for struct in ():

    R = struct(Filter(tf=tfM))
    SS = R.dSS
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
            print("\n FXPF are determined using %s additional steps.") % ((additionalSteps))
            print ("Wordlengths: %s") % (wl)
            print ("MSBs: %s") % (msb)
            print("LSBs: %s") % (lsb)
            print("Errors: %s") % (error)
        except:
            pass



    for state in range(0, order):
        plt.plot(wordlengths, [err[0, state] for err in errors], label="state_" + str(state + 1))

    plt.plot(wordlengths, [err[0, -1] for err in errors], label="output", linewidth=2, marker='o')
    plt.title("Realization " + R.structureName)
    plt.legend()
    plt.show()




