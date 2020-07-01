#!/usr/bin/env python
import sys, os

interactive = True
saveplots = False

import matplotlib
if not interactive:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import tables

from emit_analyze import *

####################################################################

basedir = os.environ["SCRATCH"] + "/elens/"

dir_1M_err0_mc0_scc0 = "offdiag_magiccomp0.00"

#dir_1M_err01_mc0_scc0 = "offdiag_1M_err01_magiccomp0.00"
#dir_1M_err01_foc_mc0_scc0 = "offdiag_1M_err01_foc_magiccomp0.00"

pltnum = 0
#keepopen = [18,19,20, 21,22,23,24,25,26]
keepopen = range(10)
def saveplt():
    global keepopen
    global pltnum
    pltname = "plot_emit_err_%04d"%pltnum
    if saveplots:
        print "saving plot: ", pltname
        plt.savefig(pltname)
    else:
        print "I would save plot: ", pltname, " put plot saving is disabled"
    # if not pltnum in keepopen:
    #     plt.close()
    pltnum = pltnum + 1
    if not interactive:
        plt.close()
    return


# 1% lattice error
dir_1M_err01_mc16_scc0p0 = "offdiag_1M_err01_magiccomp16_scc0.0.00"
dir_1M_err01_mc16_scc0p4 = "offdiag_1M_err01_magiccomp16_scc0.4.00"
dir_1M_err01_mc16_scc0p8 = "offdiag_1M_err01_magiccomp16_scc0.8.00"
dir_1M_err01_mc16_scc1p2 = "offdiag_1M_err01_magiccomp16_scc1.2.00"
dir_1M_err01_mc16_scc1p6 = "offdiag_1M_err01_magiccomp16_scc1.6.00"
dir_1M_err01_mc16_scc2p0 = "offdiag_1M_err01_magiccomp16_scc2.0.00"
dir_1M_err01_mc16_scc2p4 = "offdiag_1M_err01_magiccomp16_scc2.4.00"
dir_1M_err01_mc16_scc2p8 = "offdiag_1M_err01_magiccomp16_scc2.8.00"
dir_1M_err01_mc16_scc3p2 = "offdiag_1M_err01_magiccomp16_scc3.2.00"
dir_1M_err01_mc16_scc3p6 = "offdiag_1M_err01_magiccomp16_scc3.6.00"
dir_1M_err01_mc16_scc4p0 = "offdiag_1M_err01_magiccomp16_scc4.0.00"
dir_1M_err01_mc16_scc4p4 = "offdiag_1M_err01_magiccomp16_scc4.4.00"
dir_1M_err01_mc16_scc4p8 = "offdiag_1M_err01_magiccomp16_scc4.8.00"
dir_1M_err01_mc16_scc5p2 = "offdiag_1M_err01_magiccomp16_scc5.2.00"
dir_1M_err01_mc16_scc5p6 = "offdiag_1M_err01_magiccomp16_scc5.6.00"
dir_1M_err01_mc16_scc6p0 = "offdiag_1M_err01_magiccomp16_scc6.0.00"

onepct_emit_err00 = {}
onepct_emit_2M_err00 = {}
onepct_emit_4M_err00 = {}
onepct_emit_8M_err00 = {}
onepct_emit_16M_err00 = {}

onepct_emit_err01 = {}
onepct_emit_errm01 = {}
onepct_emit_err02 = {}
onepct_emit_err03 = {}
onepct_emit_err01_mc13 = {}
onepct_emit_err01_mc112 = {}
onepct_emit_err02_mc112 = {}
onepct_emit_4M_err01_mc16 = {}
onepct_emit_16M_err01_mc16 = {}
onepct_emit_err01_mc118 = {}

onepct_emit_err00[0] = np.load(basedir + "offdiag_magiccomp0.00/emittances.npy")[()]
onepct_emit_2M_err00[0] = np.load(basedir + "offdiag_2M_magiccomp0.00/emittances.npy")[()]
onepct_emit_4M_err00[0] = np.load(basedir + "offdiag_4M_magiccomp0.00/emittances.npy")[()]
onepct_emit_8M_err00[0] = np.load(basedir + "offdiag_8M_magiccomp0.00/emittances.npy")[()]
onepct_emit_16M_err00[0] = np.load(basedir + "offdiag_16M_magiccomp0.00/emittances.npy")[()]

for sccten in range(0, 61, 4):
    scc = sccten/10.0
    dir1M = "offdiag_1M_err01_magiccomp16_scc%3.1f.00"%scc
    onepct_emit_err01[sccten] = np.load(basedir + dir1M +
                                        "/emittances.npy")[()]

    dir1M = "offdiag_1M_errm01_magiccomp16_scc%3.1f.00"%scc
    onepct_emit_errm01[sccten] = np.load(basedir + dir1M +
                                        "/emittances.npy")[()]

    dir1M = "offdiag_1M_err02_magiccomp16_scc%3.1f.00"%scc
    onepct_emit_err02[sccten] = np.load(basedir + dir1M +
                                        "/emittances.npy")[()]

    dir1M = "offdiag_1M_err03_magiccomp16_scc%3.1f.00"%scc
    onepct_emit_err03[sccten] = np.load(basedir + dir1M +
                                        "/emittances.npy")[()]

    dir1M = "offdiag_1M_err01_magiccomp13_scc%3.1f.00"%scc
    onepct_emit_err01_mc13[sccten] = np.load(basedir + dir1M +
                                        "/emittances.npy")[()]

    dir1M = "offdiag_1M_err01_magiccomp112_scc%3.1f.00"%scc
    onepct_emit_err01_mc112[sccten] = np.load(basedir + dir1M +
                                        "/emittances.npy")[()]

    dir1M = "offdiag_1M_err02_magiccomp112_scc%3.1f.00"%scc
    onepct_emit_err02_mc112[sccten] = np.load(basedir + dir1M +
                                        "/emittances.npy")[()]

    dir1M = "offdiag_4M_err01_magiccomp16_scc%3.1f.00"%scc
    if os.path.isdir(basedir + dir1M):
        onepct_emit_4M_err01_mc16[sccten] = np.load(basedir + dir1M +
                                        "/emittances.npy")[()]

    dir1M = "offdiag_16M_err01_magiccomp16_scc%3.1f.00"%scc
    if os.path.isdir(basedir + dir1M):
        onepct_emit_16M_err01_mc16[sccten] = np.load(basedir + dir1M +
                                        "/emittances.npy")[()]

    dir1M = "offdiag_1M_err01_magiccomp118_scc%3.1f.00"%scc
    onepct_emit_err01_mc118[sccten] = np.load(basedir + dir1M +
                                        "/emittances.npy")[()]

plotgroup = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
#plotgroup = [1, 8, 9]
#plotgroup = [10]

################################
if 0 in plotgroup:
    plt.figure()
    plt.title("x RMS emittance growth err 0, 0.01, 0.02, 0.03, scc 0")
    emit0 = onepct_emit_err01[0].xemitrms[0]

    plt.plot(onepct_emit_err00[0].turns, onepct_emit_err00[0].xemitrms/emit0,
             label="lattice error 0.00", lw=2)
    plt.plot(onepct_emit_err01[0].turns, onepct_emit_err01[0].xemitrms/emit0,
             label="lattice error 0.01", lw=2)
    plt.plot(onepct_emit_err02[0].turns, onepct_emit_err02[0].xemitrms/emit0,
             label="lattice error 0.02", lw=2)
    plt.plot(onepct_emit_err03[0].turns, onepct_emit_err03[0].xemitrms/emit0,
             label="lattice error 0.03", lw=2)

    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()
    print "x RMS emittance growth for different lattice errors"
    print 0, onepct_emit_err00[0].xemitrms[-1]/emit0
    print 0.01, onepct_emit_err01[0].xemitrms[-1]/emit0
    print 0.02, onepct_emit_err02[0].xemitrms[-1]/emit0
    print 0.03, onepct_emit_err03[0].xemitrms[-1]/emit0
    print

    emit0 = onepct_emit_err01[0].xemit999[0]
    print "x 99.9% emittance for different lattice errors"
    print 0, onepct_emit_err00[0].xemit999[-1]/emit0
    print 0.01, onepct_emit_err01[0].xemit999[-1]/emit0
    print 0.02, onepct_emit_err02[0].xemit999[-1]/emit0
    print 0.03, onepct_emit_err03[0].xemit999[-1]/emit0
    print

    ####
    emit0 = onepct_emit_err01[0].xemitrms[0]
    plt.figure()
    plt.title("x RMS emittance growth err0, 1M, 2M, 4M, 8M, 16M scc 0")
    plt.plot(onepct_emit_err00[0].turns, onepct_emit_err00[0].xemitrms/emit0,
             label="1M particles")
    plt.plot(onepct_emit_2M_err00[0].turns, onepct_emit_2M_err00[0].xemitrms/emit0,
             label="2M particles")
    plt.plot(onepct_emit_4M_err00[0].turns, onepct_emit_4M_err00[0].xemitrms/emit0,
             label="4M particles")
    plt.plot(onepct_emit_8M_err00[0].turns, onepct_emit_8M_err00[0].xemitrms/emit0,
             label="8M particles")
    plt.plot(onepct_emit_16M_err00[0].turns, onepct_emit_16M_err00[0].xemitrms/emit0,
             label="16M particles")

    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()
    print "x RMS emittance growth for different numbers of MP"
    print "1M: ",onepct_emit_err00[0].xemitrms[-1]/emit0
    print "2M: ",onepct_emit_2M_err00[0].xemitrms[-1]/emit0
    print "4M: ",onepct_emit_4M_err00[0].xemitrms[-1]/emit0
    print "8M: ",onepct_emit_8M_err00[0].xemitrms[-1]/emit0
    print "16M: ",onepct_emit_16M_err00[0].xemitrms[-1]/emit0
    print

################################

if 1 in plotgroup:
    plt.figure()
    plt.title("x RMS emittance growth err 0.01 scc 0-6")
    for k in range(0, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01[k].xemitrms[0]
        plt.plot(onepct_emit_err01[k].turns,
                 onepct_emit_err01[k].xemitrms/emit0, label="err 0.01 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.01 scc 2.4-4.4")
    for k in [24, 28, 32, 36, 40, 44]:
        scc = k/10.0
        emit0 = onepct_emit_err01[k].xemitrms[0]
        plt.plot(onepct_emit_err01[k].turns,
                 onepct_emit_err01[k].xemitrms/emit0, label="err 0.01 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.01 scc 3.6-6")
    for k in [36, 40, 44, 48, 52, 56, 60]:
        scc = k/10.0
        emit0 = onepct_emit_err01[k].xemitrms[0]
        plt.plot(onepct_emit_err01[k].turns,
                 onepct_emit_err01[k].xemitrms/emit0, label="err 0.01 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()
    print "x RMS emittance growth err 01 scc4.4"
    print onepct_emit_err01[44].xemitrms[-1]/emit0
    print
    
################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.01 scc 0-6")
    for k in range(0, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01[k].xemit999[0]
        plt.plot(onepct_emit_err01[k].turns,
                 onepct_emit_err01[k].xemit999/emit0, label="err 0.01 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.01 scc 2.4-4.4")
    for k in range(24, 45, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01[k].xemit999[0]
        plt.plot(onepct_emit_err01[k].turns,
                 onepct_emit_err01[k].xemit999/emit0, label="err 0.01 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.01 scc 3.2-6")
    for k in range(32, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01[k].xemit999[0]
        plt.plot(onepct_emit_err01[k].turns,
                 onepct_emit_err01[k].xemit999/emit0, label="err 0.01 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()
    print "x 99.9% emittance growth err 01 scc3.2"
    print onepct_emit_err01[32].xemit999[-1]/emit0
    print

################################

if 2 in plotgroup:
    plt.figure()
    plt.title("x RMS emittance growth err 0.01 vs. err -0.01")
    for k in [0, 12, 24, 36]:
        scc = k/10.0
        emit0 = onepct_emit_err01[k].xemitrms[0]
        plt.plot(onepct_emit_err01[k].turns,
                 onepct_emit_err01[k].xemitrms/emit0, label="err 0.01 scc%3.1f"%scc, lw=2)
        plt.plot(onepct_emit_errm01[k].turns,
                 onepct_emit_errm01[k].xemitrms/emit0, label="err -0.01 scc%3.1f"%scc, lw=2)

    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.01 vs. err -0.01")
    for k in [0, 12, 24, 36]:
        scc = k/10.0
        emit0 = onepct_emit_err01[k].xemit999[0]
        plt.plot(onepct_emit_err01[k].turns,
                 onepct_emit_err01[k].xemit999/emit0, label="err 0.01 scc%3.1f"%scc, lw=2)
        plt.plot(onepct_emit_errm01[k].turns,
                 onepct_emit_errm01[k].xemit999/emit0, label="err -0.01 scc%3.1f"%scc, lw=2)

    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

if 3 in plotgroup:
    plt.figure()
    plt.title("x RMS emittance growth err 0.02 scc 0-6")
    for k in range(0, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err02[k].xemitrms[0]
        plt.plot(onepct_emit_err02[k].turns,
                 onepct_emit_err02[k].xemitrms/emit0, label="err 0.02 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.02 scc 2.4-4.4")
    for k in [24, 28, 32, 36, 40, 44]:
        scc = k/10.0
        emit0 = onepct_emit_err02[k].xemitrms[0]
        plt.plot(onepct_emit_err02[k].turns,
                 onepct_emit_err02[k].xemitrms/emit0, label="err 0.02 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.02 scc 3.6-6")
    for k in [36, 40, 44, 48, 52, 56, 60]:
        scc = k/10.0
        emit0 = onepct_emit_err02[k].xemitrms[0]
        plt.plot(onepct_emit_err02[k].turns,
                 onepct_emit_err02[k].xemitrms/emit0, label="err 0.02 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()
    print "x RMS emittance growth err 02 scc4.4"
    print onepct_emit_err02[44].xemitrms[-1]/emit0
    print

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.02 scc 0-6")
    for k in range(0, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err02[k].xemit999[0]
        plt.plot(onepct_emit_err02[k].turns,
                 onepct_emit_err02[k].xemit999/emit0, label="err 0.02 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.02 scc 2.4-4.4")
    for k in range(24, 45, 4):
        scc = k/10.0
        emit0 = onepct_emit_err02[k].xemit999[0]
        plt.plot(onepct_emit_err02[k].turns,
                 onepct_emit_err02[k].xemit999/emit0, label="err 0.02 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.02 scc 3.2-6")
    for k in range(32, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err02[k].xemit999[0]
        plt.plot(onepct_emit_err02[k].turns,
                 onepct_emit_err02[k].xemit999/emit0, label="err 0.02 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()
    print "x 99.9% emittance growth err 02 scc3.2"
    print onepct_emit_err02[32].xemit999[-1]/emit0
    print

################################

if 4 in plotgroup:
    plt.figure()
    plt.title("x RMS emittance growth err 0.03 scc 0-6")
    for k in range(0, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err03[k].xemitrms[0]
        plt.plot(onepct_emit_err03[k].turns,
                 onepct_emit_err03[k].xemitrms/emit0, label="err 0.03 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.03 scc 2.4-4.4")
    for k in [24, 28, 32, 36, 40, 44]:
        scc = k/10.0
        emit0 = onepct_emit_err03[k].xemitrms[0]
        plt.plot(onepct_emit_err03[k].turns,
                 onepct_emit_err03[k].xemitrms/emit0, label="err 0.03 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.03 scc 3.6-6")
    for k in [36, 40, 44, 48, 52, 56, 60]:
        scc = k/10.0
        emit0 = onepct_emit_err03[k].xemitrms[0]
        plt.plot(onepct_emit_err03[k].turns,
                 onepct_emit_err03[k].xemitrms/emit0, label="err 0.03 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()
    print "x RMS emittance growth err 03 scc4.4"
    print onepct_emit_err03[44].xemitrms[-1]/emit0
    print

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.03 scc 0-6")
    for k in range(0, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err03[k].xemit999[0]
        plt.plot(onepct_emit_err03[k].turns,
                 onepct_emit_err03[k].xemit999/emit0, label="err 0.03 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.03 scc 2.4-4.4")
    for k in range(24, 45, 4):
        scc = k/10.0
        emit0 = onepct_emit_err03[k].xemit999[0]
        plt.plot(onepct_emit_err03[k].turns,
                 onepct_emit_err03[k].xemit999/emit0, label="err 0.03 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.03 scc 3.2-6")
    for k in range(32, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err03[k].xemit999[0]
        plt.plot(onepct_emit_err03[k].turns,
                 onepct_emit_err03[k].xemit999/emit0, label="err 0.03 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()
    print "x 99.9% emittance growth err 03 scc3.2"
    print onepct_emit_err03[32].xemit999[-1]/emit0
    print

################################

if 5 in plotgroup:
    plt.figure()
    plt.title("x RMS emittance growth err 0.01 1/3 scc 0-6")
    for k in range(0, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01_mc13[k].xemitrms[0]
        plt.plot(onepct_emit_err01_mc13[k].turns,
                 onepct_emit_err01_mc13[k].xemitrms/emit0, label="err 0.01 1/3 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.01 1/3 scc 1.6-3.2")
    for k in [16, 20, 24, 28, 32]:
        scc = k/10.0
        emit0 = onepct_emit_err01_mc13[k].xemitrms[0]
        plt.plot(onepct_emit_err01_mc13[k].turns,
                 onepct_emit_err01_mc13[k].xemitrms/emit0, label="err 0.01 1/3 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.01 1/3 scc 3.2-6")
    for k in range(32, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01_mc13[k].xemitrms[0]
        plt.plot(onepct_emit_err01_mc13[k].turns,
                 onepct_emit_err01_mc13[k].xemitrms/emit0, label="err 0.01 1/3 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()
    print "x RMS emittance growth err 01 1/3 scc2.8"
    print onepct_emit_err01_mc13[28].xemitrms[-1]/emit0
    print

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.01 1/3 scc 0-6")
    for k in range(0, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01_mc13[k].xemit999[0]
        plt.plot(onepct_emit_err01_mc13[k].turns,
                 onepct_emit_err01_mc13[k].xemit999/emit0, label="err 0.01 1/3 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.01 1/3 scc 1.6-3.2")
    for k in range(16, 33, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01_mc13[k].xemit999[0]
        plt.plot(onepct_emit_err01_mc13[k].turns,
                 onepct_emit_err01_mc13[k].xemit999/emit0, label="err 0.01 1/3 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.01 1/3 scc 3.2-6")
    for k in range(32, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01_mc13[k].xemit999[0]
        plt.plot(onepct_emit_err01_mc13[k].turns,
                 onepct_emit_err01_mc13[k].xemit999/emit0, label="err 0.01 1/3 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()
    print "x 99.9% emittance growth err 01 1/3 scc2.4"
    print onepct_emit_err01_mc13[24].xemit999[-1]/emit0
    print

################################

if 6 in plotgroup:
    plt.figure()
    plt.title("x RMS emittance growth err 0.01 1/12 scc 0-6")
    for k in range(0, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01_mc112[k].xemitrms[0]
        plt.plot(onepct_emit_err01_mc112[k].turns,
                 onepct_emit_err01_mc112[k].xemitrms/emit0, label="err 0.01 1/12 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.01 1/12 scc 0-1.6")
    for k in range(0, 17, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01_mc112[k].xemitrms[0]
        plt.plot(onepct_emit_err01_mc112[k].turns,
                 onepct_emit_err01_mc112[k].xemitrms/emit0, label="err 0.01 1/12 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.01 1/12 scc 1.6-3.2")
    for k in [16, 20, 24, 28, 32]:
        scc = k/10.0
        emit0 = onepct_emit_err01_mc112[k].xemitrms[0]
        plt.plot(onepct_emit_err01_mc112[k].turns,
                 onepct_emit_err01_mc112[k].xemitrms/emit0, label="err 0.01 1/12 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.01 1/12 scc 3.2-6")
    for k in range(32, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01_mc112[k].xemitrms[0]
        plt.plot(onepct_emit_err01_mc112[k].turns,
                 onepct_emit_err01_mc112[k].xemitrms/emit0, label="err 0.01 1/12 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()
    print "x RMS emittance growth err 01 1/12 scc 0.8"
    print onepct_emit_err01_mc112[8].xemitrms[-1]/emit0
    print

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.01 1/12 scc 0-6")
    for k in range(0, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01_mc112[k].xemit999[0]
        plt.plot(onepct_emit_err01_mc112[k].turns,
                 onepct_emit_err01_mc112[k].xemit999/emit0, label="err 0.01 1/12 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.01 1/12 scc 0-1.6")
    for k in range(0, 17, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01_mc112[k].xemit999[0]
        plt.plot(onepct_emit_err01_mc112[k].turns,
                 onepct_emit_err01_mc112[k].xemit999/emit0, label="err 0.01 1/12 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.01 1/12 scc 1.6-3.2")
    for k in range(16, 33, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01_mc112[k].xemit999[0]
        plt.plot(onepct_emit_err01_mc112[k].turns,
                 onepct_emit_err01_mc112[k].xemit999/emit0, label="err 0.01 1/12 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.01 1/12 scc 3.2-6")
    for k in range(32, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01_mc112[k].xemit999[0]
        plt.plot(onepct_emit_err01_mc112[k].turns,
                 onepct_emit_err01_mc112[k].xemit999/emit0, label="err 0.01 1/12 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()
    print "x 99.9% emittance growth err 01 1/12 scc 1.2"
    print onepct_emit_err01_mc112[12].xemit999[-1]/emit0
    print

################################

if 7 in plotgroup:
    plt.figure()
    plt.title("x RMS emittance growth err 0.02 1/12 scc 0-6")
    for k in range(0, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err02_mc112[k].xemitrms[0]
        plt.plot(onepct_emit_err02_mc112[k].turns,
                 onepct_emit_err02_mc112[k].xemitrms/emit0, label="err 0.02 1/12 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.02 1/12 scc 0-1.6")
    for k in range(0, 17, 4):
        scc = k/10.0
        emit0 = onepct_emit_err02_mc112[k].xemitrms[0]
        plt.plot(onepct_emit_err02_mc112[k].turns,
                 onepct_emit_err02_mc112[k].xemitrms/emit0, label="err 0.02 1/12 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.02 1/12 scc 1.6-3.2")
    for k in [16, 20, 24, 28, 32]:
        scc = k/10.0
        emit0 = onepct_emit_err02_mc112[k].xemitrms[0]
        plt.plot(onepct_emit_err02_mc112[k].turns,
                 onepct_emit_err02_mc112[k].xemitrms/emit0, label="err 0.02 1/12 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.02 1/12 scc 3.2-6")
    for k in range(32, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err02_mc112[k].xemitrms[0]
        plt.plot(onepct_emit_err02_mc112[k].turns,
                 onepct_emit_err02_mc112[k].xemitrms/emit0, label="err 0.02 1/12 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()
    print "x RMS emittance growth err 02 1/12 scc 0.8"
    print onepct_emit_err02_mc112[8].xemitrms[-1]/emit0
    print

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.02 1/12 scc 0-6")
    for k in range(0, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err02_mc112[k].xemit999[0]
        plt.plot(onepct_emit_err02_mc112[k].turns,
                 onepct_emit_err02_mc112[k].xemit999/emit0, label="err 0.02 1/12 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.02 1/12 scc 0-1.6")
    for k in range(0, 17, 4):
        scc = k/10.0
        emit0 = onepct_emit_err02_mc112[k].xemit999[0]
        plt.plot(onepct_emit_err02_mc112[k].turns,
                 onepct_emit_err02_mc112[k].xemit999/emit0, label="err 0.02 1/12 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.02 1/12 scc 1.6-3.2")
    for k in range(16, 33, 4):
        scc = k/10.0
        emit0 = onepct_emit_err02_mc112[k].xemit999[0]
        plt.plot(onepct_emit_err02_mc112[k].turns,
                 onepct_emit_err02_mc112[k].xemit999/emit0, label="err 0.02 1/12 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.02 1/12 scc 3.2-6")
    for k in range(32, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err02_mc112[k].xemit999[0]
        plt.plot(onepct_emit_err02_mc112[k].turns,
                 onepct_emit_err02_mc112[k].xemit999/emit0, label="err 0.02 1/12 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()
    print "x 99.9% emittance growth err 02 1/12 scc 1.6"
    print onepct_emit_err02_mc112[16].xemit999[-1]/emit0
    print

################################

if 8 in plotgroup:
    plt.figure()
    plt.title("x RMS emittance growth 4M err 0.01 1/6 scc 0, 3.6-6")
    for k in [0, 32, 36, 40, 44, 48, 52, 56, 60]:
        scc = k/10.0
        emit0 = onepct_emit_4M_err01_mc16[k].xemitrms[0]
        plt.plot(onepct_emit_4M_err01_mc16[k].turns,
                 onepct_emit_4M_err01_mc16[k].xemitrms/emit0,
                 label="err 0.01 4M 1/6 scc%3.1f"%scc, lw=2)
        plt.xlabel("turns")
        plt.ylabel("emittance growth")
        plt.legend(loc='best')
    saveplt()
    print "x RMS emittance growth 4M err 01 scc4.4"
    print onepct_emit_4M_err01_mc16[44].xemitrms[-1]/emit0
    print

    plt.figure()
    plt.title("x 99.9% emittance growth 4M err 0.01 1/6 scc 0, 3.6-6")
    for k in [0, 32, 36, 40, 44, 48, 52, 56, 60]:
        scc = k/10.0
        emit0 = onepct_emit_4M_err01_mc16[k].xemit999[0]
        plt.plot(onepct_emit_4M_err01_mc16[k].turns,
                 onepct_emit_4M_err01_mc16[k].xemit999/emit0,
                 label="err 0.01 4M 1/6 scc%3.1f"%scc, lw=2)
        plt.xlabel("turns")
        plt.ylabel("emittance growth")
        plt.legend(loc='best')
    saveplt()
    print "x 99.9% emittance growth 4M err 01 1/6 scc 3.6"
    print onepct_emit_4M_err01_mc16[36].xemit999[-1]/emit0
    print

################################

if 9 in plotgroup:
    plt.figure()
    plt.title("x RMS emittance growth 16M err 0.01 1/6 scc 0, 3.6-6")
    for k in [0, 32, 36, 40, 44, 48, 52, 56, 60]:
        scc = k/10.0
        emit0 = onepct_emit_16M_err01_mc16[k].xemitrms[0]
        plt.plot(onepct_emit_16M_err01_mc16[k].turns,
                 onepct_emit_16M_err01_mc16[k].xemitrms/emit0,
                 label="err 0.01 16M 1/6 scc%3.1f"%scc, lw=2)
        plt.xlabel("turns")
        plt.ylabel("emittance growth")
        plt.legend(loc='best')
    saveplt()
    print "x RMS emittance growth 16M err 01 scc4.4"
    print onepct_emit_16M_err01_mc16[44].xemitrms[-1]/emit0
    print

    plt.figure()
    plt.title("x 99.9% emittance growth 16M err 0.01 1/6 scc 0, 3.6-6")
    for k in [0, 32, 36, 40, 44, 48, 52, 56, 60]:
        scc = k/10.0
        emit0 = onepct_emit_16M_err01_mc16[k].xemit999[0]
        plt.plot(onepct_emit_16M_err01_mc16[k].turns,
                 onepct_emit_16M_err01_mc16[k].xemit999/emit0,
                 label="err 0.01 16M 1/6 scc%3.1f"%scc, lw=2)
        plt.xlabel("turns")
        plt.ylabel("emittance growth")
        plt.legend(loc='best')
    saveplt()
        
    print "x 99.9% emittance growth 16M err 01 1/6 scc 3.2"
    print onepct_emit_16M_err01_mc16[32].xemit999[-1]/emit0
    print

################################

if 10 in plotgroup:
    plt.figure()
    plt.title("x RMS emittance growth err 0.01 1/18 scc 0-6")
    for k in range(0, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01_mc118[k].xemitrms[0]
        plt.plot(onepct_emit_err01_mc118[k].turns,
                 onepct_emit_err01_mc118[k].xemitrms/emit0, label="err 0.01 1/18 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.01 1/18 scc 0-1.6")
    for k in range(0, 17, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01_mc118[k].xemitrms[0]
        plt.plot(onepct_emit_err01_mc118[k].turns,
                 onepct_emit_err01_mc118[k].xemitrms/emit0, label="err 0.01 1/18 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

    print "x RMS emittance growth err 01 1/18 scc 1.6"
    print onepct_emit_err01_mc118[16].xemitrms[-1]/emit0
    print

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.01 1/18 scc 1.6-3.2")
    for k in [16, 20, 24, 28, 32]:
        scc = k/10.0
        emit0 = onepct_emit_err01_mc118[k].xemitrms[0]
        plt.plot(onepct_emit_err01_mc118[k].turns,
                 onepct_emit_err01_mc118[k].xemitrms/emit0, label="err 0.01 1/18 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.01 1/18 scc 3.2-6")
    for k in range(32, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01_mc118[k].xemitrms[0]
        plt.plot(onepct_emit_err01_mc118[k].turns,
                 onepct_emit_err01_mc118[k].xemitrms/emit0, label="err 0.01 1/18 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()
    print "x RMS emittance growth err 01 1/18 scc 1.6"
    print onepct_emit_err01_mc118[16].xemitrms[-1]/emit0
    print

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.01 1/18 scc 0-6")
    for k in range(0, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01_mc118[k].xemit999[0]
        plt.plot(onepct_emit_err01_mc118[k].turns,
                 onepct_emit_err01_mc118[k].xemit999/emit0, label="err 0.01 1/18 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

    print "x 99.9% emittance growth err 01 1/18 scc 0.8"
    print onepct_emit_err01_mc118[8].xemit999[-1]/emit0
    print

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.01 1/18 scc 0-1.6")
    for k in range(0, 17, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01_mc118[k].xemit999[0]
        plt.plot(onepct_emit_err01_mc118[k].turns,
                 onepct_emit_err01_mc118[k].xemit999/emit0, label="err 0.01 1/18 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.01 1/18 scc 1.6-3.2")
    for k in range(16, 33, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01_mc118[k].xemit999[0]
        plt.plot(onepct_emit_err01_mc118[k].turns,
                 onepct_emit_err01_mc118[k].xemit999/emit0, label="err 0.01 1/18 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.01 1/18 scc 3.2-6")
    for k in range(32, 61, 4):
        scc = k/10.0
        emit0 = onepct_emit_err01_mc118[k].xemit999[0]
        plt.plot(onepct_emit_err01_mc118[k].turns,
                 onepct_emit_err01_mc118[k].xemit999/emit0, label="err 0.01 1/18 scc%3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()
    print "x 99.9% emittance growth err 01 1/18 scc 1.2"
    print onepct_emit_err01_mc118[12].xemit999[-1]/emit0
    print

################################

plt.show()

sys.exit(0)


# -1% lattice error
dir_1M_errm01_mc16_scc0p0 = "offdiag_1M_errm01_magiccomp16_scc0.0.00"
dir_1M_errm01_mc16_scc0p4 = "offdiag_1M_errm01_magiccomp16_scc0.4.00"
dir_1M_errm01_mc16_scc0p8 = "offdiag_1M_errm01_magiccomp16_scc0.8.00"
dir_1M_errm01_mc16_scc1p2 = "offdiag_1M_errm01_magiccomp16_scc1.2.00"
dir_1M_errm01_mc16_scc1p6 = "offdiag_1M_errm01_magiccomp16_scc1.6.00"
dir_1M_errm01_mc16_scc2p0 = "offdiag_1M_errm01_magiccomp16_scc2.0.00"
dir_1M_errm01_mc16_scc2p4 = "offdiag_1M_errm01_magiccomp16_scc2.4.00"
dir_1M_errm01_mc16_scc2p8 = "offdiag_1M_errm01_magiccomp16_scc2.8.00"
dir_1M_errm01_mc16_scc3p2 = "offdiag_1M_errm01_magiccomp16_scc3.2.00"
dir_1M_errm01_mc16_scc3p6 = "offdiag_1M_errm01_magiccomp16_scc3.6.00"
dir_1M_errm01_mc16_scc4p0 = "offdiag_1M_errm01_magiccomp16_scc4.0.00"
dir_1M_errm01_mc16_scc4p4 = "offdiag_1M_errm01_magiccomp16_scc4.4.00"
dir_1M_errm01_mc16_scc4p8 = "offdiag_1M_errm01_magiccomp16_scc4.8.00"
dir_1M_errm01_mc16_scc5p2 = "offdiag_1M_errm01_magiccomp16_scc5.2.00"
dir_1M_errm01_mc16_scc5p6 = "offdiag_1M_errm01_magiccomp16_scc5.5.00"
dir_1M_errm01_mc16_scc6p0 = "offdiag_1M_errm01_magiccomp16_scc6.0.00"

# % 2% lattice error
dir_1M_err02_mc16_scc0p0 = "offdiag_1M_err02_magiccomp16_scc0.0.00"
dir_1M_err02_mc16_scc0p4 = "offdiag_1M_err02_magiccomp16_scc0.4.00"
dir_1M_err02_mc16_scc0p8 = "offdiag_1M_err02_magiccomp16_scc0.8.00"
dir_1M_err02_mc16_scc1p2 = "offdiag_1M_err02_magiccomp16_scc1.2.00"
dir_1M_err02_mc16_scc1p6 = "offdiag_1M_err02_magiccomp16_scc1.6.00"
dir_1M_err02_mc16_scc2p0 = "offdiag_1M_err02_magiccomp16_scc2.0.00"
dir_1M_err02_mc16_scc2p4 = "offdiag_1M_err02_magiccomp16_scc2.4.00"
dir_1M_err02_mc16_scc2p8 = "offdiag_1M_err02_magiccomp16_scc2.8.00"
dir_1M_err02_mc16_scc3p2 = "offdiag_1M_err02_magiccomp16_scc3.2.00"
dir_1M_err02_mc16_scc3p6 = "offdiag_1M_err02_magiccomp16_scc3.6.00"
dir_1M_err02_mc16_scc4p0 = "offdiag_1M_err02_magiccomp16_scc4.0.00"
dir_1M_err02_mc16_scc4p4 = "offdiag_1M_err02_magiccomp16_scc4.4.00"
dir_1M_err02_mc16_scc4p8 = "offdiag_1M_err02_magiccomp16_scc4.8.00"
dir_1M_err02_mc16_scc5p2 = "offdiag_1M_err02_magiccomp16_scc5.2.00"
dir_1M_err02_mc16_scc5p6 = "offdiag_1M_err02_magiccomp16_scc5.6.00"
dir_1M_err02_mc16_scc6p0 = "offdiag_1M_err02_magiccomp16_scc6.0.00"

#  3% lattice error
dir_1M_err03_mc16_scc0p0 = "offdiag_1M_err03_magiccomp16_scc0.0.00"
dir_1M_err03_mc16_scc0p4 = "offdiag_1M_err03_magiccomp16_scc0.4.00"
dir_1M_err03_mc16_scc0p8 = "offdiag_1M_err03_magiccomp16_scc0.8.00"
dir_1M_err03_mc16_scc1p2 = "offdiag_1M_err03_magiccomp16_scc1.2.00"
dir_1M_err03_mc16_scc1p6 = "offdiag_1M_err03_magiccomp16_scc1.6.00"
dir_1M_err03_mc16_scc2p0 = "offdiag_1M_err03_magiccomp16_scc2.0.00"
dir_1M_err03_mc16_scc2p4 = "offdiag_1M_err03_magiccomp16_scc2.4.00"
dir_1M_err03_mc16_scc2p8 = "offdiag_1M_err03_magiccomp16_scc2.8.00"
dir_1M_err03_mc16_scc3p2 = "offdiag_1M_err03_magiccomp16_scc3.2.00"
dir_1M_err03_mc16_scc3p6 = "offdiag_1M_err03_magiccomp16_scc3.6.00"
dir_1M_err03_mc16_scc4p0 = "offdiag_1M_err03_magiccomp16_scc4.0.00"
dir_1M_err03_mc16_scc4p4 = "offdiag_1M_err03_magiccomp16_scc4.4.00"
dir_1M_err03_mc16_scc4p8 = "offdiag_1M_err03_magiccomp16_scc4.8.00"
dir_1M_err03_mc16_scc5p2 = "offdiag_1M_err03_magiccomp16_scc5.2.00"
dir_1M_err03_mc16_scc5p6 = "offdiag_1M_err03_magiccomp16_scc5.6.00"
dir_1M_err03_mc16_scc6p0 = "offdiag_1M_err03_magiccomp16_scc6.0.00"



emit_err01_mc0_scc0 = np.load(basedir + dir_1M_err01_mc0_scc0 + "/emittances.npy")[()]

emit_err01_mc16_scc0p8 = np.load(basedir + dir_1M_err01_mc16_scc0p8 + "/emittances.npy")[()]
#emit_err01_mc16_scc1p6 = np.load(basedir + dir_1M_err01_mc16_scc1p6 + "/emittances.npy")[()]
emit_err01_mc16_scc2p4 = np.load(basedir + dir_1M_err01_mc16_scc2p4 + "/emittances.npy")[()]
emit_err01_mc16_scc3p2 = np.load(basedir + dir_1M_err01_mc16_scc3p2 + "/emittances.npy")[()]
emit_err01_mc16_scc4p0 = np.load(basedir + dir_1M_err01_mc16_scc4p0 + "/emittances.npy")[()]

emit_errm01_mc16_scc0 = np.load(basedir + dir_1M_errm01_mc16_scc0 + "/emittances.npy")[()]
emit_errm01_mc16_scc0p8 = np.load(basedir + dir_1M_errm01_mc16_scc0p8 + "/emittances.npy")[()]
emit_errm01_mc16_scc1p6 = np.load(basedir + dir_1M_errm01_mc16_scc1p6 + "/emittances.npy")[()]
emit_errm01_mc16_scc2p4 = np.load(basedir + dir_1M_errm01_mc16_scc2p4 + "/emittances.npy")[()]
emit_errm01_mc16_scc3p2 = np.load(basedir + dir_1M_errm01_mc16_scc3p2 + "/emittances.npy")[()]
emit_errm01_mc16_scc4p0 = np.load(basedir + dir_1M_errm01_mc16_scc4p0 + "/emittances.npy")[()]


plt.figure()
plt.title("compare err .01 to -.01")

emit0 = emit_err01_mc0_scc0.xemitrms[0]
emit999 = emit_err01_mc0_scc0.xemit999[0]

plt.plot(emit_err01_mc0_scc0.turns, emit_err01_mc0_scc0.xemitrms/emit0, label="err 0.01 scc0")
plt.plot(emit_errm01_mc16_scc0.turns, emit_errm01_mc16_scc0.xemitrms/emit0, label="err -0.01 scc0")

plt.plot(emit_err01_mc16_scc0p8.turns, emit_err01_mc16_scc0p8.xemitrms/emit0, label="err 0.01 scc0.8")
plt.plot(emit_errm01_mc16_scc0p8.turns, emit_errm01_mc16_scc0p8.xemitrms/emit0, label="err -0.01 scc0.8")

#plt.plot(emit_err01_mc16_scc1p6.turns, emit_err01_mc16_scc1p6.xemitrms/emit0, label="err 0.01 scc1.6")
#plt.plot(emit_errm01_mc16_scc1p6.turns, emit_errm01_mc16_scc1p6.xemitrms/emit0, label="err -0.01 scc1.6")

plt.plot(emit_err01_mc16_scc2p4.turns, emit_err01_mc16_scc2p4.xemitrms/emit0, label="err 0.01 scc2.4")
plt.plot(emit_errm01_mc16_scc2p4.turns, emit_errm01_mc16_scc2p4.xemitrms/emit0, label="err -0.01 scc2.4")


plt.xlabel("turns")
plt.ylabel("emittance growth")
plt.legend(loc='best')
plt.savefig("emit_growth_err01_errm01.png")

plt.show()

