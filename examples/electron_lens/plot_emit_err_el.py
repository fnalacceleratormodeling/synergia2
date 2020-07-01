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
    pltname = "plot_emit_err_el_%04d"%pltnum
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

onepct_emit_err01_el = {}
onepct_emit_err01_elgauss = {}
onepct_emit_err01_adaptive_el = {}
onepct_emit_err01_adaptive_elgauss = {}

for sccten in (0, 20, 40, 60, 80, 90, 100, 120, 140, 160, 180, 200):
    print "reading emittances for ", sccten
    scc = sccten/10.0
    dir1M = "offdiag_1M_err01_el%3.1f.00"%scc
    onepct_emit_err01_el[sccten] = np.load(basedir + dir1M +
                                        "/emittances.npy")[()]

for sccten in (0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100):
    print "reading emittances for ", sccten
    scc = sccten/10.0
    dir1M = "offdiag_1M_err01_elgauss%3.1f.00"%scc
    onepct_emit_err01_elgauss[sccten] = np.load(basedir + dir1M +
                                        "/emittances.npy")[()]

for sccten in (0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100):
    print "reading emittances for ", sccten
    scc = sccten/10.0
    dir1M = "offdiag_1M_err01_adaptive_el%3.1f.00"%scc
    onepct_emit_err01_adaptive_el[sccten] = np.load(basedir + dir1M +
                                        "/emittances.npy")[()]

for sccten in (0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100):
    print "reading emittances for ", sccten
    scc = sccten/10.0
    dir1M = "offdiag_1M_err01_adaptive_elgauss%3.1f.00"%scc
    onepct_emit_err01_adaptive_elgauss[sccten] = np.load(basedir + dir1M +
                                        "/emittances.npy")[()]

#plotgroup = [1, 2, 3]
plotgroup = [3, 4]
#plotgroup = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
#plotgroup = [1, 8, 9]
#plotgroup = [10]

################################

if 1 in plotgroup:
    plt.figure()
    plt.title("x RMS emittance growth err 0.01 el 0-9.0 A-m/cell")
    #for k in (0, 20, 40, 60, 80, 90, 100, 120, 140, 160, 180, 200):
    for k in (0, 20, 40, 60, 80, 90):
        print "plotting emittance growth ", k
        print onepct_emit_err01_el[k]
        scc = k/10.0
        if not hasattr(onepct_emit_err01_el[k], "xemitrms"):
            print "no data for ", scc
            continue
        emit0 = onepct_emit_err01_el[k].xemitrms[0]
        plt.plot(onepct_emit_err01_el[k].turns,
                 onepct_emit_err01_el[k].xemitrms/emit0, label="err 0.01 el %3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.01 el 0, 10.0-20.0 A-m/cell")
    #for k in (0, 10, 20, 30, 40, 45, 50, 60, 70, 80, 90, 100):
    for k in (0, 90, 100, 120, 140, 160, 180, 200):
        print "plotting emittance growth ", k
        print onepct_emit_err01_el[k]
        scc = k/10.0
        if not hasattr(onepct_emit_err01_el[k], "xemitrms"):
            print "no data for ", scc
            continue
        emit0 = onepct_emit_err01_el[k].xemitrms[0]
        plt.plot(onepct_emit_err01_el[k].turns,
                 onepct_emit_err01_el[k].xemitrms/emit0, label="err 0.01 el %3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.01 el 0-20.0 A-m/cell")
    for k in (0, 20, 40, 60, 80, 90, 100, 120, 140, 160, 180, 200):
        scc = k/10.0
        emit0 = onepct_emit_err01_el[k].xemit999[0]
        plt.plot(onepct_emit_err01_el[k].turns,
                 onepct_emit_err01_el[k].xemit999/emit0, label="err 0.01 el %3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

if 2 in plotgroup:
    plt.figure()
    plt.title("x RMS emittance growth err 0.01 elgauss 0-10.0 A-m/cell")
    for k in (0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100):
        print "plotting emittance growth ", k
        print onepct_emit_err01_elgauss[k]
        scc = k/10.0
        if not hasattr(onepct_emit_err01_elgauss[k], "xemitrms"):
            print "no data for ", scc
            continue
        emit0 = onepct_emit_err01_elgauss[k].xemitrms[0]
        plt.plot(onepct_emit_err01_elgauss[k].turns,
                 onepct_emit_err01_elgauss[k].xemitrms/emit0, label="err 0.01 el %3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.01 elgauss 0-6.0 A-m/cell")
    for k in (0, 10, 20, 30, 40, 50, 60):
        print "plotting emittance growth ", k
        print onepct_emit_err01_elgauss[k]
        scc = k/10.0
        if not hasattr(onepct_emit_err01_elgauss[k], "xemitrms"):
            print "no data for ", scc
            continue
        emit0 = onepct_emit_err01_elgauss[k].xemitrms[0]
        plt.plot(onepct_emit_err01_elgauss[k].turns,
                 onepct_emit_err01_elgauss[k].xemitrms/emit0, label="err 0.01 el %3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.01 elgauss 6.0-10.0 A-m/cell")
    for k in (60, 70, 80, 90, 100):
        print "plotting emittance growth ", k
        print onepct_emit_err01_elgauss[k]
        scc = k/10.0
        if not hasattr(onepct_emit_err01_elgauss[k], "xemitrms"):
            print "no data for ", scc
            continue
        emit0 = onepct_emit_err01_elgauss[k].xemitrms[0]
        plt.plot(onepct_emit_err01_elgauss[k].turns,
                 onepct_emit_err01_elgauss[k].xemitrms/emit0, label="err 0.01 el %3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.01 elgauss 0-10.0 A-m/cell")
    for k in (0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100):
        scc = k/10.0
        emit0 = onepct_emit_err01_elgauss[k].xemit999[0]
        plt.plot(onepct_emit_err01_elgauss[k].turns,
                 onepct_emit_err01_elgauss[k].xemit999/emit0, label="err 0.01 el %3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

if 3 in plotgroup:
    plt.figure()
    plt.title("x RMS emittance growth err 0.01 el adaptive 0-10.0 A-m/cell")
    for k in (0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100):
        print "plotting emittance growth ", k
        print onepct_emit_err01_adaptive_el[k]
        scc = k/10.0
        if not hasattr(onepct_emit_err01_adaptive_el[k], "xemitrms"):
            print "no data for ", scc
            continue
        emit0 = onepct_emit_err01_adaptive_el[k].xemitrms[0]
        plt.plot(onepct_emit_err01_adaptive_el[k].turns,
                 onepct_emit_err01_adaptive_el[k].xemitrms/emit0, label="err 0.01 adaptive el %3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.01 adaptive el 0-6.0 A-m/cell")
    for k in (0, 10, 20, 30, 40, 50, 60):
        print "plotting emittance growth ", k
        print onepct_emit_err01_adaptive_el[k]
        scc = k/10.0
        if not hasattr(onepct_emit_err01_adaptive_el[k], "xemitrms"):
            print "no data for ", scc
            continue
        emit0 = onepct_emit_err01_adaptive_el[k].xemitrms[0]
        plt.plot(onepct_emit_err01_adaptive_el[k].turns,
                 onepct_emit_err01_adaptive_el[k].xemitrms/emit0, label="err 0.01 adaptive el %3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.01 adaptive el 6.0-10.0 A-m/cell")
    for k in (60, 70, 80, 90, 100):
        print "plotting emittance growth ", k
        print onepct_emit_err01_adaptive_el[k]
        scc = k/10.0
        if not hasattr(onepct_emit_err01_adaptive_el[k], "xemitrms"):
            print "no data for ", scc
            continue
        emit0 = onepct_emit_err01_adaptive_el[k].xemitrms[0]
        plt.plot(onepct_emit_err01_adaptive_el[k].turns,
                 onepct_emit_err01_adaptive_el[k].xemitrms/emit0, label="err 0.01 adaptive el %3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.01 adaptive el 0-10.0 A-m/cell")
    for k in (0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100):
        scc = k/10.0
        emit0 = onepct_emit_err01_adaptive_el[k].xemit999[0]
        plt.plot(onepct_emit_err01_adaptive_el[k].turns,
                 onepct_emit_err01_adaptive_el[k].xemit999/emit0, label="err 0.01 adaptive el %3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################
################################

if 4 in plotgroup:
    plt.figure()
    plt.title("x RMS emittance growth err 0.01 elgauss adaptive 0-10.0 A-m/cell")
    for k in (0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100):
        print "plotting emittance growth ", k
        print onepct_emit_err01_adaptive_elgauss[k]
        scc = k/10.0
        if not hasattr(onepct_emit_err01_adaptive_elgauss[k], "xemitrms"):
            print "no data for ", scc
            continue
        emit0 = onepct_emit_err01_adaptive_elgauss[k].xemitrms[0]
        plt.plot(onepct_emit_err01_adaptive_elgauss[k].turns,
                 onepct_emit_err01_adaptive_elgauss[k].xemitrms/emit0, label="err 0.01 adaptive el %3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.01 adaptive elgauss 0-6.0 A-m/cell")
    for k in (0, 10, 20, 30, 40, 50, 60):
        print "plotting emittance growth ", k
        print onepct_emit_err01_adaptive_elgauss[k]
        scc = k/10.0
        if not hasattr(onepct_emit_err01_adaptive_elgauss[k], "xemitrms"):
            print "no data for ", scc
            continue
        emit0 = onepct_emit_err01_adaptive_elgauss[k].xemitrms[0]
        plt.plot(onepct_emit_err01_adaptive_elgauss[k].turns,
                 onepct_emit_err01_adaptive_elgauss[k].xemitrms/emit0, label="err 0.01 adaptive elgauss %3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x RMS emittance growth err 0.01 adaptive elgauss 6.0-10.0 A-m/cell")
    for k in (60, 70, 80, 90, 100):
        print "plotting emittance growth ", k
        print onepct_emit_err01_adaptive_elgauss[k]
        scc = k/10.0
        if not hasattr(onepct_emit_err01_adaptive_elgauss[k], "xemitrms"):
            print "no data for ", scc
            continue
        emit0 = onepct_emit_err01_adaptive_elgauss[k].xemitrms[0]
        plt.plot(onepct_emit_err01_adaptive_elgauss[k].turns,
                 onepct_emit_err01_adaptive_elgauss[k].xemitrms/emit0, label="err 0.01 adaptive elgauss %3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

    plt.figure()
    plt.title("x 99.9% emittance growth err 0.01 adaptive elgauss 0-10.0 A-m/cell")
    for k in (0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100):
        scc = k/10.0
        emit0 = onepct_emit_err01_adaptive_elgauss[k].xemit999[0]
        plt.plot(onepct_emit_err01_adaptive_elgauss[k].turns,
                 onepct_emit_err01_adaptive_elgauss[k].xemit999/emit0, label="err 0.01 adaptive el %3.1f"%scc, lw=2)
    plt.xlabel("turns")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    saveplt()

################################

plt.show()
