#!/usr/bin/env python
import sys, os
import matplotlib.pyplot as plt
import tables

from emit_analyze import *

######################################################################

# return xemitrms from diagnostics at every 25 turns

def xemitrms(d):
    h5 = tables.open_file(d + "/full2_0.h5")
    xemitrms = h5.root.emitx[::25]
    h5.close()
    return xemitrms

######################################################################

# returns a0, a1, a2 which are arrays that fit emittance
#   emit = a0 + a1 * (1/N) + a2 * (1/N)**2

def fit_noise(dir_1M, dir_2M, dir_4M):
    h5_1M = tables.open_file(dir_1M + "/full2_0.h5")
    emit_1M = h5_1M.root.emitx[::25]
    h5_1M.close()
    h5_2M = tables.open_file(dir_2M + "/full2_0.h5")
    emit_2M = h5_2M.root.emitx[::25]
    h5_2M.close()
    h5_4M = tables.open_file(dir_4M + "/full2_0.h5")
    emit_4M = h5_4M.root.emitx[::25]
    h5_4M.close()

    n = len(emit_1M)
    x = np.array([1/1.0e6, 1/2.0e6, 1/4.0e6])

    a0 = np.zeros(n, dtype='d')
    a1 = np.zeros(n, dtype='d')
    a2 = np.zeros(n, dtype='d')

    for i in range(n):
        # solve A * x = y
        # where x is is the coefficients a0, a1, a2, y is the emittance for
        # 1M, 2M, 4M particles and A is the array with rows:
        #  1 1/N 1/N**2

        A = np.zeros((3,3), dtype='d')
        for j, overn in enumerate(x):
            A[j, 0] = 1.0
            A[j, 1] = overn
            A[j, 2] = overn**2

        y = np.array([emit_1M[i], emit_2M[i], emit_4M[i]])
        ass = np.linalg.solve(A, y)
        a0[i] = ass[0]
        a1[i] = ass[1]
        a2[i] = ass[2]

    return a0, a1, a2

####################################################################

basedir = os.environ["SCRATCH"] + "/elens/"

dir_1M_err0_mc0_scc0 = "offdiag_magiccomp0.00"
dir_1M_err01_mc0_scc0 = "offdiag_1M_err01_magiccomp0.00"
dir_1M_err01_foc_mc0_scc0 = "offdiag_1M_err01_foc_magiccomp0.00"
dir_1M_err01_mc16_scc0p4 = "offdiag_1M_err01_magiccomp16_scc0.4.00"
dir_1M_err01_mc16_scc0p8 = "offdiag_1M_err01_magiccomp16_scc0.8.00"
dir_1M_err01_mc16_scc1p0 = "offdiag_1M_err01_magiccomp16_scc1.0.00"
dir_1M_err01_mc16_scc1p2 = "offdiag_1M_err01_magiccomp16_scc1.2.00"
dir_1M_err01_mc16_scc1p4 = "offdiag_1M_err01_magiccomp16_scc1.4.00"
dir_1M_err01_mc16_scc1p8 = "offdiag_1M_err01_magiccomp16_scc1.8.00"
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

dir_2M_err0_mc0_scc0 = "offdiag_2M_magiccomp0.00"
dir_2M_err01_mc0_scc0 = "offdiag_2M_err01_magiccomp0.00"
dir_2M_err01_mc16_scc0p4 = "offdiag_2M_err01_magiccomp16_scc0.4.00"
dir_2M_err01_mc16_scc0p8 = "offdiag_2M_err01_magiccomp16_scc0.8.00"
dir_2M_err01_mc16_scc1p0 = "offdiag_2M_err01_magiccomp16_scc1.0.00"
dir_2M_err01_mc16_scc1p2 = "offdiag_2M_err01_magiccomp16_scc1.2.00"
dir_2M_err01_mc16_scc1p4 = "offdiag_2M_err01_magiccomp16_scc1.4.00"
dir_2M_err01_mc16_scc1p8 = "offdiag_2M_err01_magiccomp16_scc1.8.00"
dir_2M_err01_mc16_scc2p0 = "offdiag_2M_err01_magiccomp16_scc2.0.00"

dir_4M_err0_mc0_scc0 = "offdiag_4M_magiccomp0.00"
dir_4M_err01_mc0_scc0 = "offdiag_4M_err01_magiccomp0.00"
dir_4M_err01_mc16_scc0p4 = "offdiag_4M_err01_magiccomp16_scc0.4.00"
dir_4M_err01_mc16_scc0p8 = "offdiag_4M_err01_magiccomp16_scc0.8.00"
dir_4M_err01_mc16_scc1p0 = "offdiag_4M_err01_magiccomp16_scc1.0.00"
dir_4M_err01_mc16_scc1p2 = "offdiag_4M_err01_magiccomp16_scc1.2.00"
dir_4M_err01_mc16_scc1p4 = "offdiag_4M_err01_magiccomp16_scc1.4.00"
dir_4M_err01_mc16_scc1p8 = "offdiag_4M_err01_magiccomp16_scc1.8.00"
dir_4M_err01_mc16_scc2p0 = "offdiag_4M_err01_magiccomp16_scc2.0.00"

xemit_mc0_d = xemitrms(basedir + dir_1M_err01_mc0_scc0)
xemit_mc0_f = xemitrms(basedir + dir_1M_err01_foc_mc0_scc0)
t = np.arange(0.0, 1001, 25.0)
plt.figure()
plt.title("err01 f quad vs. d quad scc0")
plt.plot(t, xemit_mc0_d/xemit_mc0_d[0], label="0.01 D quad error")
plt.plot(t, xemit_mc0_f/xemit_mc0_f[0], label="0.01 F quad error")
plt.xlabel("turn")
plt.ylabel("emittance growth")
plt.legend(loc="best")

emit_err0_mc0_nonoise = fit_noise(basedir + dir_1M_err0_mc0_scc0,
                                  basedir + dir_2M_err0_mc0_scc0,
                                  basedir + dir_4M_err0_mc0_scc0)




emit_err0_mc0_scc0 = np.load(basedir + dir_1M_err0_mc0_scc0 + "/emittances.npy")[()]

plt.figure()
plt.title("1M, 2M, 4M, inf no error scc0")
plt.plot(np.arange(0, 1001., 25), emit_err0_mc0_nonoise[0]/emit_err0_mc0_nonoise[0][0], label='no noise emit RMS')
plt.plot(emit_err0_mc0_scc0.turns, emit_err0_mc0_scc0.xemitrms/emit_err0_mc0_scc0.xemitrms[0], label="1M emit RMS")
emit2m = xemitrms(basedir + dir_2M_err0_mc0_scc0)
plt.plot(emit_err0_mc0_scc0.turns, emit2m/emit2m[0], label="2M emit RMS")
emit4m = xemitrms(basedir + dir_4M_err0_mc0_scc0)
plt.plot(emit_err0_mc0_scc0.turns, emit4m/emit4m[0], label="4M emit RMS")
plt.xlabel("turns")
plt.ylabel("emittance growth")
plt.legend(loc='best')

emit_err01_mc0_scc0 = np.load(basedir + dir_1M_err01_mc0_scc0 + "/emittances.npy")[()]
emit_err01_mc16_scc0p4 = np.load(basedir + dir_1M_err01_mc16_scc0p4 + "/emittances.npy")[()]

emit_err01_mc16_nonoise = fit_noise(basedir + dir_1M_err01_mc16_scc0p4,
                                    basedir + dir_2M_err01_mc16_scc0p4,
                                    basedir + dir_4M_err01_mc16_scc0p4)

emit1m = xemitrms(basedir + dir_1M_err01_mc16_scc0p4)
emit2m = xemitrms(basedir + dir_2M_err01_mc16_scc0p4)
emit4m = xemitrms(basedir + dir_4M_err01_mc16_scc0p4)

plt.figure()
plt.title("1M, 2M, 4M, inf err01 scc0.4")
#plt.plot(np.arange(0, 1001., 25), emit_err01_mc16_nonoise[0]/emit_err01_mc16_nonoise[0][0], label='no noise emit RMS')
plt.plot(emit_err01_mc16_scc0p4.turns, emit1m/emit1m[0], label="1M emit RMS")
plt.plot(emit_err01_mc16_scc0p4.turns, emit2m/emit2m[0], label="2M emit RMS")
plt.plot(emit_err01_mc16_scc0p4.turns, emit4m/emit4m[0], label="4M emit RMS")
plt.xlabel("turns")
plt.ylabel("emittance growth")
plt.legend(loc='best')

emit_err01_mc16_scc0p8 = np.load(basedir + dir_1M_err01_mc16_scc0p8 + "/emittances.npy")[()]
emit_err01_mc16_scc1p0 = np.load(basedir + dir_1M_err01_mc16_scc1p0 + "/emittances.npy")[()]

emit_err01_mc16_nonoise = fit_noise(basedir + dir_1M_err01_mc16_scc1p0,
                                  basedir + dir_2M_err01_mc16_scc1p0,
                                  basedir + dir_4M_err01_mc16_scc1p0)

emit1m = xemitrms(basedir + dir_1M_err01_mc16_scc1p0)
emit2m = xemitrms(basedir + dir_2M_err01_mc16_scc1p0)
emit4m = xemitrms(basedir + dir_4M_err01_mc16_scc1p0)

plt.figure()
plt.title("1M, 2M, 4M, inf err01 scc1.0")
#plt.plot(np.arange(0, 1001., 25), emit_err01_mc16_nonoise[0]/emit_err01_mc16_nonoise[0][0], label='no noise emit RMS')
plt.plot(emit_err01_mc16_scc1p0.turns, emit1m/emit1m[0], label="1M emit RMS")
plt.plot(emit_err01_mc16_scc1p0.turns, emit2m/emit2m[0], label="2M emit RMS")
plt.plot(emit_err0_mc0_scc0.turns, emit4m/emit4m[0], label="4M emit RMS")
plt.xlabel("turns")
plt.ylabel("emittance growth")
plt.legend(loc='best')

emit_err01_mc16_scc1p2 = np.load(basedir + dir_1M_err01_mc16_scc1p2 + "/emittances.npy")[()]
emit_err01_mc16_scc1p4 = np.load(basedir + dir_1M_err01_mc16_scc1p4 + "/emittances.npy")[()]
emit_err01_mc16_scc1p8 = np.load(basedir + dir_1M_err01_mc16_scc1p8 + "/emittances.npy")[()]
emit_err01_mc16_scc2p0 = np.load(basedir + dir_1M_err01_mc16_scc2p0 + "/emittances.npy")[()]

emit_err01_mc16_scc2p4 = np.load(basedir + dir_1M_err01_mc16_scc2p4 + "/emittances.npy")[()]
emit_err01_mc16_scc2p8 = np.load(basedir + dir_1M_err01_mc16_scc2p8 + "/emittances.npy")[()]
emit_err01_mc16_scc3p2 = np.load(basedir + dir_1M_err01_mc16_scc3p2 + "/emittances.npy")[()]
emit_err01_mc16_scc3p6 = np.load(basedir + dir_1M_err01_mc16_scc3p6 + "/emittances.npy")[()]

emit_err01_mc16_scc4p0 = np.load(basedir + dir_1M_err01_mc16_scc4p0 + "/emittances.npy")[()]
emit_err01_mc16_scc4p4 = np.load(basedir + dir_1M_err01_mc16_scc4p4 + "/emittances.npy")[()]
emit_err01_mc16_scc4p8 = np.load(basedir + dir_1M_err01_mc16_scc4p8 + "/emittances.npy")[()]

emit_err01_mc16_scc5p2 = np.load(basedir + dir_1M_err01_mc16_scc5p2 + "/emittances.npy")[()]
emit_err01_mc16_scc5p6 = np.load(basedir + dir_1M_err01_mc16_scc5p6 + "/emittances.npy")[()]

emit_err01_mc16_scc6p0 = np.load(basedir + dir_1M_err01_mc16_scc6p0 + "/emittances.npy")[()]

#####
#####  plot rms emittance scc0-scc2.0
plt.figure()
plt.title("x RMS emittance growth")

emit0 = emit_err0_mc0_scc0.xemitrms[0]
emit999 = emit_err0_mc0_scc0.xemit999[0]

plt.plot(emit_err0_mc0_scc0.turns, emit_err0_mc0_scc0.xemitrms/emit0, label="err0 scc0")
plt.plot(emit_err01_mc0_scc0.turns, emit_err01_mc0_scc0.xemitrms/emit0, label="err01 scc0")
plt.plot(emit_err01_mc16_scc0p4.turns, emit_err01_mc16_scc0p4.xemitrms/emit0, label="err01 scc0.4")
plt.plot(emit_err01_mc16_scc0p8.turns, emit_err01_mc16_scc0p8.xemitrms/emit0, label="err01 scc0.8")
plt.plot(emit_err01_mc16_scc1p0.turns, emit_err01_mc16_scc1p0.xemitrms/emit0, label="err01 scc1.0")
plt.plot(emit_err01_mc16_scc1p2.turns, emit_err01_mc16_scc1p2.xemitrms/emit0, label="err01 scc1.2")
plt.plot(emit_err01_mc16_scc1p4.turns, emit_err01_mc16_scc1p4.xemitrms/emit0, label="err01 scc1.4")
plt.plot(emit_err01_mc16_scc1p8.turns, emit_err01_mc16_scc1p8.xemitrms/emit0, label="err01 scc1.8")
plt.plot(emit_err01_mc16_scc2p0.turns, emit_err01_mc16_scc2p0.xemitrms/emit0, label="err01 scc2.0")
# plt.plot(emit_err01_mc16_scc2p4.turns, emit_err01_mc16_scc2p4.xemitrms/emit0, label="err01 scc2.4")
# plt.plot(emit_err01_mc16_scc2p8.turns, emit_err01_mc16_scc2p8.xemitrms/emit0, label="err01 scc2.8")
# plt.plot(emit_err01_mc16_scc3p2.turns, emit_err01_mc16_scc3p2.xemitrms/emit0, label="err01 scc3.2")
# plt.plot(emit_err01_mc16_scc3p6.turns, emit_err01_mc16_scc3p6.xemitrms/emit0, label="err01 scc3.6")
# plt.plot(emit_err01_mc16_scc4p0.turns, emit_err01_mc16_scc4p0.xemitrms/emit0, label="err01 scc4.0")
# plt.plot(emit_err01_mc16_scc4p4.turns, emit_err01_mc16_scc4p4.xemitrms/emit0, label="err01 scc4.4")
# plt.plot(emit_err01_mc16_scc4p8.turns, emit_err01_mc16_scc4p8.xemitrms/emit0, label="err01 scc4.8")
# plt.plot(emit_err01_mc16_scc5p2.turns, emit_err01_mc16_scc5p2.xemitrms/emit0, label="err01 scc5.2")
# plt.plot(emit_err01_mc16_scc5p6.turns, emit_err01_mc16_scc5p6.xemitrms/emit0, label="err01 scc5.6")
# plt.plot(emit_err01_mc16_scc6p0.turns, emit_err01_mc16_scc6p0.xemitrms/emit0, label="err01 scc6.0")

plt.xlabel("turns")
plt.ylabel("x RMS emit growth")
plt.legend(loc='best')

plt.savefig("x_rms_emittance_scc0-2.png")

####
#### plot rms emittance scc0, 1, 2-4

plt.figure()
plt.title("x RMS emittance growth")
plt.plot(emit_err0_mc0_scc0.turns, emit_err0_mc0_scc0.xemitrms/emit0, label="err0 scc0")
plt.plot(emit_err01_mc0_scc0.turns, emit_err01_mc0_scc0.xemitrms/emit0, label="err01 scc0")
#plt.plot(emit_err01_mc16_scc0p4.turns, emit_err01_mc16_scc0p4.xemitrms/emit0, label="err01 scc0.4")
#plt.plot(emit_err01_mc16_scc0p8.turns, emit_err01_mc16_scc0p8.xemitrms/emit0, label="err01 scc0.8")
plt.plot(emit_err01_mc16_scc1p0.turns, emit_err01_mc16_scc1p0.xemitrms/emit0, label="err01 scc1.0")
#plt.plot(emit_err01_mc16_scc1p2.turns, emit_err01_mc16_scc1p2.xemitrms/emit0, label="err01 scc1.2")
#plt.plot(emit_err01_mc16_scc1p4.turns, emit_err01_mc16_scc1p4.xemitrms/emit0, label="err01 scc1.4")
#plt.plot(emit_err01_mc16_scc1p8.turns, emit_err01_mc16_scc1p8.xemitrms/emit0, label="err01 scc1.8")
plt.plot(emit_err01_mc16_scc2p0.turns, emit_err01_mc16_scc2p0.xemitrms/emit0, label="err01 scc2.0")
plt.plot(emit_err01_mc16_scc2p4.turns, emit_err01_mc16_scc2p4.xemitrms/emit0, label="err01 scc2.4")
plt.plot(emit_err01_mc16_scc2p8.turns, emit_err01_mc16_scc2p8.xemitrms/emit0, label="err01 scc2.8")
plt.plot(emit_err01_mc16_scc3p2.turns, emit_err01_mc16_scc3p2.xemitrms/emit0, label="err01 scc3.2")
plt.plot(emit_err01_mc16_scc3p6.turns, emit_err01_mc16_scc3p6.xemitrms/emit0, label="err01 scc3.6")
plt.plot(emit_err01_mc16_scc4p0.turns, emit_err01_mc16_scc4p0.xemitrms/emit0, label="err01 scc4.0")
# plt.plot(emit_err01_mc16_scc4p4.turns, emit_err01_mc16_scc4p4.xemitrms/emit0, label="err01 scc4.4")
# plt.plot(emit_err01_mc16_scc4p8.turns, emit_err01_mc16_scc4p8.xemitrms/emit0, label="err01 scc4.8")
# plt.plot(emit_err01_mc16_scc5p2.turns, emit_err01_mc16_scc5p2.xemitrms/emit0, label="err01 scc5.2")
# plt.plot(emit_err01_mc16_scc5p6.turns, emit_err01_mc16_scc5p6.xemitrms/emit0, label="err01 scc5.6")
# plt.plot(emit_err01_mc16_scc6p0.turns, emit_err01_mc16_scc6p0.xemitrms/emit0, label="err01 scc6.0")

plt.xlabel("turns")
plt.ylabel("x RMS emit growth")
plt.legend(loc='best')

####
#### plot rms emittance scc0, 1, 2, 3.6, 4.0

plt.figure()
plt.title("x RMS emittance growth")
plt.plot(emit_err0_mc0_scc0.turns, emit_err0_mc0_scc0.xemitrms/emit0, label="err0 scc0")
plt.plot(emit_err01_mc0_scc0.turns, emit_err01_mc0_scc0.xemitrms/emit0, label="err01 scc0")
plt.plot(emit_err01_mc16_scc2p0.turns, emit_err01_mc16_scc2p0.xemitrms/emit0, label="err01 scc2.0")

plt.plot(emit_err01_mc16_scc3p6.turns, emit_err01_mc16_scc3p6.xemitrms/emit0, label="err01 scc3.6")
plt.plot(emit_err01_mc16_scc4p0.turns, emit_err01_mc16_scc4p0.xemitrms/emit0, label="err01 scc4.0")
plt.plot(emit_err01_mc16_scc4p4.turns, emit_err01_mc16_scc4p4.xemitrms/emit0, label="err01 scc4.4")
plt.plot(emit_err01_mc16_scc4p8.turns, emit_err01_mc16_scc4p8.xemitrms/emit0, label="err01 scc4.8")
plt.plot(emit_err01_mc16_scc5p2.turns, emit_err01_mc16_scc5p2.xemitrms/emit0, label="err01 scc5.2")
plt.plot(emit_err01_mc16_scc5p6.turns, emit_err01_mc16_scc5p6.xemitrms/emit0, label="err01 scc5.6")
plt.plot(emit_err01_mc16_scc6p0.turns, emit_err01_mc16_scc6p0.xemitrms/emit0, label="err01 scc6.0")

plt.xlabel("turns")
plt.ylabel("x RMS emit growth")
plt.legend(loc='best')

plt.savefig("x_rms_emittance_scc2-4.png")

####
#### plot rms emittance scc0, 1, 2, 3.2, 3.6

plt.figure()
plt.title("x RMS emittance growth")
plt.plot(emit_err0_mc0_scc0.turns, emit_err0_mc0_scc0.xemitrms/emit0, label="err0 scc0")
plt.plot(emit_err01_mc0_scc0.turns, emit_err01_mc0_scc0.xemitrms/emit0, label="err01 scc0")
plt.plot(emit_err01_mc16_scc2p0.turns, emit_err01_mc16_scc2p0.xemitrms/emit0, label="err01 scc2.0")

plt.plot(emit_err01_mc16_scc3p2.turns, emit_err01_mc16_scc3p2.xemitrms/emit0, label="err01 scc3.2")
plt.plot(emit_err01_mc16_scc3p6.turns, emit_err01_mc16_scc3p6.xemitrms/emit0, label="err01 scc3.6")

plt.xlabel("turns")
plt.ylabel("x RMS emit growth")
plt.legend(loc='best')

####
#### plot rms emittance scc0, 4-6

plt.figure()
plt.title("x RMS emittance growth")
plt.plot(emit_err0_mc0_scc0.turns, emit_err0_mc0_scc0.xemitrms/emit0, label="err0 scc0")
plt.plot(emit_err01_mc0_scc0.turns, emit_err01_mc0_scc0.xemitrms/emit0, label="err01 scc0")
#plt.plot(emit_err01_mc16_scc0p4.turns, emit_err01_mc16_scc0p4.xemitrms/emit0, label="err01 scc0.4")
#plt.plot(emit_err01_mc16_scc0p8.turns, emit_err01_mc16_scc0p8.xemitrms/emit0, label="err01 scc0.8")
#plt.plot(emit_err01_mc16_scc1p0.turns, emit_err01_mc16_scc1p0.xemitrms/emit0, label="err01 scc1.0")
#plt.plot(emit_err01_mc16_scc1p2.turns, emit_err01_mc16_scc1p2.xemitrms/emit0, label="err01 scc1.2")
#plt.plot(emit_err01_mc16_scc1p4.turns, emit_err01_mc16_scc1p4.xemitrms/emit0, label="err01 scc1.4")
#plt.plot(emit_err01_mc16_scc1p8.turns, emit_err01_mc16_scc1p8.xemitrms/emit0, label="err01 scc1.8")
plt.plot(emit_err01_mc16_scc2p0.turns, emit_err01_mc16_scc2p0.xemitrms/emit0, label="err01 scc2.0")
plt.plot(emit_err01_mc16_scc2p4.turns, emit_err01_mc16_scc2p4.xemitrms/emit0, label="err01 scc2.4")
plt.plot(emit_err01_mc16_scc2p8.turns, emit_err01_mc16_scc2p8.xemitrms/emit0, label="err01 scc2.8")
plt.plot(emit_err01_mc16_scc3p2.turns, emit_err01_mc16_scc3p2.xemitrms/emit0, label="err01 scc3.2")
plt.plot(emit_err01_mc16_scc3p6.turns, emit_err01_mc16_scc3p6.xemitrms/emit0, label="err01 scc3.6")
plt.plot(emit_err01_mc16_scc4p0.turns, emit_err01_mc16_scc4p0.xemitrms/emit0, label="err01 scc4.0")
plt.plot(emit_err01_mc16_scc4p4.turns, emit_err01_mc16_scc4p4.xemitrms/emit0, label="err01 scc4.4")
plt.plot(emit_err01_mc16_scc4p8.turns, emit_err01_mc16_scc4p8.xemitrms/emit0, label="err01 scc4.8")
plt.plot(emit_err01_mc16_scc5p2.turns, emit_err01_mc16_scc5p2.xemitrms/emit0, label="err01 scc5.2")
plt.plot(emit_err01_mc16_scc5p6.turns, emit_err01_mc16_scc5p6.xemitrms/emit0, label="err01 scc5.6")
plt.plot(emit_err01_mc16_scc6p0.turns, emit_err01_mc16_scc6p0.xemitrms/emit0, label="err01 scc6.0")

plt.xlabel("turns")
plt.ylabel("x RMS emit growth")
plt.legend(loc='best')

plt.savefig("x_rms_emittance_scc_2-6.png")

####
#### plot 99.9% emittance 0-2

plt.figure()
plt.title("99.9% emittance growth")
plt.plot(emit_err0_mc0_scc0.turns, emit_err0_mc0_scc0.xemit999/emit999, label="err0 scc0")
plt.plot(emit_err01_mc0_scc0.turns, emit_err01_mc0_scc0.xemit999/emit999, label="err01 scc0")
plt.plot(emit_err01_mc16_scc0p4.turns, emit_err01_mc16_scc0p4.xemit999/emit999, label="err01 scc0.4")
plt.plot(emit_err01_mc16_scc0p8.turns, emit_err01_mc16_scc0p8.xemit999/emit999, label="err01 scc0.8")
plt.plot(emit_err01_mc16_scc1p0.turns, emit_err01_mc16_scc1p0.xemit999/emit999, label="err01 scc1.0")
plt.plot(emit_err01_mc16_scc1p2.turns, emit_err01_mc16_scc1p2.xemit999/emit999, label="err01 scc1.2")
plt.plot(emit_err01_mc16_scc1p4.turns, emit_err01_mc16_scc1p4.xemit999/emit999, label="err01 scc1.4")
plt.plot(emit_err01_mc16_scc1p8.turns, emit_err01_mc16_scc1p8.xemit999/emit999, label="err01 scc1.8")
plt.plot(emit_err01_mc16_scc2p0.turns, emit_err01_mc16_scc2p0.xemit999/emit999, label="err01 scc2.0")


plt.xlabel("turns")
plt.ylabel("x 99.9% emit growth")
plt.legend(loc='best')

plt.savefig("x_999_emittance_scc0-2.png")

####
#### plot 99.9% emittance 0, 2-4

plt.figure()
plt.title("99.9% emittance growth")
plt.plot(emit_err0_mc0_scc0.turns, emit_err0_mc0_scc0.xemit999/emit999, label="err0 scc0")
plt.plot(emit_err01_mc0_scc0.turns, emit_err01_mc0_scc0.xemit999/emit999, label="err01 scc0")
plt.plot(emit_err01_mc16_scc2p0.turns, emit_err01_mc16_scc2p0.xemit999/emit999, label="err01 scc2.0")
plt.plot(emit_err01_mc16_scc2p4.turns, emit_err01_mc16_scc2p4.xemit999/emit999, label="err01 scc2.4")
plt.plot(emit_err01_mc16_scc2p8.turns, emit_err01_mc16_scc2p8.xemit999/emit999, label="err01 scc2.8")
plt.plot(emit_err01_mc16_scc3p2.turns, emit_err01_mc16_scc3p2.xemit999/emit999, label="err01 scc3.2")
plt.plot(emit_err01_mc16_scc3p6.turns, emit_err01_mc16_scc3p6.xemit999/emit999, label="err01 scc3.6")
plt.plot(emit_err01_mc16_scc4p0.turns, emit_err01_mc16_scc4p0.xemit999/emit999, label="err01 scc4.0")


plt.xlabel("turns")
plt.ylabel("x 99.9% emit growth")
plt.legend(loc='best')

plt.savefig("x_999_emittance_scc2-4.png")

####
####
#### plot 99.9% emittance 0, 2-6

plt.figure()
plt.title("99.9% emittance growth")
plt.plot(emit_err0_mc0_scc0.turns, emit_err0_mc0_scc0.xemit999/emit999, label="err0 scc0")
plt.plot(emit_err01_mc0_scc0.turns, emit_err01_mc0_scc0.xemit999/emit999, label="err01 scc0")
plt.plot(emit_err01_mc16_scc2p0.turns, emit_err01_mc16_scc2p0.xemit999/emit999, label="err01 scc2.0")
plt.plot(emit_err01_mc16_scc2p4.turns, emit_err01_mc16_scc2p4.xemit999/emit999, label="err01 scc2.4")
plt.plot(emit_err01_mc16_scc2p8.turns, emit_err01_mc16_scc2p8.xemit999/emit999, label="err01 scc2.8")
plt.plot(emit_err01_mc16_scc3p2.turns, emit_err01_mc16_scc3p2.xemit999/emit999, label="err01 scc3.2")
plt.plot(emit_err01_mc16_scc3p6.turns, emit_err01_mc16_scc3p6.xemit999/emit999, label="err01 scc3.6")
plt.plot(emit_err01_mc16_scc4p0.turns, emit_err01_mc16_scc4p0.xemit999/emit999, label="err01 scc4.0")
plt.plot(emit_err01_mc16_scc4p4.turns, emit_err01_mc16_scc4p4.xemit999/emit999, label="err01 scc4.4")
plt.plot(emit_err01_mc16_scc4p8.turns, emit_err01_mc16_scc4p8.xemit999/emit999, label="err01 scc4.8")
plt.plot(emit_err01_mc16_scc5p2.turns, emit_err01_mc16_scc5p2.xemit999/emit999, label="err01 scc5.2")
plt.plot(emit_err01_mc16_scc5p6.turns, emit_err01_mc16_scc5p6.xemit999/emit999, label="err01 scc5.6")
plt.plot(emit_err01_mc16_scc6p0.turns, emit_err01_mc16_scc6p0.xemit999/emit999, label="err01 scc6.0")


plt.xlabel("turns")
plt.ylabel("x 99.9% emit growth")
plt.legend(loc='best')

plt.savefig("x_999_emittance_scc2-6.png")

####
#### plot kurtosis

plt.figure()
plt.plot(emit_err0_mc0_scc0.turns, emit_err0_mc0_scc0.xkurt, label="x kurtosis err0 scc0")
plt.plot(emit_err01_mc0_scc0.turns, emit_err01_mc0_scc0.xkurt, label="x kurtosis err01 scc0")
plt.plot(emit_err01_mc16_scc0p4.turns, emit_err01_mc16_scc0p4.xkurt, label="x kurtosis scc0.4")
plt.plot(emit_err01_mc16_scc0p8.turns, emit_err01_mc16_scc0p8.xkurt, label="x kurtosis err01 scc0.8")
plt.plot(emit_err01_mc16_scc1p0.turns, emit_err01_mc16_scc1p0.xkurt, label="x kurtosis err01 scc1.0")
plt.plot(emit_err01_mc16_scc1p2.turns, emit_err01_mc16_scc1p2.xkurt, label="x kurtosis err01 scc1.2")
plt.plot(emit_err01_mc16_scc1p4.turns, emit_err01_mc16_scc1p4.xkurt, label="x kurtosis err01 scc1.4")
plt.plot(emit_err01_mc16_scc1p8.turns, emit_err01_mc16_scc1p8.xkurt, label="x kurtosis err01 scc1.8")
plt.plot(emit_err01_mc16_scc2p0.turns, emit_err01_mc16_scc2p0.xkurt, label="x kurtosis err01 scc2.0")

plt.xlabel("turns")
plt.ylabel("x kurtosis")
plt.legend(loc='best')

####
#### plot x rms

plt.figure()
plt.plot(emit_err0_mc0_scc0.turns, emit_err0_mc0_scc0.xrms, label="x RMS err0 scc0")
plt.plot(emit_err01_mc0_scc0.turns, emit_err01_mc0_scc0.xrms, label="x RMS err01 scc0")
plt.plot(emit_err01_mc16_scc0p4.turns, emit_err01_mc16_scc0p4.xrms, label="x RMS scc0.4")
plt.plot(emit_err01_mc16_scc0p8.turns, emit_err01_mc16_scc0p8.xrms, label="x RMS err01 scc0.8")
plt.plot(emit_err01_mc16_scc1p0.turns, emit_err01_mc16_scc1p0.xrms, label="x RMS err01 scc1.0")
plt.plot(emit_err01_mc16_scc1p2.turns, emit_err01_mc16_scc1p2.xrms, label="x RMS err01 scc1.2")
plt.plot(emit_err01_mc16_scc1p4.turns, emit_err01_mc16_scc1p4.xrms, label="x RMS err01 scc1.4")
plt.plot(emit_err01_mc16_scc1p8.turns, emit_err01_mc16_scc1p8.xrms, label="x RMS err01 scc1.8")
plt.plot(emit_err01_mc16_scc2p0.turns, emit_err01_mc16_scc2p0.xrms, label="x RMS err01 scc2.0")

plt.xlabel("turns")
plt.ylabel("x RMS")
plt.legend(loc='best')

plt.show()
