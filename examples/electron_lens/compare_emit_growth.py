#!/usr/bin/env python
import sys, os
import matplotlib.pyplot as plt

from emit_analyze import *

basedir = os.environ["SCRATCH"] + "/elens/"

dir_scc0 = "offdiag_magiccomp0.00"
dir_scc0p2 = "offdiag_magiccomp16_scc0.2.00"
dir_scc0p4 = "offdiag_magiccomp16_scc0.4.00"
dir_scc0p6 = "offdiag_magiccomp16_scc0.6.00"
dir_scc0p8 = "offdiag_magiccomp16_scc0.8.00"
dir_scc1p0 = "offdiag_magiccomp16.00"
dir_scc1p6 = "offdiag_magiccomp16_scc1.6.00"
dir_scc2p0 = "offdiag_magiccomp16_scc2.0.00"

emit_scc0 = dir_emittances(basedir + dir_scc0)
emit_scc0p2 = dir_emittances(basedir + dir_scc0p2)
emit_scc0p4 = dir_emittances(basedir + dir_scc0p4)
emit_scc0p6 = dir_emittances(basedir + dir_scc0p6)
emit_scc0p8 = dir_emittances(basedir + dir_scc0p8)
emit_scc1p0 = dir_emittances(basedir + dir_scc1p0)
emit_scc1p6 = dir_emittances(basedir + dir_scc1p6)
emit_scc2p0 = dir_emittances(basedir + dir_scc2p0)

h5scc0 = tables.open_file(basedir + dir_scc0 + "/full2_0.h5")
xrms_scc0 = h5scc0.root.std[0,:]
h5scc0.close()

plt.figure()
plt.subplot(211)
plt.plot(xrms_scc0, label="x rms no compensation")
plt.ylabel("x RMS")
plt.legend(loc='best')
plt.subplot(212)
plt.plot(emit_scc0.turns, emit_scc0.betas_x, label='beta x')
plt.xlabel("turn")
plt.ylabel("beta x")
plt.legend(loc='best')

plt.figure()
plt.title("x RMS and sqrt(beta*emit999) no compensation")
plt.plot(xrms_scc0, label="x rms", lw=2)
plt.plot(emit_scc0.turns, np.sqrt(emit_scc0.betas_x * emit_scc0.xemit999), label="sqrt(betax * xemit999)", lw=2)
plt.xlabel("turn")
plt.ylabel("RMS [m]")
plt.legend(loc='best')


plt.figure()
plt.title("RMS emittance growth")
xemitrms0 = emit_scc0.xemitrms[0]
plt.plot(emit_scc0.turns, emit_scc0.xemitrms/xemitrms0, label="no compensation",lw=2)
plt.plot(emit_scc0p2.turns, emit_scc0p2.xemitrms/xemitrms0, label="0.2 compensation",lw=2)
plt.plot(emit_scc0p4.turns, emit_scc0p4.xemitrms/xemitrms0, label="0.4 compensation",lw=2)
plt.plot(emit_scc0p6.turns, emit_scc0p6.xemitrms/xemitrms0, label="0.6 compensation",lw=2)
plt.plot(emit_scc0p8.turns, emit_scc0p8.xemitrms/xemitrms0, label="0.8 compensation",lw=2)
plt.plot(emit_scc1p0.turns, emit_scc1p0.xemitrms/xemitrms0, label="1.0 compensation",lw=2)
plt.plot(emit_scc1p6.turns, emit_scc1p6.xemitrms/xemitrms0, label="1.6 compensation",lw=2)
plt.plot(emit_scc2p0.turns, emit_scc2p0.xemitrms/xemitrms0, label="2.0 compensation",lw=2)
plt.legend(loc='best')
plt.xlabel("turn")
plt.ylabel("x RMS emittance growth")

plt.figure()
plt.title("99.9% emittance growth")

xemit9990 = emit_scc0.xemit999[0]
plt.plot(emit_scc0.turns, emit_scc0.xemit999/xemit9990, label="no compensation",lw=2)
plt.plot(emit_scc0p2.turns, emit_scc0p2.xemit999/xemit9990, label="0.2 compensation",lw=2)
plt.plot(emit_scc0p4.turns, emit_scc0p4.xemit999/xemit9990, label="0.4 compensation",lw=2)
plt.plot(emit_scc0p6.turns, emit_scc0p6.xemit999/xemit9990, label="0.6 compensation",lw=2)
plt.plot(emit_scc0p8.turns, emit_scc0p8.xemit999/xemit9990, label="0.8 compensation",lw=2)
plt.plot(emit_scc1p0.turns, emit_scc1p0.xemit999/xemit9990, label="1.0 compensation",lw=2)
plt.plot(emit_scc1p6.turns, emit_scc1p6.xemit999/xemit9990, label="1.6 compensation",lw=2)
plt.plot(emit_scc2p0.turns, emit_scc2p0.xemit999/xemit9990, label="2.0 compensation",lw=2)
plt.legend(loc='best')
plt.xlabel("turn")
plt.ylabel("x 99.9% emittance growth")

plt.show()
