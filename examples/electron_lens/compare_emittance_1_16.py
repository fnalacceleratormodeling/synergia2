#!/usr/bin/env python
import sys, os
import numpy as np
import matplotlib.pyplot as plt

from emit_analyze import *

basedir = os.environ["SCRATCH"] + "/elens/"

emit_1M = np.load(basedir + "offdiag_magiccomp0.00/emittances.npy")[()]
emit_2M = np.load(basedir + "offdiag_2M_magiccomp0.00/emittances.npy")[()]
emit_4M = np.load(basedir + "offdiag_4M_magiccomp0.00/emittances.npy")[()]
emit_8M = np.load(basedir + "offdiag_8M_magiccomp0.00/emittances.npy")[()]
emit_16M = np.load(basedir + "offdiag_16M_magiccomp0.00/emittances.npy")[()]

plt.plot(emit_1M.turns, emit_1M.xemitrms/emit_1M.xemitrms[0], label="x RMS emittance growth 1M particles")
plt.plot(emit_2M.turns, emit_2M.xemitrms/emit_2M.xemitrms[0], label="x RMS emittance growth 2M particles")
plt.plot(emit_4M.turns, emit_4M.xemitrms/emit_4M.xemitrms[0], label="x RMS emittance growth 4M particles")
plt.plot(emit_8M.turns, emit_8M.xemitrms/emit_8M.xemitrms[0], label="x RMS emittance growth 8M particles")
plt.plot(emit_16M.turns, emit_16M.xemitrms/emit_16M.xemitrms[0], label="x RMS emittance growth 16M particles")
plt.xlabel("turn")
plt.ylabel("emittance growth")
plt.legend(loc='best')

plt.figure()
plt.plot(emit_1M.turns, emit_1M.xemit999/emit_1M.xemit999[0], label="99.9% emittance growth 1M particles")
plt.plot(emit_2M.turns, emit_2M.xemit999/emit_2M.xemit999[0], label="99.9% emittance growth 2M particles")
plt.plot(emit_4M.turns, emit_4M.xemit999/emit_4M.xemit999[0], label="99.9% emittance growth 4M particles")
plt.plot(emit_8M.turns, emit_8M.xemit999/emit_8M.xemit999[0], label="99.9% emittance growth 8M particles")
plt.plot(emit_16M.turns, emit_16M.xemit999/emit_16M.xemit999[0], label="99.9% emittance growth 16M particles")
plt.xlabel("turn")
plt.ylabel("emittance growth")
plt.legend(loc='best')

emittances = np.vstack((emit_1M.turns, emit_1M.xemitrms, emit_2M.xemitrms, emit_4M.xemitrms, emit_8M.xemitrms, emit_16M.xemitrms)).transpose()
np.savetxt("emittances_by_N.txt", emittances)

plt.show()
