#!/usr/bin/env python
import sys, os

interactive = False
saveplots = True

import numpy as np
import matplotlib
if not interactive:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import tables

##########################
pltnum = 0
def saveplt():
    global keepopen
    global pltnum
    pltname = "plot_cut_particles_%04d"%pltnum
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
##########################

basedir = os.environ["SCRATCH"] + "/elens/"

h5_4sig_scc00 = "offdiag_16M_err01_magiccomp16_scc0.0_4sig.00"
h5_4sig_scc44 = "offdiag_16M_err01_magiccomp16_scc4.4_4sig.00"

h5_5sig_scc00 = "offdiag_16M_err01_magiccomp16_scc0.0_5sig.00"
h5_5sig_scc44 = "offdiag_16M_err01_magiccomp16_scc4.4_5sig.00"

el8_4sig_dir=basedir + "offdiag_16M_err01_el8.0_4sig.00"
el8_5sig_dir=basedir + "offdiag_16M_err01_el8.0_5sig.00"
elgauss6_4sig_dir=basedir + "offdiag_16M_err01_elgauss6.0_4sig.00"
elgauss6_5sig_dir=basedir + "offdiag_16M_err01_elgauss6.0_5sig.00"
eladapt5_4sig_dir=basedir + "offdiag_16M_err01_adaptive_el5.0_4sig.00"
eladapt5_5sig_dir=basedir + "offdiag_16M_err01_adaptive_el5.0_5sig.00"
elgaussadapt6_4sig_dir=basedir + "offdiag_16M_err01_adaptive_elgauss6.0_4sig.00"
elgaussadapt6_5sig_dir=basedir + "offdiag_16M_err01_adaptive_elgauss6.0_5sig.00"

elgauss6_0p32_0p44_dir = basedir + "offdiag_16M_err01_elgauss6.0_0p32_0p44_4sig.00"
elgauss6_0p42_0p54_dir = basedir + "offdiag_16M_err01_elgauss6.0_0p42_0p54_4sig.00"
elgauss6_0p52_0p64_dir = basedir + "offdiag_16M_err01_elgauss6.0_0p52_0p64_4sig.00"
elgauss6_0p62_0p74_dir = basedir + "offdiag_16M_err01_elgauss6.0_0p62_0p74_4sig.00"
elgauss6_0p72_0p84_dir = elgauss6_4sig_dir
elgauss6_0p82_0p94_dir = basedir + "offdiag_16M_err01_elgauss6.0_0p82_0p94_4sig.00"

err05_4sig_scc0_dir = basedir + "offdiag_16M_err05_magiccomp16_scc0.0_4sig.00"
err05_elgauss6_4sig_dir = basedir + "offdiag_16M_err05_elgauss6.0_4sig.00"

scc00_4sig_diag = tables.open_file(basedir + h5_4sig_scc00 + "/full2_0.h5")
scc44_4sig_diag = tables.open_file(basedir + h5_4sig_scc44 + "/full2_0.h5")

scc00_5sig_diag = tables.open_file(basedir + h5_5sig_scc00 + "/full2_0.h5")
scc44_5sig_diag = tables.open_file(basedir + h5_5sig_scc44 + "/full2_0.h5")

turns = np.arange(1001.0)
part_4sig_00 = scc00_4sig_diag.root.num_particles.read()
part_4sig_44 = scc44_4sig_diag.root.num_particles.read()

scc00_4sig_diag.close()
scc44_4sig_diag.close()

part_5sig_00 = scc00_5sig_diag.root.num_particles.read()
part_5sig_44 = scc44_5sig_diag.root.num_particles.read()

scc00_5sig_diag.close()
scc44_5sig_diag.close()



norig = part_4sig_00[0]

plt.figure()
plt.title("particles lost to 4 sigma aperture cut")
plt.subplot(211)
plt.plot(turns, norig-part_4sig_00, label="no compensation")
plt.legend(loc='best')
plt.subplot(212)
plt.plot(turns, norig-part_4sig_44, label="best compensation")
plt.xlabel("turn")
plt.ylabel("lost particles")
plt.legend(loc='best')
saveplt()

plt.figure()
plt.title("particles lost to 5 sigma aperture cut")
plt.subplot(211)
plt.plot(turns, norig-part_5sig_00, label="no compensation")
plt.legend(loc='best')
plt.subplot(212)
plt.plot(turns, norig-part_5sig_44, label="best compensation")
plt.xlabel("turn")
plt.ylabel("lost particles")
plt.legend(loc='best')
saveplt()

plt.figure()
plt.title("particles lost to 4 sigma aperture cut")
plt.plot(turns, norig-part_4sig_00, label="no compensation")
plt.plot(turns, norig-part_4sig_44, label="best compensation")
plt.xlabel("turn")
plt.ylabel("lost particles")
plt.legend(loc='best')
saveplt()

plt.figure()
plt.title("particles lost to 5 sigma aperture cut")
plt.plot(turns, norig-part_5sig_00, label="no compensation")
plt.plot(turns, norig-part_5sig_44, label="best compensation")
plt.xlabel("turn")
plt.ylabel("lost particles")
plt.legend(loc='best')
saveplt()

plt.figure()
plt.title("aperture cut particles")
plt.plot(turns, norig-part_4sig_00, label="4 sigma no compensation")
plt.plot(turns, norig-part_4sig_44, label="4 sigma best compensation")
plt.plot(turns, norig-part_5sig_00, label="5 sigma no compensation")
plt.plot(turns, norig-part_5sig_44, label="5 sigma best compensation")
plt.xlabel("turn")
plt.ylabel("lost particles")
plt.legend(loc='best')
saveplt()

print "losses 4 sigma no compensation: ", float(norig-part_4sig_00[-1])/norig
print "losses 4 sigma optimal compensation: ", float(norig-part_4sig_44[-1])/norig
print "losses 5 sigma no compensation: ", float(norig-part_5sig_00[-1])/norig
print "losses 5 sigma optimal compensation: ", float(norig-part_5sig_44[-1])/norig

print "part_4sig_00.shape: ", part_4sig_00.shape
print "part_4sig_44.shape: ", part_4sig_44.shape

nturns = part_4sig_00.shape[0]



#cut = np.vstack((np.arange(nturns), norig-part00, norig-part44)).transpose()
#np.savetxt("cut_particles.txt", cut, header="#turn no-compensation-cut-particles optimal-compensation-cut-particles")

# el8 4 and 5 sigma, no long or transverse adapt
h5 = tables.open_file(el8_4sig_dir + "/full2_0.h5")
el8_4sig_particles = h5.root.num_particles.read()
el8_4sig_xemit = h5.root.emitx.read()
el8_4sig_xstd = h5.root.std[0]
h5.close()
h5 = tables.open_file(el8_5sig_dir + "/full2_0.h5")
el8_5sig_particles = h5.root.num_particles.read()
el8_5sig_xemit = h5.root.emitx.read()
el8_5sig_xstd = h5.root.std[0]
h5.close()

# elgauss6.0 4 and 5 sigma, longitudinal shaping, fixed transverse
h5 = tables.open_file(elgauss6_4sig_dir + "/full2_0.h5")
elgauss6_4sig_particles = h5.root.num_particles.read()
elgauss6_4sig_xemit = h5.root.emitx.read()
h5.close()
h5 = tables.open_file(elgauss6_5sig_dir + "/full2_0.h5")
elgauss6_5sig_particles = h5.root.num_particles.read()
elgauss6_5sig_xemit = h5.root.emitx.read()
h5.close()

# el5_adaptive 4 and 5 sigma, adaptive transverse, DC longitudinal
h5 = tables.open_file(eladapt5_4sig_dir + "/full2_0.h5")
eladapt5_4sig_particles = h5.root.num_particles.read()
eladapt5_4sig_xemit = h5.root.emitx.read()
h5.close()
h5 = tables.open_file(eladapt5_5sig_dir + "/full2_0.h5")
eladapt5_5sig_particles = h5.root.num_particles.read()
eladapt5_5sig_xemit = h5.root.emitx.read()
h5.close()

# elgauss6 4 and 5 sigma adaptive, adaptive transverse, longitudinal pulsed
h5 = tables.open_file(elgaussadapt6_4sig_dir + "/full2_0.h5")
elgaussadapt6_4sig_particles = h5.root.num_particles.read()
elgaussadapt6_4sig_xemit = h5.root.emitx.read()
h5.close()
h5 = tables.open_file(elgaussadapt6_5sig_dir + "/full2_0.h5")
elgaussadapt6_5sig_particles = h5.root.num_particles.read()
elgaussadapt6_5sig_xemit = h5.root.emitx.read()
h5.close()

plt.figure()
plt.title("aperture cut particles el 8 A-m, fixed transverse, DC longitudinal")
plt.plot(norig-el8_4sig_particles, label="el 8.0 A-m 4 sigma")
plt.plot(norig-el8_5sig_particles, label="el 8.0 A-m 5 sigma")
plt.xlabel("turns")
plt.ylabel("cut particles")
plt.legend(loc='best')
saveplt()
plt.figure()
plt.title("x RMS emittance 4 and 5 sigma cut el 8.0 A-m")
plt.plot(el8_4sig_xemit, label="4 sigma")
plt.plot(el8_5sig_xemit, label="5 sigma")
plt.xlabel("turns")
plt.ylabel("x RMS emittance")
plt.legend(loc='best')
saveplt()
plt.figure()
plt.title("x std 4 and 5 sigma cut el 8.0 A-m")
plt.plot(el8_4sig_xstd, label="4 sigma")
plt.plot(el8_5sig_xstd, label="5 sigma")
plt.xlabel("turns")
plt.ylabel("x std")
plt.legend(loc='best')
saveplt()
print "el 8.0 A-m 4sigma cut losses: ", float(norig-el8_4sig_particles[-1])/norig
print "el 8.0 A-m 5sigma cut losses: ", float(norig-el8_5sig_particles[-1])/norig

plt.figure()
plt.title("aperture cut particles elgauss 6.0 A-m, fixed tranverse, pulsed el")
plt.plot(norig-elgauss6_4sig_particles, label="4 sigma")
plt.plot(norig-elgauss6_5sig_particles, label="5 sigma")
plt.xlabel("turns")
plt.ylabel("cut particles")
plt.legend(loc='best')
saveplt()
plt.figure()
plt.title("x RMS emittance 4 and 5 sigma cut 6 A-m pulse current, fixed transverse")
plt.plot(elgauss6_4sig_xemit, label="4 sigma")
plt.plot(elgauss6_5sig_xemit, label="5 sigma")
plt.xlabel("turns")
plt.ylabel("x RMS emittance")
plt.legend(loc='best')
saveplt()
print "elgauss 6.0 A-m 4sigma cut losses: ", float(norig-elgauss6_4sig_particles[-1])/norig
print "elgauss 6.0 A-m 5sigma cut losses: ", float(norig-elgauss6_5sig_particles[-1])/norig

plt.figure()
plt.title("aperture cut particles eladapt 5.0 A-m, adaptive tranverse, DC el")
plt.plot(norig-eladapt5_4sig_particles, label="4 sigma")
plt.plot(norig-eladapt5_5sig_particles, label="5 sigma")
plt.xlabel("turns")
plt.ylabel("cut particles")
plt.legend(loc='best')
saveplt()
plt.figure("x RMS emittance 4 and 5 sigma cut eladapt 5 A-m, adapt transverse, DC el")
plt.plot(eladapt5_4sig_xemit, label="4 sigma")
plt.plot(eladapt5_5sig_xemit, label="5 sigma")
plt.xlabel("turns")
plt.ylabel("x RMS emittance")
plt.legend(loc='best')
saveplt()
print "eladapt 5.0 A-m 4sigma cut losses: ", float(norig-eladapt5_4sig_particles[-1])/norig
print "eladapt 5.0 A-m 5sigma cut losses: ", float(norig-eladapt5_5sig_particles[-1])/norig

plt.figure()
plt.title("aperture cut particles elgaussadapt 6.0 A-m, adapt tranverse, pulsed el")
plt.plot(norig-elgaussadapt6_4sig_particles, label="4 sigma")
plt.plot(norig-elgaussadapt6_5sig_particles, label="5 sigma")
plt.xlabel("turns")
plt.ylabel("cut particles")
plt.legend(loc='best')
saveplt()
plt.figure()
plt.title("x RMS emittance 4 and 5 sigma cut 6 A-m pulse current, adapt transverse")
plt.plot(elgaussadapt6_4sig_xemit, label="4 sigma")
plt.plot(elgaussadapt6_5sig_xemit, label="5 sigma")
plt.xlabel("turns")
plt.ylabel("x RMS emittance")
plt.legend(loc='best')
saveplt()
print "elgaussadapt 6.0 A-m 4sigma cut losses: ", float(norig-elgaussadapt6_4sig_particles[-1])/norig
print "elgaussadapt 6.0 A-m 5sigma cut losses: ", float(norig-elgaussadapt6_5sig_particles[-1])/norig

############################

h5 = tables.open_file(elgauss6_0p32_0p44_dir + "/full2_0.h5")
elgauss6_0p32_0p44_4sig_particles = h5.root.num_particles.read()
elgauss6_0p32_0p44_4sig_xemit = h5.root.emitx.read()
h5.close()
h5 = tables.open_file(elgauss6_0p42_0p54_dir + "/full2_0.h5")
elgauss6_0p42_0p54_4sig_particles = h5.root.num_particles.read()
elgauss6_0p42_0p54_4sig_xemit = h5.root.emitx.read()
h5.close()
h5 = tables.open_file(elgauss6_0p52_0p64_dir + "/full2_0.h5")
elgauss6_0p52_0p64_4sig_particles = h5.root.num_particles.read()
elgauss6_0p52_0p64_4sig_xemit = h5.root.emitx.read()
h5.close()
h5 = tables.open_file(elgauss6_0p62_0p74_dir + "/full2_0.h5")
elgauss6_0p62_0p74_4sig_particles = h5.root.num_particles.read()
elgauss6_0p62_0p74_4sig_xemit = h5.root.emitx.read()
h5.close()
elgauss6_0p72_0p84_4sig_particles = elgauss6_4sig_particles
elgauss6_0p72_0p84_4sig_xemit = elgauss6_4sig_xemit
h5 = tables.open_file(elgauss6_0p82_0p94_dir + "/full2_0.h5")
elgauss6_0p82_0p94_4sig_particles = h5.root.num_particles.read()
elgauss6_0p82_0p94_4sig_xemit = h5.root.emitx.read()
h5.close()

plt.figure()
plt.title("Tune scan particle losses, 4 sigma aperture, 6.0 A-m fixed transverse, pulsed current lens")
plt.plot(norig-elgauss6_0p32_0p44_4sig_particles, label="(0.32,0.44) 4 sigma")
plt.plot(norig-elgauss6_0p42_0p54_4sig_particles, label="(0.42,0.54) 4 sigma")
plt.plot(norig-elgauss6_0p52_0p64_4sig_particles, label="(0.52,0.64) 4 sigma")
plt.plot(norig-elgauss6_0p62_0p74_4sig_particles, label="(0.62,0.74) 4 sigma")
plt.plot(norig-elgauss6_0p72_0p84_4sig_particles, label="(0.72,0.84) 4 sigma")
plt.plot(norig-elgauss6_0p82_0p94_4sig_particles, label="(0.82,0.94) 4 sigma")
plt.xlabel("turns")
plt.ylabel("cut particles")
plt.legend(loc='best')
saveplt()

print "tune scan cut losses"
print "(0.32, 0.44): ", 1.0 - float(elgauss6_0p32_0p44_4sig_particles[-1])/norig
print "(0.42, 0.54): ", 1.0 - float(elgauss6_0p42_0p54_4sig_particles[-1])/norig
print "(0.52, 0.64): ", 1.0 - float(elgauss6_0p52_0p64_4sig_particles[-1])/norig
print "(0.62, 0.74): ", 1.0 - float(elgauss6_0p62_0p74_4sig_particles[-1])/norig
print "(0.72, 0.84): ", 1.0 - float(elgauss6_0p72_0p84_4sig_particles[-1])/norig
print "(0.82, 0.94): ", 1.0 - float(elgauss6_0p82_0p94_4sig_particles[-1])/norig

h5 = tables.open_file(err05_4sig_scc0_dir + "/full2_0.h5")
err05_4sig_particles = h5.root.num_particles.read()
err05_4sig_xemit = h5.root.emitx.read()
err05_4sig_xstd = h5.root.std[0,:]

h5.close()
h5 = tables.open_file(err05_elgauss6_4sig_dir + "/full2_0.h5")
err05_elgauss6_4sig_particles = h5.root.num_particles.read()
h5.close()

plt.figure()
plt.title("particles cut 4 sigma aperture, no compensation")
plt.plot(turns, norig-part_4sig_00, label="1% lattice error")
plt.plot(turns, norig-err05_4sig_particles, label="5% lattice error")
plt.xlabel("turns")
plt.ylabel("particles")
plt.legend(loc='best')
saveplt()

plt.figure()
plt.title("err05 4 sigma cut x RMS emittance")
plt.plot(turns, err05_4sig_xemit)
plt.xlabel("turn")
plt.ylabel("x emittance")
saveplt()

plt.figure()
plt.title("err05 4 sigma cut x RMS")
plt.plot(turns, err05_4sig_xstd)
plt.xlabel("turns")
plt.ylabel("x RMS")
saveplt()

plt.figure()
#plt.title("err05 4 sigma cut 6.0 A-m fixed transverse, pulsed longitudinal EL")
plt.title("Particles cut 4 sigma aperture, 5% lattice error")
plt.plot(turns, norig-err05_4sig_particles, label="no compensation")
plt.plot(turns, norig-err05_elgauss6_4sig_particles, label="6.0 A-m fixed transverse, pulse longitudinal")
plt.xlabel("turns")
plt.ylabel("particles")
plt.legend(loc='best')
saveplt()

print
print "5% lattice error losses with 4 sigma aperture"
print "no compensation: ", 1.0 - float(err05_4sig_particles[-1])/norig
print "6.0 A-m fixed transverse, pulsed longitudinal: ", 1.0 - float(err05_elgauss6_4sig_particles[-1])/norig

#########################3

plt.figure()
plt.title("Integrated losses")
plt.plot(turns, (norig-part_4sig_00)/float(norig), lw=2, label="no compensation")
plt.plot(turns, (norig-elgauss6_4sig_particles)/float(norig), lw=2, label="66% compensation, gaussian pulsed lens")
plt.xlabel("turns", fontsize='large')
plt.ylabel("fraction lost", fontsize='large')
plt.legend(loc='best', fontsize='x-large')
saveplt()

print
print "1% lattice error losses with 4 sigma aperture"
print "no compensation: ", 1.0 - float(part_4sig_00[-1])/norig
print "6.0 A-m fixed transverse, pulsed longitudinal: ", 1.0 - float(elgauss6_4sig_particles[-1])/norig


plt.show()
