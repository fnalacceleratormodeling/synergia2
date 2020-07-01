#!/usr/bin/env python
import sys,os
import numpy as np
import tables
import matplotlib.pyplot as plt

def calculate_actions(coords):
    x = coords[0]
    xp = x*bpm_alpha_x + coords[1]*bpm_beta_x
    xaction = (x**2 + xp**2)/(2.0*bpm_beta_x)
    y = coords[2]
    yp = y*bpm_alpha_y + coords[3]*bpm_beta_y
    yaction = (y**2 + yp**2)/(2.0*bpm_beta_y)
    return (xaction, yaction)

# p is the particle array [npart, coord]
def calculate_emittances(p):

    bpm_beta_x = 4.70792560035
    bpm_alpha_x = -0.449832076233
    bpm_beta_y = 46.1594967296
    bpm_alpha_y = 3.37300523129

    beta_x = 17.3259852015
    alpha_x = 1.85063532729
    beta_y = 17.2246414528
    alpha_y = -1.90226103118

    xactions = (p[:,0]**2 + (alpha_x*p[:,0]+beta_x*p[:,1])**2)/(2.0*beta_x)
    yactions = (p[:,2]**2 + (alpha_y*p[:,2]+beta_y*p[:,3])**2)/(2.0*beta_y)

    npart = p.shape[0]

    xactions.sort()
    yactions.sort()
    i95 = int(npart * 0.95 + 1)
    i99 = int(npart * 0.99 + 1)
    i999 = int(npart * 0.999 + 1)
    xemitrms = xactions.mean()
    xemit95 = xactions[i95]
    xemit99 = xactions[i99]
    xemit999 = xactions[i999]
    yemitrms = yactions.mean()
    yemit95 = yactions[i95]
    yemit99 = yactions[i99]
    yemit999 = yactions[i999]
    return (xemitrms, xemit95, xemit99, xemit999, yemitrms, yemit95, yemit99, yemit999)

if __name__ == "__main__":
    os.environ["HDF5_DISABLE_VERSION_CHECK"] = "2"

    filelist = sys.argv[1:]
    xemitrms = []
    xemit95 = []
    xemit99 = []
    xemit999 = []
    yemitrms = []
    yemit95 = []
    yemit99 = []
    yemit999 = []
    for fn in filelist:
        print "reading ", fn
        h5 = tables.open_file(fn)
        p = h5.root.particles.read()
        npart = p.shape[0]
        h5.close()
        (xem1, xem2, xem3, xem4, yem1, yem2, yem3, yem4) = calculate_emittances(p)
        xemitrms.append(xem1)
        xemit95.append(xem2)
        xemit99.append(xem3)
        xemit999.append(xem4)
        yemitrms.append(yem1)
        yemit95.append(yem2)
        yemit99.append(yem3)
        yemit999.append(yem4)
        del(p)

    turns = np.arange(0.0, 1001.0, 25.0)

    plt.figure()
    plt.title("x emittances")
    plt.plot(turns, xemitrms,'-o', label="x rms")
    plt.plot(turns, xemit95, '-o', label="x 95%")
    plt.plot(turns, xemit99, '-o', label="x 99%")
    plt.plot(turns, xemit999, '-o', label="x 99.9%")
    plt.xlabel("turn")
    plt.ylabel("emittance")
    plt.legend(loc='upper right')

    plt.figure()
    plt.title("x emittances/emitx0")
    plt.plot(xemitrms/xemitrms[0],'-o', label="x rms")
    plt.plot(xemit95/xemitrms[0], '-o', label="x 95%")
    plt.plot(xemit99/xemitrms[0], '-o', label="x 99%")
    plt.plot(xemit999/xemitrms[0], '-o', label="x 99.9%")
    plt.xlabel("turn")
    plt.ylabel("emittance growth")
    plt.legend(loc='upper right')
    plt.savefig("x_emittances_growth.png")

    plt.figure()
    plt.title("x emittance growth RMS and 99.9%")
    ax1 = plt.subplot(211)
    ax1.set_title("x RMS emittance growth")
    ax1.plot(turns, xemitrms/xemitrms[0], '-o')
    plt.setp(ax1.get_xticklabels(), visible=False)
    #ax1.set_xticklabels(ax1.get_xticklabels(), visible=False)
    #ax1.legend(loc='best')

    ax2 = plt.subplot(212, sharex=ax1)
    ax2.set_title("x 99.9% emittance growth")
    ax2.plot(turns,xemit999/xemit999[0], '-o')
    #ax2.legend(loc='best')

    ax2.set_xlabel("turn")
    ax2.set_ylabel("emittance growth")
    plt.savefig("xemit_growth_rms_999.png")

    plt.figure()
    plt.title("y emittances")
    plt.plot(turns, yemitrms, '-o', label="y rms")
    plt.plot(turns, yemit95,'-o', label="y 95%")
    plt.plot(turns, yemit99, '-o', label="y 99%")
    plt.plot(turns, yemit999, '-o', label="y 99.9%")
    plt.xlabel("turn")
    plt.ylabel("emittance")
    plt.legend(loc='upper right')

    plt.show()
    # h5 = tables.openFile(sys.argv[1])
    # p = h5.root.particles.read()
    # h5.close()
    # (xemitrms, xemit95, xemit99, yemitrms, yemit95, yemit99) = calculate_emittances(p)
    # print "xemitrms: ", xemitrms
    # print "xemit95: ", xemit95
    # print "xemit99: ", xemit99
    # print "yemitrms: ", yemitrms
    # print "yemit95: ", yemit95
    # print "yemit99: ", yemit99
