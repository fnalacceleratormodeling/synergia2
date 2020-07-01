#!/usr/bin/env python
import sys,os
import numpy as np
import tables
import matplotlib.pyplot as plt

plt.rcParams['xtick.labelsize'] = 'large'
plt.rcParams['ytick.labelsize'] = 'large'
plt.rcParams['axes.labelsize'] = 'x-large'

def calculate_actions(coords):
    x = coords[0]
    xp = x*alpha_x + coords[1]*beta_x
    xaction = (x**2 + xp**2)/(2.0*beta_x)
    y = coords[2]
    yp = y*alpha_y + coords[3]*beta_y
    yaction = (y**2 + yp**2)/(2.0*beta_y)
    return (xaction, yaction)

# p is the particle array [npart, coord]
def calculate_emittances(p):
    #  these lattice functions are the unequal beta lattice
    # beta_x = 17.3259852015
    # alpha_x = 1.85063532729
    # beta_y = 17.2246414528
    # alpha_y = -1.90226103118

    # these are for the equal beta function lattice with tunes 3.74
    beta_x = 17.2336287
    alpha_x = 1.84393013
    beta_y = 17.2336287
    alpha_y = -1.84393013

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

    if len(sys.argv) > 2:
        plttitle = sys.argv[2]
    else:
        plttitle = ""

    h5 = tables.open_file(sys.argv[1])
    npart = h5.root.track_coords.shape[0]
    nsteps = h5.root.track_coords.shape[2]
    xemitrms = []
    xemit95 = []
    xemit99 = []
    xemit999 = []
    yemitrms = []
    yemit95 = []
    yemit99 = []
    yemit999 = []
    for stp in range(nsteps):
    #for stp in range(20):
        print "step: ", stp
        p = h5.root.track_coords[:,:,stp]
        (xem1, xem2, xem3, xem4, yem1, yem2, yem3, yem4) = calculate_emittances(p)
        xemitrms.append(xem1)
        xemit95.append(xem2)
        xemit99.append(xem3)
        xemit999.append(xem4)
        yemitrms.append(yem1)
        yemit95.append(yem2)
        yemit99.append(yem3)
        yemit999.append(yem4)

    plt.figure()
    plt.title("x emittances "+plttitle)
    plt.plot(xemitrms, label="x rms")
    plt.plot(xemit95, label="x 95%")
    plt.plot(xemit99, label="x 99%")
    plt.plot(xemit999, label="x 99.9%")
    plt.xlabel("cell")
    plt.ylabel("emittance [pi m-rad]")
    plt.legend(loc='best')

    plt.figure()
    plt.title("x emittances/initial rms emittance "+plttitle)
    plt.plot(xemitrms/xemitrms[0], label="x rms")
    plt.plot(xemit95/xemitrms[0], label="x 95%")
    plt.plot(xemit99/xemitrms[0], label="x 99%")
    plt.plot(xemit999/xemitrms[0], label="x 99.9%")
    plt.xlabel("cell")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')
    plt.savefig("xemittances_growth_"+plttitle+".png")

    plt.figure()
    plt.title("x emittance growth " + plttitle)
    ax1 = plt.subplot(211)
    ax1.set_title("x RMS emittance growth " + plttitle)
    ax1.plot(xemitrms/xemitrms[0])
    plt.setp(ax1.get_xticklabels(), visible=False)
    #ax1.set_xticklabels(ax1.get_xticklabels(), visible=False)
    #ax1.legend(loc='best')

    ax2 = plt.subplot(212, sharex=ax1)
    ax2.set_title("x 99.9% emittance growth " + plttitle)
    ax2.plot(xemit999/xemit999[0])
    #ax2.legend(loc='best')

    ax2.set_xlabel("cell")
    ax2.set_ylabel("emittance growth")
    plt.savefig("xemit_rms_99.9_growth_"+plttitle+".png")

    # plt.figure()
    # plt.title("x emittances relative growth")
    # plt.plot(xemitrms/xemitrms[0], label="x emit rms growth")
    # plt.plot(xemit95/xemit95[0], label="x 95% growth")
    # plt.plot(xemit99/xemit99[0], label="x 99% growth")
    # plt.plot(xemit999/xemit999[0], label="x 99.9% growth")
    # plt.xlabel("cell")
    # plt.ylabel("emittance growth")
    # plt.legend(loc='best')

    plt.figure()
    plt.title("y emittances/initial rms emittance")
    plt.plot(yemitrms/yemitrms[0], label="y rms")
    plt.plot(yemit95/yemitrms[0], label="y 95%")
    plt.plot(yemit99/yemitrms[0], label="y 99%")
    plt.plot(yemit999/yemitrms[0], label="y 99.9%")
    plt.xlabel("cell")
    plt.ylabel("emittance growth")
    plt.legend(loc='best')

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
