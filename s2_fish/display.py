#!/usr/bin/env python

import pylab
import matplotlib
import math
import numpy

def display888(sf):
    index = 1
    max = sf.max()
    for i in range(0,4):
        for j in range(0,2):
            pylab.subplot(2,4,index)
            pylab.pcolor(sf[:,:,index-1],
                         norm=matplotlib.colors.normalize(vmax=max,vmin=0),
                         vmax=max,vmin=0.0)
            index +=1

def display888_trans(sf):
    index = 1
    max = sf.max()
    for i in range(0,4):
        for j in range(0,2):
            pylab.subplot(2,4,index)
            pylab.pcolor(numpy.transpose(sf[index-1,:,:]),
                         norm=matplotlib.colors.normalize(vmax=max,vmin=0),
                         vmax=max,vmin=0.0)
            index +=1

def display161616(sf):
    index = 1
    max = sf.max()
    for i in range(0,4):
        for j in range(0,4):
            pylab.subplot(4,4,index)
            pylab.pcolor(sf[:,:,index-1],
                         norm=matplotlib.colors.normalize(vmax=max,vmin=0),
                         vmax=max,vmin=0.0)
            index +=1

def display16916(sf):
    index = 1
    for i in range(0,4):
        for j in range(0,4):
            pylab.subplot(4,4,index)
            pylab.pcolor(sf[:,:,index-1])
            index +=1

def display_axes(sf,index):
    pylab.subplot(2,2,1)
    pylab.plot(sf[index[2],index[1],:])
    pylab.subplot(2,2,2)
    pylab.plot(sf[index[2],:,index[0]])
    pylab.subplot(2,2,3)
    pylab.plot(sf[:,index[1],index[0]])

def display_axes_transpose(sf,index):
    pylab.subplot(2,2,1)
    pylab.plot(sf[:,index[1],index[2]])
    pylab.subplot(2,2,2)
    pylab.plot(sf[index[0],:,index[2]])
    pylab.subplot(2,2,3)
    pylab.plot(sf[index[0],index[1],:])

def print_ratio(ours, theirs,show_pass=0):
    zerotol = 1.0e-10
    tol = 1.0e-5
    passed = 0
    for i in range(0,8):
        for j in range(0,8):
            for k in range (0,8):
                if theirs[k,j,i] != 0.0:
                    ratio = ours[i,j,k]/theirs[k,j,i]
                    if abs(ratio-1.0)>tol:
                        try:
                            print i,j,k,ours[i,j,k],theirs[k,j,i], \
                                  "ratio=",ratio
                        except:
                            print "funky in",i,j,k
                    else:
                        if show_pass:
                            print "passed",i,j,k
                        passed += 1
                else:
                    diff = ours[i,j,k]-theirs[k,j,i]
                    if abs(diff) > zerotol:
                        print i,j,k,ours[i,j,k],theirs[k,j,i],"diff=",diff
                    else:
                        if show_pass:
                            print "passed",i,j,k,ours[i,j,k],theirs[k,j,i],"diff=",diff
                        passed += 1
    print passed,"passed"
                
if __name__ == "__main__":
    import arrayio

    grn = arrayio.readfort("../fort.50",[9,9,9])
    G2 = arrayio.readpetscvec("G2_petsc")
    print "G2 compared to grn"
    print_ratio(G2.real,grn)
#     display161616(G2.real)
#     pylab.figure()
#     display161616(G2.imag)
#     display_axes_transpose(grn,[0,0,0])
#     display_axes(G2,[0,0,0])
    
    pylab.figure(1)
    impact_rho = arrayio.readfort("../fort.40",[8,8,8])
    display888_trans(impact_rho)

    pylab.figure(2)
    rho = arrayio.readtensor("rho")
    display888(rho)

    pylab.figure(5)
    rho_in = arrayio.readscalarfield("rho_in")
    display888(rho_in)

#     print "rho compared"
#     print_ratio(rho,impact_rho)

#     pylab.figure()
#     phi2 = arrayio.readpetscvec("phi2_petsc")
#     display161616(phi2.real)

    pylab.figure(3)
    impact_phi = arrayio.readfort("../fort.60",[8,8,8])
    display888(impact_phi)

    pylab.figure(4)
    phi = arrayio.readtensor("phi")
    display888(phi)

#     print "phi compared"
#     print_ratio(phi,impact_phi)

#     pylab.figure()
#     rho_hat_real_impact = arrayio.readfort("../fort.41",[16,9,16])
#     rho_hat_imag_impact = arrayio.readfort("../fort.42",[16,9,16])
# #    display888_trans(rho_hat_real_impact[0:8,0:8,0:8])
#     display888_trans(rho_hat_imag_impact[0:8,0:8,0:8])

#     pylab.figure()
#     rho_hat = arrayio.readpetscvec("complex_rho_hat2_petsc")
#     display888(rho_hat.real[0:8,0:8,0:8])
#     pylab.figure()
# #    display161616(rho_hat.real)
#     display161616(rho_hat.imag)

#     print "compare Re(rho hat)"
#     print_ratio(rho_hat_real_impact[0:8,0:8,0:8],rho_hat.real[0:8,0:8,0:8],
#                 1)

#     print "compare Im(rho hat)"
#     print_ratio(rho_hat_imag_impact[0:8,0:8,0:8],rho_hat.imag[0:8,0:8,0:8],
#                 1)

#     pylab.figure()
#     phi2 = arrayio.readpetscvec("phi2_petsc")
#     display161616(phi2.real)
#     pylab.figure()
#     reGhat = arrayio.readfort("../fort.45",[16,9,16])
#     display_axes(reGhat,[8,5,8])
#     pylab.figure()
#     display16916(reGhat)

#     pylab.figure()
#     imGhat = arrayio.readfort("../fort.46",[16,9,16])
#     display_axes(imGhat,[8,5,8])
#     pylab.figure()
#     display16916(imGhat)

    
    pylab.show()
