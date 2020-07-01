#!/usr/bin/env python
import sys, os
import matplotlib.pyplot as plt

# calculate lattice functions based on particle distribution.
# this presupposes that the planes are uncoupled, the normal form exists
# and the bunch distribution is in equilibrium.
def calculate_alpha_beta(p):
    m = np.cov(p[:, 0:4], rowvar=False)

    mxx  = m[0,0]
    mxxp = m[0,1]
    mxpxp= m[1,1]

    myy  = m[2,2]
    myyp = m[2,3]
    mypyp= m[3,3]

    betaxsq = mxx**2/(mxx * mxpxp - mxxp**2)
    betax = np.sqrt(betaxsq)
    alphax = -betax * mxxp/mxx

    betaysq = myy**2/(myy * mypyp - myyp**2)
    betay = np.sqrt(betaysq)
    alphay = -betay * myyp/myy

    return (alphax, betax, alphay, betay)

#
# calculate the turn by turn lattice functions for a run
def run_lf(diag_file):
    h5diag = tables.open_file(diag_file)
    nsteps = h5diag.root.mom2.shape[2]
    b_x = []
    a_x = []
    b_y = []
    a_y = []
    for i in range(nsteps):
        m = h5diag.root.mom2[0:4, 0:4, i]
        ax, bx, ay, by = calculate_alpha_beta(m)
        b_x.append(bx)
        b_y.append(by)
        a_x.append(ax)
        a_y.append(ay)

    h5diag.close()
    return np.array(a_x), np.array(b_x), np.array(a_y), np.array(b_y)

# calculate emittance growths for turns 50 to 1000 for the bulk track file
# return array of xemit growth
def emit_growths(bulk_file, alphas_x, betas_x, i0, i1):
    h5bulk = tables.open_file(bulk_file)
    npart = h5bulk.root.track_coords.shape[0]
    egrowth = np.zeros(npart)
    a0 = alphas_x[i0]
    a1 = alphas_x[i1]
    b0 = betas_x[i0]
    b1 = betas_x[i1]
    for p in range(npart):
        p0 = h5bulk.root.track_coords[p, :, i0]
        p1 = h5bulk.root.track_coords[p, :, i1]
        x0 = p0[0]
        x1 = p1[0]
        xp0 = p0[1]
        zp1 = p1[1]
        e0 = (x0**2 + (a0*x0 + b0*xp0)**2)/(2*b0)
        e1 = (x1**2 + (a1*x1 + b1*xp1)**2)/(2*b1)
        egrowth[p] = e1/e0

    h5bulk.close()
    return egrowth
