#!/usr/bin/env python

import sys
import os
import numpy as np

# eig_match_tolerance is how close to reciprocal do eigenvalues need to be
# to signal a match.  This may need to be tuned for different maps.
eig_match_tolerance = 1.0e-6
# conj_match_tolerance is how close to conjugates matching eigenvectors
# should be
conj_match_tolerance = 1.0e-8

#-------------------------------------------------------------

# given an input map of dimension 2n return the matrix that
# transforms the coordinates into normal form coordinates and the
# inverse of that matrix which does the reverse.  The matrix is
# "sensibly-ordered" by the Michelotti definition.
def normal_form(map):
    if len(map.shape) != 2:
        raise TypeError, "input map should be a 2-dimension array!"

    twondim = map.shape[0]
    if map.shape[1] != twondim:
        raise ValueError, "input matrix is not square!"

    if twondim%2 != 0:
        raise ValueError, "input matrix must have even size"

    ndim = twondim/2

    # generate "J" matrix based on input dimension
    bigJ = np.zeros((twondim,twondim), 'd')
    bigJ[0:ndim,ndim:twondim] = np.eye(ndim, dtype='d')
    bigJ[ndim:twondim, 0:ndim] = -np.eye(ndim, dtype='d')

    # get eigenvalues and eigenvectors
    lmap, vmap = np.linalg.eig(map)

    # reorder eigenvectors
    vreorder = np.zeros((twondim,twondim),dtype='D') # this has to be complex
    lreorder = np.zeros((twondim,twondim),dtype='D')

    unused = [1]*twondim
    filling_coord=0
    while unused.count(1) != 0:
        working_on = unused.index(1)
        eig1 = lmap[working_on]
        unused[working_on] = 0
        for find_pair in range(working_on+1,twondim):
            if unused[find_pair] and abs(eig1*lmap[find_pair] - 1.0) < eig_match_tolerance:
                eig2 = lmap[find_pair]
                #found it
                unused[find_pair] = 0
                lreorder[filling_coord] = eig1
                lreorder[filling_coord+ndim] = eig2
                vreorder[:,filling_coord] = vmap[:,working_on]
                vreorder[:,filling_coord+ndim] = vmap[:,find_pair]
                filling_coord = filling_coord+1
                break
        else:
            print [lmap[jj] for jj in range(twondim) if unused[jj]!=0]
            raise RuntimeError, "Didn't find matching eigenvalue pair for eigenvalue: %s!!"%repr(eig1)

    # make sure that corresponding eigenvectors are also complex conjugates
    for j in range(ndim):
        diff1 = vreorder[:,j] - vreorder[:,j+ndim].conj()   
        if np.vdot(diff1, diff1) > conj_match_tolerance:
            raise RuntimeError, "matched eigenvectors %d and %d are not complex conjugates"%(j,j+ndim)

    # Construct i*vreorder.T * J * vreorder * J (should be diagonal)

    VTJVJ = 1.0j*np.dot(np.dot(np.dot(vreorder.T, bigJ), vreorder), bigJ)

    # if a diagonal element of VTJVJ is negative, switch it with it's
    # conjugate and switch eigenvalue array also (although I won't use it)
    for j in range(ndim):
        if VTJVJ[j,j].real < 0:
            tcolumn = vreorder[:,j].copy()
            teigv = lreorder[j].copy()
            vreorder[:,j] = vreorder[:,j+ndim]
            lreorder[j] = lreorder[j+ndim]
            vreorder[:,j+ndim] = tcolumn
            lreorder[j+ndim] =teigv


    # calculate VTJVJ again to normalize eigenvalues
    VTJVJ = 1.0j*np.dot(np.dot(np.dot(vreorder.T, bigJ), vreorder), bigJ)

    # normalize eigenvectors
    for j in range(ndim):
        nfact = np.sqrt(VTJVJ[j,j].real)
        vreorder[:,j] = vreorder[:,j]/nfact
        vreorder[:,j+ndim] = vreorder[:,j+ndim]/nfact

    # Calculate vreorder^{-1} = i J vreorder.T J
    vinv = 1.0j * np.dot(np.dot(bigJ, vreorder.T), bigJ)

    return (vreorder, vinv)
              
#-------------------------------------------------------------

# Given a eigenvector matrix which converts to normal form, and an
#   list of standard deviations, return array of mean actions for
#   each plane such that particles generated with those mean
#   actions will have the correct standard deviations. nfmatrix will
#   be complex in general.  This is solving equation (16) for I_k given
#   the matrix B and C_ii in the Fermilab-FN-0826-CD linear normal form memo.

def get_mean_actions(nfmatrix, stds):

    # check input parameters.  nfmatrix must be square and have
    # dimensionality of 2Nx2N.  stds must have length N.

    if len(nfmatrix.shape) != 2:
        raise TypeError, "normal form matrix should be a 2-dimension array!"

    twondim = nfmatrix.shape[0]
    if nfmatrix.shape[1] != twondim:
        raise ValueError, "input matrix is not square!"

    if twondim%2 != 0:
        raise ValueError, "input matrix must have even size"

    ndim = twondim/2

    if len(stds) != ndim:
        raise ValueError, "size of the list of stds(%d) does not match size of normal form transformation matrix(%d)"%(len(stds),ndim)

    cmat = np.zeros((ndim,ndim),'d')
    for i in range(ndim):
        for j in range(ndim):
            cmat[i,j] = 2.0*abs(nfmatrix[i,j])**2

    Ik = np.linalg.solve(cmat, stds)
    return Ik

#-------------------------------------------------------------

# generate npart gaussian distribution of particles in normal coordinate space
# with given mean actions and normal form transformation matrix.  Transform
# the particles to normal space.

# The returned particle array has shape (ndim, npart) where ndim
# is the number of normal form dimensions determined by the shape
# of the normal form matrix.

# the normal form coordinates are generated following equation (13) in the
# Fermilab-FN-0826-CD memo.

def generate_particles(npart, actions, nfmatrix):
    # check input parameters.  nfmatrix must be square and have
    # dimensionality of 2Nx2N.  stds must have length N.

    if len(nfmatrix.shape) != 2:
        raise TypeError, "normal form matrix should be a 2-dimension array!"

    twondim = nfmatrix.shape[0]
    if nfmatrix.shape[1] != twondim:
        raise ValueError, "input matrix is not square!"

    if twondim%2 != 0:
        raise ValueError, "input matrix must have even size"

    ndim = twondim/2

    if len(actions) != ndim:
        raise ValueError, "size of the list of actions does not match size of normal form transformation matrix"

    nfpart = np.zeros((twondim, npart),'D')

    # generate particles in normal form space
    for pln in range(ndim):
        # generate accotding to eqn. 13
        # action distributed exponentially with mean actions[pln], with
        # (complex) phase uniformly 0 to 2*pi.
        nfpart[pln,:] = 1.0j * np.sqrt(np.random.exponential(scale=actions[pln],size=npart)) * np.exp(-2.0j*np.pi*np.random.uniform(size=npart))
        nfpart[pln+ndim,:] = nfpart[pln,:].conj()

    # convert to regular space
    particles = np.dot(nfmatrix, nfpart)
    # I expect no large imaginary components
    impart = np.amax(abs(particles.imag))
    if impart > 1.0e-8:
        print "?? error, large imaginary component"

    return particles.real
                                       
#-------------------------------------------------------------

# remaps matrix map from x-x'-y-y' ordering to "sensibly ordered"
def sensible_order(map):
    if len(map.shape) != 2:
        raise TypeError, "input map should be a 2-dimension array!"

    twondim = map.shape[0]
    if map.shape[1] != twondim:
        raise ValueError, "input matrix is not square!"

    if twondim%2 != 0:
        raise ValueError, "input matrix must have even size"

    ndim = twondim/2

    shufflemap = [0]*twondim

    todim=0
    for j in range(ndim):
        fromdim = 2*j
        shufflemap[fromdim] = todim
        shufflemap[fromdim+1] = todim+ndim
        todim = todim+1

    remap = np.zeros((twondim,twondim),type(map[0,0]))

    for i in range(twondim):
        newi = shufflemap[i]
        for j in range(twondim):
            newj = shufflemap[j]
            remap[newi,newj] = map[i,j]

    return remap

#-------------------------------------------------------------

# remaps matrix map from "sensibly-ordered" to x-x'-y-y' ordering.
def usual_order(map):
    if len(map.shape) != 2:
        raise TypeError, "input map should be a 2-dimension array!"

    twondim = map.shape[0]
    if map.shape[1] != twondim:
        raise ValueError, "input matrix is not square!"

    if twondim%2 != 0:
        raise ValueError, "input matrix must have even size"

    ndim = twondim/2

    shufflemap = [0]*twondim

    for j in range(ndim):
        todim = 2*j
        shufflemap[j] = todim
        shufflemap[j+ndim] = todim+1

    remap = np.zeros((twondim,twondim),type(map[0,0]))

    for i in range(twondim):
        newi = shufflemap[i]
        for j in range(twondim):
            newj = shufflemap[j]
            remap[newi,newj] = map[i,j]

    return remap
