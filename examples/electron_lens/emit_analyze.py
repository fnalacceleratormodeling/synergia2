import numpy as np
import tables
import os
from glob import glob
import scipy
import scipy.stats

class emittances:
    def __init__(self, turn, alphax, betax, alphay, betay, xrms, yrms, xemitrms, xemit95, xemit99, xemit999, xkurt, yemitrms, yemit95, yemit99, yemit999, ykurt):
        self.turn = turn
        self.alpha_x = alphax
        self.beta_x = betax
        self.alpha_y = alphay
        self.beta_y = betay
        self.xrms = xrms
        self.yrms = yrms
        self.xemitrms = xemitrms
        self.xemit95 = xemit95
        self.xemit99 = xemit99
        self.xemit999 = xemit999
        self.xkurt = xkurt
        self.yemitrms = yemitrms
        self.yemit95 = yemit95
        self.yemit99 = yemit99
        self.yemit999 = yemit999
        self.ykurt = ykurt

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

def file_emittances(filename):
    h5 = tables.open_file(filename)
    turn = h5.root.rep[()]
    # for some reason, rep is not incremented for the first
    # set of files, so it is 25 more than what the file shows
    # for any file on _0000.h5.
    if not filename.endswith("_0_0000.h5"):
        turn += 25

    p = h5.root.particles

    alpha_x, beta_x, alpha_y, beta_y = calculate_alpha_beta(p)

    xactions = (p[:,0]**2 + (alpha_x*p[:,0]+beta_x*p[:,1])**2)/(2.0*beta_x)
    yactions = (p[:,2]**2 + (alpha_y*p[:,2]+beta_y*p[:,3])**2)/(2.0*beta_y)

    xrms = p[:,0].std()
    yrms = p[:,2].std()

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
    xkurt = scipy.stats.kurtosis(p[:,0])
    yemitrms = yactions.mean()
    yemit95 = yactions[i95]
    yemit99 = yactions[i99]
    yemit999 = yactions[i999]
    ykurt = scipy.stats.kurtosis(p[:,2])

    h5.close()

    return emittances(turn, alpha_x, beta_x, alpha_y, beta_y, xrms, yrms, xemitrms, xemit95, xemit99, xemit999, xkurt, yemitrms, yemit95, yemit99, yemit999, ykurt)

class dir_emittances:
    def __init__(self, dir):
        pfile_list = glob(dir+"/turn_particles_0_*.h5")
        pfile_list.sort()
        #print "pfile_list: ", pfile_list
        e_list = []
        for f in pfile_list:
            #print "reading file ", f
            e_list.append(file_emittances(f))
        self.turns = np.array([e.turn for e in e_list])
        self.alphas_x = np.array([e.alpha_x for e in e_list])
        self.betas_x = np.array([e.beta_x for e in e_list])
        self.alphas_y = np.array([e.alpha_y for e in e_list])
        self.betas_y = np.array([e.beta_y for e in e_list])
        self.xrms = np.array([e.xrms for e in e_list])
        self.yrms = np.array([e.yrms for e in e_list])
        self.xemitrms = np.array([e.xemitrms for e in e_list])
        self.xemit999 = np.array([e.xemit999 for e in e_list])
        self.xkurt = np.array([e.xkurt for e in e_list])
        self.yemitrms = np.array([e.yemitrms for e in e_list])
        self.yemit999 = np.array([e.yemit999 for e in e_list])
        self.ykurt = np.array([e.ykurt for e in e_list])

# if called as the main program, calculate the emittances for the current
# directory and save them as a numpy file

if __name__ == "__main__":
    emittances = dir_emittances(".")
    np.save("emittances.npy", emittances)
