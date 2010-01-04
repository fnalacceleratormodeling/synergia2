#!/usr/bin/env python

import local_paths
import sys

from bmlfactory import *
from basic_toolkit import *
from mxyzptlk import *
from beamline import *
from physics_toolkit import *
from physics_constants import *

import loadfile

order = 1
JetParticle.createStandardEnvironments(order)
file_name = 'ertml_filecalls.xsif'
lattice_name = 'ERTML'
xfactory = XSIFFactory(file_name)
energy = 5.0
positron = Positron(energy)
brho = positron.ReferenceBRho()
beamline_orig = xfactory.create_beamline(lattice_name,brho)
#~ beamline_orig.flatten()
bmln = DriftsToSlots(beamline_orig)
context = BeamlineContext(positron,bmln)

BETX = 1.6479
ALFX =  0.4982
BETY = 8.8630
ALFY = -2.3771

initialCovariance = CovFunc()
initialCovariance.arcLength = 0.0
initialCovariance.beta.hor  = BETX
initialCovariance.beta.ver  = BETY
initialCovariance.alpha.hor = ALFX
initialCovariance.alpha.ver = ALFY
errcode = BmlUtil.makeCovariance( initialCovariance, positron )
context.setInitialCovariance( initialCovariance )

initiallattFunc = lattFunc()
initiallattFunc.arcLength = 0.0
initiallattFunc.beta.hor  = BETX
initiallattFunc.beta.ver  = BETY
initiallattFunc.alpha.hor = ALFX
initiallattFunc.alpha.ver = ALFY
initiallattFunc.dispersion.hor = 0.0
initiallattFunc.dispersion.ver = 0.0
initiallattFunc.dPrime.hor = 0.0
initiallattFunc.dPrime.ver = 0.0
initiallattFunc.psi.hor = 0.0
initiallattFunc.psi.ver= 0.0
context.setInitialTwiss(initiallattFunc)

cov_array = context.getCovarianceArray()
#~ for cov in cov_array:
    #~ print cov.arcLength, cov.beta.hor, cov.beta.ver

initialdis = DispFunc()
initialdis.arcLength = 0.0
initialdis.closedOrbit.hor = 0.0
initialdis.closedOrbit.ver= 0.0
initialdis.dispersion.hor = 0.0
initialdis.dispersion.ver= 0.0
initialdis.dPrime.hor = 0.0
initialdis.dPrime.ver= 0.0

context.setInitialDispersion( initialdis )
dis_array = context.getDispersionArray()
#~ for dis in dis_array:
    #~ print dis.arcLength, dis.dispersion.hor, dis.dispersion.ver

twiss_array = context.getTwissArray()
#~ for twiss in twiss_array:
    #~ print twiss.arcLength,twiss.beta.hor,twiss.beta.ver,\
        #~ twiss.dispersion.hor,twiss.dispersion.ver

#~ for i in range(0,len(cov_array)):
    #~ cov = cov_array[i]
    #~ twiss = twiss_array[i]
    #~ dis = dis_array[i]
    #~ print cov.arcLength, cov.beta.hor, cov.beta.ver, \
        #~ dis.dispersion.hor, dis.dispersion.ver, \
        #~ twiss.beta.hor,twiss.beta.ver

pt = loadfile.loadfile("pt.dat");
index_chef = 0
index_pt = 0
f = open("test6_out.dat",'w')
for index in range(0,len(pt)):
    s_chef = cov_array[index_chef].arcLength
    s_pt = pt[index_pt,0]
    if abs(s_chef - s_pt) < 1.0e-5:
        f.write("%g %g %g %g %g" % tuple(pt[index_pt,:]))
        f.write(" %g %g" % (cov_array[index_chef].beta.hor,\
            cov_array[index_chef].beta.ver))
        f.write(" %g %g" % (twiss_array[index_chef].beta.hor,\
            twiss_array[index_chef].beta.ver))
        f.write(" %g %g" % (dis_array[index_chef].dispersion.hor,\
            dis_array[index_chef].dPrime.hor))        
        f.write("\n");
        index_chef += 1
        index_pt +=1
    elif s_chef > s_pt:
        index_pt +=1
    elif s_pt > s_chef:
        index_chef +=1
    else:
        print "didn't expect to get here"
f.close()

print "success!"
