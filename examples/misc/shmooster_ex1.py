#!/usr/bin/env python

gourmet_path = "/home3/amundson/work/fnal/branches/jfa1/python-bindings/src/.libs"
impact_path = "/home3/amundson/work/Layer-kruger_0_1/Forthon_Interfaces"
###mad_path = "/home3/amundson/work/Layer-head"
mad_path = "/home3/amundson/work/fnal/branches/jfa1/python-bindings/src"

import sys
sys.path.append(gourmet_path)
sys.path.append(impact_path)

from basic_toolkit import *
from mxyzptlk       import *
from beamline       import *
from physics_toolkit  import *
from bmlfactory      import *

PH_NORM_mp  = 0.938               #  [GeV]  Proton Rest Mass 
energy = 0.4 + PH_NORM_mp

#
# define a first order environment in 6D phase space
#
order = 1
Jet.BeginEnvironment(order)
x    = coord(0.0)
y    = coord(0.0)
ct   = coord(0.0)
npx  = coord(0.0)
npy  = coord(0.0)
np   = coord(0.0)
JetC.setLastEnv(JetC.CreateEnvFrom(Jet.EndEnvironment()))

bfact    = bmlfactory("%s/shmooster.mad"%mad_path)


#
# instantiate a beamline
#
pr = Proton( energy )
jp = JetProton(energy)

print "Reference BRho = ", pr.ReferenceBRho()

bcel01 = bfact.create_beamline("bcel01", pr.ReferenceBRho() )
bit = DeepBeamlineIterator(bcel01)

### walk through beamline
# i=0
# be = bit.reset()
# be = bit.next()
# while   be:
# 	print i, be.Type()
# 	be = bit.next()
# 	i = i+1   

### Calculate line length
be = bit.reset()
be = bit.next()
line_length = 0
while   be:
        line_length += be.OrbitLength(pr)
	be = bit.next()
print "total length =", line_length

### Lattice functions
bcel01.propagateJetParticle(jp)
lfsage = LattFuncSage(bcel01,0)

foo = lfsage.TuneCalc(jp, 1)
lfsage.Disp_Calc(jp)

lfsage.Fast_CS_Calc(jp)

be = bit.reset()
be = bit.next()
while   be:
	lf  = be.dataHook.find("Twiss",1)
        if lf:
            print lf.arcLength, \
                  lf.beta.hor, \
                  lf.beta.ver, \
                  lf.alpha.hor, \
                  lf.alpha.ver, \
                  lf.dispersion.hor, \
                  lf.dispersion.ver, \
                  lf.dPrime.hor, \
                  lf.dPrime.ver
        else:
            print "No lf for", be.Type()

	be = bit.next()


### element maps
maps = []
names = []
types = []
lengths = []
be = bit.reset()
be = bit.next()
while be:
    jpr = JetProton(energy)
    print be.Type()
    be.propagateJetParticle(jpr)
###    print jpr.State().Jacobian()
    maps.append(jpr.State().Jacobian())
    types.append(be.Type())
    names.append(be.Name())
    lengths.append(be.OrbitLength(pr))
    if be.OrbitLength(pr) == 0.0:
        print "element",be.Name(),"of type",be.Type(),"has zero length"
                  
    be = bit.next()

print types
print names
print lengths

print "I guess I'll just hang now..."

