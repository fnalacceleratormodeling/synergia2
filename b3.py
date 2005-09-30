#!/usr/bin/env python
import sys
sys.path.append(".libs")
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

bfact    = bmlfactory("shmooster.mad")


#
# instantiate a beamline
#
pr = Proton( energy )
jp = JetProton(energy)

print "Reference BRho = ", pr.ReferenceBRho()

bcel01 = bfact.create_beamline("hack", pr.ReferenceBRho() )

### walk through beamline
i=0
bit = DeepBeamlineIterator(bcel01)
be = bit.reset()
be = bit.next()
while   be:
	print i, be.Type()
#        print dir(be)
	be = bit.next()
	i = i+1   

### Calculate line length
be = bit.reset()
be = bit.next()
line_length = 0
while   be:
        line_length += be.OrbitLength(pr)
	be = bit.next()
print "total length =", line_length

### Lattice functions
# bcel01.propagateJetParticle(jp)
# lfsage = LattFuncSage(bcel01,0)

# foo = lfsage.TuneCalc(jp, 1)
# lfsage.Disp_Calc(jp)

# lfsage.Fast_CS_Calc(jp)

# be = bit.reset()
# be = bit.next()
# while   be:
# 	lf  = be.dataHook.find("Twiss",1)
#         if lf:
#             print lf.arcLength, \
#                   lf.beta.hor, \
#                   lf.beta.ver, \
#                   lf.alpha.hor, \
#                   lf.alpha.ver, \
#                   lf.dispersion.hor, \
#                   lf.dispersion.ver, \
#                   lf.dPrime.hor, \
#                   lf.dPrime.ver
#         else:
#             print "No lf for", be.Type()

# 	be = bit.next()


### element maps
be = bit.reset()
be = bit.next()
while be:
    jpr = JetProton(energy)
    print be.Type()
    be.propagateJetParticle(jpr)
    ###print dir(jpr)
    print jpr.State().Jacobian()
    be = bit.next()


print "yup"
