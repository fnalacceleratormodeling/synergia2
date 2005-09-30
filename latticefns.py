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
from physics_constants import *

import math

accuracy_marker = marker()
pacifier = drift("pacifier",0.0)

def insert_markers(bmline, num_markers_per_element, energy, momentum):
	pr = Proton(energy)
	deep_beamline_iterator = DeepBeamlineIterator(bmline);
	master_insertion_point = 0.0
	insertion_list = InsertionList(momentum)
	element = deep_beamline_iterator.next()
	ile_list = []
	while element:
		if element.OrbitLength(pr) > 0:
			marker_interval = element.OrbitLength(pr)/ \
			(num_markers_per_element + 1.0)
			insertion_point = master_insertion_point
			for i in range(0,num_markers_per_element):
				insertion_point += marker_interval
				ile = InsertionListElement(insertion_point,accuracy_marker)
				ile_list.append(ile)
				insertion_list.Append(ile)
		master_insertion_point += element.OrbitLength(pr)
		element = deep_beamline_iterator.next()
	removed_elements = slist()
	s_0 = 0.0
	print "about to insert"
	bmline.InsertElementsFromList(s_0, insertion_list, removed_elements)
	print "inserted"
	bmline.append(accuracy_marker)
	
#mad_file = "%s/shmooster.mad" % mad_path
#line_name = "hack"

mad_file = sys.argv[1]
line_name = sys.argv[2]

kinetic_energy = 0.4;
mass = PH_NORM_mp;
energy   = kinetic_energy + mass;
momentum = math.sqrt(energy*energy - mass*mass);

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

pr = Proton( energy )
jp = JetProton(energy)

bml_fact    = bmlfactory(mad_file)
bmline_orig = bml_fact.create_beamline(line_name, pr.ReferenceBRho())
#bmline_orig.flatten()
#bmline_orig.insert(pacifier)
#bmline_orig.append(pacifier)
#bmline = DriftsToSlots(bmline_orig)
#bmline = bmline_orig

#insert_markers(bmline, 20, energy, momentum)

bmln_context = BeamlineContext(0,bmline)
print "isRing =",bmln_context.isRing()
print "isTreatedAsRing =",bmln_context.isTreatedAsRing()
print "modifying..."
bmln_context.handleAsRing()
print "isRing =",bmln_context.isRing()
print "isTreatedAsRing =",bmln_context.isTreatedAsRing()


bit = DeepBeamlineIterator(bmline)

# ### Lattice functions
# bmline.propagateJetParticle(jp)
# lfsage = LattFuncSage(bmline,0)

# foo = lfsage.TuneCalc(jp, 1)
# lfsage.Disp_Calc(jp)

# lfsage.Fast_CS_Calc(jp)

be = bit.reset()
be = bit.next()
s = []
betax = []
betay = []
i = 0
while   be:
###	lf  = be.dataHook.find("Twiss",1)
	print i,be.Name(),be.Type()
	lf = bmln_context.getLattFuncPtr(i)
        if lf:
              s.append(lf.arcLength)
	      betax.append(lf.beta.hor)
	      betay.append(lf.beta.ver)
        else:
            print "No lf for", be.Type()
	be = bit.next()
	i+=1


import pylab


bx = pylab.plot(s,betax,'-ob',linewidth=1.5)
by = pylab.plot(s,betay,'-og',linewidth=1.5)
pylab.legend((bx,by),('beta x','beta y'))
pylab.xlabel('s (m)')
pylab.ylabel('beta (m)')

pylab.show()


print "I guess I'll just hang now..."

