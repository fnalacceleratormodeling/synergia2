import GaussSC
import numpy


foo = numpy.array([[1.e-7,0.0,5.e-9,-0.001,0.0,0.000]])
print "Send particles =",foo

numpart = 1

sigmaX = 3.6706e-05
sigmaY = 4.2563e-06 
sigmaZ = 0.1

gamma = 1000.
LengthScale = 1.0
PartPersigmaZ = 2.e10

tau=1.

iret = GaussSC.apply_BasErs_kick(foo, numpart, sigmaX, sigmaY, sigmaZ, gamma, tau, PartPersigmaZ, LengthScale)

print "GaussSC returns",iret
print "foo particles after step", foo

print "success is questionable"
