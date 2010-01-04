#!/usr/bin/env python

import arrayio

def write_vector(file,t):
    for i in range(0,len(t)):
        file.write(str(t[i]))
        if i<len(t)-1:
            file.write(" ")
    file.write("\n")
h = map(float,open("../fort.51").readline().split())
print h
n = map(int,open("../fort.39").readline().split())
print n
physical_size=[h[0]*(n[0]-1),h[1]*(n[1]-1),h[2]*(n[2]-1)]
print physical_size
physical_offset = [0,0,0]

fort_points = arrayio.readfort("../fort.40",n)

f = open("rho_in","w")
write_vector(f,physical_size)
write_vector(f,physical_offset)
write_vector(f,n)
for i in range(0,n[0]):
    for j in range(0,n[1]):
        for k in range(0,n[2]):
            f.write("%g " % fort_points[k,j,i])
f.write("\n")
f.close()
