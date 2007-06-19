#!/usr/bin/env python

import sys
import string

def doit(dirname):
    f = open(dirname +"/command")
    parts = f.readline().split()
    for part in parts:
        subparts = part.split("=")
        if subparts[0] == "numproc":
            numproc = int(subparts[1])
    print numproc
    f.close()
    
    etime = 0
    results = {}
    f = open(dirname + "/synergia.out")
    for line in f.readlines():
        parts = line.split(" ")
        if parts[0][0:4] == "bcel":
            line = string.join(parts[1:])
        parts = line.split(":")
        if len(parts) == 2:
            time = float(parts[1])
            if results.has_key(parts[0]):
                results[parts[0]] += time
            else:
                results[parts[0]] = time
        else:
            eparts = line.split("=")
            if eparts[0] == "elapsed time ":
                etime = float(eparts[1])
    #~ print results
    #~ print etime
    return numproc, etime, results

keys = None
bigresults = {}
etimes = {}
for dirname in sys.argv[2:]:
    ok = 1
    try:
        numproc,etime,results = doit(dirname)
    except:
        ok = 0
        numproc = 0
        etime = {}
        results = {}
    if keys:
        tmp = results.keys()
        tmp.sort()
        if tmp != keys:
            print "keys do not match:",dirname
            print keys
            print tmp
            ok = 0
    else:
        keys = results.keys()
        keys.sort()
    if ok:
        bigresults[numproc] = results
        etimes[numproc] = etime
f = open(sys.argv[1] + ".dat", "w")
nps = bigresults.keys()
nps.sort()
print nps
for np in nps:
    f.write("%d %g " % (np,etimes[np]))
    for key in keys:
        f.write("%g " % bigresults[np][key])
    f.write("\n")
f.close()

f = open(sys.argv[1]  + "-key.dat","w")
for key in keys:
    f.write("%s\n" % key)
f.close()
