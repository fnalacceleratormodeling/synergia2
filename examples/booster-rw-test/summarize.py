#!/usr/bin/env python

import sys

def addiflonger(hash,key,value):
    if hash.has_key(key):
        if len(hash[key]) < len(value):
            hash[key] = value
    else:
        hash[key] = value

dummyhash = {}
results = {}
exclude = ["summarizedir","jobdir","createjob",
    "trackfraction","saveperiod","submit","track","space_charge",
    "numproc"]
for dir in sys.argv[1:]:
    addiflonger(dummyhash,"summarizedir",dir)
    results[dir] = {}
    f = open("%s/command"%dir)
    items = f.readline().split()
    for item in items:
        pair = item.split('=')
        if len(pair)==2:
            results[dir][pair[0]] = pair[1]
            addiflonger(dummyhash,pair[0],pair[0])
            addiflonger(dummyhash,pair[0],pair[1])
    f.close()

dirformat = "%%%ds" % len(dummyhash["summarizedir"])
print dirformat%"",
for key in dummyhash.keys():
    if exclude.count(key) == 0:
        format = "%%%ds" % len(dummyhash[key])
        print format % key,
print

for dir in sys.argv[1:]:
    print dirformat % dir,
    for key in dummyhash.keys():
        if exclude.count(key) == 0:
            format = "%%%ds" % len(dummyhash[key])
            val = ""
            if results[dir].has_key(key):
                val = results[dir][key]
            print format % val,
    print
