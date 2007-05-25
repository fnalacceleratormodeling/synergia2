#!/usr/bin/env python

import time

t0 = 0.0

def mytimer(message):
    global t0
    if t0 == 0.0:
        t0 = time.time()
    t1 = time.time()
    print "%s: %g"%(message,t1-t0)
    t0 = time.time()
