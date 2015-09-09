#!/usr/bin/env python
import sys
import os

# memory related contents of /proc/self/status from http://locklessinc.com/articles/memory_usage/
# VmPeak:	Peak virtual memory usage
# VmSize:	Current virtual memory usage
# VmLck:	Current mlocked memory
# VmHWM:	Peak resident set size
# VmRSS:	Resident set size
# VmData:	Size of "data" segment
# VmStk:	Size of stack
# VmExe:	Size of "text" segment
# VmLib:	Shared library usage
# VmPTE:	Pagetable entries size
# VmSwap:	Swap space used

# VmSwap not present in SL5

vm_names = [
    "VmPeak:", "VmSize:", "VmLck:", "VmHWM:", "VmRSS:", "VmData:",
    "VmStk:", "VmExe:", "VmLib:", "VmPTE:", "VmSwap:"]

# convert units to MB
vm_units = {
    "kB": 1.0/1024.0, "KB": 1.0/1024.0, "mB": 1.0, "MB": 1.0 }

# memusage() return a dictionary containing the current memory usage of the process
def memusage():
    f = open("/proc/self/status")
    statuslines = f.readlines()
    f.close()

    memdict = {}
    
    for sl in statuslines:
        wds = sl.split(None, 3)
        if len(wds) != 3:
            pass # not right format
        if wds[0] in vm_names:
            # this is a vm line
            memsize = float(wds[1])*vm_units[wds[2]]
            memdict[wds[0]] = memsize
    return memdict

def memusage_as_string():
    vmu = memusage()
    return"VmSize: %g, VmRSS: %g"%(vmu["VmSize:"],vmu["VmRSS:"])

def print_memusage(ostream=sys.stdout):
    vmu = memusage()
    ostream.write(memusage_as_string()+"\n")


