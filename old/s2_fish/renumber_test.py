#!/usr/bin/env python

import sys,os,shutil
import re

filename = sys.argv[1]
shutil.copy(filename,filename+".bak")
tmpfilename = filename+".tmp"
outfile = open(tmpfilename,"w")
infile = open(filename,"r")
test_counter = None
for line in infile.readlines():
    outline = line
    if re.search("class.*TestCase\):",line):
        test_counter = 0
    if re.search("def test_[0-9]+_",line):
        test_counter += 1
        outline = re.sub("def test_[0-9]+_","def test_%02d_" % test_counter,
                         line)
        print outline.rstrip()
    outfile.write(outline)
infile.close()
outfile.close()
os.rename(tmpfilename,filename)


