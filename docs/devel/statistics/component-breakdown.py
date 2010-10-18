#!/usr/bin/python

from pylab import *
from commands import getoutput

def count(suffix, test):
    if test:
        yesorno = ""
    else:
        yesorno = " -v "
    out = getoutput("find src -name \*." + suffix +
                    " | grep " + yesorno +
                    " tests| grep -v eigen2 | xargs wc -l | tail -n 1").split()
    print suffix, "( test =", test, ")", out
    return int(out[0])

# make a square figure and axes
figure(1, figsize=(6, 6))
ax = axes([0.1, 0.1, 0.8, 0.8])

cc = count("cc", False)
testcc = count("cc", True)
h = count("h", False)
testh = count("h", True)
py = count("py", False)
testpy = count("py", True)

labels = ["C++ source\n" + str(cc + h), "C++ test\n" + str(testcc + testh),
           "Python source\n" + str(py), "Python test\n" + str(testpy)]
fracs = [cc + h, testcc + testh, py, testpy]
explode = (0.0, 0.05, 0.05, 0.05)
pie(fracs, labels=labels, shadow=True, explode=explode)
title("Total lines: " + str(cc + h + testcc + testh + py + testpy))
savefig("language_breakdown.png", dpi=100)
show()

figure(1, figsize=(6, 6))
ax = axes([0.1, 0.1, 0.8, 0.8])
labels = ["Source", "Test"]
fracs = [cc + h + py, testcc + testh + testpy]
explode = (0, 0.05)
pie(fracs, labels=labels, autopct='%1.1f%%', shadow=True, explode=explode)
savefig("test_breakdown.png", dpi=100)
show()
