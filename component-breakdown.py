#!/usr/bin/python

from pylab import *
from commands import getoutput

def count(suffix, test):
    if test:
        yesorno = ""
    else:
        yesorno = " -v "
    out1 = getoutput("find components -name \*." + suffix +
                    " | grep " + yesorno +
                    " tests| xargs wc -l | tail -n 1").split()
    out2 = getoutput("find utils -name \*." + suffix +
                     "| grep -v eigen2 " +
                     " | grep " + yesorno +
                     " tests| xargs wc -l | tail -n 1").split()
    print suffix, "(test =", test, ")", out1, out2
    return int(out1[0]) + int(out2[0])

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
