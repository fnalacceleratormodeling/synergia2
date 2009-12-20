#!/usr/bin/env python

from python_test_options import opts

if __name__ == "__main__":
    print "hello world from python_test.py"
    print "the value of x is",opts.get("x")
