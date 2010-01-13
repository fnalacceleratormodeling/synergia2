#!/usr/bin/env python

from boost_python_test_options import opts
import bp_hello

if __name__ == "__main__":
    print "about to call C++ function:"
    retval = bp_hello.hello(opts.get("x"))
    print "function returned",retval
    print "boost_python_test.py finished successfully"
