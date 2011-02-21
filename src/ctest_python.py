#!/usr/bin/env python

import os
import re
import string
import time
import commands
import sys
import shutil

def find_ctesttestfiles(dir):
    matches = []
    for root, dirnames, filenames in os.walk(dir):
        for filename in filenames:
            if filename == 'CTestTestfile.cmake':
                matches.append(os.path.join(root, filename))
    return matches

def extract_command_args(concat_command):
    retval = []
    in_string = False
    string = ''
    for ch in concat_command:
        if ch == '"':
            if in_string:
                retval.append(string)
                string = ''
                in_string = False
            else:
                in_string = True
        else:
            if in_string:
                string += ch
    return retval

class Test:
    def __init__(self, name, command, dir):
        self.name = name
        self.command = "PATH=.:$PATH; export PATH;" + command
        self.dir = dir

    def run(self, index, max_index):
        print "      Start:", self.name

        here = os.getcwd()
        os.chdir(self.dir)
        t0 = time.time()
        status, output = commands.getstatusoutput(self.command)
        t1 = time.time()
        os.chdir(here)
        beginning = "%2d/%2d" % (index, max_index) + ' '
        beginning += "Test #%d:" % index + ' '
        beginning += self.name
        end = ''
        if status:
            end += "***Failed   "
        else:
            end += "   Passed   "
        end += "%0.2f sec" % (t1 - t0)
        dots = '...................................................................'
        num_dots = 78 - len(beginning) - len(end)
        print "%s %s %s" % (beginning, dots[0:num_dots], end)
        return status, output

def extract_tests(ctesttestfile):
    f = open(ctesttestfile, 'r')
    tests = []
    names = []
    commands = []
    for line in f.readlines():
        match = re.match('ADD_TEST\((.*)\)', line)
        if match:
            split = match.group(1).split(' ')
            name = split[0]
            concat_command = string.join(split[1:])
            command = extract_command_args(concat_command)
            tests.append(Test(name, concat_command, os.path.dirname(ctesttestfile)))
    return tests

if __name__ == '__main__':
    cttfiles = find_ctesttestfiles('.')
    all_tests = []
    for cttfile in cttfiles:
        tests = extract_tests(cttfile)
        all_tests.extend(tests)
    logdir = os.path.abspath('./ctest_python_logs')
    if os.path.exists(logdir):
        shutil.rmtree(logdir)
    os.mkdir(logdir)
    count = 0
    errors = []
    error_names = []
    t0 = time.time()
    for test in all_tests:
        count += 1
        status, output = test.run(count, len(all_tests))
        if status:
            errors.append(count)
            error_names.append(test.name)
        open(os.path.join(logdir, '%d.log' % count), 'w').write(output + '\n')
    t1 = time.time()

    print
    print "%d%%" % (100 * (1.0 - len(errors) / (1.0 * len(all_tests)))),
    print "tests passed,",
    print len(errors), "tests failed out of", len(all_tests)
    print
    print "Total Test time (real) = %0.2f" % (t1 - t0) , "sec"
    print
    if len(errors) > 0:
        print "The following tests FAILED:"
        for error, error_name in zip(errors, error_names):
            print "%8d" % error, '-', error_name, '(Failed)'
        print 'Errors while running ctest_python.py'
        sys.exit(1)
