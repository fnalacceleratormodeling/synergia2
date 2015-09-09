#!/usr/bin/env python

import numpy
import re
import sys

def extract_times(filename):
    f = open(filename, 'r')
    times_detail = {}
    propagate_time = None
    for line in f.readlines():
        match = re.match('simple_timer:(.*):([0-9].*)', line)
        if match:
            label = match.group(1)
            time = float(match.group(2))
            if times_detail.has_key(label):
                times_detail[label].append(time)
            else:
                times_detail[label] = [time]
        match = re.match('propagate time = (.*)', line)
        if match:
            propagate_time = float(match.group(1))
    f.close()
    times = {}
    for key in times_detail.keys():
        n = len(times_detail[key])
        times[key] = numpy.array(times_detail[key]).sum()
        times[key+"_std"] = numpy.array(times_detail[key]).std() * numpy.sqrt(n)
    total = 0.0
    for key in times.keys():
        total += times[key]
    times['total'] = total
    if propagate_time:
        times['propagate'] = propagate_time
    return times_detail, times

def write_summary_file(file, times, csv):
    keys = times.keys()
    keys.sort()
    keys.remove('total')
    keys.append('total')
    if csv:
        for key in keys:
            file.write(key + ',' + str(times[key]) + '\n')
    else:
        maxlen = 0
        for key in keys:
            if len(key) > maxlen:
                maxlen = len(key)
        format = "%%%ds  %%g\n" % maxlen
        for key in keys:
            file.write(format % (key, times[key]))

def write_detail_file(file, times, csv):
    keys = times.keys()
    keys.sort()
    if csv:
        for key in keys:
            file.write(key)
            for time in times[key]:
                file.write(',' + str(time))
            file.write('\n')
    else:
        maxlen = 0
        for key in keys:
            if len(key) > maxlen:
                maxlen = len(key)
        format = "%%%ds " % maxlen
        for key in keys:
            file.write(format % key)
            for time in times[key]:
                file.write(' ')
                file.write(str(time))
            file.write('\n')

def script_usage(retval):
    print "usage:"
    print "simple_timer.py [--csv] filename [summary_filename] [detail_filename]"
    sys.exit(retval)

if __name__ == "__main__":
    args = []
    csv = False
    for arg in sys.argv[1:]:
        if arg == '--help':
            script_usage(0)
        if arg == '--csv':
            csv = True
        else:
            args.append(arg)
    if len(args) < 1:
        script_usage(1)
    filename = args[0]

    if len(args) > 1:
        summary_filename = args[1]
    else:
        summary_filename = None

    if len(args) > 2:
        detail_filename = args[2]
    else:
        detail_filename = None

    times_detail, times = extract_times(filename)

    if summary_filename:
        summary_file = open(summary_filename, 'w')
    else:
        summary_file = sys.stdout
    write_summary_file(summary_file, times, csv)
    if summary_filename:
        summary_file.close()

    if detail_filename:
        detail_file = open(detail_filename, 'w')
        write_detail_file(detail_file, times_detail, csv)
        detail_file.close()

