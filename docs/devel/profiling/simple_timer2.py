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
        times[key] = numpy.array(times_detail[key]).sum()
    if propagate_time:
        times['propagate'] = propagate_time
    return times_detail, times

def recursive_print(file, format, level0, level0s, items, times):
    sum = 0.0
    maxlen = 0
    nexts = []
    total = level0 + '-total'
    for key in items[level0]:
        if key != total:
            file.write(format % (key, times[key]) + '\n')
            sum += times[key]
            split = key.split('-')
            if len(split) > 1:
                next = split[1]
                if not next in nexts:
                    nexts.append(next)
    if total in items[level0]:
        file.write(format % (total, times[total]) + '\n')
    if len(items[level0])>1:
        file.write(format % ('** ' + level0 + ' sum', sum) + ' **\n')
    file.write('\n')
    level0s.remove(level0)
    for next in nexts:
        if next in level0s:
            recursive_print(file, format, next, level0s, items, times)
    if len(level0s) > 0:
        next = level0s[0]
    else:
        next = None
    return next
        
def write_summary_file(file, times, csv):
    keys = times.keys()
    keys.sort()
    level0s = []
    items = {}
    misc = 'misc'
    items[misc] = []
    for key in keys:
        split = key.split('-')
        if len(split) > 1:
            level0 = key.split('-')[0]
            if not level0 in level0s:
                level0s.append(level0)
                items[level0] = [key]
            else:
                items[level0].append(key)
        else:
            items[misc].append(key)
    maxlen = 0
    for key in keys:
        if len(key) > maxlen:
            maxlen = len(key)
    format = "%%%ds  %%12.5f" % maxlen
    next = 'propagate'
    while len(level0s) > 0:
        next = recursive_print(file, format, next, level0s, items, times)
    recursive_print(file, format, misc, [misc], items, times)

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

