#!/usr/bin/env python

import matplotlib
from matplotlib import pyplot
import numpy
import re
import sys

def get_time(dir):
    f = open(dir + "/synergia.out", 'r')
    times = {}
    propagate_time = None
    for line in f.readlines():
        match = re.match('([^ ]*):([0-9].*)', line)
        if match:
            label = match.group(1)
            time = float(match.group(2))
            if times.has_key(label):
                times[label].append(time)
            else:
                times[label] = [time]
        match = re.match('propagate time = (.*)', line)
        if match:
            propagate_time = float(match.group(1))
    f.close()
    time = {}
    for key in times.keys():
        time[key] = numpy.array(times[key]).sum()
    total = 0.0
    for key in time.keys():
        total += time[key]
    time['total'] = total
    time['propagate'] = propagate_time
    return time

def get_node_proc(dir):
    f = open(dir + "/job", 'r')
    node = None
    proc = None
    for line in f.readlines():
        match = re.match('.*nodes=(.*)', line)
        if match:
            node = int(match.group(1))
        match = re.match('.*-np ([0-9]+) ', line)
        if match:
            proc = int(match.group(1))
    return node, proc

def get_summary(dirs, procs_per_node):
    procs = []
    times = {}
    for dir in dirs:
        node, proc = get_node_proc(dir)
        time = get_time(dir)
        if (procs_per_node == 0) or (node * procs_per_node == proc):
            procs.append(proc)
            for key in time.keys():
                if times.has_key(key):
                    times[key].append(time[key])
                else:
                    times[key] = [time[key]]
    return procs, times

def plot_total(procs8, times8, procs12, times12):
    pyplot.loglog(procs8, times8['total'], 'o-', label='8 cores/node', basex=2)
    pyplot.loglog(procs12, times12['total'], 'o-', label='12 cores/node', basex=2)
    pyplot.legend()
    pyplot.xlabel('cores')
    pyplot.ylabel('time [s]')
    pyplot.axis([6, 64, 8, 15])

def plot_biggest(procs, times, minfrac):
    length = len(procs)
    end = length - 1
    other = numpy.zeros([length],numpy.float64)
    pyplot.loglog(procs, times['total'], 'o-', basex=2, label='total')
    for key in times.keys():
        if (key != 'total') and (key != 'propagate'):
            if times[key][end] / times['total'][end] > minfrac:
                pyplot.loglog(procs, times[key], basex=2, label=key)
            else:
                for i in range(0,length):
                    other[i] += times[key][i]
    pyplot.loglog(procs, other, 'x-', basex=2, label='other')
    pyplot.legend(loc='lower left', ncol=3,
                  prop=matplotlib.font_manager.FontProperties(size='x-small'))
    pyplot.xlabel('cores')
    pyplot.ylabel('time [s]')
    pyplot.axis([6, 48, .1, 15])

def plot_sc(procs, times):
    length = len(procs)
    end = length - 1
    total = numpy.zeros([length],numpy.float64)
    for key in times.keys():
        if re.match('sc-.*', key):
            pyplot.loglog(procs, times[key], basex=2, label=key)
            for i in range(0,length):
                total[i] += times[key][i]
    pyplot.loglog(procs, total, 'o-', basex=2, label='sc total')
    pyplot.legend(loc='lower left', ncol=3,
                  prop=matplotlib.font_manager.FontProperties(size='x-small'))
    pyplot.xlabel('cores')
    pyplot.ylabel('time [s]')
    pyplot.axis([6, 48, .1, 15])

def plot_sc_nonsc(procs, times):
    length = len(procs)
    end = length - 1
    sc = numpy.zeros([length],numpy.float64)
    nonsc = numpy.zeros([length],numpy.float64)
    other = numpy.zeros([length],numpy.float64)
    pyplot.loglog(procs, times['total'], 'o-', basex=2, label='total')
    for key in times.keys():
        if (key != 'total') and (key != 'propagate'):
            if re.match('sc-.*', key):
                for i in range(0,length):
                    sc[i] += times[key][i]
            else:
                for i in range(0,length):
                    nonsc[i] += times[key][i]
    pyplot.loglog(procs, nonsc, basex=2, label='non-sc')
    pyplot.loglog(procs, sc, basex=2, label='sc')
    pyplot.legend(loc='lower left', ncol=3,
                  prop=matplotlib.font_manager.FontProperties(size='x-small'))
    pyplot.xlabel('cores')
    pyplot.ylabel('time [s]')
    pyplot.axis([6, 48, .1, 15])

def plot_all(procs, times):
    index = 0
    pyplot.subplots_adjust(hspace=0.5)
    for key in times.keys():
        index += 1
        pyplot.subplot(5,3,index)
        pyplot.loglog(procs, times[key], basex=2)
        pyplot.axis([6, 48, .01, 15])
        pyplot.title(key,
                  fontsize='x-small')
        ax = pyplot.gca()
        fontsize = 'xx-small'
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)


dirs = sys.argv[1:]
procs8, times8 = get_summary(dirs, 8)
procs12, times12 = get_summary(dirs, 12)

pyplot.figure(1)
plot_total(procs8, times8, procs12, times12)
pyplot.savefig('total.pdf')

pyplot.figure(2)
plot_biggest(procs8, times8, 0.05)
pyplot.savefig('biggest8.pdf')

pyplot.figure(3)
plot_biggest(procs12, times12, 0.05)
pyplot.savefig('biggest8.pdf')

pyplot.figure(4)
plot_sc(procs8, times8)
pyplot.savefig('sc8.pdf')

pyplot.figure(5)
plot_sc(procs12, times12)
pyplot.savefig('sc12.pdf')

pyplot.figure(6)
plot_sc_nonsc(procs8, times8)
pyplot.savefig('scnonsc8.pdf')

pyplot.figure(7)
plot_sc_nonsc(procs12, times12)
pyplot.savefig('scnonsc12.pdf')
pyplot.figure(8)
plot_all(procs8, times8)
pyplot.savefig('all8.pdf')

pyplot.figure(9)
plot_all(procs12, times12)
pyplot.savefig('all12.pdf')

pyplot.show()
