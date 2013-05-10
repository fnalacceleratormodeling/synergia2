#!/usr/bin/env python

import glob
import re
import sys
import os.path

def start_file(filename):
    f = open(filename, 'w')
    s = os.path.basename(filename).split('.')
    header = 'sphinx-module-headers/' + s[0] + '-header.' + s[1]
    if os.path.exists(header):
        hf = open(header,'r')
        for line in hf.readlines():
            f.write(line)
        hf.close()
    else:
        sys.stderr.write(sys.argv[0] + 
                         ":failed to open " + header + "\n")
        sys.exit(1)
    return f

def print_classes(f, classes):
    f.write('''-------
Classes
-------

''')
    classes.sort()
    for class_ in classes:
        f.write('.. doxygenclass:: %s\n' % class_)
# jfa: DANGER DANGER
#      hardwired to workaround a bug in Sphinx by not
#      documenting members in Logger.
        if class_ != 'Logger':
            f.write('    :members:\n')
        f.write('\n')

def print_typedefs(f, typedefs):
    if len(typedefs) == 0:
        return
    f.write('''--------
Typedefs
--------

''')
    typedefs.sort()
    for typedef in typedefs:
        f.write('''.. doxygentypedef:: %s
    :project: synergia

''' % typedef)

def get_classes_typedefs(dir, verbose):
    classes = []
    typedefs = []
    for fname in glob.glob(os.path.join(dir,'*.h')):
        f = open(fname,'r')
        in_multiline_typedef = False
        for line in f.readlines():
            if in_multiline_typedef:
                m2 = re.search('([^ ]+);', line)
                if m2:
                    in_multiline_typedef = False
                    m3 = re.search('; *// *syndoc:include', line)
                    typedef = m2.group(1)
                    if m3:
                        typedefs.append(typedef)
                    else:
                        if verbose:
                            print 'skipping typedef', typedef,'in', fname
            else:
                m = re.search('^[ \t]*class (.*)', line)
                if m and (m.group(1).find(';') == -1):
                    classname = m.group(1)
                    colon_loc = classname.find(':')
                    if colon_loc > 0:
                        classname = classname[0:colon_loc]
                    semicolon_loc = classname.find(';')
                    if semicolon_loc > 0:
                        classname = classname[0:semicolon_loc]
                    classes.append(classname)

                m = re.search('^[ \t]*typedef (.*)', line)
                if m:
                    m2 = re.search('([^ ]+);', m.group(1))
                    if m2:
                        m3 = re.search('; *// *syndoc:include', m.group(1))
                        typedef = m2.group(1)
                        if m3:
                            typedefs.append(typedef)
                        else:
                            if verbose:
                                print 'skipping typedef', typedef,'in', fname
                    else:
                        in_multiline_typedef = True
    return classes, typedefs

def main(argv):
    verbose = False
    if argv[1] == '--verbose':
        verbose = True
        argv = argv[1:]
    if len(argv) != 3:
        sys.stderr.write('usage: ' + sys.argv[0] + ' [--verbose] <dir> <output_filename.rst>\n')
        sys.exit(1)
    classes, typedefs = get_classes_typedefs(argv[1], verbose)
    f = start_file(argv[2])
    print_classes(f, classes)
    print_typedefs(f, typedefs)

if __name__ == '__main__':
    main(sys.argv)