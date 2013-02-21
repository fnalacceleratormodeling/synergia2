#!/usr/bin/env python

import sys, os.path
from synergia.utils import Commxx
from synergia.simulation import Propagator, Resume

class Options:
    def __init__(self):
        self.unspecified_str = ""
        self.unspecified_int = -1
        self.directory = Propagator.default_checkpoint_dir
        self.new_checkpoint_directory = self.unspecified_str
        self.checkpoint_period = self.unspecified_int
        self.max_turns = self.unspecified_int
        self.verbosity = self.unspecified_int

def do_error(message):
    sys.stderr.write(message + '\n')
    sys.exit(1)

def do_help():
    print "usage: synergia-pyresume [options] [checkpoint directory]";
    print "  options:";
    print "    --help: this message";
    print "    --new-dir=<dir>: directory name to use for subsequent checkpointing";
    print "    --period=<period>: period for subsequent checkpointing";
    print "    --max=<num>: maximum number of turns for this run";
    print "    --verbosity=<num>: verbosity for this run";
    sys.exit(0)

def handle_args(args):
    options = Options()
    for arg in args:
        if arg[0] == '-':
            if arg == '--help':
                do_help()
            elif arg.find('--new-dir') == 0:
                options.new_checkpoint_directory = arg.split('=')[1]
            elif arg.find('--period') == 0:
                options.checkpoint_period = int(arg.split('=')[1])
            elif arg.find('--max') == 0:
                options.max_turns = int(arg.split('=')[1])
            elif arg.find('--verbosity') == 0:
                options.verbosity = int(arg.split('=')[1])
            else:
                do_error('Unknown argument "%s"' % arg)
        else:
            options.directory = arg
    return options

def check_size(options):
    f = open(os.path.join(options.directory,
                          Propagator.description_file_name), 'r')
    for line in f.readlines():
        s = line.split('=')
        if s[0] == 'mpi_size':
            original_size = int(s[1])
    current_size = Commxx().get_size()
    if (original_size != current_size):
        sys.stderr.write('''synergia-pyresume: Error. Number of MPI processes (%d) must be equal to
                   number of MPI processes in original job (%d).\n''' % \
            (current_size, original_size))
        sys.exit(1)

def run(options):
    check_size(options)
    resume = Resume(options.directory)
    if options.checkpoint_period != options.unspecified_int:
        resume.set_checkpoint_period(options.checkpoint_period)
    if options.new_checkpoint_directory != options.unspecified_str:
        resume.set_new_checkpoint_dir(options.new_checkpoint_directory)
    new_max_turns = (options.max_turns != options.unspecified_int)
    new_verbosity = (options.verbosity != options.unspecified_int)
    resume.propagate(new_max_turns, options.max_turns, new_verbosity,
                     options.verbosity)

if __name__ == '__main__':
    options = handle_args(sys.argv[1:])
    run(options)
