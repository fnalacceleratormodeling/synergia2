#!/usr/bin/env python

import sys
import os, os.path
import shutil
import time
import string
import commands
import re
import options

job_mgr_opts = options.Options("job_manager")
job_mgr_opts.add("createjob", False, "Whether to create new job directory", bool)
#job_mgr_opts.add("resumejob", 0, "Whether to resume a previously checkpointed job", int)
job_mgr_opts.add("jobdir", "run", "Job directory", str)
#job_mgr_opts.add("resumedir", "run", "Directory containing checkpointed files", str)
job_mgr_opts.add("numproc", 1, "Number of processes", int)
job_mgr_opts.add("procspernode", 1, "Number of processes per node", int)
job_mgr_opts.add("submit", False, "Whether to immediately submit job", bool)
job_mgr_opts.add("run", False, "Whether to immediately run job", bool)
job_mgr_opts.add("overwrite", False, "Whether to overwrite existing job directory", bool)
#job_mgr_opts.add("walltime", "00:30:00", "Limit job to given wall time", str)
job_mgr_opts.add("synergia_executable", "synergia", "Name or path of synergia executable", str)

#job_mgr_opts.add("checkpoint", 0, "enable generation of checkpoint", int)
#job_mgr_opts.add("checkpoint_freq", 100, "frequency of checkpoint generation", int)

def get_synergia_directory(die_on_failure=1):
    if os.environ.has_key("SYNERGIA2DIR"):
        synergia_dir = os.environ["SYNERGIA2DIR"]
    else:
        if os.path.isfile("./local_paths.py"):
            synergia_dir = "."
        else:
            if os.path.isfile("../local_paths.py"):
                synergia_dir = ".."
            else:
                print "Unable to determine Synergia directory. Defaults"
                print "are \".\" and \"..\". Set the environment variable"
                print "SYNERGIA2DIR to point to the directory containing"
                print "local_paths.py."
                if die_on_failure:
                    sys.exit(1)
                else:
                    synergia_dir = None
    if synergia_dir:
        synergia_dir = os.path.abspath(synergia_dir)
    return synergia_dir

def get_script_templates_dir():
    if os.environ.has_key('SYNERGIA2TEMPLATES'):
        return os.environ['SYNERGIA2TEMPLATES']
    else:
        return os.path.join(get_synergia_directory(), 'synergia-script-templates')

def get_default_script_templates_dir():
    return os.path.join(get_synergia_directory(), 'synergia-script-templates')

def add_local_opts():
    found_local_options = False
    local_options_location = None
    if os.path.exists(os.path.join(get_script_templates_dir(), 'local_opts.py')):
        found_local_options = True
        local_opts_location = get_script_templates_dir()
    else:
        if os.path.exists(os.path.join(get_default_script_templates_dir(), 'local_opts.py')):
            found_local_options = True
            local_opts_location = get_default_script_templates_dir()
            if get_script_templates_dir() != get_default_script_templates_dir():
                print "Note: using local_opts.py from", get_default_script_templates_dir()
    if found_local_options:
        sys.path.insert(0, local_opts_location)
        import local_opts
        if hasattr(local_opts, 'opts'):
            job_mgr_opts.add_suboptions(local_opts.opts)
        else:
            print 'Warning: local_opts.py found in'
            print "%s," % local_opts_location
            print "but no opts object was found there."

class Job_manager:
    def __init__(self, script, opts, extra_files=None, extra_dirs=None,
                 extra_opt_files=None, extra_opt_dirs=["lattice_cache"],
                 standalone=False, argv=sys.argv):
        add_local_opts()
        self.directory = None
        self.real_script = os.path.abspath(script)
        options_file = os.path.splitext(script)[0] + '_options.py'
        if os.path.exists(options_file):
            self.real_options_file = os.path.abspath(options_file)
        else:
            self.real_options_file = None
        self.argv = argv
        self.standalone = standalone
        if self.standalone:
            job_mgr_opts.add("executable", script, "Name or path of standalone executable", str)
        self.opts = opts
        self.opts.add_suboptions(job_mgr_opts)
        if len(self.argv) > 2:
            if (self.argv[2] == "--create-cxx-options-source"):
                self.create_cxx_options_source()
                sys.exit(0)
        if len(self.argv) > 1:
            if (self.argv[1] == "--create-cxx-options-source"):
                self.create_cxx_options_source()
                sys.exit(0)
            if self.argv[1] == "--strip-options-file":
                self.argv = self.argv[2:]
        self.synergia_dir = get_synergia_directory()

        self.opts.parse_argv(self.argv)

        # if we are resuming a checkpointed job, get the absolute
        # path of the resumedir to pass to the resumed job.  I don't check
        # whether the resumedir exists because I might be submitting a
        # series of dependent jobs to run one after another.
        if self.opts.get("resumejob"):
            real_resumedir = os.path.abspath(self.opts.get("resumedir"))
            # edit the absolute resumedir path into the resumedir= argument
            foundresumedir = False
            for argidx in range(len(self.argv)):
                splitarg = string.split(self.argv[argidx], "=")
                if len(splitarg) > 1:
                    if splitarg[0] == "resumedir":
                        foundresumedir = True
                        self.argv[argidx] = "resumedir=" + real_resumedir
            # if I've made it all the way through the list without finding
            # a resumedir argument, add it now
            if not foundresumedir:
                self.argv.append("resumedir=" + real_resumedir)

        if self.opts.get("createjob"):
            self.create_job(self.opts.get("jobdir"), extra_files, extra_dirs,
                            extra_opt_files, extra_opt_dirs)
            if not opts.job_manager.run:
                print "created",
                if opts.get("submit"):
                    print "and submitted",
                print "job in directory"
            sys.exit(0)

    def _args_to_string(self, args, strip=[None]):
        retval = ""
        for arg in args:
	    argout = arg
	    splitarg = string.split(arg, "=")
	    if len(splitarg) > 1:
                if splitarg[0] in strip:
                    argout = None
		elif(string.count(splitarg[1], " ")) > 0:
		    argout = splitarg[0] + '="'
		    splitarg.pop(0)
		    argout = argout + string.join(splitarg, "=") + '"'
            if argout:
                if retval != "":
                    retval += " "
                retval += argout
        return retval

    def create_script(self, template, name, directory, subs):
        script_templates_dir = get_script_templates_dir()
        default_script_templates_dir = get_default_script_templates_dir()
        template_path = os.path.join(script_templates_dir, template)
        if ((not os.path.exists(template_path)) and
            (not script_templates_dir == default_script_templates_dir)):
            template_path = os.path.join(default_script_templates_dir,
                                         template)
            print "Note: taking", template, "from default path", default_script_templates_dir
        if not os.path.exists(template_path):
            template_example = template + "_example"
            print "Warning: using", template_example, "for", template, \
                "template."
            print "You should create a template for your system in"
            print os.path.join(script_templates_dir, template)
            template_path = os.path.join(default_script_templates_dir,
                                        template_example)
        output_path = os.path.join(directory, name)
        process_template(template_path, output_path, subs)

    def create_job(self, directory, extra_files, extra_dirs,
                   extra_opt_files, extra_opt_dirs):
###        real_script = os.path.abspath(self.argv[0])
        old_cwd = os.getcwd()
        overwrite = self.opts.get("overwrite")
        directory = create_new_directory(directory, 0, overwrite)
        self.directory = directory
        os.chdir(directory)
        shutil.copy(self.real_script, ".")
        if self.real_options_file:
            shutil.copy(self.real_options_file, ".")
        commandfile = open("command", "w")
        commandfile.write("%s\n" % self._args_to_string(self.argv))
        commandfile.close()
        os.chdir(old_cwd)
        subs = {}
        for sub in job_mgr_opts.options():
            subs[sub] = job_mgr_opts.get(sub)
        subs["numproc"] = self.opts.get("numproc")
        numnode = (self.opts.get("numproc") + self.opts.get("procspernode") - 1) / \
                  self.opts.get("procspernode")
        subs["procspernode"] = self.opts.get("procspernode")
        subs["numnode"] = numnode
        subs["synergia2dir"] = self.synergia_dir
        subs["args"] = self._args_to_string(self.argv[1:], ["createjob"])
        subs["jobdir"] = os.path.abspath(self.directory)
        if self.standalone:
            subs["synergia_executable"] = self.real_script
            subs["script"] = ""
        else:
            subs["script"] = os.path.basename(self.real_script)
        job_name = os.path.basename(directory) + "_job"
        self.create_script("job", job_name, directory, subs)
        self.create_script("cleanup", "cleanup", directory, subs)
        if extra_files:
            self.copy_extra_files(extra_files)
        if extra_dirs:
            self.copy_extra_dirs(extra_dirs)
        if extra_opt_files:
            self.copy_extra_files(extra_opt_files, optional=True)
        if extra_opt_dirs:
            self.copy_extra_dirs(extra_opt_dirs, optional=True)
        if self.opts.get("submit"):
            os.chdir(directory)
            os.system("qsub %s" % job_name)
            os.chdir(old_cwd)
        if self.opts.get("run"):
            os.chdir(directory)
            os.system("./" + job_name)
            os.chdir(old_cwd)

    def copy_extra_files(self, files, optional=False):
        for file in files:
            if os.path.exists(file):
                shutil.copy(file, self.directory)
            else:
                if not optional:
                    raise RuntimeError("Job_manager: required file " + file + " not found")

    def copy_extra_dirs(self, dirs, optional=False):
        for dir in dirs:
            if os.path.exists(dir):
                shutil.copytree(dir, self.directory)
            else:
                if not optional:
                    raise RuntimeError("Job_manager: required directory " + dir + " not found")

    def create_cxx_options_source(self):
        filename_base = self.real_script + "_options"
        classname = os.path.basename(self.real_script) + "_options"
        classname = classname[0].upper() + classname[1:]
        includeguard = classname.upper() + "_H_"
        disclaimer = '// this file was automatically generated by the command\n' + \
            '//     synergia ' + sys.argv[0] + ' ' + sys.argv[-1] + '\n' + \
            '// DO NOT EDIT\n'

        header_filename = filename_base + ".h"
        header = open(header_filename, "w")
        header.write('#ifndef ' + includeguard + '\n')
        header.write('#define ' + includeguard + '\n')
        header.write('#include <string>\n')
        header.write('\n')
        header.write(disclaimer)
        header.write('\n')
        header.write('struct ' + classname + '\n')
        header.write('{\n')
        header.write('    ' + classname + '(int argc, char **argv);\n')
        for optname in self.opts.dict.keys():
            opt = self.opts.dict[optname]
            header.write('    ' + cxx_typename(opt.val_type) + ' ' + optname + ';\n')
        header.write('};\n')
        header.write('#endif /* ' + includeguard + ' */\n')
        header.close()

        source = open(filename_base + ".cc", "w")
        source.write('#include "' + os.path.basename(header_filename) + '"\n')
        source.write('#include "synergia/utils/command_line_arg.h"\n')
        source.write('\n')
        source.write(disclaimer)
        source.write('\n')
        source.write(classname + '::' + classname + '(int argc, char **argv) :\n')
        count = 0
        for optname in self.opts.dict.keys():
            if (count > 0):
                source.write(',\n')
            source.write('    ' + optname + '(' + \
                         cxx_source_value(self.opts.get(optname)) + ')')
            count += 1
        source.write('\n')
        source.write('{\n')
        source.write('    for (int i = 1; i < argc; ++i) {\n')
        source.write('        Command_line_arg arg(argv[i]);\n')
        source.write('        if (arg.is_equal_pair()) {\n')
        count = 0
        for optname in self.opts.dict.keys():
            opt = self.opts.dict[optname]
            if (count > 0):
                source.write(' else ')
            else:
                source.write('            ')
            source.write('if (arg.get_lhs() == "' + optname + '") {\n')
            source.write('                ' + optname + ' = arg.extract_value<' +\
                         cxx_typename(opt.val_type) + ' >();\n')
            source.write('            }')
            count += 1
        source.write(''' else if (arg.get_lhs() == "synergia_executable") {
                // ignore
            } else if (arg.get_lhs() == "run") {
                // ignore
            } else if (arg.get_lhs() == "submit") {
                // ignore
            } else if (arg.get_lhs() == "overwrite") {
                // ignore
            } else if (arg.get_lhs() == "numproc") {
                // ignore
            } else if (arg.get_lhs() == "procspernode") {
                // ignore
            } else {
                throw std::runtime_error("Unknown argument " + arg.get_lhs());
            }
        } else {
            throw std::runtime_error("Bad argument " + std::string(argv[i]));
        }
    }
}
''')
        source.close()

def cxx_typename(type):
    if type == int:
        retval = "int"
    elif type == float:
        retval = "double"
    elif type == bool:
        retval = "bool"
    elif type == str:
        retval = "std::string"
    else:
        raise RuntimeError("Unable to map type " + str(type) + " to C++ type")
    return retval

def cxx_source_value(val):
    retval = str(val)
    if type(val) == bool:
        if val:
            retval = "true"
        else:
            retval = "false"
    if type(val) == str:
        retval = '"' + val + '"'
    return retval

def process_template(template_name, output_name, subs):
    template = open(template_name, "r")
    output = open(output_name, "w")
    unknown_vars = []
    for line in template.readlines():
        match = re.search("@@[A-z0-9]+@@", line)
        while match:
            var = string.replace(match.group(), "@@", "")
            if subs.has_key(var):
                replacement = str(subs[var])
            else:
                replacement = ""
                unknown_vars.append(var)
            original = match.group()
            line = string.replace(line, original, replacement)
            match = re.search("@@[A-z0-9]+@@", line)
        match = re.search("__([A-z0-9]+){{(.*)}}{{(.*)}}__", line)
        while match:
            var = match.group(1)
            if subs.has_key(var):
                replacement = match.group(2)
            else:
                replacement = match.group(3)
                if unknown_vars.count(var) > 0:
                    unknown_vars.remove(var)
            original = "__%s{{%s}}{{%s}}__" % (match.group(1), match.group(2),
                                               match.group(3))
            line = string.replace(line, original, replacement)
            match = re.search("@@[A-z0-9]+@@", line)
        output.write(line)
    for var in unknown_vars:
        print "process_template warning: variable \"%s\" unkown." % var
    template.close()
    output.close()
    os.chmod(output_name, 0755)

def create_new_directory(directory, version, overwrite):
    if version > 9999:
        print "Sanity check failure: attempt to create directory version %d."\
              % version
        print "Maximum is 9999."
        sys.exit(1)
    if overwrite:
        created_directory = directory
    else:
        created_directory = "%s.%02d" % (directory, version)
    if os.path.isdir(created_directory):
        if overwrite:
            shutil.rmtree(created_directory)
            os.mkdir(created_directory)
            print "created directory", created_directory
        else:
            created_directory = create_new_directory(directory,
                                                     version + 1,
                                                     overwrite)
    else:
        os.mkdir(created_directory)
        print "created directory", created_directory
    return created_directory

if __name__ == "__main__":
    subs = {}
    subs["NODES"] = str(72)
    subs["PROCESSES"] = str(144)
    subs["WALLTIME"] = "01:00:00"
    #subs["MAIL"] = None
    process_template("generic_pbs_template.sh", "synergia-pbs.sh", subs)
