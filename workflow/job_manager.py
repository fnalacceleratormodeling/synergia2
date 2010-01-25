#!/usr/bin/env python

import sys
import os, os.path
import shutil
import time
import string
import commands
import re
import options

job_mgr_opts = options.Options("Job Manager")
job_mgr_opts.add("createjob", 0, "Whether to create new job directory", int)
job_mgr_opts.add("jobdir", "run", "Job directory", str)
job_mgr_opts.add("numproc", 1, "Number of processors", int)
job_mgr_opts.add("submit", 0, "Whether to immediately submit job", int)
job_mgr_opts.add("overwrite", 0, "Whether to overwrite existing job directory", int)
job_mgr_opts.add("walltime", None, "Limit job to given wall time", str)

class Job_manager:
    def __init__(self, script, opts, extra_files=None, argv=sys.argv):
        self.directory = None
        self.real_script = os.path.abspath(script)
        options_file = os.path.splitext(script)[0] + '_options.py'
        if os.path.exists(options_file):
            self.real_options_file = os.path.abspath(options_file)
        else:
            self.real_options_file = None
        self.argv = argv
        if len(self.argv) > 1:
            if self.argv[1] == "--strip-options-file":
                self.argv = self.argv[2:]
        self.synergia_dir = get_synergia_directory()
        self.opts = opts
        self.opts.add_suboptions(job_mgr_opts)
        self.opts.parse_argv(self.argv)
        
        if self.opts.get("createjob"):
            self.create_job(self.opts.get("jobdir"))
            if extra_files:
                self.copy_extra_files(extra_files)
            print "created",
            if opts.get("submit"):
                print "and submitted",
            print "job"
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
 
    def create_script(self,template,name,directory,subs):
        script_templates_dir = get_script_templates_dir()
        default_script_templates_dir = get_default_script_templates_dir()
        template_path = os.path.join(script_templates_dir,template)
        if ((not os.path.exists(template_path)) and 
            (not script_templates_dir == default_script_templates_dir)):
            template_path = os.path.join(default_script_templates_dir,
                                         template)
            print "Note: taking",template,"from default path",default_script_templates_dir
        if not os.path.exists(template_path):
            template_example = template + "_example"
            print "Warning: using", template_example, "for", template, \
                "template."
            print "You should create a template for your system in"
            print os.path.join(script_templates_dir,template)
            template_path = os.path.join(default_script_templates_dir,
                                        template_example)
        output_path = os.path.join(directory,name)
        process_template(template_path,output_path,subs)

    def create_job(self, directory):
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
        subs["numproc"] = self.opts.get("numproc")
        subs["synergia2dir"] = self.synergia_dir
        subs["args"] = self._args_to_string(self.argv[1:], ["createjob"])
        subs["jobdir"] = os.path.abspath(self.directory)
        subs["script"] = os.path.basename(self.real_script)
        job_name = directory + "_job"
        self.create_script("job", job_name, directory, subs)
        self.create_script("cleanup", "cleanup", directory, subs)
        if self.opts.get("submit"):
            os.chdir(directory)
            os.system("qsub %s" % job_name)
            os.chdir(old_cwd)

    def copy_extra_files(self, files):
        for file in files:
            shutil.copy(file, self.directory)

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
    if version > 499:
        print "Sanity check failure: attempt to create directory version %d."\
              % version
        print "Maximum is 499."
        sys.exit(1)
    if version == 0:
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
    if os.environ.has_key('SYNERGIA_SCRIPT_TEMPLATES'):
        return os.environ['SYNERGIA_SCRIPT_TEMPLATES']
    else:
        return os.path.join(get_synergia_directory(),'script-templates')

def get_default_script_templates_dir():
    return os.path.join(get_synergia_directory(),'script-templates')


if __name__ == "__main__":
    subs = {}
    subs["NODES"] = str(72)
    subs["PROCESSES"] = str(144)
    subs["WALLTIME"] = "01:00:00"
    #subs["MAIL"] = None
    process_template("generic_pbs_template.sh", "synergia-pbs.sh", subs)
