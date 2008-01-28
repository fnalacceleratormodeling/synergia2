#!/usr/bin/env python

import sys
import os, os.path
import shutil
import time
import string
import commands
import re
import options

opts = options.Options("Job Manager")
opts.add("createjob",0,"Whether to create new job directory",int)
opts.add("jobdir","run","Job directory",str)
opts.add("numproc",1,"(single!!!) Number of processors",int)
opts.add("submit",0,"Whether to immediately submit job",int)
###options.add("block",0,"Whether to block until job has finished",int)
###options.add("polltime",60.0,"Time delay in seconds between polls of queue when blocking"
###            ,float)
opts.add("overwrite",0,"Whether to overwrite existing job directory",int)
opts.add("walltime",None,"Limit job to given wall time",str)
###options.add("remote",0,"Prepare job for remote host",int)


class Job_manager:
    def __init__(self, arguments, opts, extra_files = None):
        self.directory = None
        self.arguments = arguments
        self.synergia_dir = get_synergia_directory()
        self.opts = opts
        if self.opts.get("createjob"):
            self.create_job(self.opts.get("jobdir"))
            if extra_files:
                self.copy_extra_files(extra_files)
            print "created",
            if opts.get("submit"):
                print "and submitted",
            print "job"
            sys.exit(0)

    def _args_to_string(self,args,strip=[None]):
        retval = ""
        for arg in args:
	    argout = arg
	    splitarg = string.split(arg,"=")
	    if len(splitarg)>1:
                if splitarg[0] in strip:
                    argout = None
		elif(string.count(splitarg[1]," "))>0:
		    argout = splitarg[0] + '="'
		    splitarg.pop(0)
		    argout = argout + string.join(splitarg,"=") + '"'
            if argout:
                if retval != "":
                    retval += " "
                retval += argout
        return retval
 
    def create_script(self,template,name,directory,subs):
        template_path = os.path.join(self.synergia_dir,"script-templates",
                                     template)
        output_path = os.path.join(directory,name)
        process_template(template_path,output_path,subs)

    def create_job(self, directory):
        real_script = os.path.abspath(self.arguments[0])
        old_cwd = os.getcwd()
        overwrite = self.opts.get("overwrite")
        directory = create_new_directory(directory,0,overwrite)
        self.directory = directory 
        os.chdir(directory)
        shutil.copy(real_script,".")
        commandfile = open("command","w")
        commandfile.write("%s\n" % self._args_to_string(self.arguments))
        commandfile.close()
        os.chdir(old_cwd)
        subs = {}
        subs["numproc"] = opts.get("numproc")
        subs["synergia2dir"] = self.synergia_dir
        subs["args"] = self._args_to_string(self.arguments[1:],["createjob"])
        subs["jobdir"] = os.path.abspath(self.directory)
        subs["script"] = os.path.basename(real_script)
        job_name = directory + "_job"
        self.create_script("job",job_name,directory,subs)
        self.create_script("cleanup","cleanup",directory,subs)
        if self.opts.get("submit"):
            os.chdir(directory)
            os.system("qsub %s" % job_name)
            os.chdir(old_cwd)

    def copy_extra_files(self,files):
        for file in files:
            shutil.copy(file,self.directory)

def process_template(template_name,output_name,subs):
    template = open(template_name,"r")
    output = open(output_name,"w")
    unknown_vars = []
    for line in template.readlines():
        match = re.search("@@[A-z0-9]+@@",line)
        while match:
            var = string.replace(match.group(),"@@","")
            if subs.has_key(var):
                replacement = str(subs[var])
            else:
                replacement = ""
                unknown_vars.append(var)
            original = match.group()
            line = string.replace(line,original,replacement)
            match = re.search("@@[A-z0-9]+@@",line)
        match = re.search("__([A-z0-9]+){{(.*)}}{{(.*)}}__",line)
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
            line = string.replace(line,original,replacement)
            match = re.search("@@[A-z0-9]+@@",line)
        output.write(line)
    for var in unknown_vars:
        print "process_template warning: variable \"%s\" unkown." % var
    template.close()
    output.close()
    os.chmod(output_name,0755)


def create_new_directory(directory, version, overwrite):
    if version > 499:
        print "Sanity check failure: attempt to create directory version %d."\
              %version
        print "Maximum is 499."
        sys.exit(1)
    if version == 0:
        created_directory = directory
    else:
        created_directory = "%s.%02d" % (directory,version)
    if os.path.isdir(created_directory):
        if overwrite:
            shutil.rmtree(created_directory)
            os.mkdir(created_directory)
            print "created directory",created_directory
        else:
            created_directory = create_new_directory(directory,
                                                     version + 1,
                                                     overwrite)
    else:
        os.mkdir(created_directory)
        print "created directory",created_directory
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


if __name__ == "__main__":
    subs = {}
    subs["NODES"] = str(72)
    subs["PROCESSES"] = str(144)
    subs["WALLTIME"] = "01:00:00"
    #subs["MAIL"] = None
    process_template("generic_pbs_template.sh","synergia-pbs.sh",subs)
