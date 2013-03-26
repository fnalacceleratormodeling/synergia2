#!/usr/bin/env python
import string
import sys

class _option:
    def __init__(self, name, default_value, doc_string, val_type, valid_values):
        self.name = name
        self.value = default_value
        self.doc_string = doc_string
        self.val_type = val_type
        if type(default_value) == type([]):
            self.length = len(default_value)
            if self.val_type == None:
                self.val_type = type(default_value[0])
        else:
            self.length = 1
            if self.val_type == None:
                self.val_type = type(default_value)
        self.valid_values = valid_values

    def _apply_val_type(self, val):
        val = str(val)
        have_multiple = False
        if len(val) > 0:
            if val[0] == '_':
                have_multiple = True
        if have_multiple:
            retval = val
        else:
            if self.val_type == bool:
                if (val == "True" or val == "true" or val == "t" or val == "1"):
                    retval = True
                elif (val == "False" or val == "false" or val == "f" or val == "0"
                      or val == "nil"):
                    retval = False
                else:
                    raise RuntimeError('Cannot convert "' + val + '" to boolean')
            else:
                retval = self.val_type(val)
        return retval

    def get(self):
        if self.value == None:
            return None
        else:
            if self.val_type:
                if self.length == 1:
                    return self._apply_val_type(self.value)
                else:
                    return map(self._apply_val_type, self.value)
            else:
                return self.value

    def set(self, value):
        if self.length == 1:
            self.value = value
        else:
            self.value = string.split(value, ",")
            if len(self.value) != self.length:
                raise RuntimeError("Options: expected a % d - tuple" % self.length)
        if self.valid_values:
            if self.valid_values.count(value) == 0:
                raise RuntimeError("Options: " + value \
                                   + " is not a valid value for " + self.name \
                                   + "\n" \
                                   +"valid values are " + str(self.valid_values) )

class Options:
    '''Define a set of command-line options.
    Hierarchical sets of option may be created with add_suboptions'''
    def __init__(self, name):
        self.name = name
        self.dict = {}
        self.suboptions = []
        self.is_options = True

    def options_name(self):
        return self.name

    def add(self, option, default_value, doc_string, val_type=None, valid_values=None):
        '''Add a new option definition'''
        if hasattr(self, option):
            raise RuntimeError('Options: option name "' + option +
                              '" already in use')
        self.dict[option] = _option(option, default_value, doc_string, val_type, valid_values)
        setattr(self, option, self.get(option))

    def get(self, option):
        if self.dict.has_key(option):
            return self.dict[option].get()
        else:
            for suboption in self.suboptions:
                if suboption.has_option(option):
                    return suboption.get(option)

    def set(self, option, value):
        if self.dict.has_key(option):
            self.dict[option].set(value)
            setattr(self, option, self.get(option))
        else:
            found = 0
            for suboption in self.suboptions:
                if suboption.has_option(option):
                    suboption.set(option, value)
                    found = 1
            if not found:
                print "Error: option", option, "not found."

    def has_option(self, option):
        if self.dict.has_key(option):
            return 1
        else:
            for suboption in self.suboptions:
                if suboption.has_option(option):
                    return 1
        return 0

    def options(self, include_suboptions=True):
        '''Returns a list of all options, including suboptions'''
        list = self.dict.keys()
        if include_suboptions:
            for suboption in self.suboptions:
                list = list + suboption.options()
        return list

    def add_suboptions(self, suboptions):
        name = suboptions.options_name()
        if hasattr(self, name):
            raise RuntimeError('Options: option name "' + name +
                              '" already in use')
        setattr(self, name, suboptions)
        self.suboptions.append(suboptions)

    def override(self, overrides):
        for name in dir(overrides):
            if (name[0] != '_') and (name != "is_override"):
                if self.has_option(name):
                    self.set(name, getattr(overrides, name))
                else:
                    print "warning: override", name, "not found in existing options"


    def usage(self):
        '''Print usage message to stdout'''
        for suboption in self.suboptions:
            suboption.usage()
        print "%s options:" % self.name
        all_options = self.options(include_suboptions=0)
        all_options.sort()
        for option in all_options:
            sys.stdout.write(" %s=" % option)
            val_type = self.dict[option].val_type
            if val_type == type(1):
                typename = "int"
            elif val_type == type(1.0):
                typename = "float"
            elif val_type == type(""):
                typename = "str"
            elif val_type == type(True):
                typename = "bool"
            else:
                typename = "x"
            if self.dict[option].length == 1:
                print "<%s> " % typename,
            else:
                sys.stdout.write("<")
                for i in range(1, self.dict[option].length):
                    sys.stdout.write("%s, " % typename)
                print "%s> " % typename,
            print "%s," % self.dict[option].doc_string,
            print "default=",
            for item in range(0, self.dict[option].length):
                if self.dict[option].length == 1:
                    val = self.dict[option].get()
                else:
                    val = self.dict[option].get()[item]
                if val_type == type(1.0):
                    if val != None:
                        sys.stdout.write("%g" % val)
                    else:
                        sys.stdout.write("None")
                else:
                    sys.stdout.write(str(val))
                if item + 1 < self.dict[option].length:
                    print ", ",
            if self.dict[option].valid_values:
                print ", valid values: ", self.dict[option].valid_values,
            print

        print

    def parse_argv(self, argv):
        '''Parse command-line arguments from argv'''
        for arg in argv[1:]:
            if arg == "--help" or arg == "help":
                self.usage()
                sys.exit(0)
            pair = string.split(arg, "=")
            if len(pair) < 2:
                self.usage_error(arg)
            if self.has_option(pair[0]):
		first = pair.pop(0)
                self.set(first, string.join(pair, "="))
            else:
                self.usage_error(arg)

    def usage_error(self, unknown_argument):
        self.usage()
        sys.stderr.write("Error: unknown argument \"%s\"\n\n" % unknown_argument)
        sys.exit(1)

class Override:
    def __init__(self):
        self.is_override = True

if __name__ == "__main__":
    stupid = Options("stupid")
    stupid.add("fred", 1, "fred option")
    stupid.add("barney", 2, "barney option", float)
    stupid.add("wilma", [1, 2, 3], "wilma vector", float)
#    stupid.add("add", 1, "no worky", int)

    really_stupid = Options("really_stupid")
    really_stupid.add("daffy", "duck", "daffy's species", str)
    really_stupid.add("bugs", "male", "bugs's gender", str)
    stupid.add_suboptions(really_stupid)

    stupid.parse_argv(sys.argv)
    for option in stupid.options():
        print option, stupid.get(option)
    print "stupid.fred =", stupid.fred
    print "stupid.really_stupid.bugs =", stupid.really_stupid.bugs
