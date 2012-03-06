#!/usr/bin/env python
import os, stat, sys, hashlib

from synergia.lattice import Mad8_parser
from synergia.lattice import Lattice_element, Element_adaptor_map, Lattice
from synergia.lattice import binary_save_lattice, binary_load_lattice
from synergia.foundation import pconstants, Four_momentum, Reference_particle
from mpi4py import MPI

def islist(x):
    return hasattr(x, "__iter__")

class Cache_entry:
    def __init__(self, hash_hex, index):
        self.hash_hex = hash_hex
        self.index = index

class Lattice_cache:
    def __init__(self, file_name, line_name):
        self.file_name = file_name
        self.real_file = (file_name != None)
        if self.real_file:
            self.line_name = line_name
            self.cache_dir = os.path.join(os.path.dirname(file_name),
                                           "lattice_cache")
            cache_file_basename = os.path.basename(file_name)
            cache_file_basename = cache_file_basename.replace('.','_') + "_cache"
            self.cache_file_name = os.path.join(self.cache_dir,
                                                cache_file_basename)
            self.cache = None
            self.max_index = 0

    def _get_hash_hex(self):
        hash = hashlib.sha1(open(self.file_name, 'r').read())
        return hash.hexdigest()

    def _read_cache(self):
        self.cache = {}
        if os.path.exists(self.cache_file_name):
            cache_file = open(self.cache_file_name, 'r')
            line_name = cache_file.readline().rstrip()
            while line_name:
                hash_hex = cache_file.readline().rstrip()
                index = int(cache_file.readline().rstrip())
                if (index > self.max_index):
                    self.max_index = index
                self.cache[line_name] = Cache_entry(hash_hex, index)
                line_name = cache_file.readline().rstrip()
            cache_file.close()

    def is_readable(self):
        if not self.real_file:
            return False
        if not self.cache:
            self._read_cache()
        readable = False
        exists = self.cache.has_key(self.line_name)
        if exists:
            readable = self.cache[self.line_name].hash_hex == self._get_hash_hex()
        return readable

    def _binary_file_name(self, index):
        return self.cache_file_name + ("_%d.bina" % index)

    def read(self):
        retval = None
        if self.is_readable():
            retval = Lattice()
            index = self.cache[self.line_name].index
            try:
                binary_load_lattice(retval, self._binary_file_name(index))
            except RuntimeError, e:
                retval = None
        return retval

    def write(self, lattice):
        if not self.real_file:
            return
        if not self.cache:
            self._read_cache()
        hash_hex = self._get_hash_hex()
        index = None
        need_to_write = True
        if self.cache.has_key(self.line_name):
            index = self.cache[self.line_name].index
            if self.cache[self.line_name].hash_hex == hash_hex:
                need_to_write = False
        if index == None:
            self.max_index += 1
            index = self.max_index
        if need_to_write:
            self.cache[self.line_name] = Cache_entry(hash_hex, index)
            if not os.path.exists(self.cache_dir):
                os.mkdir(self.cache_dir)
            cache_file = open(self.cache_file_name, 'w')
            for line_name in self.cache.keys():
                cache_file.write(line_name + '\n')
                cache_file.write('%s\n' % self.cache[line_name].hash_hex)
                cache_file.write('%d\n' % self.cache[line_name].index)
            cache_file.close()
            binary_save_lattice(lattice, self._binary_file_name(index))

class Mad8_reader:
    def __init__(self, element_adaptor_map=None):
        if element_adaptor_map:
            self.element_adaptor_map = element_adaptor_map
        else:
            self.element_adaptor_map = Element_adaptor_map()
        self.parser = None

    def get_element_adaptor_map(self):
        return self.element_adaptor_map

    def parse_string(self, string):
        self.parser = Mad8_parser()
        self.parser.parse(string)

    def parse(self, filename):
        try:
            the_file = open(filename, 'r')
        except IOError, e:
            os.system("ls")
            msg = "Mad8_reader: unable to open: " + filename + "\n"
            msg += "current working directory: " + os.getcwd() + "\n"
            raise RuntimeError(msg)
        self.parse_string(the_file.read())

    def _parser_check(self, filename, method):
        if not self.parser:
            if not filename:
                raise RuntimeError, "Mad8_reader." + method + \
            ": filename must be specified if no file or string has been parsed"
            self.parse(filename)

    def get_lines(self, filename=None):
        self._parser_check(filename, "get_lines")
        return self.parser.lines.keys()

    def _expand_type(self, short):
        retval = None
        for long in self.element_adaptor_map.get_adaptor_names():
            if long.find(short) == 0:
                if retval:
                    raise RuntimeError, 'Mad8_reader: ambiguous abbreviated type "'\
                        + short + '": matches "' + retval + '" and "' + long
                else:
                    retval = long
        return retval

    def get_lattice_element(self, label, filename=None):
        self._parser_check(filename, "get_lattice_element")
        type = self._expand_type(self.parser.labels[label].name)
        attributes = self.parser.labels[label].attributes
        element = Lattice_element(type, label)
        for attribute in attributes:
            if attributes[attribute] == None:
                element.set_string_attribute(attribute, '')
            else:
                try:
                    double_value = float(attributes[attribute])
                    element.set_double_attribute(attribute, double_value)
                except ValueError:
                    string_value = attributes[attribute]
                    element.set_string_attribute(attribute, string_value)
        return element

    def _extract_element(self, label, ancestors, lattice):
        element = self.get_lattice_element(label)
        for ancestor in ancestors:
            element.add_ancestor(ancestor)
        lattice.append(element)

    def _reverse(self, entries):
        retval = []
        if islist(entries):
            index = 0
            while index < len(entries):
                entry = entries[index]
                if islist(entry):
                    retval.insert(0, self._reverse(entries[index]))
                    index += 1
                elif (entry[len(entry) - 1] == '*'):
                    retval.insert(0, self._reverse(entries[index + 1]))
                    retval.insert(0, entry)
                    index += 2
                else:
                    retval.insert(0, entries[index])
                    retval.insert(0, '-')
                    index += 1
        else:
            retval = entries
        return retval

    def _extract_elements(self, entries, ancestors, lattice, reverse=False):
        if not islist(entries):
            entries = [entries]
        i = 0
        while i < len(entries):
            entry = entries[i]
            if islist(entry):
                if reverse:
                    entry = self._reverse(entry)
                self._extract_elements(entry, ancestors, lattice, reverse)
                i += 1
            elif entry in self.parser.lines:
                ancestors.append(entry)
                subentry = self.parser.lines[entry]
                if reverse:
                    subentry = self._reverse(subentry)
                self._extract_elements(subentry, ancestors, lattice, reverse)
                ancestors.pop()
                i += 1
            else:
                if entry == '-':
                    subentry = entries[i + 1]
                    self._extract_elements(self._reverse(subentry), ancestors,
                                            lattice, True)
                    i += 2
                elif entry[len(entry) - 1] == '*':
                    repeat_num = int(entry[0:(len(entry) - 1)])
                    subentry = entries[i + 1]
                    for repeat in range(0, repeat_num):
                        self._extract_elements(subentry, ancestors, lattice, reverse)
                    i += 2
                else:
                    self._extract_element(entry, ancestors, lattice)
                    i += 1

    def _extract_reference_particle(self, lattice):
        for command in self.parser.commands:
            particle = None
            mass = None
            charge = None
            pc = None
            energy = None
            gamma = None
            found_beam = False
            if command.name == 'beam':
                found_beam = True
                for attribute in command.attributes:
                    if attribute == "particle":
                        if command.attributes[attribute] == "proton":
                            mass = pconstants.mp
                            charge = pconstants.proton_charge
                        elif command.attributes[attribute] == "antiproton":
                            mass = pconstants.mp
                            charge = pconstants.antiproton_charge
                        elif command.attributes[attribute] == "electron":
                            mass = pconstants.me
                            charge = pconstants.electron_charge
                        elif command.attributes[attribute] == "positron":
                            mass = pconstants.me
                            charge = pconstants.positron_charge
                        elif command.attributes[attribute] == "muon":
                            mass = pconstants.mmu
                            charge = pconstants.muon_charge
                        elif command.attributes[attribute] == "antimuon":
                            mass = pconstants.mmu
                            charge = pconstants.antimuon_charge
                    elif attribute == "mass":
                        mass = float(command.attributes[attribute])
                    elif attribute == "charge":
                        charge = float(command.attributes[attribute])
                    elif attribute == "pc":
                        pc = float(command.attributes[attribute])
                    elif attribute == "energy":
                        energy = float(command.attributes[attribute])
                    elif attribute == "gamma":
                        gamma = float(command.attributes[attribute])
            if found_beam:
                four_momentum = Four_momentum(mass)
                if pc:
                    four_momentum.set_momentum(pc)
                if energy:
                    four_momentum.set_total_energy(energy)
                if gamma:
                    four_momentum.set_gamma(gamma)
                reference_particle = Reference_particle(charge, four_momentum)
                lattice.set_reference_particle(reference_particle)

    def get_lattice(self, line_name, filename=None, enable_cache_write=True,
                    enable_cache_read=True):
        lattice_cache = Lattice_cache(filename, line_name)
        lattice = None
        if enable_cache_read and lattice_cache.is_readable():
            lattice = lattice_cache.read()
        if not lattice:
            self._parser_check(filename, "get_lattice")
            lattice = Lattice(line_name)
            self._extract_elements(self.parser.lines[line_name], [line_name], lattice)
            self._extract_reference_particle(lattice)
            if enable_cache_write:
                MPI.COMM_WORLD.Barrier()
                if MPI.COMM_WORLD.Get_rank() == 0:
                    lattice_cache.write(lattice)
        return lattice

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 2:
        line = sys.argv[1]
        filename = sys.argv[2]
    else:
        line = 'fodo'
        filename = 'fodo.lat'

    lattice = Mad8_reader().get_lattice(line, filename, False, False)
    print "has reference particle:", lattice.has_reference_particle()
    lattice.print_()
