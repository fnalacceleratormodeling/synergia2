#!/usr/bin/env python

from mad8_parser import Mad8_parser
from pylattice import Lattice_element, Element_adaptor_map, Lattice
from pyfoundation import pconstants, Four_momentum, Reference_particle

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
        self.parse_string(open(filename, 'r').read())

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

    def _extract_elements(self, entries, ancestors, lattice, reverse=False):
        repeat_num = 1
        if type(entries) <> type([]):
            entries = [entries]
        for entry in entries:
            if type(entry) == type([]):
                if reverse:
                    entry.reverse()
                for i in range(0, repeat_num):
                    self._extract_elements(entry, ancestors, lattice, reverse)
                if reverse:
                    entry.reverse()
                reverse = False
                repeat_num = 1
            elif entry in self.parser.lines:
                ancestors.append(entry)
                subentry = self.parser.lines[entry]
                if reverse:
                    subentry.reverse()
                for i in range(0, repeat_num):
                    self._extract_elements(subentry, ancestors, lattice, reverse)
                if reverse:
                    subentry.reverse()
                reverse = False
                repeat_num = 1
                ancestors.pop()
            else:
                if entry[len(entry) - 1] == '*':
                    repeat_num = int(entry[0:(len(entry) - 1)])
                elif entry == '-':
                    reverse = True
                else:
                    for i in range(0, repeat_num):
                        self._extract_element(entry, ancestors, lattice)
                    repeat_num = 1

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

    def get_lattice(self, line_name, filename=None):
        self._parser_check(filename, "get_lattice")
        lattice = Lattice(line_name)
        self._extract_elements(self.parser.lines[line_name], [line_name], lattice)
        self._extract_reference_particle(lattice)
        return lattice

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 2:
        line = sys.argv[1]
        filename = sys.argv[2]
    else:
        line = 'fodo'
        filename = 'fodo.lat'

    lattice = Mad8_reader().get_lattice(line, filename)
    print "has reference particle:", lattice.has_reference_particle()
    lattice.print_()
