#!/usr/bin/env python

import sys
import os
import synergia
import numpy as np
import tables
import pygsl.errno as pygslerrno
import pygsl
from pygsl import multiroots

# class that calculates corrector settings to create a 3 kick local orbit bump
class Three_bump:

    ##################################################

    # lattice is the lattice in which to create the bump
    # start_name is the name of the starting element for the bump
    #    (use a unique marker if necessary)
    # end_name is the name of the ending element for the bump
    #    (use a unique marker if necessary)
    # hcorr_names is a length=3 sequence of names of horizontal corrector
    #    elements that have a "kick=" attribute
    # vcorr_names is a length=3 sequence of names of vertical corrector
    #    elements that have a "kick=" attribute
    # target name is the name of the element at which the bump offset will
    #    be specified
    # optional argument coords is a two component sequence specifying
    #     (i, j) where i,j are in [0,1,2,3,4,5] meaning
    #     [x,xp,y,yp,cdt,dpop] are the coordinates that
    #     will be aimed for at the target location.  The default (0, 2)
    #     specifies (x,y) or (0,1) specifies (x, xp)
    # verbose = False/True on whether the module is chatty

    def __init__(self, lattice, start_name, end_name, hcorr_names, vcorr_names, target_name, coords=(0,2), verbose=False):
        self.lattice = lattice
        # I keep the elements separately so I can adjust them at the end when I know
        # what the settings are
        self.lattice_elements = lattice.get_elements()
        self.lattice_elements_idx = range(len(self.lattice_elements))
        self.start_name = start_name
        self.end_name = end_name
        self.hcorr_names = hcorr_names
        self.vcorr_names = vcorr_names
        self.target_name = target_name
        self.coords = coords
        self.verbose = verbose

        self._construct()

    def _construct(self):
        # until the Lattice::get_element_adaptor_map() method gets wrapped this routine will be
        # forced to assume the element adaptor type is MadX_adaptor_map.  Apparently, setting the
        # bump with one adaptor element type and using the main lattice of another element adaptor type
        # doesn't work.
        #self.bump_lattice = synergia.lattice.Lattice("bump", self.lattice.get_element_adaptor_map())
        self.bump_lattice = synergia.lattice.Lattice("bump", synergia.lattice.MadX_adaptor_map())
        self.bump_lattice.set_reference_particle(self.lattice.get_reference_particle())

        elem_names = [e.get_name() for e in self.lattice_elements]
        try:
            start_idx = elem_names.index(self.start_name)
        except:
            raise RuntimeError, "Three_bump: start_name: %s not found"%self.start_name
        try:
            end_idx = elem_names.index(self.end_name)
        except:
            raise RuntimeError, "Three_bump: end_name: %s not found"%self.end_name
    
        if start_idx < end_idx:
            for elem in self.lattice_elements[start_idx:end_idx+1]:
                self.bump_lattice.append(elem)
            self.bump_idx = self.lattice_elements_idx[start_idx:end_idx+1]
        else:
            for elem in self.lattice_elements[start_idx:] + self.lattice_elements[:end_idx+1]:
                self.bump_lattice.append(elem)
            self.bump_idx = self.lattice_elements_idx[start_idx:] + self.lattice_elements_idx[:end_idx+1]

        # self.bump_idx[] is the index pointing to the original elements in self.lattice_elements

        for elem in self.bump_lattice.get_elements():
            elem.set_string_attribute("extractor_type", "chef_propagate")

        bump_elements = self.bump_lattice.get_elements()
        bump_enames = [e.get_name() for e in bump_elements]

        # get the corrector elements and keep track of their original positions in the lattice

        # get horizontal corrector elements
        self.hcorr_elements = []
        self.hcorr_idx = []
        for i in range(3):
            try:
                hc_idx = bump_enames.index(self.hcorr_names[i])
            except:
                raise RuntimeError, "Three_bump: hcorr_name: %s not found"%self.hcorr_names[i]

            self.hcorr_elements.append(bump_elements[hc_idx])
            self.hcorr_idx.append(self.bump_idx[hc_idx])

        # get vertical corrector elements
        self.vcorr_elements = []
        self.vcorr_idx = []
        for i in range(3):
            try:
                vc_idx = bump_enames.index(self.vcorr_names[i])
            except:
                raise RuntimeError, "Three_bump: vcorr_name: %s not found"%self.vcorr_names[i]

            self.vcorr_elements.append(bump_elements[vc_idx])
            self.vcorr_idx.append(self.bump_idx[vc_idx])

        try:
            target_idx = bump_enames.index(self.target_name)
        except:
            raise RuntimeError, "Three_bump: target_name: %s not found"%self.target_name

        self.target_elem = bump_elements[target_idx]
        self.target_elem.set_string_attribute("force_diagnostics", "true")
        bump_elements[-1].set_string_attribute("force_diagnostics", "true")

    ##################################################

    # print out information about the bump settings

    def information(self):
        print "bump_lattice: ", len(self.bump_lattice.get_elements()), " elements, length: ", self.bump_lattice.get_length()
        print "horizontal correctors: "
        for i in range(3):
            print "\t%s, kick = %g"%(self.hcorr_elements[i].get_name(), self.hcorr_elements[i].get_double_attribute("kick"))
        print "vertical correctors: "
        for i in range(3):
            print "\t%s, kick = %g"%(self.vcorr_elements[i].get_name(), self.vcorr_elements[i].get_double_attribute("kick"))
        print "target element name: ", self.target_name


    ##################################################

    # Set the horizontal and vertical corrector values in preparation
    #    for the propagation of particles.
    # hcorr_values and vcorr_values are a length 3 sequence

    def set_corrector_elements(self, hcorr_values, vcorr_values):
        if len(self.hcorr_elements) != len(hcorr_values):
            raise RuntimeError, "set_corrector_elements: len(hcorr_elements) != len(hcorr_values)"

        if len(self.vcorr_elements) != len(vcorr_values):
            raise RuntimeError, "set_corrector_elements: len(vcorr_elements) != len(vcorr_values)"

        for i in range(len(hcorr_values)):
            self.hcorr_elements[i].set_double_attribute("kick", hcorr_values[i])
        for i in range(len(vcorr_values)):
            self.vcorr_elements[i].set_double_attribute("kick", vcorr_values[i])

    ##################################################

    # Propagate particle at 0,0,0,0,0,0 through the bump section saving
    #     diagnostics and returning the final position as array
    #     mean[6,2].  mean[:,0] is the position at the target location
    #     mean[:,1] is the position at the end of the bump section

    def propagate_zero(self):
        if self.verbose:
            verbosity = 1
        else:
            verbosity = 0
        comm = synergia.utils.Commxx()
        refpart = self.bump_lattice.get_reference_particle()
        stepper = synergia.simulation.Independent_stepper(self.bump_lattice, 1, 1)
        # 3 particles is the minimum so that the diagnostics don't crash
        bunch = synergia.bunch.Bunch(refpart, 3, 1.0e10, comm)
        bunch.get_local_particles()[:,0:6] = 0.0
        bunch_simulator = synergia.simulation.Bunch_simulator(bunch)
        bunch_simulator.add_per_forced_diagnostics_step(synergia.bunch.Diagnostics_basic("bump_basic.h5"))
        #bunch_simulator.add_per_step(synergia.bunch.Diagnostics_basic("step_basic.h5"))
        propagator = synergia.simulation.Propagator(stepper)
        propagator.propagate(bunch_simulator, 1, 1, verbosity)

        del propagator
        del stepper
        del bunch_simulator
        del bunch

        h5 = tables.openFile("bump_basic.h5")
        mean = h5.root.mean.read()
        h5.close()
        return mean

    ##################################################

    # utility function for fitting the bump.  given corrector settings
    # returns the positions and momenta of the 0 particle

    # x are the corrector settings (hc1, hc2, hc3, vc1, vc2, vc3). params are the desired (x,y) position at the midpoint.
    def bump_f(self, x, params):
        hc = (x[0], x[1], x[2])
        vc = (x[3], x[4], x[5])
        self.set_corrector_elements(hc, vc)
        mean = self.propagate_zero()
        return (mean[self.coords[0],0]-params[0], mean[self.coords[1],0]-params[1], mean[0, 1], mean[1,1], mean[2,1], mean[3,1])

    ##################################################

    # adjust the correctors to achieve an orbit bump.  Optional argument

    #  Returns the
    #     final values for the correctors as array
    #     [hc1, hc2, hc3, vc1, vc2, vc3]

    def set_bump(self, desired_position, coords=(0, 2)):
        params = desired_position
        mysys = multiroots.gsl_multiroot_function(self.bump_f, params, 6)
        solver = multiroots.hybrids(mysys, 6)

        tmp = np.zeros(6)
        solver.set(tmp)

        if self.verbose:
            print "  bump solver residuals:"
            #print "  %5s %9s %9s %9s %9s  %9s  %9s" %("iter", "hc1", "hc2", "hc3", "vc1", "vc2", "vc3")

        for iter in range(100):
            status = solver.iterate()
            r = solver.root()
            x = solver.getx()
            f = solver.getf()
            if self.verbose:
                print "  %5d % .7g % .7g % .7g % .7g  % .7g  % .7g" %(iter, f[0], f[1], f[2], f[3], f[4], f[5])
 
            status = multiroots.test_residual(f, 1.0e-13)
            if status == pygslerrno.GSL_SUCCESS:
                if self.verbose: print "Converged!!"
                break
        else:
            raise ValueError, "too many iterations"

        # set the correctors in the original lattice
        for i in range(3):
            self.lattice_elements[self.hcorr_idx[i]].set_double_attribute("kick", x[i])
        for i in range(3):
            self.lattice_elements[self.vcorr_idx[i]].set_double_attribute("kick", x[i+3])

        return x
        
    ##################################################
    ##################################################
#  just a little tester for the class
if __name__ == "__main__":
    lattice = synergia.lattice.MadX_reader().get_lattice("model", "foborodobo128.madx")
    print "read lattice: ", len(lattice.get_elements()), " elements, length = ", lattice.get_length()
    hcorr_names = ('hc1', 'hc2', 'hc3')
    vcorr_names = ('vc1', 'vc2', 'vc3')
    three_bump = Three_bump(lattice, 'm1', 'm2', hcorr_names, vcorr_names, 'm3', (0,2), True)
    three_bump.information()

    bump_settings = three_bump.set_bump((0.001, -0.0005))
    print "bump_settings: ", bump_settings[0], bump_settings[1], bump_settings[2], bump_settings[3], bump_settings[4], bump_settings[5]

    # propagate the whole lattice now
    comm = synergia.utils.Commxx()
    refpart = lattice.get_reference_particle()
    stepper = synergia.simulation.Independent_stepper(lattice, 1, 1)
    # 3 particles is the minimum so that the diagnostics don't crash
    bunch = synergia.bunch.Bunch(refpart, 3, 1.0e10, comm)
    bunch.get_local_particles()[:,0:6] = 0.0
    bunch_simulator = synergia.simulation.Bunch_simulator(bunch)
    bunch_simulator.add_per_step(synergia.bunch.Diagnostics_basic("step_basic.h5"))
    propagator = synergia.simulation.Propagator(stepper)
    propagator.propagate(bunch_simulator, 1, 1, 1)
    print "final coordinates: ", np.array2string(bunch.get_local_particles()[0, 0:6])
