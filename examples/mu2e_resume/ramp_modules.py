#!/usr/bin/env python

import numpy
import synergia

class Pickle_helper:
    __getstate_manages_dict__ = 1
    def __init__(self, *args):
        self.args = args
    def __getinitargs__(self):
        return self.args
    def __getstate__(self):
        return self.__dict__
    def __setstate__(self, state):
        self.__dict__ = state

class Ramp_actions(synergia.simulation.Propagate_actions, Pickle_helper): 
    def __init__(self, ramp_turns, turns_to_extract, initial_k1, final_k1, final_k2l, rfko_kicker):
        synergia.simulation.Propagate_actions.__init__(self)
        Pickle_helper.__init__(self, ramp_turns, turns_to_extract, initial_k1, final_k1, final_k2l, rfko_kicker)
        self.ramp_turns = ramp_turns
        self.turns_to_extract = turns_to_extract
        self.initial_k1 = initial_k1
        self.final_k1 = final_k1
        self.final_k2l = final_k2l
        self.rfko_kicker = rfko_kicker
    def turn_end_action(self, stepper, bunch, next_turn):
        synergia_elements = stepper.get_lattice_simulator().get_lattice().get_elements()
        next_turn += 1
        myrank = synergia.utils.Commxx().get_rank()
        global avg_rate, old_intensity, kp
        k0 = 1.0
        #k1 = 0.28
        #if (next_turn <= 1500):
        #    k0 = 0.025
        #elif ((next_turn > 1500) and (next_turn <= 2500)):
        #    k0 = 0.1
        #elif ((next_turn > 2500) and (next_turn <= 3500)):
        #    k0 = 0.17
        #else:
        #    k0 = k1

###############################################################################
#           sextupole ramping...
###############################################################################
        if next_turn <= self.ramp_turns:
            index = 0
            for element in synergia_elements:
                if element.get_type() == "multipole":
                    new_k2l = self.final_k2l[index]*next_turn/self.ramp_turns
                    element.set_double_attribute("k2l", new_k2l)
                    index += 1
                    #if myrank == 0:
                    #    print
                    #    print "    turn                             :", 
                    #    print next_turn
                    #    print "    updated multipole                :", 
                    #    print element.get_name()
                    #    print "    final k2l                        :", 
                    #    print self.final_k2l[index - 1], "1/m^3"
                    #    print "    new k2l                          :", 
                    #    print new_k2l, "1/m^3"
                    #    print "    new k2l (real)                   :",
                    #    print element.get_double_attribute("k2l"), "1/m^3"
            if next_turn == 1:
                old_intensity = bunch.get_total_num()
                n0 = old_intensity
                avg_rate = n0 / float(self.turns_to_extract - self.ramp_turns)
            if next_turn == self.ramp_turns:
                kp = k0

        if next_turn > self.ramp_turns and next_turn <= self.turns_to_extract:
###############################################################################
#           tune quads ramping...
###############################################################################
            total_turn = self.turns_to_extract - self.ramp_turns
            epsilon = 1. * (next_turn - self.ramp_turns) / total_turn

            #~a*sqrt(t) + b
            #fe = numpy.sqrt(epsilon)
            #~piece-wise linear
            #f1 = 4.15
            #factor = [12.9, 6.2, 4.15, f1, f1, f1]
            #dturn = [300, 300, 300, 300, 300]
            nturn = next_turn - self.ramp_turns
            f1 = 0.42
            factor = [12.92, 9.77, 5.70, 4.40, 3.56, 2.91, 2.37, 2.04, 1.81, \
                       1.62, 1.48, 1.36, 1.19, 1.08, 1.01, 0.86, 0.77, 0.67, \
                       0.60, 0.52, 0.50, 0.49, 0.45, 0.48, 0.43, 0.42, 0.54, \
                       0.42,   f1,   f1,   f1,   f1,   f1,   f1,   f1,   f1]
            dturn  = [  200,  200,  200,  200,  200,  300,  300,  300,  300, \
                        300,  300,  400,  500,  600,  600, 1000, 1200, 1400, \
                       1600, 1800, 3000, 2000, 2000, 1200,  200,  200,  200, \
                        200, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000]
            #           300   500   700   900  1100  1400  1700  2000  2300
            #          2600  2900  3300  3800  4400  5000  6000  7200  8600
            #         10200 12000 15000 17000 19000 20200 20400 20600 20800
            #         21000
            aturn = 0
            id = 0
            while (nturn > aturn):
                turni = aturn
                aturn += dturn[id]
                turnf = aturn
                id += 1
            id -= 1
            if (nturn > turni) and (nturn <= turnf):
                csum = 0.0
                for j in range(0, id + 1):
                    depsilon = 1. * dturn[j] / total_turn
                    csum += (factor[j] - factor[id]) * depsilon
            fe = factor[id] * epsilon + csum
            #~linear
            #fe = epsilon
            #~no squeeze (constant separatrix)
            #fe = 0.0

            index = 0
            for element in synergia_elements:
                if element.get_type() == "quadrupole":
                    # quadratic
                    #coeff = 0.0 #1.0e-12
                    #new_k1 = (coeff * total_turn * total_turn * epsilon \
                    #        +(self.final_k1[index] - self.initial_k1[index])) \
                    #        * (epsilon - 1.0) + self.final_k1[index]

                    # piecewise-linear
                    new_k1 = (self.final_k1[index] - self.initial_k1[index]) \
                            * fe + self.initial_k1[index]

                    element.set_double_attribute("k1", new_k1)
                    index += 1
                    #if myrank == 0:
                    #    print
                    #    print "    turn                             :",
                    #    print next_turn
                    #    print "    updated quadrupole               :", 
                    #    print element.get_name()
                    #    print "    delta k1                         :",
                    #    print self.final_k1[index-1] - self.initial_k1[index-1]
                    #    print "    epsilon                          :", epsilon
                    #    print "    initial k1                       :", 
                    #    print self.initial_k1[index - 1], "1/m^2"
                    #    print "    final k1                         :", 
                    #    print self.final_k1[index - 1], "1/m^2"
                    #    print "    new k1                           :", 
                    #    print new_k1, "1/m^2"
                    #    print "    new k1 (real)                    :", 
                    #    print element.get_double_attribute("k1"), "1/m^2"
###############################################################################
#           rfko kicker strength setting...
###############################################################################
            id = next_turn - self.ramp_turns - 1
            feedback_loop = 200
            if (id % feedback_loop) == 0:
                new_intensity = bunch.get_total_num()
                if (id == 0):
                    spill_rate = float(old_intensity - new_intensity) / self.ramp_turns
                else:
                    spill_rate = float(old_intensity - new_intensity) / feedback_loop
                kp = k0 #* avg_rate / spill_rate
                if myrank == 0:
                    print
                    print "    turn                             :", next_turn
                    print "    intensity                        :",
                    print next_turn, old_intensity, new_intensity
                    print "    spill rate                       :",
                    print next_turn, avg_rate, spill_rate, kp
                old_intensity = new_intensity


            for element in synergia_elements:
                if element.get_name() == "rfko_kicker":
                    new_kicker_strength = kp * self.rfko_kicker[id]
                    element.set_double_attribute("kick", new_kicker_strength)
                    #if myrank == 0:
                    #    print
                    #    print "    turn                             :",
                    #    print next_turn
                    #    print "    rfko kicker strength             :",
                    #    print new_kicker_strength * 1e6, "urad"
                    #    print "    kp                               :", kp
                    #    print "    k0                               :", k0

        if (next_turn % 100) == 0:
            #if myrank == 0:
            #    print
            #    print "    horizontal chromaticity          :",
            #    print stepper.get_lattice_simulator().get_horizontal_chromaticity()
            #    print "    vertical chromaticity            :",
            #    print stepper.get_lattice_simulator().get_vertical_chromaticity()
            for element in synergia_elements:
                if element.get_name() == "mpsa":
                    lattice_functions = stepper.get_lattice_simulator().get_lattice_functions(element)
                    ax = lattice_functions.alpha_x
                    bx = lattice_functions.beta_x
                    ay = lattice_functions.alpha_y
                    by = lattice_functions.beta_y
                    ex = element.get_double_attribute("phasespace_aperture_emittance_x")
                    element.set_double_attribute("phasespace_aperture_alpha_x", ax)
                    element.set_double_attribute("phasespace_aperture_beta_x", bx)
                    ey = element.get_double_attribute("phasespace_aperture_emittance_y")
                    element.set_double_attribute("phasespace_aperture_alpha_y", ay)
                    element.set_double_attribute("phasespace_aperture_beta_y", by)
                    if myrank == 0:
                        print
                        print "....Phasespace aperture parameters"
                        print "    turn                             :", 
                        print next_turn
                        print "    emitx                            :",
                        print ex, "mm-mrad"
                        print "    emity                            :",
                        print ey, "mm-mrad"
                        print "    alpha_x                          :",
                        print ax, "m^{1/2}"
                        print "    alpha_y                          :",
                        print ay, "m^{1/2}"
                        print "    beta_x                           :", bx, "m"
                        print "    beta_y                           :", by, "m"
        stepper.get_lattice_simulator().update()

        # print tune-quad strengths
        #tqf_id = 0
        #for element in synergia_elements:
        #    if element.get_name() == 'tqf':
        #        if myrank == 0 and tqf_id == 0:
        #            print "tqf strength :",
        #            print next_turn, element.get_double_attribute("k1")
        #        tqf_id += 1

        # print twiss parameters
        #if next_turn == self.ramp_turns:
        #    index = 0
        #    length = 0
        #    for element in synergia_elements:
        #        lattice_functions = stepper.get_lattice_simulator().get_lattice_functions(element)
        #        ax = lattice_functions.alpha_x
        #        bx = lattice_functions.beta_x
        #        rx = (1.0 + ax * ax) / bx
        #        ay = lattice_functions.alpha_y
        #        by = lattice_functions.beta_y
        #        ry = (1.0 + ay * ay) / by
        #        length += element.get_length()
        #        if myrank == 0:
        #            print "%5d %11s %11s %9.5f %9.5f %12.8f %12.8f %12.8f \
        #            %12.8f %12.8f %12.8f" % (index, element.get_name(), \
        #            element.get_type(), element.get_length(), length, ax, bx, \
        #            rx, ay, by, ry)
        #        index += 1
