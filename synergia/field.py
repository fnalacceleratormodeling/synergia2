#!/usr/bin/env python

from UberPkgpy import *
import physics_constants
import math

class Field:
    def __init__(self, beam_parameters, processor_grid, computational_grid,
                 pipe_radius):
        self.pipe_radius = pipe_radius
        self.period_length = 2.0*math.pi*physics_constants.PH_MKS_c/\
                             beam_parameters.get_omega() *\
                             beam_parameters.get_beta() *\
                             beam_parameters.get_gamma()
        self.compdom = CompDom()
        construct_CompDom_external(self.compdom,
                                   beam_parameters.get_distparam(),
                                   beam_parameters.get_nparam(),
                                   beam_parameters.get_dist_type(),
                                   computational_grid.get_nx(),
                                   computational_grid.get_ny(),
                                   computational_grid.get_nz(),
                                   processor_grid.get_pgrid2d(),
                                   processor_grid.get_nrow(),
                                   processor_grid.get_ncol(),
                                   computational_grid.get_bc_num(),
                                   self.get_pipe_radius(),
                                   self.get_pipe_radius(),
                                   self.get_period_length())
        self.fieldquant = FieldQuant()
        construct_FieldQuant_external(self.fieldquant,
                                      computational_grid.get_nx(),
                                      computational_grid.get_ny(),
                                      computational_grid.get_nz(),
                                      self.compdom,
                                      processor_grid.get_pgrid2d())

    def get_fieldquant(self):
        return self.fieldquant
    def get_compdom(self):
        return self.compdom
    def get_period_length(self):
        return self.period_length
    def get_pipe_radius(self):
        return self.pipe_radius
