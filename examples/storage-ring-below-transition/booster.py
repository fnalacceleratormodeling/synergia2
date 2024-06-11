import sys, os
import numpy as np
import h5py
import synergia
from mpi4py import MPI

from booster_options import opts

from booster_rf_ramp import *
#from booster_momentum import e_vs_t


DEBUG = False

macroparticles = opts.macroparticles
real_particles = opts.real_particles
emitx = opts.emitx
emity = opts.emity
stddpop = opts.stddpop

accel_turns = opts.accel_turns
turns = opts.turns

gridx = opts.gridx
gridy = opts.gridy
gridz = opts.gridz
steps = opts.steps

commsize = opts.comm_group_size

#-------------------------------------------------------------------------------

# quick and dirty twiss parameter calculator from 2x2 courant-snyder map array
def map2twiss(csmap):
    cosmu = 0.5 * (csmap[0,0]+csmap[1,1])
    asinmu = 0.5*(csmap[0,0]-csmap[1,1])

    if abs(cosmu) > 1.0:
        raise RuntimeError("map is unstable")

    mu =np.arccos(cosmu)

    # beta is positive
    if csmap[0,1] < 0.0:
        mu = 2.0 * np.pi - mu

    beta = csmap[0,1]/np.sin(mu)
    alpha = asinmu/np.sin(mu)
    tune = mu/(2.0*np.pi)

    return (alpha, beta, tune)

#-------------------------------------------------------------------------------

def set_adjust_markers(lattice):
    for elem in lattice.get_elements():
        # focussing quads are in the cps with name qsxx
        # defocussing quads in  the cpl with name qlxx
        # focussing sextupoles are in the cps with name sxsxx
        # defocussing quads in  the cpl with name sxlxx
        if elem.get_name() == "qsxx":
            elem.set_marker(synergia.lattice.marker_type.h_tunes_corrector)
        elif elem.get_name() == "qlxx":
            elem.set_marker(synergia.lattice.marker_type.v_tunes_corrector)
        elif elem.get_name() == "sxsxx":
            elem.set_marker(synergia.lattice.marker_type.h_chrom_corrector)
        elif elem.get_name() == "sxlxx":
            elem.set_marker(synergia.lattice.marker_type.v_chrom_corrector)

#-------------------------------------------------------------------------------

def set_tunes_and_chromaticities(lattice, xtune, ytune, xchrom, ychrom, logger):

    set_adjust_markers(lattice)
    
    print(f'Setting  xtune: {xtune}, ytune: {ytune}', file=logger)
    synergia.simulation.Lattice_simulator.adjust_tunes(lattice, xtune, ytune, 1.0e-6)

    print(f'Setting xchrom {xchrom}, ychrom: {ychrom}', file=logger)
    synergia.simulation.Lattice_simulator.adjust_chromaticities(lattice, xchrom, ychrom)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def print_bunch_statistics(bunch, logger=sys.stdout):

    parts = bunch.get_particles_numpy()
    print(parts.shape,  ", ", parts.size , file=logger)
    print("shape: {0}, {1}".format(parts.shape[0], parts.shape[1]), file=logger)

    mean = synergia.bunch.Core_diagnostics.calculate_mean(bunch)
    std = synergia.bunch.Core_diagnostics.calculate_std(bunch, mean)
    print("mean = {}".format(mean), file=logger)
    print("std = {}".format(std), file=logger)

#-------------------------------------------------------------------------------

# Read the lattice from a previously created json file
def get_lattice():
    with open(opts.lattice_file, 'r') as f:
        jsonlattice = f.read()
    lattice = synergia.lattice.Lattice.load_from_json(jsonlattice)

    return lattice

#-------------------------------------------------------------------------------

# set the voltage and tune the lattice
# voltage is total voltage in GV

def set_rf(lattice, voltage, harmno, phase, above_transition=False):

    # above transition, phase needs to be (pi - phase) for longitudinal stability
    if above_transition:
        phase_set = np.pi - phase
    else:
        phase_set = phase

    if DEBUG:
        print('setrf: lattice: ', id(lattice))
        print('set_rf: voltage=', voltage, ', harmno = ', harmno, ', phase = ', phase_set, end='')
    # count RF cavities
    cavities = 0
    for elem in lattice.get_elements():
        if elem.get_type() == synergia.lattice.element_type.rfcavity:
            cavities = cavities + 1
    if DEBUG:
        print(' for ', cavities, ' cavities')

    # Set the RF cavity voltage
    for elem in lattice.get_elements():
        if elem.get_type() == synergia.lattice.element_type.rfcavity:
            elem.set_double_attribute('volt', 1000*voltage/cavities)
            elem.set_double_attribute('lag', phase_set/(2*np.pi))
            elem.set_double_attribute('harmon', harmno)

    synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice)

    for elem in lattice.get_elements():
        if DEBUG and elem.get_type() == synergia.lattice.element_type.rfcavity:
            print('set_rf: ', elem)
            break

    return lattice


#-------------------------------------------------------------------------------

# set apertures around the Booster for short-straight, long-straight
# and RF cavity elements.

# Set a 2.25" diameter circular aperture everywhere.

def set_apertures(lattice):

    inches_to_m = 0.0254
    short_straight_aperture = 4.5
    long_straight_aperture = 3.125
    rf_aperture = 2.25

    fmag_vertices = [(3.74,0.506), (-2.16, +1.09), (-2.16, -1.09), (3.74, -0.506)]
    dmag_vertices = [(3.50,1.52), (-2.40, +0.901), (-2.40, -0.901), (3.50, -1.52)]

    for elem in lattice.get_elements():

        if elem.has_string_attribute('type'):
            etype = elem.get_string_attribute('type')

            if etype == "shortstraight":
                # short straight aperture 4.5"
                elem.set_string_attribute("aperture_type", "circular")
                elem.set_double_attribute("circular_aperture_radius",
                                          short_straight_aperture * inches_to_m/2)
            

            elif etype == "longstraight":
                # long straight aperture 3.125"
                elem.set_string_attribute("aperture_type", "circular")
                elem.set_double_attribute("circular_aperture_radius",
                                          long_straight_aperture * inches_to_m/2)

                
            elif etype == "fmag":

                # focussing magnet
                elem.set_string_attribute("aperture_type", "polygon")
                elem.set_double_attribute("the_number_of_vertices", 4)
                # vertex coordinates in inches
                vertices = fmag_vertices
                for i, vtx in enumerate(vertices):
                    elem.set_double_attribute("pax%d"%(i+1), vtx[0]*inches_to_m)
                    elem.set_double_attribute("pay%d"%(i+1), vtx[1]*inches_to_m)
                    # find distance from center to aperture edge
                    # need the equation for the between (-2.54, 0.925) to (3.46, 0.925)
                    # line that goes between (x1, y1) and (x2, y2)
                    #  A = (y2-y1), B = -(x2-x1), C = -x1*(y2-y1) + y1*(x2-x1)
                    # min_radius2 = C**2/(A**2 + B**2)
                    A = vertices[1][1] - vertices[0][1]
                    B = -(vertices[1][0] - vertices[0][0])
                    C = -vertices[0][0]*(vertices[1][1] - vertices[0][1]) + vertices[0][1]*(vertices[1][0] - vertices[0][0])
                    elem.set_double_attribute("min_radius2", C**2/(A**2 + B**2))

            elif etype == "dmag":
                # defocussing magnet
                elem.set_string_attribute("aperture_type", "polygon")
                elem.set_double_attribute("the_number_of_vertices", 4)
                # vertex coordinates in inches
                vertices = dmag_vertices
                for i, vtx in enumerate(vertices):
                    elem.set_double_attribute("pax%d"%(i+1), vtx[0]*inches_to_m)
                    elem.set_double_attribute("pay%d"%(i+1), vtx[1]*inches_to_m)
                    # find distance from center to aperture edge
                    # need the equation for the between (-2.54, 0.925) to (3.46, 0.925)
                    # line that goes between (x1, y1) and (x2, y2)
                    #  A = (y2-y1), B = -(x2-x1), C = -x1*(y2-y1) + y1*(x2-x1)
                    # min_radius2 = C**2/(A**2 + B**2)
                    A = vertices[1][1] - vertices[0][1]
                    B = -(vertices[1][0] - vertices[0][0])
                    C = -vertices[0][0]*(vertices[1][1] - vertices[0][1]) + vertices[0][1]*(vertices[1][0] - vertices[0][0])
                    elem.set_double_attribute("min_radius2", C**2/(A**2 + B**2))

            elif etype == "rfaperture":
                elem.set_string_attribute("aperture_type", "circular")
                elem.set_double_attribute("circular_aperture_radius",
                                          rf_aperture * inches_to_m/2)


#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

def lattice_functions(lattice):
    synergia.simulation.Lattice_simulator.CourantSnyderLatticeFunctions(lattice)
    synergia.simulation.Lattice_simulator.calc_dispersions(lattice)
    elements = lattice.get_elements()
    return elements[-1].lf

#-------------------------------------------------------------------------------

# return stdx, stdy
def bunch_parameters(lattice, emitx, emity, stddpop):
    # emittances are geometric RMS
    ex = emitx
    ey = emity
    lf = lattice_functions(lattice)
    betax = lf.beta.hor
    dispx = lf.dispersion.hor
    betay = lf.beta.ver
    stdx = np.sqrt(ex * betax + dispx**2 * stddpop**2)
    stdy = np.sqrt(ey * betay)
    return  (stdx, stdy)    

#-------------------------------------------------------------------------------

# Use the lattice from from the pip-ii injection point and the flat booster
# start to transfer pip-ii particles into this lattice.

# particles are the particle coordinates with shape [:, 0:6]

# The transformed particles are placed back into the input
# particles array.

def transfer_map_from_pip2(particles, logger):

    # derived from reading the last lattice
    # postinjection.02/orbump_rampdown_0099.lsx
    # and extracting lattice functions with Synergia2.
    alphax1 = 0.10909911380287352
    betax1 = 6.281951552320663
    alphay1 = 0.0355584253729662
    betay1 = 19.07484635252504
    betaz1 = 972.5834349469636

    # these parameters were derived from the beam distribution
    # after postinjection.02
    # betax1 =  7.074348593176607
    # alphax1 =  0.09300611292219932
    # betay1 =  22.15477233177723
    # alphay1 = 0.03179490512610093

    # these lattice functions are calculated with Synergia3
    # from the sbbooster.madx file.
    alphax2 = -4.49357088517076e-07
    betax2 = 33.53638957494266
    alphay2 = 7.703299144983451e-15
    betay2 = 5.348053006973944
    betaz2 = 974.5243838758502

    map = np.zeros((6,6), dtype='d')
    map[0,0] = np.sqrt(betax2/betax1)
    #map[0,1] = 0 because it is multiplied by sin(phase)
    map[1,1] = np.sqrt(betax1/betax2)
    map[1,0] = (alphax1 - alphax2)/np.sqrt(betax1*betax2)
    map[2,2] = np.sqrt(betay2/betay1)
    # map[2,3] = 0 because it is multiplied by sin(phase)
    map[3,3] = np.sqrt(betay1/betay2)
    map[3,2] = (alphay1 - alphay2)/np.sqrt(betay1*betay2)
    map[4,4] = np.sqrt(betaz2/betaz1)
    map[5,5] = np.sqrt(betaz1/betaz2)

    mapped = np.dot(map, particles[:,0:6].transpose())
    
    particles[:, 0:6] = mapped.transpose()


#-------------------------------------------------------------------------------

def populate_bunch_emit(lattice, bunch, stdx, stdy, stddpop):
    map = synergia.simulation.Lattice_simulator.get_linear_one_turn_map(lattice)
    beta = lattice.get_reference_particle().get_beta()
    covars = synergia.bunch.get_correlation_matrix(map, stdx, stdy, stddpop, beta, (0,2,5))

    means = np.zeros(6, dtype='d')
    dist = synergia.foundation.Random_distribution(opts.seed, MPI.COMM_WORLD.rank)
    synergia.bunch.populate_6d(dist, bunch, means, covars)

    #-------------------------------------------------------------------------------

def populate_bunch_covar(bunch, covar):
    covars = eval(opts.covar)
    means = np.zeros(6, dtype='d')
    dist = synergia.foundation.Random_distribution(opts.seed, MPI.COMM_WORLD.rank)
    synergia.bunch.populate_6d(dist, bunch, means, covars)


#-------------------------------------------------------------------------------

def create_simulator(refpart, num_particles, real_particles):
    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(
        refpart, num_particles, real_particles)
    return sim


#-------------------------------------------------------------------------------

def register_diagnostics(sim):
    # diagnostics
    diag_full2 = synergia.bunch.Diagnostics_full2("diag.h5")
    sim.reg_diag_per_turn(diag_full2)

    diag_bt = synergia.bunch.Diagnostics_bulk_track("tracks.h5", opts.tracks, 0)
    sim.reg_diag_per_turn(diag_bt)

    if opts.save_particles:
        diag_part = synergia.bunch.Diagnostics_particles("particles.h5")
        sim.reg_diag_per_turn(diag_part, opts.particles_period)

    if opts.step_tracks:
        diag_step_track = synergia.bunch.Diagnostics_bulk_track("step_tracks.h5", opts.step_tracks, 0)
        sim.reg_diag_per_step(diag_step_track)
    
    if opts.step_diag:
        diag_stp = synergia.bunch.Diagnostics_full2("diag_step.h5")
        sim.reg_diag_per_step(diag_stp)

#-------------------------------------------------------------------------------

def get_impedance_op(lattice_length):
    wn = opts.wave_number[:]
    ImpedStr = synergia.collective.Impedance_options(wake_file=opts.wakefile_f, wake_type=opts.waketype, z_grid=opts.imp_zgrid)
    ImpedStr.orbit_length = lattice_length
    ImpedStr.nstored_turns = opts.registered_turns
    ImpedStr.mwf_xlead = opts.wf_xlead
    ImpedStr.mwf_xtrail = opts.wf_xtrail
    ImpedStr.mwf_ylead = opts.wf_ylead
    ImpedStr.mwf_ytrail = opts.wf_ytrail
    ImpedStr.mwf_zwake = opts.wf_zwake
    return ImpedStr

#-------------------------------------------------------------------------------


def get_propagator(lattice):
    if  DEBUG:
        print('get_propagator operating on lattice: ', id(lattice))

    if not opts.collective:
        sc_ops = synergia.collective.Dummy_CO_options()
    elif opts.collective == "rectangular":
        grid = [opts.gridx, opts.gridy, opts.gridz]
        pipesize = [opts.pipesizex, opts.pipesizey, opts.pipesizez]
        sc_ops = synergia.collective.Space_charge_rectangular_options(grid, pipesize)
    else:
        # unknown collective operator not defined
        raise RuntimeError(f'unhandled collective operator specified: {opts.collective}')

    if opts.stepper == "elements":
        if opts.collective:
            raise RuntimeError('may not use elements stepper with a collective operator')

        stepper = synergia.simulation.Independent_stepper_elements(steps)
    else:
        stepper = synergia.simulation.Split_operator_stepper(sc_ops, steps)

    if opts.impedance:
        imp_coll_op = get_impedance_op(lattice.get_length())
        stepper.append_collective_op(imp_coll_op)

    propagator = synergia.simulation.Propagator(lattice, stepper)

    if DEBUG:
        print('lattice from propagator: ', id(propagator.get_lattice()))

    return propagator

#-------------------------------------------------------------------------------

def print_map_and_tune(lattice, outfs):
    map = synergia.simulation.Lattice_simulator.get_linear_one_turn_map(lattice)
    print('map: ', file=outfs)
    print(np.array2string(map, max_line_width=200), file=outfs)
    ax,bx,qx = map2twiss(map[0:2, 0:2])
    print('beta_x: ', bx, ', alpha_x: ', ax, ', qx: ', qx, file=outfs)
    ay,by,qy = map2twiss(map[2:4, 2:4])
    print('beta_y: ', by, ', alpha_y: ', ay, ', qy: ', qy, file=outfs)
    longitudinal_map = map[4:6, 4:6]
    # correct sign of cdt components
    longitudinal_map[0, 1] = -longitudinal_map[0, 1]
    longitudinal_map[1, 0] = -longitudinal_map[1, 0]
    a, b, n = map2twiss(longitudinal_map)
    print('longitudinal beta: ', b, file=outfs)
    print('synchrotron tune: ', n, file=outfs)
 
#-------------------------------------------------------------------------------

def save_json_lattice(lattice, jlfile='init_lattice.json'):
    # save the lattice
    if MPI.COMM_WORLD.rank == 0:
        f = open(jlfile, 'w')
        print(lattice.as_json(), file=f)
        f.close()

#-------------------------------------------------------------------------------

# ---------------- get lattice RF parameters


def get_ringparm(in_turn, inlattice):
    if DEBUG:
        print('get_ringparm: inlattice: ', id(inlattice))
    class p:
        turn = in_turn
    # get the RF voltage, and phase
    p.lattice_E = inlattice.get_lattice_energy()

    rffreq = None
    rfvolt = None
    rflag = None
    ncavities = 0
    for elem in inlattice.get_elements():
        if elem.get_type() == synergia.lattice.element_type.rfcavity:
            ncavities = ncavities + 1
            if (not rfvolt) or (not rflag):
                rfvolt = elem.get_double_attribute("volt")
                rflag = elem.get_double_attribute("lag")
                rffreq = elem.get_double_attribute("freq")
    if not ncavities:
        raise RuntimeError('Did not locate any RF cavities in get_ringparm')
    p.rfvolt = rfvolt
    p.phase = 2*np.pi*rflag
    p.freq = rffreq*1.0e6
    p.dE = ncavities * p.rfvolt * np.sin(p.phase)
    return p

# ---------------- end get lattice RF parameters




def main():
    logger = synergia.utils.Logger() # logs to screen with verbosity INFO

    lattice = get_lattice()
    print('Read lattice, {} elements, length: {}'.format(len(lattice.get_elements()), lattice.get_length()), file=logger)
    # print('before lattice changes map and tunes', file=logger)
    # print_map_and_tune(lattice, logger)

    # turn on transverse aperture maybe
    if opts.trans_aperture:
        set_apertures(lattice)

    # print('after aperturess map and tunes', file=logger)
    # print_map_and_tune(lattice, logger)

    if DEBUG:
        print('lattice: ', id(lattice))

    harmonic_number = 84
    set_rf(lattice, opts.init_voltage, harmonic_number, 0.0)
    bucket_length = lattice.get_length()/harmonic_number
    
    refpart = lattice.get_reference_particle()
    energy = refpart.get_total_energy()
    momentum = refpart.get_momentum()
    gamma = refpart.get_gamma()
    beta = refpart.get_beta()

    print('energy: ', energy, file=logger)
    print('momentum: ', momentum, file=logger)
    print('gamma: ', gamma, file=logger)
    print('beta: ', beta, file=logger)

    synergia.simulation.Lattice_simulator.set_closed_orbit_tolerance(1.0e-11)

    # print('after get_closed_orbit map and tunes', file=logger)
    # print_map_and_tune(lattice, logger)

    target_xtune, target_ytune, orig_cdt = synergia.simulation.Lattice_simulator.calculate_tune_and_cdt(lattice)

    if opts.set_xtune:
        target_xtune = opts.set_xtune
    if opts.set_ytune:
        target_ytune = opts.set_ytune

    orig_chrom = synergia.simulation.Lattice_simulator.get_chromaticities(lattice)
    target_xchrom = orig_chrom.horizontal_chromaticity
    target_ychrom = orig_chrom.vertical_chromaticity
    if opts.set_xchrom:
        target_xchrom = opts.set_xchrom
    if opts.set_ychrom:
        target_ychrom = opts.set_ychrom

    # print('before set_tunes_and_chromaticitiess map and tunes', file=logger)
    # print_map_and_tune(lattice, logger)

    set_tunes_and_chromaticities(lattice, target_xtune, target_ytune, target_xchrom, target_ychrom, logger)

    # print the map and twiss parameters
    print_map_and_tune(lattice, logger)

    chrom = synergia.simulation.Lattice_simulator.get_chromaticities(lattice)
    print('horizontal chromaticity: ', chrom.horizontal_chromaticity, ', hchrom_prime: ', chrom.horizontal_chromaticity_prime, file=logger)
    print('vertical chromaticity: ', chrom.vertical_chromaticity, ', vchrom_prime: ', chrom.vertical_chromaticity_prime, file=logger)
    print('momentum compaction: ', chrom.momentum_compaction, file=logger)
    transition_energy = synergia.foundation.pconstants.mp * np.sqrt(1/chrom.momentum_compaction)

    if opts.generate_bunch:
        num_particles = macroparticles
        print('num_particles: ', num_particles, file=logger)
        sim = create_simulator(refpart, num_particles, real_particles)
        # set behavior of longitudinal boundary
        bdy = getattr(synergia.bunch.LongitudinalBoundary, opts.z_boundary)
        sim.get_bunch(0, 0).set_longitudinal_boundary(bdy, bucket_length)
        lb, value = sim.get_bunch(0,0).get_longitudinal_boundary()
        print(f'Longitudinal boundary: {lb}, size: {value}', file=logger)

        if opts.matching == "emit":
            stdx, stdy = bunch_parameters(lattice, emitx, emity, stddpop)

            populate_bunch_emit(lattice, sim.get_bunch(0, 0), stdx, stdy, stddpop)
        elif opts.matching == "covar":
            covars = eval(opts.covar)
            populate_bunch_covar(sim.get_bunch(0, 0), covars)
        else:
            raise RuntimeError(f'invalid matching option: "{opts.matching}"')            
    else:
        # create the bunch that will receive the particles by first
        # reading the bunch parameters from the file itself.

        if opts.read_openpmd:
            series = io.Series(opts.particles_file,
                               io.Access_Type.read_only)
            if 'mass' in series.attributes:
                file_mass = series.get_attribute('mass')
            else:
                file_mass = opts.openpmd_mass
            if 'pz' in series.attributes:
                file_pz = series.get_attribute('pz')
            else:
                file_pz = opts.openpmd_pz

            file_energy = np.sqrt(file_pz**2 + file_mass**2)
            num_particles = series.iterations[opts.openpmd_iter].particles["bunch_particles"].to_df().shape[0]
            series.close()
        else: # read bunch parameters from HDF5 file
            h5 = h5py.File(opts.particles_file, 'r')
            num_particles = h5.get('particles').shape[0]
            # Have to read the momentum from the particles file which might be
            # different than the energy/momentum read from the lattice if we're
            # restarting an acceleration simulation. Set the reference particle
            # to be consistent.
            file_pz = h5.get('pz')[()]
            file_mass = h5.get('mass')[()]
            file_energy = np.sqrt(file_pz**2 + file_mass**2)
            h5.close()
        refpart.set_total_energy(file_energy)
        lattice.set_lattice_energy(file_energy)
        print('Populating bunch with {} macroparticles, real charge {} from file {}'.format(num_particles, real_particles, opts.particles_file), file=logger)
        print('Beam energy from file: ', file_energy, file=logger)
        sim = create_simulator(refpart, num_particles, real_particles)
        bunch = sim.get_bunch(0, 0)
        # set behavior of longitudinal boundary
        bdy = getattr(synergia.bunch.LongitudinalBoundary, opts.z_boundary)
        bunch.set_longitudinal_boundary(bdy, bucket_length)
        lb, value = sim.get_bunch(0,0).get_longitudinal_boundary()
        print(f'Longitudinal boundary: {lb}, size: {value}', file=logger)

        # Now read the particle data from the files
        if opts.read_openpmd:
            bunch.read_openpmd_file(opts.particles_file)
        else:
            bunch.read_file_legacy(opts.particles_file)

        particles = bunch.get_particles_numpy()[:,0:6]
        if opts.correct_longitudinal_offset:
            # Need to determine the offsets to be corrected
            bunchmeans = synergia.bunch.Core_diagnostics.calculate_mean(bunch)
            print(f'correcting longitudinal offset: cdt: {bunchmeans[4]:.3f}, dp/p: {bunchmeans[5]:.5g}', file=logger)
            particles[:, 4] = particles[:, 4] - bunchmeans[4]
            particles[:, 5] = particles[:, 5] - bunchmeans[5]

        transfer_map_from_pip2(particles, logger)
        bunch.checkin_particles()

    # Bunch now created, set the start time
    sim.get_bunch(0,0).get_reference_particle().set_bunch_abs_time(opts.tstart)
    print('Starting simulation at t0 = ', opts.tstart, file=logger)

    # print statistics of populsted bunch
    print_bunch_statistics(sim.get_bunch(0, 0), logger)

    register_diagnostics(sim)

    # logger for simulation
    simlog = synergia.utils.parallel_utils.Logger(0, 
            synergia.utils.parallel_utils.LoggerV.INFO_TURN)
            #synergia.utils.parallel_utils.LoggerV.INFO)
            #synergia.utils.parallel_utils.LoggerV.INFO_STEP)

    # save lattice as json
    save_json_lattice(lattice)

    # ---------------- booster_ramp_context saves state for turn_end_action

    # save state for Booster ramp
    class booster_ramp_context:
        pass

    #    end booster_ramp_context -------------------------------




    # ---------------- Initialize booster_ramp_context 

    booster_ramp_context.transition_energy = transition_energy
    booster_ramp_context.pip2ramp = PIP2ramp()
    # start at vstep 1 which assumes injection has finished and
    # beginning the squeeze, so start determining vstep starting
    # at 1

    booster_ramp_context.logger = logger

    # The target ring energy at the end of the _next_ turn
    # we will initialize to the current energy because we are starting without
    # acceleration. This will have to be address if there acceleration at the beginning
    # of the run
    booster_ramp_context.next_turn_energy = lattice.get_reference_particle().get_total_energy()

    # smooth the bunch mean cdt
    n_cdt_smooth = 10
    booster_ramp_context.smooth_cdt = np.zeros(n_cdt_smooth)

    # ---------------- End initialize booster_ramp_context 

    # ---------------- define the turn end action
    # jiggle the reference particle at the ends of turns so that if the bunch
    # momentum has changed the lattice reference particle keeps up with it.

    fh = open('rf_history.txt', 'w')
    print("#turn time lattice-energy bunch-design-energy bunch-energy rf-voltage rf-phase rf-frequency dE", file=fh, flush=True)

    def turn_end_action(sim, inlattice, turn):
        if DEBUG:
            print('turn_end_action turn: ', turn, ', lattice: ', id(inlattice), flush=True)

        if id(inlattice) != booster_ramp_context.lattice_address:
            raise RuntimeError("!!!!!!!!!!!! wrong lattice address passed to turn_end_action!!!")

        # save lattice this turn
        save_json_lattice(inlattice, f'booster-{turn:05d}.json')

        bunch = sim.get_bunch()
        bunch_design_E = bunch.get_design_reference_particle().get_total_energy()
        bunch_E = bunch.get_reference_particle().get_total_energy()
        lattice_E = inlattice.get_lattice_energy()

        # the reference particle energy may have increased because of
        # accelerating phase which may have been set to actually accelerator and/or
        # to add energy to compensate impedance energy losses.
        bunch.get_design_reference_particle().set_total_energy(bunch_E)

        # Renormalize the bunch reference energy and the design reference energy to match what
        # the machine energy should
        # be as set in the ramp_context.
        bunch_E = booster_ramp_context.next_turn_energy
        bunch.adjust_bunch_particles_reference_energy(bunch_E)
        bunch.get_design_reference_particle().set_total_energy(bunch_E)
        bunch.get_reference_particle().set_total_energy(bunch_E)

        if DEBUG:
            print('    Entry bunch design refpart energy: ', bunch_design_E, flush=True)
            print('    Entry bunch refpart energy: ', bunch_E, flush=True)
            print('    Entry lattice energy: ', lattice_E, flush=True)
            print('    Entry context next turn energy: ', booster_ramp_context.next_turn_energy)

        # need to set the lattice energy and tune the RF cavities for the
        # new energy
        inlattice.set_lattice_energy(bunch_E)

        if DEBUG:
            print('Turn: ', turn, flush=True)
            new_bunch_design_E = bunch.get_design_reference_particle().get_total_energy()
            new_lattice_E = inlattice.get_lattice_energy()
            print('    exit bunch design refpart energy: ', new_bunch_design_E, flush=True)
            print('    exit lattice energy: ', new_lattice_E, flush=True)
            #assert (new_lattice_E == bunch_E), "why didn't the lattice energy change?"

        synergia.simulation.Lattice_simulator.tune_circular_lattice(inlattice)

        p = get_ringparm(turn, inlattice)

        current_time = bunch.get_reference_particle().get_bunch_abs_time()

        print(turn, current_time, lattice_E, bunch_design_E, bunch_E, p.rfvolt, p.phase, p.freq, p.dE, file=fh, flush=True)

        # Storage ring is constant energy so this should not change
        current_energy = lattice.get_lattice_energy()
        if DEBUG:
            print('current Booster energy: ', current_energy, flush=True)
        # what will the energy be at the end of the current turn (assuming velocity constant?)
        velocity = bunch.get_reference_particle().get_beta()*synergia.foundation.pconstants.c
        turn_time = inlattice.get_length()/velocity
        end_of_turn_time = current_time + turn_time
        end_of_turn_energy = current_energy
        # Set it in the context so it's available next turn
        booster_ramp_context.next_turn_energy = end_of_turn_energy
        if DEBUG:
            print('end-of-current-turn Booster energy: ', end_of_turn_energy, flush=True)

        # calculate the lattice parameters for correcting the energy
        # here after the lattice energy has been updated for acceleration.

        # I'm going to need the longitudinal lattice parameters so get the
        # map of this lattice
        turnmap = synergia.simulation.Lattice_simulator.get_linear_one_turn_map(inlattice)
        longitudinal_map = np.array(turnmap[4:6, 4:6]) # create new object
        # correct sign of cdt components
        longitudinal_map[0, 1] = -longitudinal_map[0, 1]
        longitudinal_map[1, 0] = -longitudinal_map[1, 0]
        a, beta_z, q_s = synergia.simulation.Lattice_simulator.map_to_twiss(longitudinal_map)
        print(f'longitudinal beta turn {turn}: {beta_z}, Qs: {q_s}, phase_advance/turn: {2.0*np.pi*q_s}', file=booster_ramp_context.logger, flush=True)
        print(f'map components:\n{turnmap[4, 0]} {turnmap[4, 1]} {turnmap[4,4]} {-turnmap[4,5]}', file=booster_ramp_context.logger, flush=True)
        print(f'{turnmap[5, 0]} {turnmap[5, 1]} {-turnmap[5,4]} {turnmap[5,5]}', file=booster_ramp_context.logger, flush=True)
        sin_phase_advance = np.sin(2.0*np.pi*q_s)
        # dE = dp/p * p * beta
        p_beta = bunch.get_reference_particle().get_beta() * bunch.get_reference_particle().get_momentum()

        # what is the actual beam energy currently?
        beam_energy = bunch.get_reference_particle().get_total_energy()

        print('beam energy: ', beam_energy, ', p_beta: ', p_beta, file=booster_ramp_context.logger, flush=True)

        # with longitudinal impedance, the beam energy is reduced but
        # it only shows up as a negative mean dp/p. If compensating is
        # active, adjust the beam_energy to account for this.
        if opts.impedance and opts.imp_compensate:
            # calculate mean dp/p
            bunch.checkout_particles()
            bunch_means = synergia.bunch.Core_diagnostics.calculate_mean(bunch)
            mean_cdtpos = bunch_means[4]
            mean_dpop = bunch_means[5]
            mean_x = bunch_means[0]
            mean_xp = bunch_means[1]
            print(f'bunch_mean cdt: {mean_cdtpos}, dpop: {mean_dpop}, mean x: {mean_x}, mean xp: {mean_xp}', file=booster_ramp_context.logger, flush=True)
            # dE = beta * dp


        # We need to possibly add energy to the beam over the next turn
        # What is the current voltage?
        # Storage ring RF voltage is (should be) static so it's the same as initial
        current_rf_voltage = opts.init_voltage
        if DEBUG:
            print('current RF voltage: ', current_rf_voltage, flush=True)

        # this is the nominal energy error without including impedance
        delta_E = end_of_turn_energy - beam_energy
        print('turn energy gain: ', delta_E, file=booster_ramp_context.logger, flush=True)

        if opts.impedance and (opts.imp_compensate == 1):
            # add/subtract energy to compensate for impedance energy
            # loss. Assume impedance energy loss results in <dp/p> != 0
            # so calculate energy gain to restore it to 0 in one turn.

            # dE = p beta * dp/p
            dE_comp_kick = -mean_dpop*p_beta
            print('apparent required energy booster: ', dE_comp_kick, file=booster_ramp_context.logger, flush=True)

            delta_E = delta_E + dE_comp_kick
            print('total delta_E: ', delta_E, file=booster_ramp_context.logger, flush=True)

        elif opts.impedance and (opts.imp_compensate == 2)

            # Add/subtract additional energy boost if beam centroid is
            # slipping longitudinally.
            # in these units, +dp/p leads to -cdt. Add energy with the same sign
            # as the cdt to compensate drift.
            dpop_comp_kick = -opts.imp_comp_gain*mean_cdtpos*turnmap[4,4]/turnmap[4,5]
            print('comp_dpop_kick: ', dpop_comp_kick, file=booster_ramp_context.logger, flush=True)
            # dE_comp_kick = comp_dpop_kick * p_beta
            # print('comp_dE_kick: ', dE_comp_kick, file=booster_ramp_context.logger, flush=True)

            # delta_E = delta_E + dE_comp_kick
            # print('total delta_E: ', delta_E, file=booster_ramp_context.logger, flush=True)

        if DEBUG:
            print('Need to add ', delta_E, ' to beam energy', flush=True)
        # set phase to:
        if abs(delta_E/current_rf_voltage) > 1.0:
            if MPI.COMM_WORLD.rank == 0:
                # Oh no! can't keep up with ramp
                print("requested larger energy jump than possible with current RF voltage: requested: ", delta_E, ", voltage: ", current_rf_voltage, ", turn: ", turn, flush=True)
            # use max voltage I have
            desired_phase = np.pi/2
        else:
            desired_phase = np.arcsin(delta_E/current_rf_voltage)

        print('RF voltage: ', current_rf_voltage, ' desired_phase: ', desired_phase, file=booster_ramp_context.logger, flush=True)


        if DEBUG:
            print('Set voltage to ', current_rf_voltage, ' ,set phase to ', desired_phase, flush=True)

        # set RF including whether we're above transition
        set_rf(inlattice, current_rf_voltage, 84, desired_phase,
               (beam_energy > booster_ramp_context.transition_energy))

        if DEBUG:
            # check RF parameters
            for elem in propagator.get_lattice().get_elements():
                if elem.get_type() == synergia.lattice.element_type.rfcavity:
                    assert (abs(22*elem.get_double_attribute('volt')*0.001/current_rf_voltage-1) < 1.0e-6), "cavity voltage: {} should be {}".format(elem.get_double_attribute('volt')*0.022, current_rf_voltage)
                    assert (abs( (elem.get_double_attribute('lag')/(2*np.pi))/desired_phase - 1) < 1.0e-6) 
                break
        

        # end of turn-end-action

    #---------------------------------------------------------------------------


    #---------------------------------------------------------------------------
    propagator = get_propagator(lattice)

    booster_ramp_context.lattice = propagator.get_lattice()
    booster_ramp_context.lattice_address = id(propagator.get_lattice())

    if DEBUG:
        print('propagator created containing lattice: ', id(propagator.get_lattice()))

    sim.reg_prop_action_turn_end(turn_end_action)

    # Initial conditions for propagation
    current_rf_voltage = opts.init_voltage
    print('setting Voltage to: ', current_rf_voltage, file=logger)
    set_rf(propagator.get_lattice(), current_rf_voltage, 84, 0)  # activate this for real simulations

    # set_rf(propagator.get_lattice(), 0.001, 84, 0.01)  # this is for testing only

    propagator.propagate(sim, simlog, opts.turns)


#-------------------------------------------------------------------------------

if __name__ == '__main__':
    main()

