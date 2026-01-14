#!/usr/bin/env python
# coding: utf-8

# In[93]:


import sys, os
import re
import numpy as np
import synergia
import matplotlib.pyplot as plt
import h5py


# In[32]:


lattice_txt = """
f1 = 15.0;
f2 = 14.0;
fu: multipole, knl={0, 1/f1};
fd: multipole, knl={0, 1/f1};
du: multipole, knl={0, -1/f2};
dd: multipole, knl={0, -1/f2};

r: rfcavity, l=0.0, volt=%volt, harmon=5;

mp: multipole, knl={0.0, 0.0, %sxn, %ocn}, ksl={0.0, %qks, %sxs, %ocs};

bcell: sequence, l=20.0;
r, at=0.0;
fd, at=2.0;
du, at=6.0;
mp, at = 10.0;
dd, at=14.0;
fu, at=18.0;
endsequence;

beam, particle=proton, energy=pmass+0.8;
"""


# In[28]:


def get_lattice(volt=0.2, qks=0.0, sxn=0.0, sxs=0.0, ocn=0.0, ocs=0.0):
    latstring = lattice_txt
    latstring = re.sub('%volt', repr(volt), latstring)
    latstring = re.sub('%qks', repr(qks), latstring)
    latstring = re.sub('%sxn', repr(sxn), latstring)
    latstring = re.sub('%sxs', repr(sxs), latstring)
    latstring = re.sub('%ocn', repr(ocn), latstring)
    latstring = re.sub('%ocs', repr(ocs), latstring)

    reader = synergia.lattice.MadX_reader()
    reader.parse(latstring)
    
    lattice = reader.get_lattice('bcell')
    synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice)
    return lattice


# In[92]:


def lattice_info(lattice):

    print('lattice length = ', lattice.get_length())
    print(lattice)
    
    refpart = lattice.get_reference_particle()
    energy = refpart.get_total_energy()
    momentum = refpart.get_momentum()
    gamma = refpart.get_gamma()
    beta = refpart.get_beta()

    print('energy: ', energy)
    print('momentum: ', momentum)
    print('gamma: ', gamma)
    print('beta: ', beta)
    
    xtune, ytune, cdt = synergia.simulation.Lattice_simulator.calculate_tune_and_cdt(lattice)
    print('xtune: ', xtune, 'ytune: ', ytune, 'cdT: ', cdt)

    chrom = synergia.simulation.Lattice_simulator.get_chromaticities(lattice)
    print('h chromaticity: ', chrom.horizontal_chromaticity, 'v chromaticity: ', chrom.vertical_chromaticity)
    print('momentum compaction: ', chrom.momentum_compaction)
    print('slip factor: ', chrom.slip_factor)
    print('bucket length: ', synergia.simulation.Lattice_simulator.get_bucket_length(lattice))
    
    one_turn_map = synergia.simulation.Lattice_simulator.get_linear_one_turn_map(lattice)
    
    # print('one turn map')
    # print(one_turn_map)
    v, l = np.linalg.eig(one_turn_map)
    for z in v:
        print('eigenvalue: ', z, 'tune: ', np.log(z).imag/(2*np.pi))

        synergia.simulation.Lattice_simulator.CourantSnyderLatticeFunctions(lattice)
    synergia.simulation.Lattice_simulator.calc_dispersions(lattice)

    elem = lattice.get_elements()[0]
    beta_x = elem.lf.beta.hor
    alpha_x = elem.lf.alpha.hor
    beta_y = elem.lf.beta.ver
    alpha_y = elem.lf.alpha.ver

    print('beta_x: ', beta_x, 'alpha_x: ', alpha_x)
    print('beta_y: ', beta_y, 'alpha_y: ', alpha_y)




# In[67]:


def create_propagator(lattice):
    stepper = synergia.simulation.Independent_stepper_elements(1)
    propagator = synergia.simulation.Propagator(lattice, stepper)
    return propagator


# In[123]:


def create_sim(lattice):
    # need the lattice for the reference particle
    refpart = lattice.get_reference_particle()
    # 51 test particles in X
    # 50 test particles in Y
    # 50 test particles in s
    NX = 51
    NY = 50
    NS = 50
    num_particles = NX + NY + NS
    real_particles = 0.5e9
    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(
        refpart, num_particles, real_particles)
    
    # populate the bunch with particles aout to  5 transverse sigma and full bucket
    
    synergia.simulation.Lattice_simulator.CourantSnyderLatticeFunctions(lattice)
    synergia.simulation.Lattice_simulator.calc_dispersions(lattice)
    
    elem = lattice.get_elements()[0]
    beta_x = elem.lf.beta.hor
    alpha_x = elem.lf.alpha.hor
    beta_y = elem.lf.beta.ver
    alpha_y = elem.lf.alpha.ver
    
    betagamma = refpart.get_beta() * refpart.get_gamma()
    emitx = 20.0e-6/betagamma
    emity = 20.0e-6/betagamma
    print('emitx: ', emitx)
    print('emity: ', emity)
    sigx = np.sqrt(emitx * beta_x)
    sigy = np.sqrt(emity * beta_y)
    print('sigx: ', sigx)
    print('sigy: ', sigy)
    bucketlength = synergia.simulation.Lattice_simulator.get_bucket_length(lattice)
    
    bunch = sim.get_bunch(0, 0)
    bunch.checkout_particles()
    localpart = bunch.get_particles_numpy()
    dx = sigx/10.0
    dy = sigy/10.0
    ds = bucketlength/(2*50.0*refpart.get_beta())
    localpart[:, 0:6] = 0.0
    localpart[0:NX, 0] = dx * np.arange(0, NX)
    localpart[NX:NX+NY, 2] = dy * np.arange(1, NY+1)
    localpart[NX+NY:NX+NY+NS, 4] = ds * np.arange(0, NS)
    bunch.checkin_particles()
    
    return sim


# In[70]:


def register_diagnostics(sim):
    npart = sim.get_bunch(0, 0).get_total_num()

    diag = synergia.bunch.Diagnostics_bulk_track('tracks.h5', npart, 0)
    sim.reg_diag_per_turn(diag)
    return


# In[113]:


def run(volt=0.2, qks=0.0, sxn=0.0, sxs=0.0, ocn=0.0, ocs=0.0):
    
    print('Running with options')
    print(f'\tvolt: {volt} (RF voltage [MV])')
    print(f'\tqks: {qks} (skew quadrupole)')
    print(f'\tsxn: {sxn} (normal sextupole)')
    print(f'\tsxs: {sxs} (skew sextupole)')
    print(f'\tocn: {ocn} (normal octupole)')
    print(f'\tocs: {ocs} (skew octupole)')

    lattice = get_lattice(volt, qks, sxn, sxs, ocn, ocs)
    
    sim = create_sim(lattice)
    
    # normal form
    nf = synergia.simulation.Lattice_simulator.calculate_normal_form_o3(lattice)

    register_diagnostics(sim)

    fh = open("nf.dat", "w")

    def action(sim, lattice, turn):
        bunch = sim.get_bunch()
        bunch.checkout_particles()
        hparts = bunch.get_particles_numpy()

        for p in range(bunch.size()):
            hp = np.zeros(6, dtype='d')
            hp[:] = hparts[p, 0:6]
            tm = nf.convert_xyz_to_normal(hp)

            print('{0} {1} {2} {3} {4} {5} '.format(
                tm[0].real, tm[0].imag,
                tm[1].real, tm[1].imag,
                tm[2].real, tm[2].imag), file=fh, end="")

        print(file=fh)
        fh.flush()

    # register turn end action for writing particles out
    sim.reg_prop_action_turn_end(action)

    propagator = create_propagator(lattice)
    simlog = synergia.utils.parallel_utils.Logger(0, 
            synergia.utils.parallel_utils.LoggerV.INFO_TURN)

    propagator.propagate(sim, simlog, 1000)

    fh.close()

    print('run finished')


# In[112]:


def plot_traj(track_file, nf_file):
    h5 = h5py.File(track_file, 'r')
    tracks = h5.get('track_coords')
    # get mass and pz to convert cdT to z
    mass = h5.get('mass')[()]
    pz = h5.get('pz')[()]
    betagamma = pz/mass
    gamma = np.sqrt(betagamma**2 + 1)
    beta = betagamma/gamma

    nfdata = np.loadtxt(nf_file)
    nturn = nfdata.shape[0]
    npart = nfdata.shape[1]//6
    nftracks = np.zeros((nturn, npart, 6))
    for i in range(nturn):
        nftracks[i, :, :] = nfdata[i, :].reshape(npart, 6)
    del nfdata

    plt.figure()
    plt.title('x vs xp')
    for part in range(1, 51, 10):
        plt.plot(tracks[:, part, 0], tracks[:, part, 1], '.', label=f'particle {part}')
    plt.xlabel('x [m]')
    plt.ylabel('xp')
    plt.legend(loc='best')
    
    plt.figure()
    plt.title('real X vs. imag X')
    for part in range(1,51, 10):
        plt.plot(nftracks[:, part, 0], nftracks[:, part, 1], '.', label=f'particle {part}')
    plt.xlabel('real X')
    plt.ylabel('imag X')
    plt.legend(loc='best')

    plt.figure()
    plt.title('y vs. yp')
    for part in range(51, 101, 10):
        plt.plot(tracks[:, part, 2], tracks[:, part, 3], '.', label=f'particle {part}')
    plt.xlabel('y [m]')
    plt.ylabel('yp')
    plt.legend(loc='best')

    plt.figure()
    plt.title('real Y vs. imag Y')
    for part in range(51,101, 10):
        plt.plot(nftracks[:, part, 2], nftracks[:, part, 3], '.', label=f'particle {part}')
    plt.xlabel('real Y')
    plt.ylabel('imag Y')
    plt.legend(loc='best')

    plt.figure()
    plt.title('z vs. dp/p')
    for part in list(range(101, 151, 10)) + [150]:
        plt.plot(-tracks[:, part, 4]*beta, tracks[:, part, 5], '.', label=f'particle {part}')
    plt.xlabel('z [m]')
    plt.ylabel('dp/p')
    plt.legend(loc='best')

    plt.figure()
    plt.title('real S vs. imag S')
    for part in list(range(101, 151, 10)) + [150]:
        plt.plot(nftracks[:, part, 4], nftracks[:, part, 5], '.', label=f'particle {part}')
    plt.xlabel('real S')
    plt.ylabel('imag S')
    plt.legend(loc='best')
    
    plt.show()

    h5.close()

if __name__ == "__main__":
    if sys.argv[1] == "run":
        run()
    elif sys.argv[1] == "plot":
        plot_traj('tracks.h5', 'nf.dat')
    else:
        pass


