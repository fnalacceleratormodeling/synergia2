#!/usr/bin/env python

from math import pi
import synergia_workflow

opts = synergia_workflow.Options("booster")

opts.add("seed", 12345791, "Pseudorandom number generator seed", int)

opts.add('lattice_file', 'booster-06800.json', 'JSON format lattice object')

opts.add("enable_rf", True, "turn on the RF cavities")
opts.add('init_voltage', 0.042727272727272725*22e-3, 'initial RF voltage [GV]')

opts.add('tstart',0.011779574949254912, 'Time start for simulation [s]')

opts.add('matching', 'emit', 'matching (either emit or covar)')
opts.add('emitx', 6.737134323979276e-07, 'geometric x emittance')
opts.add('emity', 6.647081159425567e-07, 'geometric y emittance')
opts.add('stddpop', 1.04873e-3, "RMS dp/p")

#opts.add('covar', '\nnp.array([[ 2.91636681e-05,  8.94577986e-09,  5.99375367e-09,\n         1.90809797e-09,  6.97308230e-05,  3.54114740e-06],\n       [ 8.94577986e-09,  1.55662798e-08,  1.78012700e-10,\n        -3.30751876e-12, -8.88754563e-08,  1.19126016e-09],\n       [ 5.99375367e-09,  1.78012700e-10,  3.51790066e-06,\n        -3.59780090e-09, -1.43518523e-06, -1.93106639e-09],\n       [ 1.90809797e-09, -3.30751876e-12, -3.59780090e-09,\n         1.25600427e-07, -4.33983936e-07, -6.05312035e-10],\n       [ 6.97308230e-05, -8.88754563e-08, -1.43518523e-06,\n        -4.33983936e-07,  4.25663831e-02,  2.21651577e-05],\n       [ 3.54114740e-06,  1.19126016e-09, -1.93106639e-09,\n        -6.05312035e-10,  2.21651577e-05,  1.09983463e-06]])\n', 'covariance matrix as string')
opts.add('covar', '\nnp.array([[ 2.91636681e-05,  8.94577986e-09,  5.99375367e-09,\n         1.90809797e-09,  6.97308230e-05,  3.54114740e-06],\n       [ 8.94577986e-09,  1.55662798e-08,  1.78012700e-10,\n        -3.30751876e-12, -8.88754563e-08,  1.19126016e-09],\n       [ 5.99375367e-09,  1.78012700e-10,  3.51790066e-06,\n        -3.59780090e-09, -1.43518523e-06, -1.93106639e-09],\n       [ 1.90809797e-09, -3.30751876e-12, -3.59780090e-09,\n         1.25600427e-07, -4.33983936e-07, -6.05312035e-10],\n       [ 6.97308230e-05, -8.88754563e-08, -1.43518523e-06,\n        -4.33983936e-07,  4.25663831e-02,  2.21651577e-05],\n       [ 3.54114740e-06,  1.19126016e-09, -1.93106639e-09,\n        -6.05312035e-10,  2.21651577e-05,  1.09983463e-06]])\n', 'covariance matrix as string')

#opts.add("stdz", 0.0001, "z standard deviation")
opts.add('real_particles', 82716049382.71603, 'bunch charge = 6.7e12/81')

opts.add("macroparticles", 10000, "number of macro particles")
opts.add("init_turns", 100, "number of initial turns with no acceleration")
opts.add("squeeze_turns", 1000, 'number of turns to squeeze')
opts.add("accel_turns", 1000, "number of turns with acceleration on")
opts.add("turns", 1000, "number of turns")

opts.add('generate_bunch', True, 'Generate a new bunch')
#opts.add('particles_file', '/pscratch/sd/e/egstern/pip2/postinjection.02/pip-ii-injected.h5', 'File with particles')
opts.add('read_openpmd', False, 'Read particles file in openpmd format')
opts.add('openpmd_iter', 0, 'The iteration to read of the openpmd file')

# proton mass by default
opts.add('openpmd_mass', 0.938272046, 'mass to use for bunch if openpmd file does not contain  missing mass attribute')
# Booster injection energy KE=0.8 by default
opts.add('openpmd_pz', 1.4632960307470262, 'pz to use for bunch if openpmd file does not contain pz attribute')

#opts.add('particles_file', '/wclustre/accelsim/egstern/pip2/postinjection.02/pip-ii-injected.h5', 'File with particles')
opts.add('particles_file', '/pscratch/sd/e/egstern/pip2/postinjection.02/pip-ii-injected.h5', 'File with particles')
opts.add('correct_longitudinal_offset', True, "whether to correct the longitudinal offset in read-in particles")

opts.add("collective", None, 'collective operator [off|rectangular', str)

opts.add("gridx", 64, "x grid size")
opts.add("gridy", 64, "y grid size")
opts.add("gridz", 256, "z grid size")
opts.add("comm_group_size", 1, "Communication group size for space charge solvers (must be 1 on GPUs), probably 16 on CPU", int)

opts.add('pipesizex', 0.14986, 'pipe size x direction')
opts.add('pipesizey', 0.077216,'pipe size y direction')
opts.add('pipesizez', 474.202752/84, 'pipe size z direction')

opts.add('impedance', False, 'Enable Impedance (wakes) if True')
#opts.add('wf_xlead', 1.0, 'scale factor x leading')
#opts.add('wf_xtrail', 1.0, 'scale factor x trailing')
#opts.add('wf_ylead', 1.0, 'scale factor y leading')
#opts.add('wf_ytrail', 1.0, 'scale factor y tailing')

# disable transverse wakes
opts.add('wf_xlead', 0.0, 'scale factor x leading')
opts.add('wf_xtrail', 0.0, 'scale factor x trailing')
opts.add('wf_ylead', 0.0, 'scale factor y leading')
opts.add('wf_ytrail', 0.0, 'scale factor y tailing')
opts.add('wf_zwake', 1.0, 'scale factor z wake')

# The following come from booster_options.h in synergia2
# archived-applications/booster
opts.add('wakefile_f', 'Fwake.dat', 'File containing F magnet wakes')
opts.add('wakefile_d', 'Dwake.dat', 'File containing D magnet wakes')
opts.add('waketype', "XLXTYLYTZpp", 'Type of wake file')
opts.add('registered_turns', 15, 'Number of turns for long range wakes')
opts.add('full_machine', False, 'Full machine populated if True')
opts.add('wave_number', [0, 0, 0], 'wave number for multibunch instability')
opts.add('imp_zgrid', 1000, 'The size of the zgrid for impedance')
opts.add('imp_compensate', 2, 'RF phase compensates energy loss from impedance')
opts.add("imp_comp_gain", 1.0, "gain to apply to restoring kick")

opts.add("set_xtune", 0.71, "adjust x tune", float)
opts.add("set_ytune", 0.81, "adjust y tune", float)
#opts.add("set_xtune", None, "adjust x tune", float)
#opts.add("set_ytune", None, "adjust y tune", float)

# chromaticity adjustments
#opts.add("set_xchrom", None, "adjust x chromaticity", float)
#opts.add("set_ychrom", None, "adjust y chromaticity", float)
opts.add("set_xchrom", -8.0, "adjust x chromaticity", float)
opts.add("set_ychrom", -8.0, "adjust y chromaticity", float)

# apertures
opts.add("trans_aperture", True, "apply transverse aperture")

opts.add("z_boundary", 'periodic', "longitudinal boundary condition: {open| periodic| aperture}  bucket_barrier}")

#opts.add("stepper", "splitoperator", "which stepper to use independent|elements|splitoperator")
opts.add("stepper", "elements", "which stepper to use independent|elements|splitoperator")
#opts.add("steps", 6*24, "# steps")
opts.add("steps", 1, "# steps")

opts.add("step_basic", False, "Basic diagnostics each step")
opts.add("tracks", 100, "number of particles to track")
opts.add("step_diag", 0, "diags per step")
opts.add("step_tracks", 0, "Number of particles to track by step")
#opts.add("save_particles", False, "if True, save particles")
opts.add("save_particles", True, "if True, save particles")
#opts.add("particles_period", 1, "0: save every turn, n!=0, save particles every n turns")
opts.add("particles_period", 25, "0: save every turn, n!=0, save particles every n turns")

job_mgr = synergia_workflow.Job_manager("booster.py", opts, ["booster-06800.json", "booster_rf_ramp.py", "Fwake.dat", "Dwake.dat"])
