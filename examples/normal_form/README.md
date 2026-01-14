# Normal form and mapping examples

1. `normal_form.cc` Propagate particles through a nonlinear lattice, write out normal form coordinates.

2. `normal_form.py` Propagate particles through a nonlinear lattice, write out normal form coordinates in Python.

3. `channel_map.cc` Read in lattice file `channel.madx`. Print out nonlinear terms in the one-turn map using trigon.

4. `channel_map.py` Read in lattice file `channel.madx`. Print out nonlinear terms in the one-turn map using trigon.

5. `booster_map.py` Read in a simplified booster lattice in json form from file `booster_init_lattice.json`. Print out nonliner terms in the one-turn map using trigon.

6. `run_normal_form.py` A module containing two routines:
    - `run(qks=0.0, volt=0.2, sxn=0.0, sxs=0.0, ocn=0.0, ocs=0.0)`
    Run with skew quadrupole strength `qks`\
    Run with RF voltage `volt` [MV]\
    Run with normal sextupole strength `sxn`\
    Run with skew sextupole strengh `sxs`\
    Run with normal octupole strength `ocn`\
    Run with skew octupole strength `ocs`

    This writes out files tracks.h5 with particle trajectories in
    xyz coordinates and `nf.dat` with particle trajectories
    in normal form coordinates.

    - `plot_traj(tracks.h5, nf.dat)` Plot the trajectories of particles in
    xyz coordinates as well as  normal form coordinates.
    
    Example:

    ```
    import run_normal_form
    # run with skew quad and normal sextupole
    run_normal_form.run(qks=0.001, sxn=0.5)
    # plot results
    run_normal_form.plot_traj(tracks.h5 nf.dat)
    ```

8. `run_normal_form.ipynb` ipython notebook demonstrating same functionality
   as run_normal_form.py with example runs.
   