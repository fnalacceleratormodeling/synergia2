# IOTA example

Example using the `t1_iota_8_4.madx` lattice.

Run `iota_map.madx` file with madx. This writes out file `twiss.out` containing the one turn map.  Run the script `read_twiss_out.py` to read this file and print the map.

egstern@egstern1:~/fnal/syn2-devel3/src/synergia2/examples/Iota8-4$ python read_twiss_out.py
re: 
```
[[-3.09050448e-01  5.37734262e-01  0.00000000e+00  0.00000000e+00 -2.15266855e-06 -2.51906621e-01]
 [-1.68203362e+00 -3.09051393e-01  0.00000000e+00  0.00000000e+00 -1.08923466e-06 -3.23681972e-01]
 [ 0.00000000e+00  0.00000000e+00 -3.08994870e-01  1.30511335e+00  0.00000000e+00  0.00000000e+00]
 [ 0.00000000e+00  0.00000000e+00 -6.93060238e-01 -3.08994874e-01  0.00000000e+00  0.00000000e+00]
 [ 3.23684240e-01  2.51907619e-01  0.00000000e+00  0.00000000e+00  9.99995987e-01 -1.23484536e+00]
 [-1.93192400e-06 -1.98563244e-07  0.00000000e+00  0.00000000e+00  4.17788699e-06  9.99998432e-01]]
```

Run `print_iota_map.py` file with Synergia3.

```
egstern@egstern1:~/fnal/syn2-devel3/src/synergia2/examples/Iota8-4$ python print_iota_map.py
Kokkos::OpenMP::initialize WARNING: OMP_PROC_BIND environment variable not set
  In general, for best performance with OpenMP 4.0 or better set OMP_PROC_BIND=spread and OMP_PLACES=threads
  For best performance with OpenMP 3.1 set OMP_PROC_BIND=true
  For unit testing set OMP_PROC_BIND=false
energy:  0.15051100599422132
momentum:  0.15051013854716935
gamma:  294.542704
beta:  0.9999942366536836
map
[[-3.09503100e-01  5.37616661e-01  0.00000000e+00  0.00000000e+00  1.84770477e-06 -2.51922375e-01]
 [-1.68188094e+00 -3.09503910e-01  0.00000000e+00  0.00000000e+00  9.34298204e-07 -3.23560883e-01]
 [ 0.00000000e+00  0.00000000e+00 -3.09205480e-01  1.30511158e+00  0.00000000e+00  0.00000000e+00]
 [ 0.00000000e+00  0.00000000e+00 -6.92961418e-01 -3.09205484e-01  0.00000000e+00  0.00000000e+00]
 [-3.23564693e-01 -2.51924684e-01  0.00000000e+00  0.00000000e+00  9.99996555e-01  1.23476071e+00]
 [-1.65808361e-06 -1.70703557e-07  0.00000000e+00  0.00000000e+00 -3.58626770e-06  9.99998654e-01]]
tunex:  0.3000814973342238
tuney:  0.30003154359574474

chromaticity (x):  4.071262208198145
chromaticity (y):  1.7035287745649956

```