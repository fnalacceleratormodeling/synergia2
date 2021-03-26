synergia channel.py <options>

channel options:
----------------
RFVolt=<float>  RF Cavity voltage [MV] (try 0.005 or 10 for fun),
                default=1
octo=<float>    strength of octopole (try 0.05 for fun), default=0
offset=<float>  offset for test particles, default=0.001
order=<int>     normal form order, default=7
sext=<float>    strength of sextupole, default=0
skew=<float>    skew quad strength (try 0.005 for run), default=0
skocto=<float>  strength of skew octopole, default=0
sksext=<float>  strength of skew sextupole, default=0
turns=<int>     number of turns, default=1000

The script launches 80 particles in a FODO channel.  The particles are spaced
at uniform offset, 40 in the x direction and 40 in the y direction.  Optionally,
nonlinear sextupole, skew sextupole, octopole, skew octopole can be enabled
as well as a skew quadrupole for inter-plane mixing.  The RF cavity voltage
can also be controlled.

For each turn, the particle coordinates are converted to normal form and
written out to the file nf.dat with each containing the 6 normal form components of all 80 particles. Also bulk_track diagnostics are written to tracks.h5.

The script plot_normal_forms.py will read the nf.dat and display the particle
motion for a specific particle:

python3 plot_normal_forms.py <n>
