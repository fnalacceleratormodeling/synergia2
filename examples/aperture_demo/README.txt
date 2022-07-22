Script demonstrating aperture operation

aperture_demo options:
----------------------
aperture=<str>                       Aperture type (either 'circular',
                                     'elliptical', 'rectangular', or
                                     'polygon', default=elliptical
circular_aperture_radius=<float>     circular aperture radius,
                                     default=0.005
elliptical_horizontal_radius=<float> elliptical horizontal radius,
                                     default=0.005
elliptical_vertical_radius=<float>   elliptical vertical radius,
                                     default=0.002
hoffset=<float>                      horizontal aperture offset [m],
                                     default=0.0025
macroparticles=<int>                 Number of macro particles,
                                     default=20000
polygon_aperture_arm=<float>         length arm of star aperture,
                                     default=0.005
real_particles=<float>               Number of real particles,
                                     default=1.2e+12
rectangular_aperture_height=<float>  rectangular aperture height,
                                     default=0.004
rectangular_aperture_width=<float>   rectangular aperture width,
                                     default=0.01
seed=<int>                           Pseudorandom number generator
                                     seed, default=0
verbose=<bool>                       Verbose propagation, default=True
voffset=<float>                      vertical aperture offset [m],


Running the script transports a beam distribution through a channel
with the user specified aperture mask.  There are 4 particle dump files. particles_0000.h5 is before
the channel, particles_0001.h5 and later are after.

You can view before and after with beam_plot. After paths and environment variables are set up:

python3 aperture_demo.py aperture=polygon


python3 beam_plot.py particles_0000.h5 --bins=40 --minh=-0.024 --maxh=0.024 --minv=-0.012 --maxv=0.012 x y &
python3 beam_plot.py particles_0001.h5 --bins=40 --minh=-0.024 --maxh=0.024 --minv=-0.012 --maxv=0.012 x y &

will show before and after distributions.