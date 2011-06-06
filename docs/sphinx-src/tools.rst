Command-line tools
==================

 These are Python scripts based on the Synergia and CHEF libraries. While the user can link his own C++ code to
 these libraries, add his own analysis code and run it, it might be more convenient to use our existing scripts.  

synergia
--------

The synergia script is really just a frontend to Python. It is use for instance to convert MAD file 
to a Synergia lattice file used in a C++ code. Evidently, it can also be used to drive a complet Synergia
simulation. 

usage: synergia <synergia_script> [arguments]

	--help for help

	-i for interactive mode

	--ipython for interactive mode with ipython


synbeamplot
-----------

Plots one- and two-dimensional density projections of a beam.

usage: synbeamplot <filename> [option1] ... [optionn] <h coord> <v coord>

available options are:

    --nohist : do not show histograms (not on by default)

    --bins=<num> : number of bins in each direction

    --output=<file> : save output to file (not on by default)

    --show : show plots on screen (on by default unless --output flag is present

available coords are:

    x xp y yp z zp

syndiagplot
-----------

Plots one- or two-dimensional diagnostics vs trajectory length. Multiple plots
can be performed at once.

usage: syndiagplot <filename> [option1] ... [optionn] <plot1> ... <plotn>

available options are:

    --oneplot : put all plots on the same axis (not on by default)

    --output=<file> : save output to file (not on by default)

    --show : show plots on screen (on by default unless --output flag is present

available plots are:

    x_emit x_mean x_std x_xp_corr x_xp_mom2 x_y_corr x_y_mom2 x_yp_corr x_yp_mom2 x_z_corr x_z_mom2 x_zp_corr x_zp_mom2 xp_mean xp_std xp_y_corr xp_y_mom2 xp_yp_corr xp_yp_mom2 xp_z_corr xp_z_mom2 xp_zp_corr xp_zp_mom2 xy_emit xyz_emit y_emit y_mean y_std y_yp_corr y_yp_mom2 y_z_corr y_z_mom2 y_zp_corr y_zp_mom2 yp_mean yp_std yp_z_corr yp_z_mom2 yp_zp_corr yp_zp_mom2 z_emit z_mean z_std z_zp_corr z_zp_mom2 zp_mean zp_std

syntrackplot
------------

Plots particle tracks.

usage: syntrackplot <filename> [option1] ... [optionn] <h coord1> <v coord1> ... <h coordn> <v coordn>

available options are:

    --oneplot : put all plots on the same axis (not on by default)

    --output=<file> : save output to file (not on by default)

    --show : show plots on screen (on by default unless --output flag is present

available coords are:

    x xp y yp z zp

synpoincareplot
---------------

Performs Poincare plots of pairs of phase-space coordinates.

usage: synpoincareplot <filename1> ... <filenamen> [option1] ... [optionn] <h coord> <v coord>

available options are:

    --pointsize=<float>: size of plotted points (default=4.0)

    --output=<file> : save output to file (not on by default)

    --show : show plots on screen (on by default unless --output flag is present

available coords are:

    x xp y yp z zp

syninspecth5
------------

Inspects HDF5 files.

usage:

     syninspecth5 <hdf5_file>

         to list members, or

     syninspecth5 <hdf5_file> <member>

           to display member

synmad8toxml
------------

Converts Mad8 lattice files to Synergia XML files.

usage: synmad8toxml <mad8 file> <line name> <xml file>

    Reads line <line name> from <mad8 file> and writes to <xml file>.

