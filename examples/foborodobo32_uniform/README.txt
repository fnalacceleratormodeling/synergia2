This is a customized version of foborodobo32 that populates a uniformly
distributed beam longitudinally.  The option "matching=uniform" gives
a gaussian momentum distribution.  The  option "matching=airbag" gives
a momentum distribution that is either + or - 3*stddpop.  Tha script
analyze.py reads the output h5 files and plots some representative
quantities. Note that the RF has been disabled by setting rf_volt to
0.  Turning on the RF will destroy the uniform distribution.
