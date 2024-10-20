#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import synergia.simulation as SIM
ET = synergia.lattice.element_type
MT = synergia.lattice.marker_type
PCONST = synergia.foundation.pconstants
PGreg = synergia.bunch.ParticleGroup.regular
PGspect = synergia.bunch.ParticleGroup.spectator


#check_lattice = "BOOSTER"
check_lattice = "IOTA"
#check_lattice = "DIPEDGE"

import pytest

# Run parameters
macroparticles = 16
spectparticles = 4
#spectparticles = 0
realparticles = 285452129459.3449 # 0.1 mA
nturns = 5

synergia.simulation.Lattice_simulator.set_closed_orbit_tolerance(1.0e-12)

# prop_fixture is a propagator
#@pytest.fixture
def prop_fixture():
    # sample IOTA lattice
    iota_madx = """
kqa1r := kq01;
kq01 = -7.72652301;
kqa2r := kq02;
kq02 = 12.28401222;
kqa3r := kq03;
kq03 = -12.43016989;
kqa4r := kq04;
kq04 = 20.16074347;
kqb1r := kq05;
kq05 = -10.24365752;
kqb2r := kq06;
kq06 = 15.12808788;
kqb3r := kq07;
kq07 = -6.92311681;
kqb4r := kq08;
kq08 = -6.90057605;
kqb5r := kq09;
kq09 = 13.50655178;
kqb6r := kq10;
kq10 = -11.91343344;
kqc1r := kq11;
kq11 = -13.51948869;
kqc2r := kq12;
kq12 = 12.0339278;
kqc3r := kq13;
kq13 = -13.56878135;
kqd1r := kq14;
kq14 = -7.97007816;
kqd2r := kq15;
kq15 = 5.92322639;
kqd3r := kq16;
kq16 = -6.32915747;
kqd4r := kq17;
kq17 = 5.1636516;
kqe1r := kq18;
kq18 = -4.69712477;
kqe2r := kq19;
kq19 = 7.0326898;
kqe3 := kq20;
kq20 = -7.19881671;
kqe2l := kq19;
kqe1l := kq18;
kqd4l := kq17;
kqd3l := kq16;
kqd2l := kq15;
kqd1l := kq14;
kqc3l := kq13;
kqc2l := kq12;
kqc1l := kq11;
kqb6l := kq10;
kqb5l := kq09;
kqb4l := kq08;
kqb3l := kq07;
kqb2l := kq06;
kqb1l := kq05;
kqa4l := kq04;
kqa3l := kq03;
kqa2l := kq02;
kqa1l := kq01;
mseqari: marker;
ibpm: monitor;
ibpma1: ibpm;
qa1r: quadrupole,l:= 0.21,k1:=kqa1r ;
qa2r: quadrupole,l:= 0.21,k1:=kqa2r ;
sqa1r: quadrupole,l:= 0.1,k1s:= 0;
ibpma2r: ibpm;
qa3r: quadrupole,l:= 0.21,k1:=kqa3r ;
qa4r: quadrupole,l:= 0.21,k1:=kqa4r ;
ibpma3r: ibpm;
sqa2r: quadrupole,l:= 0.1,k1s:= 0;
mseqare: marker;
mphm1ri: marker;
dedge30: dipedge,e1:= 0,h:= 1.338646717,hgap:= 0.015042,fint:= 0.5;
m1r: sbend,l:= 0.3911403725,angle:= 0.5235987756;
mphm1re: marker;
mseqbri: marker;
sqb1r: quadrupole,l:= 0.1,k1s:= 0;
qb1r: quadrupole,l:= 0.21,k1:=kqb1r ;
qb2r: quadrupole,l:= 0.21,k1:=kqb2r ;
qb3r: quadrupole,l:= 0.21,k1:=kqb3r ;
ibpmb1r: ibpm;
nlr1: marker;
ior: marker;
nlr2: marker;
ibpmb2r: ibpm;
qb4r: quadrupole,l:= 0.21,k1:=kqb4r ;
qb5r: quadrupole,l:= 0.21,k1:=kqb5r ;
qb6r: quadrupole,l:= 0.21,k1:=kqb6r ;
sqb2r: quadrupole,l:= 0.1,k1s:= 0;
mseqbre: marker;
mphm2ri: marker;
dedge60: dipedge,e1:= 0,h:= 1.381554029,hgap:= 0.014786,fint:= 0.5;
m2r: sbend,l:= 0.757985232,angle:= 1.047197551;
mphm2re: marker;
mseqcri: marker;
sqc1r: quadrupole,l:= 0.1,k1s:= 0;
ibpmc1r: ibpm;
qc1r: quadrupole,l:= 0.21,k1:=kqc1r ;
qc2r: quadrupole,l:= 0.21,k1:=kqc2r ;
qc3r: quadrupole,l:= 0.21,k1:=kqc3r ;
ibpmc2r: ibpm;
sqc2r: quadrupole,l:= 0.1,k1s:= 0;
mseqcre: marker;
mphm3ri: marker;
m3r: sbend,l:= 0.757985232,angle:= 1.047197551;
mphm3re: marker;
mseqdri: marker;
ibpmd1r: ibpm;
sqd1r: quadrupole,l:= 0.1,k1s:= 0;
qd1r: quadrupole,l:= 0.21,k1:=kqd1r ;
qd2r: quadrupole,l:= 0.21,k1:=kqd2r ;
el1: marker;
cel: solenoid,l:= 0.7,ks:= 0;
el2: marker;
qd3r: quadrupole,l:= 0.21,k1:=kqd3r ;
sqd2r: quadrupole,l:= 0.1,k1s:= 0;
qd4r: quadrupole,l:= 0.21,k1:=kqd4r ;
ibpmd2r: ibpm;
mseqdre: marker;
mphm4ri: marker;
m4r: sbend,l:= 0.3911403725,angle:= 0.5235987756;
mphm4re: marker;
mseqei: marker;
ibpme1r: ibpm;
qe1r: quadrupole,l:= 0.21,k1:=kqe1r ;
sqe1r: quadrupole,l:= 0.1,k1s:= 0;
ibpme2r: ibpm;
qe2r: quadrupole,l:= 0.21,k1:=kqe2r ;
sqe2r: quadrupole,l:= 0.1,k1s:= 0;
qe3: quadrupole,l:= 0.21,k1:=kqe3 ;
sqe2l: quadrupole,l:= 0.1,k1s:= 0;
qe2l: quadrupole,l:= 0.21,k1:=kqe2l ;
ibpme2l: ibpm;
sqe1l: quadrupole,l:= 0.1,k1s:= 0;
qe1l: quadrupole,l:= 0.21,k1:=kqe1l ;
ibpme1l: ibpm;
mseqee: marker;
mphm4li: marker;
m4l: sbend,l:= 0.3911403725,angle:= 0.5235987756;
mphm4le: marker;
mseqdli: marker;
ibpmd2l: ibpm;
qd4l: quadrupole,l:= 0.21,k1:=kqd4l ;
sqd2l: quadrupole,l:= 0.1,k1s:= 0;
qd3l: quadrupole,l:= 0.21,k1:=kqd3l ;
rfc: rfcavity,l:= 0.05,volt:= 0.000847,lag:= 0,harmon:= 4;
qd2l: quadrupole,l:= 0.21,k1:=kqd2l ;
qd1l: quadrupole,l:= 0.21,k1:=kqd1l ;
sqd1l: quadrupole,l:= 0.1,k1s:= 0;
ibpmd1l: ibpm;
mseqdle: marker;
mphm3li: marker;
m3l: sbend,l:= 0.757985232,angle:= 1.047197551;
mphm3le: marker;
mseqcli: marker;
sqc2l: quadrupole,l:= 0.1,k1s:= 0;
ibpmc2l: ibpm;
qc3l: quadrupole,l:= 0.21,k1:=kqc3l ;
qc2l: quadrupole,l:= 0.21,k1:=kqc2l ;
qc1l: quadrupole,l:= 0.21,k1:=kqc1l ;
ibpmc1l: ibpm;
sqc1l: quadrupole,l:= 0.1,k1s:= 0;
mseqcle: marker;
mphm2li: marker;
m2l: sbend,l:= 0.757985232,angle:= 1.047197551;
mphm2le: marker;
mseqbli: marker;
sqb2l: quadrupole,l:= 0.1,k1s:= 0;
qb6l: quadrupole,l:= 0.21,k1:=kqb6l ;
qb5l: quadrupole,l:= 0.21,k1:=kqb5l ;
qb4l: quadrupole,l:= 0.21,k1:=kqb4l ;
ibpmb2l: ibpm;
nll1: marker;
nll2: marker;
ibpmb1l: ibpm;
qb3l: quadrupole,l:= 0.21,k1:=kqb3l ;
qb2l: quadrupole,l:= 0.21,k1:=kqb2l ;
qb1l: quadrupole,l:= 0.21,k1:=kqb1l ;
sqb1l: quadrupole,l:= 0.1,k1s:= 0;
mseqble: marker;
mphm1li: marker;
m1l: sbend,l:= 0.3911403725,angle:= 0.5235987756;
mphm1le: marker;
mseqali: marker;
sqa2l: quadrupole,l:= 0.1,k1s:= 0;
ibpma3l: ibpm;
qa4l: quadrupole,l:= 0.21,k1:=kqa4l ;
qa3l: quadrupole,l:= 0.21,k1:=kqa3l ;
ibpma2l: ibpm;
sqa1l: quadrupole,l:= 0.1,k1s:= 0;
qa2l: quadrupole,l:= 0.21,k1:=kqa2l ;
qa1l: quadrupole,l:= 0.21,k1:=kqa1l ;
mseqale: marker;
!iota: sequence, l = 39.95567226;
iota: sequence, l = 3.117400202;
mseqari, at = 0;
ibpma1, at = 0.02;
qa1r, at = 1.0175;
qa2r, at = 1.3625;
sqa1r, at = 1.99;
ibpma2r, at = 2.095;
qa3r, at = 2.2975;
qa4r, at = 2.6525;
ibpma3r, at = 2.865;
sqa2r, at = 2.97;
mseqare, at = 3.0405;
mphm1ri, at = 3.0405;
dedge30, at = 3.117400202;
endsequence;
"""

    foo_madx = """
m1r, at = 3.312970389;
dedge30, at = 3.508540575;
mphm1re, at = 3.585440777;
mseqbri, at = 3.585440777;
sqb1r, at = 3.715940777;
qb1r, at = 3.953440777;
qb2r, at = 4.303440777;
qb3r, at = 4.653440777;
ibpmb1r, at = 4.865940777;
nlr1, at = 4.910940777;
ior, at = 5.810940777;
nlr2, at = 6.710940777;
ibpmb2r, at = 6.755940777;
qb4r, at = 6.968440777;
qb5r, at = 7.318440777;
qb6r, at = 7.668440777;
sqb2r, at = 7.905940777;
mseqbre, at = 8.036440777;
mphm2ri, at = 8.036440777;
dedge60, at = 8.112186805;
m2r, at = 8.491179421;
dedge60, at = 8.870172037;
mphm2re, at = 8.945918065;
mseqcri, at = 8.945918065;
sqc1r, at = 9.106418065;
ibpmc1r, at = 9.211418065;
qc1r, at = 9.423918065;
qc2r, at = 9.988918065;
qc3r, at = 10.55391806;
ibpmc2r, at = 10.76641806;
sqc2r, at = 10.87141806;
mseqcre, at = 11.03191806;
mphm3ri, at = 11.03191806;
dedge60, at = 11.10766409;
m3r, at = 11.48665671;
dedge60, at = 11.86564932;
mphm3re, at = 11.94139535;
mseqdri, at = 11.94139535;
ibpmd1r, at = 12.20689535;
sqd1r, at = 12.39189535;
qd1r, at = 12.61939535;
qd2r, at = 13.24939535;
el1, at = 13.81689535;
cel, at = 14.16689535;
el2, at = 14.51689535;
qd3r, at = 15.08439535;
sqd2r, at = 15.39939535;
qd4r, at = 15.71439535;
ibpmd2r, at = 16.12689535;
mseqdre, at = 16.39239535;
mphm4ri, at = 16.39239535;
dedge30, at = 16.46929555;
m4r, at = 16.66486574;
dedge30, at = 16.86043593;
mphm4re, at = 16.93733613;
mseqei, at = 16.93733613;
ibpme1r, at = 17.20283613;
qe1r, at = 17.51533613;
sqe1r, at = 17.76283613;
ibpme2r, at = 18.74783613;
qe2r, at = 18.96033613;
sqe2r, at = 19.34783613;
qe3, at = 19.97783613;
sqe2l, at = 20.60783613;
qe2l, at = 20.99533613;
ibpme2l, at = 21.20783613;
sqe1l, at = 22.19283613;
qe1l, at = 22.44033613;
ibpme1l, at = 22.75283613;
mseqee, at = 23.01833613;
mphm4li, at = 23.01833613;
dedge30, at = 23.09523633;
m4l, at = 23.29080652;
dedge30, at = 23.4863767;
mphm4le, at = 23.56327691;
mseqdli, at = 23.56327691;
ibpmd2l, at = 23.82877691;
qd4l, at = 24.24127691;
sqd2l, at = 24.55627691;
qd3l, at = 24.87127691;
rfc, at = 25.78877691;
qd2l, at = 26.70627691;
qd1l, at = 27.33627691;
sqd1l, at = 27.56377691;
ibpmd1l, at = 27.74877691;
mseqdle, at = 28.01427691;
mphm3li, at = 28.01427691;
dedge60, at = 28.09002293;
m3l, at = 28.46901555;
dedge60, at = 28.84800817;
mphm3le, at = 28.92375419;
mseqcli, at = 28.92375419;
sqc2l, at = 29.08425419;
ibpmc2l, at = 29.18925419;
qc3l, at = 29.40175419;
qc2l, at = 29.96675419;
qc1l, at = 30.53175419;
ibpmc1l, at = 30.74425419;
sqc1l, at = 30.84925419;
mseqcle, at = 31.00975419;
mphm2li, at = 31.00975419;
dedge60, at = 31.08550022;
m2l, at = 31.46449284;
dedge60, at = 31.84348545;
mphm2le, at = 31.91923148;
mseqbli, at = 31.91923148;
sqb2l, at = 32.04973148;
qb6l, at = 32.28723148;
qb5l, at = 32.63723148;
qb4l, at = 32.98723148;
ibpmb2l, at = 33.19973148;
nll1, at = 33.24473148;
nll2, at = 35.04473148;
ibpmb1l, at = 35.08973148;
qb3l, at = 35.30223148;
qb2l, at = 35.65223148;
qb1l, at = 36.00223148;
sqb1l, at = 36.23973148;
mseqble, at = 36.37023148;
mphm1li, at = 36.37023148;
dedge30, at = 36.44713168;
m1l, at = 36.64270187;
dedge30, at = 36.83827206;
mphm1le, at = 36.91517226;
mseqali, at = 36.91517226;
sqa2l, at = 36.98567226;
ibpma3l, at = 37.09067226;
qa4l, at = 37.30317226;
qa3l, at = 37.65817226;
ibpma2l, at = 37.86067226;
sqa1l, at = 37.96567226;
qa2l, at = 38.59317226;
qa1l, at = 38.93817226;
mseqale, at = 39.95567226;
endsequence;
"""

    booster_madx = """
! Simplified Booster lattice

// From JFO 2022-12-08
// The simplified booster lattice is trivial. I have no madx lattice (you could make one
// very easily) - I just use pyorbit classes directly to instantiate a basic cell;  it is then
// replicated 24 times. I took the bending magnets lengths and 
// strengths directly from the official MADX  lattice file. 

// The basic cell is 

// d1 fmag d2 dmag d3 

// d1, d2, d3 : drifts of lengths 0.6 0.5 and 3.0 m
// fmag:  focusing      bend   L = 2.889612 m 
// dmag   defocusing bend   L = 2.889612 m

// total cell length: 19.758 m
// total ring  length = 24*19.758 = 474.20 m 

// The length, focusing strengths and curvature radius of the 
// magnets are as in the booster MADX file.  

// If you entered 1 cell correctly, you should get the periodic solution:
// bx = 33.86 m ax = 0
// by = 5.39m    ay =0 
// For 24 cells, the raw tunes are nux = 7.017 and nuy = 6.675.  You will need to tweak the nominal focusing strengths a bit to avoid resonances.


//--------- Nominal Gradient Magnet Definitions  

// EGS
// The apparent cell structure is actually:

// D1, l=0.6;
// FMAGU01;
// D2, l=0.5;
// DMAGU01;
// D3, l=6.0;
// DMAGD01;
// D4, l=0.5;
// FMAGD01;
// DR, l=0.6


! Expand structure to include corrector packages:

! Normal short straight

! drift 0.6 (end of previous cell)
! drift 0.176
! corrector package (short) l=0.168
! drift 0.256
! end of short straight

! non-RF cavity long straight
!
! drift 5.581
! correction package (long) l=0.168
! drift 0.251
!

! long straight with RF cavity
!
! drift 0.21
! rfcavity drift-cavity-drift 2.35
! drift 0.12
! rfcavity drift-cavity-drift 2.35
! 4 drifts total 0.551
! corrector packages (long) 0.168
! drift 0.251

! Corrector package:
!
! HLxx, HKICKER, l=0.024
! VLxx, VKICKER, l=0.024
! Q{S|L}xx QUADRUPOLE, l=0.024 (normal)
! MULTIPOLE Q{S|L}ERR, l=0
! MONITOR, l=0.024
! QS{S|L}xx, QUADRUPOLE, l=0.024 (skew)
! MULTIPOLE QS{S|L}ERR, l=0
! SEXTUPOLE SX{L|S}, l=0.024 (normal)
! SEXTUPOLE SS{L|S}, l=0.024 (skew)
! 

ke1 = 0.8;  !800 MeV kinetic energy at injection

rhof  :=  40.847086;   !  bending radius of focusing magnet
rhod  :=  48.034101;   !  bending radius of defocusing magnet

blength :=     2.889612;    !  arc length for both F and D magnets
blengthf :=    2.889009499; !  physical length (straight length) for F magnet
blengthd :=    2.889176299; !  physical length (straight length) for D magnet


!
! The quad field for the gradient magnet is contained in file " qsdqsf.dat" to be read in before this file !
!
! read from file at time step = 7
qsd := -57.38855012e-3;
qsf := 54.10921561e-3;


! These ssd and ssf strengths come from fitting to 01 Dec 2015 chromaticity data
! and predicts chromaticity much better than using the Drozhdin et al measurements above

ssd :=  -0.04381647074 + ke1*(0.009150934932+ ke1*(-0.0023900895  + ke1*(0.000318068028 -  ke1* 1.6353205e-05)));

ssf :=  -0.006384940088 + ke1*(0.01967542848 + ke1*( -0.006776746 + ke1*(0.00091367565 - ke1* 4.293705e-05)));

 !
 ! Gradient magnets defined by their physical length aith their bend angle
 ! being defined by the arc length/radius of curvature

!FMAG: RBEND,  L = blengthf  , ANGLE = blength/rhof, K1 = qsf  , K2 = ssf;
FMAG: SBEND,  L = blength  , ANGLE = blength/rhof, e1=blength/(2*rhof), e2=blength/(2*rhof), K1 = qsf  , K2 = ssf, type=fmag;
DMAG: SBEND,  L = blength  , ANGLE = blength/rhod, e1=blength/(2*rhod), e2=blength/(2*rhod), K1 = qsd  , K2 = ssd, type=dmag;

! drifts in the short straight section
mins: drift, l=0.5, type=shortstraight;
sc: drift, l=0.6, type=shortstraight;
sa: drift, l=0.176, type=shortstraight;
sb: drift, l=0.256, type=shortstraight;

! drifts in the long straight section
dlong: drift, l=5.581, type=longstraight; ! for no-RF cavity cells
drifta: drift, l=0.21, type=longstraight; ! start of RF cavity cell
driftb: drift, l=0.12, type=longstraight;
dmidls: drift, l=0.551, type=longstraight; ! (drift in the middle of the longstraight)
drifte: drift, l=0.251, type=longstraight; ! end of RF cavity cell

drrf: drift, l=2.35/2, type=rfaperture;
rfc: rfcavity, l=0, harmon=84, volt=0, lag=0, type=rfaperture;

! short corrector package
hsxx: drift, l=0.024, type=shortstraight;
vsxx: drift, l=0.024, type=shortstraight;
qsxx: quadrupole, k1=0.0, l=0.024, type=shortstraight;
! multipole not included
bpms: drift, l=0.024, type=shortstraight;
qssxx: drift, l=0.024, type=shortstraight; // skew quad
! multipole not included
sxsxx: sextupole, k2=0.0, l=0.024, type=shortstraight;
sssxx: drift, l=0.024, type=shortstraight; // skew sextupole

cpshort: line=(hsxx, vsxx, qsxx, bpms, qssxx, sxsxx, sssxx);

! long corrector package
hlxx: drift, l=0.024, type=longstraight;
vlxx: drift, l=0.024, type=longstraight;;
qlxx: quadrupole, k1=0.0, l=0.024, type=shortstraight;
! multipole not included
bpml: drift, l=0.024, type=shortstraight;
qlsxx: drift, l=0.024, type=longstraight; // skew quad
! multipole not included
sxlxx: sextupole, k2=0.0, l=0.024, type=longstraight;
sslxx: drift, l=0.024, type=longstraight; // skew sextupole

cplong: line=(hlxx, vlxx, qlxx, bpms, qlsxx, sxlxx, sslxx);

!!!!!!!!!!!!!!!!   beginning of ring definition

fmagu01: fmag;
fmagd01: fmag;
dmagu01: dmag;
dmagd01: dmag;

cell01 : line = (sa, cpshort, sb, fmagu01, mins, dmagu01, dlong, cplong, drifte, dmagd01, mins, fmagd01, sc);

fmagu02: fmag;
fmagd02: fmag;
dmagu02: dmag;
dmagd02: dmag;

cell02 : line = (sa, cpshort, sb, fmagu02, mins, dmagu02, dlong, cplong, drifte, dmagd02, mins, fmagd02, sc);

fmagu03: fmag;
fmagd03: fmag;
dmagu03: dmag;
dmagd03: dmag;

cell03 : line = (sa, cpshort, sb, fmagu03, mins, dmagu03, dlong, cplong, drifte, dmagd03, mins, fmagd03, sc);

fmagu04: fmag;
fmagd04: fmag;
dmagu04: dmag;
dmagd04: dmag;

cell04 : line = (sa, cpshort, sb, fmagu04, mins, dmagu04, dlong, cplong, drifte, dmagd04, mins, fmagd04, sc);

fmagu05: fmag;
fmagd05: fmag;
dmagu05: dmag;
dmagd05: dmag;

cell05 : line = (sa, cpshort, sb, fmagu05, mins, dmagu05, dlong, cplong, drifte, dmagd05, mins, fmagd05, sc);

fmagu06: fmag;
fmagd06: fmag;
dmagu06: dmag;
dmagd06: dmag;

cell06 : line = (sa, cpshort, sb, fmagu06, mins, dmagu06, dlong, cplong, drifte, dmagd06, mins, fmagd06, sc);

fmagu07: fmag;
fmagd07: fmag;
dmagu07: dmag;
dmagd07: dmag;

cell07 : line = (sa, cpshort, sb, fmagu07, mins, dmagu07, dlong, cplong, drifte, dmagd07, mins, fmagd07, sc);

fmagu08: fmag;
fmagd08: fmag;
dmagu08: dmag;
dmagd08: dmag;

cell08 : line = (sa, cpshort, sb, fmagu08, mins, dmagu08, dlong, cplong, drifte, dmagd08, mins, fmagd08, sc);
 
fmagu09: fmag;
fmagd09: fmag;
dmagu09: dmag;
dmagd09: dmag;

cell09 : line = (sa, cpshort, sb, fmagu09, mins, dmagu09, dlong, cplong, drifte, dmagd09, mins, fmagd09, sc);

fmagu10: fmag;
fmagd10: fmag;
dmagu10: dmag;
dmagd10: dmag;

cell10 : line = (sa, cpshort, sb, fmagu10, mins, dmagu10, dlong, cplong, drifte, dmagd10, mins, fmagd10, sc);

fmagu11: fmag;
fmagd11: fmag;
dmagu11: dmag;
dmagd11: dmag;

cell11 : line = (sa, cpshort, sb, fmagu11, mins, dmagu11, dlong, cplong, drifte, dmagd11, mins, fmagd11, sc);

fmagu12: fmag;
fmagd12: fmag;
dmagu12: dmag;
dmagd12: dmag;

cell12 : line = (sa, cpshort, sb, fmagu12, mins, dmagu12, dlong, cplong, drifte, dmagd12, mins, fmagd12, sc);

fmagu13: fmag;
fmagd13: fmag;
dmagu13: dmag;
dmagd13: dmag;

cell13 : line = (sa, cpshort, sb, fmagu13, mins, dmagu13, dlong, cplong, drifte, dmagd13, mins, fmagd13, sc);

fmagu14: fmag;
fmagd14: fmag;
dmagu14: dmag;
dmagd14: dmag;
rf01: line=(drrf, rfc, drrf);
rf02: line=(drrf, rfc, drrf);

cell14: line = (sa, cpshort, sb, fmagu14, mins, dmagu14, drifta, rf01, driftb, rf02, dmidls, cplong, drifte, dmagd14, mins, fmagd14, sc);


fmagu15: fmag;
fmagd15: fmag;
dmagu15: dmag;
dmagd15: dmag;
rf03: line=(drrf, rfc, drrf);
rf04: line=(drrf, rfc, drrf);

cell15: line = (sa, cpshort, sb, fmagu15, mins, dmagu15, drifta, rf03, driftb, rf04, dmidls, cplong, drifte, dmagd15, mins, fmagd15, sc);

fmagu16: fmag;
fmagd16: fmag;
dmagu16: dmag;
dmagd16: dmag;
rf05: line=(drrf, rfc, drrf);
rf06: line=(drrf, rfc, drrf);

cell16: line = (sa, cpshort, sb, fmagu16, mins, dmagu16, drifta, rf05, driftb, rf06, dmidls, cplong, drifte, dmagd16, mins, fmagd16, sc);

fmagu17: fmag;
fmagd17: fmag;
dmagu17: dmag;
dmagd17: dmag;
rf07: line=(drrf, rfc, drrf);
rf08: line=(drrf, rfc, drrf);

cell17: line = (sa, cpshort, sb, fmagu17, mins, dmagu17, drifta, rf07, driftb, rf08, dmidls, cplong, drifte, dmagd17, mins, fmagd17, sc);

fmagu18: fmag;
fmagd18: fmag;
dmagu18: dmag;
dmagd18: dmag;
rf09: line=(drrf, rfc, drrf);
rf10: line=(drrf, rfc, drrf);

cell18: line = (sa, cpshort, sb, fmagu18, mins, dmagu18, drifta, rf09, driftb, rf10, dmidls, cplong, drifte, dmagd18, mins, fmagd18, sc);

fmagu19: fmag;
fmagd19: fmag;
dmagu19: dmag;
dmagd19: dmag;
rf11: line=(drrf, rfc, drrf);
rf12: line=(drrf, rfc, drrf);

cell19: line = (sa, cpshort, sb, fmagu19, mins, dmagu19, drifta, rf11, driftb, rf12, dmidls, cplong, drifte, dmagd19, mins, fmagd19, sc);

fmagu20: fmag;
fmagd20: fmag;
dmagu20: dmag;
dmagd20: dmag;
rf13: line=(drrf, rfc, drrf);
rf14: line=(drrf, rfc, drrf);

cell20: line = (sa, cpshort, sb, fmagu20, mins, dmagu20, drifta, rf13, driftb, rf14, dmidls, cplong, drifte, dmagd20, mins, fmagd20, sc);

fmagu21: fmag;
fmagd21: fmag;
dmagu21: dmag;
dmagd21: dmag;
rf15: line=(drrf, rfc, drrf);
rf16: line=(drrf, rfc, drrf);

cell21: line = (sa, cpshort, sb, fmagu21, mins, dmagu21, drifta, rf15, driftb, rf16, dmidls, cplong, drifte, dmagd21, mins, fmagd21, sc);

fmagu22: fmag;
fmagd22: fmag;
dmagu22: dmag;
dmagd22: dmag;
rf17: line=(drrf, rfc, drrf);
rf18: line=(drrf, rfc, drrf);

cell22: line = (sa, cpshort, sb, fmagu22, mins, dmagu22, drifta, rf17, driftb, rf18, dmidls, cplong, drifte, dmagd22, mins, fmagd22, sc);

fmagu23: fmag;
fmagd23: fmag;
dmagu23: dmag;
dmagd23: dmag;
rf19: line=(drrf, rfc, drrf);
rf20: line=(drrf, rfc, drrf);

cell23: line = (sa, cpshort, sb, fmagu23, mins, dmagu23, drifta, rf19, driftb, rf20, dmidls, cplong, drifte, dmagd23, mins, fmagd23, sc);

fmagu24: fmag;
fmagd24: fmag;
dmagu24: dmag;
dmagd24: dmag;
rf21: line=(drrf, rfc, drrf);
rf22: line=(drrf, rfc, drrf);

cell24: line = (sa, cpshort, sb, fmagu24, mins, dmagu24, drifta, rf21, driftb, rf22, dmidls, cplong, drifte, dmagd24, mins, fmagd24, sc);

booster: line=(cell01, cell02, cell03, cell04, cell05, cell06, cell07, cell08, cell09, cell10, cell11, cell12, cell13, cell14, cell15, cell16, cell17, cell18, cell19, cell20, cell21, cell22, cell23, cell24);

beam, particle=proton, energy=0.8+pmass;
!beam, particle=proton, energy=1.7386206689709858;
!beam, particle=proton, energy=0.8003486229709857+pmass;
"""

    dipedge_madx = """
qa1r: quadrupole,l:= 0.21,k1:=kqa1r ;
dedge30: dipedge,e1:= 0,h:= 1.338646717,hgap:= 0.015042,fint:= 0.5;
m1r: sbend,l:= 0.3911403725,angle:= 0.5235987756;
dipedge: sequence, l=0.3911403725+0.21;
qa1r, at = 0.21/2;
dedge30, at = 0.21;
m1r, at=0.3911403725/2+0.21;
dedge30, at=0.3911403725+0.21;
endsequence;
"""

    reader = synergia.lattice.MadX_reader()
    if check_lattice == "IOTA":
        reader.parse(iota_madx)
        lattice = reader.get_lattice("iota")
        KE = 0.00250 # kinetic energy 2.5 MeV
    elif check_lattice == "DIPEDGE":
        KE = 0.00250 # kinetic energy 2.5 MeV
        reader.parse(dipedge_madx)
        lattice = reader.get_lattice("dipedge")
        KE = 0.00250 # kinetic energy 2.5 MeV
    elif check_lattice == "BOOSTER":
        KE = 0.8
        reader.parse(booster_madx)
        lattice = reader.get_lattice("booster")
    else:
        raise RuntimeError("unknown lattice")

    refpart = synergia.foundation.Reference_particle(1, PCONST.mp, PCONST.mp+KE)
    lattice.set_reference_particle(refpart)

    f = open('lattice.csv', 'w')
    s = 0.0
    for elem in lattice.get_elements():
        print(f'{s}, {elem.get_length()}, {elem.get_name()}, {elem.get_type_name()}', file=f)
        s += elem.get_length()
    f.close()

    synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice)
    stepper = synergia.simulation.Independent_stepper_elements(1)

    propagator = synergia.simulation.Propagator(lattice, stepper)

    return propagator

def create_simulator(refpart):
    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(
        refpart, macroparticles, realparticles, synergia.utils.Commxx(), spectparticles)

    bunch = sim.get_bunch()
    bunch.checkout_particles()
    lp = bunch.get_particles_numpy()
    lp[:, 0:6] = 0.0
    sp = bunch.get_particles_numpy(PGspect)
    sp[:, 0:6] = 0.0
    bunch.checkin_particles()

    return sim

def register_diagnostics(sim):
    diag_stat = synergia.bunch.Diagnostics_full2("diag.h5")
    sim.reg_diag_per_step(diag_stat)

    diag_tracks = synergia.bunch.Diagnostics_bulk_track("tracks.h5", 16, 0,
                                                        synergia.bunch.ParticleGroup.regular )

    sim.reg_diag_per_step(diag_tracks)


    if spectparticles:
        spdiag_tracks = synergia.bunch.Diagnostics_bulk_track("sptracks.h5", 4, 0,
                                                              synergia.bunch.ParticleGroup.spectator )

        sim.reg_diag_per_step(spdiag_tracks)

    #diag_part = synergia.bunch.Diagnostics_particles("particles.h5",
    #                                                     macroparticles,
    #                                                     0, spectparticles, 0)
    #sim.reg_diag_per_step(diag_part)


def test_iota_propagation(prop_fixture):
    refpart = prop_fixture.get_lattice().get_reference_particle()
    sim = create_simulator(prop_fixture.get_lattice().get_reference_particle())
    print('num particles: ', sim.get_bunch().get_total_num())
    print('spectator particle: ', sim.get_bunch().get_total_num(PGspect))

    register_diagnostics(sim)

    simlog = synergia.utils.parallel_utils.Logger(
        0, synergia.utils.parallel_utils.LoggerV.INFO_STEP, True
    )

    # turn_end_action will check coordinates of particle 0 of regular
    # particles with spectator particle 0 which should be identical

    def step_end_action(sim, lattice_in, turn, step):
        bunch =  sim.get_bunch(0, 0)
        bunch.checkout_particles()
        preg = bunch.get_particles_numpy(PGreg)
        pspect = bunch.get_particles_numpy(PGspect)
        assert preg[0, 0] == pspect[0, 0]
        assert preg[0, 1] == pspect[0, 1]
        assert preg[0, 2] == pspect[0, 2]
        assert preg[0, 3] == pspect[0, 3]
        assert preg[0, 4] == pspect[0, 4]
        assert preg[0, 5] == pspect[0, 5]
        pass # end of step_end_action

    sim.reg_prop_action_step_end(step_end_action)

    prop_fixture.propagate(sim, simlog, nturns)
    print('final num particles: ', sim.get_bunch().get_total_num())
    assert sim.get_bunch().get_total_num() != 0
    
def main():
    pf = prop_fixture()
    test_iota_propagation(pf)

if __name__ == "__main__":
    main()
