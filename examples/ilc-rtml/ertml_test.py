#!/usr/bin/env bwpython

import local_paths
import gourmet
import Numeric
import physics_constants
import bunch
import diagnostics
try:
    import pylab
except:
    print "Warning: could not import pylab"
import beam_parameters
import processor_grid
import computational_grid
import processor_grid
import matching
import time
import field
import math

import sys
import memory

import UberPkgpy #SpaceChargePkgpy

import error_eater
import options

from mpi4py import MPI

def adjust_particles(base,procs):
    retval = base
    multiple = base/(procs * 10.0)
    if not multiple == round(multiple):
        retval = round(multiple + 1) * \
                   (procs * 10)
    return int(retval)

mem0 = 0.0

if ( __name__ == '__main__'):
    t0 = time.time()

#
# Beam parameters at entrance:
# TWISS_EGETAWAY : BETA0, BETX = 1.6479, ALFX =  0.4982, &
##                          BETY = 8.8630, ALFY = -2.3771
#
##  BEAM_EGETAWAY  :  BEAM, ENERGY = 5, NPART = BunchCharge, &
##                    SIGE = 1.5E-3, SIGT = 9.E-03, &
##                    EX = 8E-6 * EMASS/5, EY = 20E-9 * EMASS/5
# where
# EMASS = 0.0005110034
# BunchCharge =  2e10
# where average Current is N/L_scale*q*c
# sigE de/e
# sigt = bunch length (sigma) in meters
# Ex/y = geometric emittance x/y



    myopts = options.Options("ertml")
    myopts.add("current",1.6e-4,"current",float)
    # transverse = 1 for longitudinally uniform beam
    myopts.add("transverse",0,"longitudinally uniform beam",int)
    myopts.add("maporder",2,"map order",int)
    myopts.add("open_all",1,"open transverse boundary conditions", int)
    myopts.add("sigmaz_m",9.E-03,"sigma z in meters",float)
    myopts.add("dpop",1.50e-3,"(delta p)/p",float)
    myopts.add("showplot",1,"show plot",int)
    myopts.add("kicksperline",1000,"kicks per line",int)
    myopts.add("xinit",0.0,"x initial",float)
    myopts.add("xpinit",0.0,"x prime initial",float)
    myopts.add("yinit",0.0,"y initial",float)
    myopts.add("ypinit",0.0,"y prime initial",float)
    myopts.add("dpinit",0.0,"dp initial",float)
    myopts.add("scgrid",[17,17,17],"Space charge grid",int)
    myopts.add("part_per_cell",1,"particlels per cell",int)
    
    myopts.parse_argv(sys.argv)
    f = open("command","w")
    for arg in sys.argv:
        f.write("%s "%arg)
    f.write("\n")
    f.close()

### start save
    scaling_frequency = 1.30e9
### end save
    pipe_radius = 0.04
    griddim = myopts.get("scgrid")
    part_per_cell = myopts.get("part_per_cell")
    num_particles = adjust_particles(griddim[0]*griddim[1]*griddim[2] *\
                                     part_per_cell,MPI.COMM_WORLD.Get_size())

##    num_particles = adjust_particles(10,MPI.COMM_WORLD.Get_size())

    ee = error_eater.Error_eater()
    ee.start()

# Need to set positrons!

    g = gourmet.Gourmet("ertml.lat","ERTML",5.0,
                        scaling_frequency,myopts.get("maporder"),particle='positron')
    g.insert_space_charge_markers(myopts.get("kicksperline"))
    units = g.get_u(g.get_initial_energy())
    print units
###    sys.exit(1)

    EMASS = 0.0005110034
    Ex = 8E-6 * EMASS/5
    Ey = 20E-9 * EMASS/5
    rtbetax=1.6479
    rtbetay=8.8630
    [sigma_x,sigma_xprime,r_x]=matching.match_twiss_emittance(Ex,0.4982,rtbetax)   
    [sigma_y,sigma_yprime,r_y]=matching.match_twiss_emittance(Ey,-2.3771,rtbetay)   

# should first do twiss matching as asked
##    [sigma_x,sigma_xprime,r_x,\
##     sigma_y,sigma_yprime,r_y] = matching.envelope_match(
##        myopts.get("emittance"),myopts.get("current"),g)


# Is it enough to set the positron mass here?
    bp = beam_parameters.Beam_parameters(physics_constants.PH_NORM_me,
                                         1.0, 5.0, 0.0, scaling_frequency,
                                         myopts.get("transverse"))
    
    pz = bp.get_gamma() * bp.get_beta() * bp.mass_GeV
    bp.x_params(sigma = sigma_x, lam = sigma_xprime * pz)
    bp.y_params(sigma = sigma_y, lam = sigma_yprime * pz)
    sigma_z_meters = myopts.get("sigmaz_m")
    bp.z_params(sigma = sigma_z_meters, lam = myopts.get("dpop")* pz)
    bp.correlation_coeffs(xpx = -r_x, ypy = -r_y)


##    bp.x_params(sigma = myopts.get("xinit"), lam = myopts.get("xpinit")*pz)
##    bp.y_params(sigma = myopts.get("yinit"), lam = myopts.get("ypinit")*pz)
##    sigma_z_meters = bp.get_beta()*physics_constants.PH_MKS_c/\
##                     scaling_frequency/math.pi * myopts.get("zfrac")
##    bp.z_params(sigma = sigma_z_meters, lam = myopts.get("dpop")* pz)
##    bp.correlation_coeffs(xpx = 0.0, ypy = 0.0)

    pgrid = processor_grid.Processor_grid(1)

    if myopts.get("open_all")==1:
        cgrid = computational_grid.Computational_grid(griddim[0],griddim[1],griddim[2],
                                                      "3d open")
    else:
        cgrid = computational_grid.Computational_grid(griddim[0],griddim[1],griddim[2],
                                                      "trans finite, long periodic round")
    piperad = 0.04
    field = field.Field(bp, pgrid, cgrid, piperad)
    b = bunch.Bunch(myopts.get("current"), bp, num_particles, pgrid)
    b.generate_particles()

    b.write_particles("initial.dat")
     
    line_length = g.orbit_length()
    tau = 0.5*line_length/myopts.get("kicksperline")
    s = 0.0
    b.write_fort(s)


    line_x = None
    line_y = None
    steps = 0
    if MPI.COMM_WORLD.Get_rank() == 0 and myopts.get("showplot"):
        pylab.ion()
        pylab.hold(0)
        xpl=[]
        xppl=[]
        ypl=[]
        yppl=[]
        xx=[]
        yy=[]
        rr=[]
##    for j in range(0,num_particles-1):
##        b.particles()[:,j]=0.0
    print "**** ", num_particles
    b.particles()[0,0]=myopts.get("xinit")*units[0]
    b.particles()[1,0]=myopts.get("xpinit")*units[1]
    b.particles()[2,0]=myopts.get("yinit")*units[2]
    b.particles()[3,0]=myopts.get("ypinit")*units[3]
    b.particles()[4,0]=0.0
    b.particles()[5,0]=myopts.get("dpinit")*units[5]
###    print " myopts ", myopts.get("xpinit")

    for kick in range(0,myopts.get("kicksperline")):
        if MPI.COMM_WORLD.Get_rank() == 0 and myopts.get("showplot"):
            if(int(b.particles()[6,0]) == 1):
                xpl.append(b.particles()[0,0]/
                           (rtbetax*units[0]))
                xppl.append(b.particles()[1,0]*
                            rtbetax/units[1])
                ypl.append(b.particles()[2,0]/
                           (rtbetay*units[2]))
                yppl.append( b.particles()[3,0]*
                             rtbetay/units[3])
                xt=b.particles()[0,0]/units[0]
                xx.append(xt)
                yt=b.particles()[2,0]/units[2]
                yy.append(yt)
                rr.append(math.sqrt(xt**2+yt**2))
            else:
                print "elapsed time =",time.time() - t0, " particle lost at turn cell ", turn, cell
                pylab.subplot(2,2,1)
                pylab.plot(xpl,xppl,'r,')
                pylab.title('x vs xp')
                pylab.subplot(2,2,2)
                pylab.plot(ypl,yppl,'g,')
                pylab.title('y vs yp')
                pylab.subplot(2,2,3)
                pylab.plot(xx,yy,'b,')
                pylab.title('x vs y')
                pylab.subplot(2,2,4)
                pylab.plot(rr,'b,')
                pylab.title('r vs s')
                pylab.show()
                sys.exit(0)

        steps += 1
###                if MPI.COMM_WORLD.Get_rank() == 0:
###                    print "turn %d, cell %d, kick %d" %\
###                          (turn,cell,kick)
        ###        wrapped_apply_map(g.maps[kick*2],b)
        g.get_fast_mapping(kick*2).apply(b.particles(), b.num_particles_local())
        sys.stdout.flush()
        UberPkgpy.Apply_SpaceCharge_external(\
                    b.get_beambunch(),
                    pgrid.get_pgrid2d(),
                    field.get_fieldquant(),
                    field.get_compdom(),
                    field.get_period_length(),
                    cgrid.get_bc_num(),
                    field.get_pipe_radius(),
                    tau, 0, scaling_frequency, 0)
        sys.stdout.flush()
        ###        wrapped_apply_map(g.maps[kick*2+1],b)
        g.get_fast_mapping(kick*2+1).apply(b.particles(), b.num_particles_local())
        s += line_length/myopts.get("kicksperline")
        b.write_fort(s)
        if kick%50 == 0:
            particle_output_str = 'section'+`kick`
            b.write_particles(particle_output_str)
            
    print "elapsed time =",time.time() - t0
    pylab.subplot(3,2,1)
    pylab.plot(xpl,xppl,'r,')
    pylab.title('x vs xp')
    pylab.subplot(3,2,2)
    pylab.plot(ypl,yppl,'g,')
    pylab.title('y vs yp')
    pylab.subplot(3,2,3)
    pylab.plot(xx,yy,'b,')
    pylab.title('x vs y')
    pylab.subplot(3,2,4)
    pylab.plot(rr,'r,')
    pylab.title('r vs s')
    pylab.subplot(3,2,5)
    pylab.plot(xx,'b,')
    pylab.title('x vs s')
    pylab.subplot(3,2,6)
    pylab.plot(yy,'b,')
    pylab.title('y vs s')
    pylab.show()
                    

