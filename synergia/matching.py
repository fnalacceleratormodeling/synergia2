#!/usr/bin/env python

from math import sqrt, sin, acos, pi
import function_cache
import numpy
import sys
from mpi4py import MPI

try:
    import octapy
except:
    pass
import os.path
from synergia import physics_constants, envelope_matching

def match_twiss_width(width,alpha,beta):
    """Calculate input parameters for a matched beam of given width
    using Courant-Snyder (Twiss) parameters. Returns
        (width_prime,r,emittance),
    where width_prime is the width in the conjugate coordinate and r
    is the correlation coefficient."""
    gamma = (1+alpha**2)/beta
    emittance = width**2/beta
    width_prime = sqrt(gamma*emittance)
    r = -alpha/sqrt(1+alpha**2)
    return (width_prime,r,emittance)

def match_twiss_emittance(emittance,alpha,beta):
    """Calculate input parameters for a matched beam of given width
    using Courant-Snyder (Twiss) parameters. Returns
        (width,width_prime,r),
    where width and width_prime are the width in the coordinate and
    its conjugate and r is the correlation coefficient."""
    gamma = (1+alpha**2)/beta
    width = sqrt(beta*emittance)
    width_prime = sqrt(gamma*emittance)
    r = -alpha/sqrt(1+alpha**2)
    return (width,width_prime,r)

def get_alpha_beta(my_gourmet):
    '''Get Courant-Snyder (Twiss) parameters from a Gourmet instance.
    Returns 
        (alpha_x, alpha_y, beta_x, beta_y).'''
    mymap = my_gourmet.get_single_linear_map()
#    print "mymap is "
#    print mymap
    u = my_gourmet.get_u(my_gourmet.get_initial_energy())
    mxx = mymap[0,0]
    mxpxp = mymap[1,1]
    mxxp = mymap[0,1]
    cos_mu = (mxx+mxpxp)/2.0
    mu = acos(cos_mu)
    # beta function is positive
    # use this to pick branch
    if mxxp/sin(mu) < 0:
        mu = 2*pi - mu	
    beta_x = mxxp/sin(mu)*u[1]/u[0]
    alpha_x = (mxx-mxpxp)/(2.0*sin(mu))

    myy = mymap[2,2]
    mypyp = mymap[3,3]
    myyp = mymap[2,3]
    cos_mu = (myy+mypyp)/2.0
    mu = acos(cos_mu)
    # beta function is positive
    # use this to pick branch
    if myyp/sin(mu) < 0:
        mu = 2*pi - mu
    beta_y = myyp/sin(mu)*u[3]/u[2]
    alpha_y = (myy-mypyp)/(2.0*sin(mu))

    return (alpha_x, alpha_y, beta_x, beta_y)
    
def  get_tunes(my_gourmet):
    mymap = my_gourmet.get_single_linear_map()
#    print "mymap is "
#    print mymap
    u = my_gourmet.get_u(my_gourmet.get_initial_energy())
    mxx = mymap[0,0]
    mxpxp = mymap[1,1]
    mxxp = mymap[0,1]
    cos_mu = (mxx+mxpxp)/2.0
    mu = acos(cos_mu) 
    if mxxp/sin(mu) < 0:
       mu = 2*pi - mu
    tune_x=mu/(2.*pi)
    
    
    myy = mymap[2,2]
    mypyp = mymap[3,3]
    myyp = mymap[2,3]
    cos_mu = (myy+mypyp)/2.0
    mu = acos(cos_mu)
    # beta function is positive
    # use this to pick branch
    if myyp/sin(mu) < 0:
        mu = 2*pi - mu
    tune_y=mu/(2.*pi) 
    
    mzz=mymap[4,4]
    mzpzp = mymap[5,5]
    mzzp = mymap[4,5]
    cos_mu = (mzz+mzpzp)/2.0
    mu = acos(cos_mu) 
    if mzzp/sin(mu) < 0:
       mu = 2*pi - mu
    tune_z=mu/(2.*pi)
    
    return (tune_x,tune_y, tune_z)

envelope_match_cache  = function_cache.Function_cache("envelope_match.cache")
def envelope_match(emitx,emity,current,g,use_cache=0,use_octave=0,do_plot=0):
    '''Match a beam with horizontal emittance emitx and vertical emittance emity to a lattice from
    a Gourmet instance g using the envelope equations.
    Returns 
        [sigma_x, sigma_xprime,r_x,sigma_y, sigma_yprime,r_y].'''
    if use_cache:
        if envelope_match_cache.in_cache(emitx,emity,current,g.get_lattice_file(),
                                         g.get_line_name(),g.get_initial_kinetic_energy()):
            return envelope_match_cache.get(emitx,emity,current,g.get_lattice_file(),
                                            g.get_line_name(),g.get_initial_kinetic_energy())
    (alpha_x, alpha_y, beta_x, beta_y) = get_alpha_beta(g)
    (s,kx,ky) = g.get_strengths()
    
    if use_octave:
      o = octapy.Octave()
      o.execute('LOADPATH="%s:";' %
                os.path.join(os.environ["SYNERGIA2DIR"],"envelope"))
      o.set_value("m",g.get_mass())
      o.set_value("alphax",alpha_x)
      o.set_value("alphay",alpha_y)
      o.set_value("betax",beta_x)
      o.set_value("betay",beta_y)
      o.set_value("s_array",s)
      o.set_value("Kx_array",kx)
      o.set_value("Ky_array",ky)
      o.set_value("emitx",emitx)
      o.set_value("emity",emity)
      o.set_value("current",current)
      o.set_value("kinetic_energy",g.get_initial_kinetic_energy())
      o.execute('[sigma_x,sigma_xprime,r_x,sigma_y,sigma_yprime,r_y] = envelope_match(m,alphax,alphay,betax,betay,s_array,Kx_array,Ky_array, kinetic_energy, current, emitx, emity, 4, 1.0e-13,0)' )
      sigma_x = o.get_value("sigma_x")
      sigma_xprime = o.get_value("sigma_xprime")
      r_x = o.get_value("r_x")
      sigma_y = o.get_value("sigma_y")
      sigma_yprime = o.get_value("sigma_yprime")
      r_y = o.get_value("r_y")
      retval = [sigma_x,sigma_xprime,r_x,sigma_y,sigma_yprime,r_y]
      if retval.count(None) == 0:
          envelope_match_cache.add(retval,emitx,emity,current,g.get_lattice_file(),
                                 g.get_line_name(),g.get_initial_kinetic_energy())
      
    else:	    
      mass=g.get_mass() 
      kinetic_energy= g.get_initial_kinetic_energy()     
      gamma =1.0+kinetic_energy/mass
      beta = sqrt(1-1/gamma**2);
    
      charge = 1.0  # electron charge in C
      A = 1.0; # atomic number
    
      lambd = current/(physics_constants.PH_MKS_e*beta*physics_constants.PH_MKS_c)	    
      xi=4.0*charge**2*physics_constants.PH_MKS_rp*lambd/(A*beta**2*gamma**3)
    
      (xwidth,xpwidth,rx) = match_twiss_emittance(emitx,alpha_x,beta_x)
      (ywidth,ypwidth,ry) = match_twiss_emittance(emity,alpha_y,beta_y)
      widths=[xwidth,xpwidth,rx,ywidth,ypwidth,ry]
      
      retval=envelope_matching.envelope_match(widths, s, kx, ky,
           xi, accuracy=1.0e-9, verbose=True, integrator=4,do_plot=do_plot)
      
#     retval= [sigma_x, sigma_xprime, r_x, sigma_y,s igma_yprime, r_y]  
                     
    if retval.count(None) == 0:
       envelope_match_cache.add(retval,emitx,emity,current,g.get_lattice_file(),
                                g.get_line_name(),g.get_initial_kinetic_energy())	     
		     
    return retval

def envelope_motion(widths_in,current,g,do_plot=0,do_match=0):
  
    
    (s,kx,ky) = g.get_strengths()        
    mass=g.get_mass() 
    kinetic_energy= g.get_initial_kinetic_energy()     
    gamma =1.0+kinetic_energy/mass
    beta = sqrt(1-1/gamma**2);    
    charge = 1.0  # electron charge in C
    A = 1.0; # atomic number
    lambd = current/(physics_constants.PH_MKS_e*beta*physics_constants.PH_MKS_c)	    
    xi=4.0*charge**2*physics_constants.PH_MKS_rp*lambd/(A*beta**2*gamma**3)
    
   
      
    retval=envelope_matching.envelope_match(widths_in, s, kx, ky,
           xi, accuracy=1.0e-9, verbose=True, integrator=4,do_map=do_match,do_plot=do_plot)
      
#     retval= [sigma_x, sigma_xprime, r_x, sigma_y,s igma_yprime, r_y]  
                     
   	     
		     
    return retval



def rms_match_3d(linear_map,beam_parameters,arms,brms,crms,rms_index,print_emittances=True):
    '''here are 3 rms input parameters,arms, brms, crms, which corresponds to  indices rms _index[0], rms _index[1], rms _index[2] 
        example: rms_index=[0,2,5]==> arms=xrms, brms=yrms, crms=Urms
         units of rms should be  X_synergia/synegia_units, i.e. [xrms]=m, [pxrms]=Gev/c, [pzrms] = Gev, 
        '''

    numpy_map = linear_map
    
    evals,evect_matrix = numpy.linalg.eig(numpy_map)
    evects = []
    for i in range(0,6):
        evects.append(evect_matrix[:,i])
    F = range(0,3)
    remaining = range(5,-1,-1)
    for i in range(0,3):
        # find complex conjugate among remaining eigenvectors
        first = remaining.pop()
        best = 1.0e30
        conj = -1
        for item in remaining:
            sum = evects[first]+evects[item]
            if abs(numpy.max(sum.imag)) < best:
                best = abs(numpy.max(sum.imag))
                conj = item
        if conj == -1:
            raise RuntimeError,"failed to find a conjugate pair in ha_match"
        remaining.remove(conj)
        tmp=numpy.outer(evects[first],
            numpy.conjugate(evects[first]))
        tmp+=numpy.outer(evects[conj],
            numpy.conjugate(evects[conj]))
        F[i]=tmp.real
    
    S=numpy.zeros((3,3),'d')
    for i in range(0,3):
        for j in range(0,3):
            S[i,j]=F[j][rms_index[i],rms_index[i]]
        
    Sinv=numpy.linalg.inv(S)   
    
    C = numpy.zeros([6,6],'d')
    Cxy, Cxpyp, Cz, Czp = beam_parameters.get_conversions()
    units=[Cxy,Cxpyp,Cxy,Cxpyp,Cz, Czp]
    cd1=arms*units[rms_index[0]]*arms*units[rms_index[0]]
    cd2=brms*units[rms_index[1]]*brms*units[rms_index[1]]
    cd3=crms*units[rms_index[2]]*crms*units[rms_index[2]]
    
    for i in range(0,3):
        C += F[i]*(Sinv[i,0]*cd1+Sinv[i,1]*cd2+Sinv[i,2]*cd3)
        
        
    if print_emittances:
        beta=beam_parameters.get_beta()
        gamma=beam_parameters.get_gamma()
        pz = gamma * beta * beam_parameters.mass_GeV
        energy=beam_parameters.get_kinetic_energy()+beam_parameters.get_mass()
        emitx=sqrt(C[0,0]*C[1,1]-C[0,1]*C[1,0])/units[0]/units[1]
        emity=sqrt(C[2,2]*C[3,3]-C[2,3]*C[3,2])/units[2]/units[3]
        emitz=sqrt(C[4,4]*C[5,5]-C[4,5]*C[5,4])/units[4]/units[5] 
        if MPI.COMM_WORLD.Get_rank() ==0: 
            print "************ BEAM MATCHED PARAMETERS *****************"
            print "*    emitx=", emitx, " meters*GeV/c   =", emitx/pz, " meters*rad =", emitx/pz/pi, " pi*meters*rad"
            print "*    emity=", emity, " meters*GeV/c   =", emity/pz, " meters*rad =", emity/pz/pi, " pi*meters*rad"
            print "*    emitz=", emitz, " meters*GeV =", emitz*1.e9/(physics_constants.PH_MKS_c*beta), " eV*s"
            print " "
            print "*    Normalized emitx=",  emitx*gamma*beta/pz, " meters*rad =", emitx*gamma*beta/pz/pi, " pi*meters*rad"
            print "*    Normalized emity=",  emity*gamma*beta/pz, " meters*rad =", emity*gamma*beta/pz/pi, " pi*meters*rad"
            print " "  
            print "*    xrms=",sqrt(C[0,0])/units[0] , " meters"
            print "*    yrms=",sqrt(C[2,2])/units[2] , " meters"
            print "*    zrms=",sqrt(C[4,4])/units[4] , " meters=",2.*pi*sqrt(C[4,4])/units[4]/beam_parameters.get_z_length(), " rad"
            print "*    pxrms=",sqrt(C[1,1])/units[1] , " GeV/c,    dpx/p=",sqrt(C[1,1])/units[1]/pz
            print "*    pyrms=",sqrt(C[3,3])/units[3] , " GeV/c,    dpy/p=",sqrt(C[3,3])/units[3]/pz
            print "*    pzrms(Erms)=",sqrt(C[5,5])/units[5] , " GeV,  deoe=",sqrt(C[5,5])/units[5]/energy ,\
             ",  dpzop=", sqrt(C[5,5])/(units[5]*energy*beta*beta) 
            print "" 
            print "*    bucket length=",beam_parameters.get_z_length(),  " meters"
            print "*    pz=",pz, "  GeV/c"
            print "*    energy=",energy,"  GeV" 
            print "****************************************************"
            
          
            
           
    return C
