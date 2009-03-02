#!/usr/bin/env python

from math import sqrt, sin, acos, pi
import function_cache
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

envelope_match_cache  = function_cache.Function_cache("envelope_match.cache")
def envelope_match(emitx,emity,current,g,use_cache=0,use_octave=0,do_plot=0):
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


