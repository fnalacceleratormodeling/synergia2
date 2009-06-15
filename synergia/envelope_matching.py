import Numeric
import numpy
import time
from math import sqrt
import os.path
import sys
from LinearAlgebra import inverse
import synergia
from pygsl import odeiv
import pygsl



def envelope_rhs(s, x, param):
    kx=param[0]
    ky=param[1]
    eps_x=param[2]
    eps_y=param[3]
    xi=param[4]
    do_map=param[5]
    
    if do_map:
	xp = numpy.zeros((20,), 'd') 
    else:
	xp = numpy.zeros((4,), 'd') 
    
    sigma_x=x[0]
    sigma_y=x[2]	

    xp[0] = x[1]
    xp[2] = x[3]	
    xp[1] = xi/(4.*(sigma_x+sigma_y)) - kx*x[0] + eps_x**2/x[0]**3	
    xp[3] = xi/(4.*(sigma_x+sigma_y)) - ky*x[2] + eps_y**2/x[2]**3



    if do_map:
	K_s = xi/(4*(sigma_x+sigma_y)**2)

	xp[0+4] = x[4+4]
	xp[1+4] = x[5+4]
	xp[2+4] = x[6+4]
	xp[3+4] = x[7+4]


	c1 = 3.*(eps_x**2)/(sigma_x**4)+K_s+kx
	xp[4+4] = -x[0+4]*c1-2.*x[8+4]*K_s
	xp[5+4] = -x[1+4]*c1-2.*x[9+4]*K_s	
	xp[6+4] = -x[2+4]*c1-2.*x[10+4]*K_s
	xp[7+4] = -x[3+4]*c1-2.*x[11+4]*K_s

	xp[8+4] = x[12+4]
	xp[9+4] = x[13+4]
	xp[10+4] = x[14+4]
	xp[11+4] = x[15+4]

	c2 = 3.*eps_y**2/sigma_y**4+K_s+ky
	xp[12+4] = -x[8+4]*c2-2.*x[0+4]*K_s
	xp[13+4] = -x[9+4]*c2-2.*x[1+4]*K_s
	xp[14+4] = -x[10+4]*c2-2.*x[2+4]*K_s
	xp[15+4] = -x[11+4]*c2-2.*x[3+4]*K_s

    return xp

#def   envelope_match(alphax, alphay, betax, betay, s_array, kx_array,\
#	ky_array,xi_in,eps_x_in, eps_y_in,\
#	accuracy=1.0e-9, verbose=True, integrator=4, do_plot=0):   
	   
def   envelope_match(widths, s_array, kx_array,ky_array,xi_in,\
	accuracy=1.0e-9, verbose=True, integrator=4, do_map=True, do_plot=0):   	   
      
      print "envelope_matching.py  called" 
   
      stepper_list = [odeiv.step_rk2,odeiv.step_rk4,odeiv.step_rkf45,\
                    odeiv.step_rkck,odeiv.step_rk8pd,odeiv.step_rk2imp,\
                    odeiv.step_rk4imp,odeiv.step_gear1,odeiv.step_gear2]
      stepper=stepper_list[integrator]
      if verbose:
        print "using integration method:  ", stepper	   
      

      epsilon = accuracy      
      eps_x = sqrt((widths[0]**2)*(widths[1]**2)*(1-widths[2]**2))
      eps_y = sqrt((widths[3]**2)*(widths[4]**2)*(1-widths[5]**2))
      xi=xi_in     #  4.0*Q**2*r0*lambd/(A*beta**2*gamma**3)
      
      param=[0.,0., eps_x,eps_y, xi, do_map] # param=[kx,ky, eps_x,eps_y, xi, True] 
   
      
      #beta_x0 = betax
      #betap_x0 = -2.*alphax
      #beta_y0 = betay
      #betap_y0 = -2.*alphay
      
      #sigma_x = sqrt(eps_x*beta_x0)
      #sigma_prime_x = betap_x0*eps_x/(2.*sigma_x)
       
      #sigma_y = sqrt(eps_y*beta_y0)
      #sigma_prime_y = betap_y0*eps_y/(2.*sigma_y)
      sigma_x = widths[0]
      sigma_prime_x=widths[1]*widths[2]

      sigma_y = widths[3]
      sigma_prime_y=widths[4]*widths[5]


      x0 = numpy.array([sigma_x,sigma_prime_x,sigma_y,sigma_prime_y], Numeric.Float)       
      delta=numpy.array([1,1,1,1])

      if verbose:
	print " " 
        print "x_initial= %e  %e  %e  %e " % (x0[0],x0[1],x0[2],x0[3])


      max_steps=100000
      it=0

      t00=time.time()
      t0=t00
      t1=time.time()-t0
      while (max(abs(delta))>epsilon):
	if verbose:      
	  print "x        = %e  %e  %e  %e   delta=  %e   ( %f s)"  %\
	   (x0[0],x0[1],x0[2],x0[3],max(abs(delta)),t1 )

	M0 = numpy.array([1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1])
        if do_map:
          elem_all=numpy.concatenate((x0,M0),1)
        else:
         elem_all=x0

	dimension = len(elem_all)
	n_elem=numpy.size(kx_array)
	if do_plot:	
	    elem_stored=numpy.zeros((n_elem, dimension), 'd')
	
	for i in range(n_elem):	
	    if i==0:
	      ss=0
	    else:	
	      ss=s_array[i-1] 
	    
	    s1=s_array[i]		
	    param[0]=kx_array[i]
	    param[1]=ky_array[i]

            step=stepper(dimension,  envelope_rhs, None, param)
	    control = odeiv.control_y_new(step, 1e-15, 1e-15)
	    evolve  = odeiv.evolve(step, control, dimension)
	    h=0.03
	  	  	
	    for j in range(max_steps):	  
	      if ss>=s1:		
		  break	
	      ss,h,elem_all=evolve.apply(ss, s1, h, elem_all)
	      if (int(pygsl.version.split('.')[1]) < 9):
	        elem_all=elem_all[-1]

	    else:
	          raise ValueError, "Maximum number of steps exceeded!"
	    if do_plot:	 
	       elem_stored[i,0:dimension]=elem_all[0:dimension]	
	
	
	if do_map:
	  M_matrix=numpy.array(elem_all[4:20])
	  M_matrix=numpy.reshape(M_matrix,(4,4))
	  Jinv=numpy.linalg.inv(numpy.identity(4)-M_matrix)
	  F = x0 - numpy.array(elem_all[0:4])
	 

          xnew = x0 - numpy.dot(Jinv, F)
	  delta=(xnew - x0)/(1.0+abs(xnew))	 
	  it=it+1
	  
	  x0 = xnew
	  t1=time.time()-t0
	  t0=time.time()
	else: 
	  break
	

      xx=numpy.array(elem_all[0:4])
      sigma_x = xx[0]
      sigma_xprime=sqrt( (eps_x/sigma_x)**2 + xx[1]**2 )
      r_x=x0[1]/sigma_xprime 
      
      sigma_y = xx[2]
      sigma_yprime=sqrt( (eps_y/sigma_y)**2 + xx[3]**2 )
      r_y=xx[3]/sigma_yprime 
    
      if verbose: 
        print "x_final  = %e  %e  %e  %e  ( %f s total time)" % (x0[0],x0[1],x0[2],x0[3],time.time()-t00)
	print " "
        print "sigma_x = %e" % sigma_x
        print "sigma_xprime= %e" % sigma_xprime
        print "r_x= %e" % r_x
        print "sigma_y = %e" % sigma_y
        print "sigma_yprime= %e " % sigma_yprime
        print "r_y= %e" % r_y
	
      if do_plot:
          import pylab
	  
	  sigma_plot=numpy.zeros(n_elem, 'd')
	  sigma_plot[0:n_elem]=elem_stored[0:n_elem, 0]
	  pylab.figure(1)
	  pylab.plot(s_array, sigma_plot,'ro',label='env eq ')   


	  for i in range(n_elem):
	     sigma_plot[i]= sqrt( (eps_x/elem_stored[i, 0])**2 + elem_stored[i, 1]**2 )
	  pylab.figure(2)
	  pylab.plot(s_array, sigma_plot,'ro',label='env eq')   	
    
      return [sigma_x, sigma_xprime,r_x,sigma_y, sigma_yprime,r_y]
 




