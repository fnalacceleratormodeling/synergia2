function [sigma_x,sigma_xprime,r_x,sigma_y,sigma_yprime,r_y] = ...
      envelope_match(alphax, alphay, betax, betay,...
		     s_array, Kx_array, Ky_array,...
		     kin_energy_GeV, current_A, eps_x_in, eps_y_in, ...
		     num_points, accuracy, verbose)

### code to solve envelope equations

  global Kx
  global Ky
  global Q
  global r0
  global lambda
  global A
  global beta
  global gamma
  global eps_x
  global eps_y
  global do_map
  do_map = true;

  if(nargin<1)
    printf("Envelope_match error: you must be a fool, usage is:\n");
    printf("[sigma_x,sigma_xprime,r_x,sigma_y,sigma_yprime,r_y] = ...\n");
    printf("envelope_match(madfile_name,lat_name,kin_energy_GeV,num_points,current_A, eps_x, eps_y, accuracy, verbose, elementcalc_path)\n");
    printf("accuracy is optional and defaults to 1.0e-9.\n");
    printf("verbose is optional and defaults to true.\n");
    printf("elementcalc_path is optional.\n");
    printf("The best matched envelope will be saved in the file envelope_match.dat\n");
    printf("with columns 1=z, 2=x, 3=x', 4=y, 5=y'\n");
    s=0;
    x=0;
    return;
  endif
  
#   if (nargin<8)
#     epsilon = 1.0e-9;
#   else
#     epsilon = accuracy;
#   endif

#   if (nargin<9)
#     verbose=true;
#   endif

  epsilon = 1.0e-9;
  verbose = true;

  eps_x = eps_x_in; eps_y = eps_y_in;


  m = 0.93827231; # mass of the proton in GeV/c^2
  gamma = 1 + kin_energy_GeV/m;
  beta = sqrt(1-1/gamma**2);
  c =  299792458; # speed of light in m/s

  r0 = 1.534698e-18; # classical radius of the proton in m
  e = 1.6021773e-19; # electron charge in C
  Q = 1.0; # charge in units of e
  A = 1.0; # atomic number
  lambda = current_A/(e*beta*c); # Line charge density in units of e/m

  beta_x0 = betax;
  betap_x0 = -2.*alphax;
  beta_y0 = betay;
  betap_y0 = -2.*alphay;

  sigma_x = sqrt(eps_x*beta_x0);
  wpx0 = 0.5*sqrt(eps_x/beta_x0)*betap_x0;
###  sigma_prime_x = sign(betap_x0)*sqrt( wpx0**2 + (eps_x/sigma_x)**2 );
  sigma_prime_x = betap_x0*eps_x/(2.*sigma_x);

  sigma_y = sqrt(eps_y*beta_y0);
  wpy0 = 0.5*sqrt(eps_y/beta_y0)*betap_y0;
###  sigma_prime_y = sign(betap_y0)*sqrt( wpy0**2 + (eps_y/sigma_y)**2 );
  sigma_prime_y = betap_y0*eps_y/(2.*sigma_y);

  x0 = [sigma_x;sigma_prime_x;sigma_y;sigma_prime_y];

  delta=[1;1;1;1];
  lsode_options("absolute tolerance",1.0e-15);


  t0 = time();
  t00 = t0;

  while (max(abs(delta))>epsilon)
    if verbose
      printf("x = %#10.6g %#10.6g %#10.6g %#10.6g, delta = %.2e (%0.2f s)\n",...
	     x0,max(abs(delta)),time()-t0);
    endif
    t0 = time();
    fflush(stdout);
    elem_x0 = x0;
    s = zeros((rows(Kx_array)-1)*num_points,1);
    x = zeros((rows(Kx_array)-1)*num_points,4);
    M0 = [1;0;0;0; 0;1;0;0; 0;0;1;0; 0;0;0;1];
    for i = 1:rows(Kx_array)-1
      Kx = Kx_array(i+1);
      Ky = Ky_array(i+1);
      length = s_array(i+1) - s_array(i);
      elem_all0=zeros(1,20);
      elem_all0(1:4) = elem_x0;
      elem_all0(5:20) = M0;

      if (length > 0.0)
	elem_s  = linspace(s_array(i),s_array(i+1),num_points);
	elem_all=lsode("envelope_rhs",elem_all0,...
		       elem_s,s_array(i+1));
	elem_x = elem_all(:,1:4);
	elem_M = elem_all(:,5:20);
	x((i-1)*num_points+1:i*num_points,:)=elem_x;
	s((i-1)*num_points+1:i*num_points,1)=elem_s';
	M0 = elem_M(rows(elem_M),:);
	elem_x0=elem_x(rows(elem_x),:);
      else
	elem_s = s_array(i);
	elem_all=elem_all0+envelope_rhs_delta(elem_all0,elem_s);
	elem_x = elem_all(1:4);
	elem_M = elem_all(5:20);
	for dumb = (i-1)*num_points+1:i*num_points
	  x(dumb,:)=elem_x;
	  s(dumb,1)=elem_s';
	endfor
	M0 = elem_M;
	elem_x0=elem_x;
      endif    
    endfor
    if verbose
      plot(s,x(:,1));
    endif
    
    Mmatrix = [M0(1),M0(2), M0(3), M0(4);  M0(5),M0(6),M0(7),M0(8);...
	       M0(9),M0(10),M0(11),M0(12); M0(13),M0(14),M0(15),M0(16)];
    Jinv = inv(eye(4) - Mmatrix);
    F = x0 - elem_x0';
    xnew = x0 - Jinv*F;
    for i=1:4
###      delta(i) = (xnew(i) - x0(i))/xnew(i);
      delta(i) = (xnew(i) - x0(i))/(1.0+abs(xnew(i)));
    endfor
    x0 = xnew;
  endwhile
  if verbose
    printf("x = %10.6g %10.6g %10.6g %10.6g, delta = %.2e (%0.2f s tot) \n",...
	   x0,max(abs(delta)),time()-t00);
    sigma_x = x0(1);
    sigma_xprime = sqrt( (eps_x/sigma_x)**2 + x0(2)**2 );
    r_x = -x0(2)/sigma_xprime;

    sigma_y = x0(3);
    sigma_yprime = sqrt( (eps_y/sigma_y)**2 + x0(4)**2 );
    r_y = -x0(4)/sigma_yprime;

  endif
#  export("envelope_match.dat",s,x);
endfunction
