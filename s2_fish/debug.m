mr = 6;
mphi = 4;
mz = 4;

rmax = 1.0;
deltar = rmax*2.0/(2.0*mr+1.0);

b = zeros(mr,mphi, mz);
r0 = 0.4;
x0 = r0/deltar + 0.5;

for i = 1:mr
  binleft = (i-1)*deltar;
  binright = i*deltar;
###  r = (i-0.5)*deltar;
  r = (binright-binleft)/2.0;
  if (binright < r0)
    b(i,:,:) = -4*pi;
  elseif (binleft < r0)
    b(i,:,:) = -4*pi*(r0**2-binleft**2)/(binright**2-binleft**2);
  endif
###  printf("i=%d, r=%f, x0=%f\n",i,r,x0);
endfor

# make non-zero on outside:
b(3,:,:) = -1.0;
b(4,:,:) = -1.0;
b(5,:,:) = -2.0;
b(6,:,:) = -2.0;
# punch an asymmetric hole:
b(2,2,2) = 100.0;

%  format long
%  b
%  format short

blm = zeros(mr, mphi, mz);
tmp = zeros(mphi,mz);
for i=1:mr
  tmp(:,:) = b(i,:,:);
  blm(i,:,:) = fft2(tmp);
endfor

# b(1,:,:)
# blm(1,:,:)

# blm

A = zeros(mr,mr);
rsoln = zeros(mr,1);
psilm = zeros(mr,mphi, mz);
psi = zeros(mr, mphi, mz);

deltar
for l=0:mphi-1
  for m=0:mz-1
    for i=1:mr
      r = (i-0.5)*deltar;
      rsoln(i) = r;
      if (i == 1)
        A(i,i+1) = 1.0/deltar**2 + 1.0/(2*deltar*r) ;
      elseif (i==mr)
        A(i,i-1) = 1.0/deltar**2 - 1.0/(2*deltar*r);
      else
        A(i,i-1) = 1.0/deltar**2 - 1.0/(2*deltar*r);
        A(i,i+1) = 1.0/deltar**2 + 1.0/(2*deltar*r);
      endif
      A(i,i) = -2.0*(1.0/deltar**2);
      A(i,i) += - l**2/r**2 - m**2;
    endfor
#     A
    psilm(:,l+1,m+1) = A\blm(:,l+1,m+1);
  endfor
endfor

# psilm

for i=1:mr
  tmp(:,:) = psilm(i,:,:);
  psi(i,:,:) = ifft2(tmp);
endfor

%  real(psi)
rho = zeros(mr,mphi, mz);
rho(:,:,1) = load("rho-cxx0.dat");
rho(:,:,2) = load("rho-cxx1.dat");
rho(:,:,3) = load("rho-cxx2.dat");
rho(:,:,4) = load("rho-cxx3.dat");

diff_rho = b - rho;
# b
# rho
%  diff_rho = diff_rho ./ b

phi = zeros(mr,mphi, mz);
phi(:,:,1) = load("phi-cxx0.dat");
phi(:,:,2) = load("phi-cxx1.dat");
phi(:,:,3) = load("phi-cxx2.dat");
phi(:,:,4) = load("phi-cxx3.dat");

psi
real(psi)
phi
diff_phi = real(psi) - phi;
diff_phi = diff_phi ./ real(psi)

fftwmz = mz/2+1;
philm = zeros(mr,mphi,fftwmz);
philm(:,:,1) = load("philm-cxx0imag.dat")*sqrt(-1);
philm(:,:,2) = load("philm-cxx1imag.dat")*sqrt(-1);
philm(:,:,3) = load("philm-cxx2imag.dat")*sqrt(-1);
philm(:,:,1) += load("philm-cxx0real.dat");
philm(:,:,2) += load("philm-cxx1real.dat");
philm(:,:,3) += load("philm-cxx2real.dat");

philm;
diff_philm = psilm(:,:,1:fftwmz) - philm;

printf("max(diff_rho) = %g\n",max(max(max(abs(diff_rho)))));
printf("max(diff_phi) = %g\n",max(max(max(abs(diff_phi)))));
printf("max(abs(diff_philm)) = %g\n",max(max(max(abs(diff_philm)))));

printf("real(psi(1,:,:))\n");
real(psi(1,:,:))
printf("phi(1,:,:))\n");
phi(1,:,:)

printf("real(psi(2,:,:))\n");
real(psi(2,:,:))
printf("phi(2,:,:))\n");
phi(2,:,:)

printf("real(psi(3,:,:))\n");
real(psi(3,:,:))
printf("phi(3,:,:))\n");
phi(3,:,:)

printf("real(psi(4,:,:))\n");
real(psi(4,:,:))
printf("phi(4,:,:))\n");
phi(4,:,:)

printf("real(psi(5,:,:))\n");
real(psi(5,:,:))
printf("phi(5,:,:))\n");
phi(5,:,:)

printf("real(psi(6,:,:))\n");
real(psi(6,:,:))
printf("phi(6,:,:))\n");
phi(6,:,:)


%  psi22 = zeros(mr,1);
%  psi22 = psi(:,1,1);
%  
%  plot(rsoln,psi22,'1-*;jim;');
%  
%  hold on
%  #x1 = 0:.05:1;
%  ### r0 = radius of charge
%  x1 = 0:.05:r0;
%  x2 = r0:.05:1;
%  #plot(x1,.25*pi*(3 - 4*x1.**2 + x1.**4),'3;exact;')
%  plot(x1,-pi*x1.*x1 + pi*r0*r0*(-2*log(r0)+1),'3;Exact Solution;')
%  plot(x2,-2*pi*r0*r0*log(x2)+2*pi*r0*r0*log(1),'3;;')
%  ###load soln.dat
%  
%  hold off
