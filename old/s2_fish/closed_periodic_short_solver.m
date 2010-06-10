function [rsoln,realpsi]=closed_periodic_short_solver(mr, mphi, mz)

rmax = 1.0;
deltar = rmax*2.0/(2.0*mr+1.0);

b = zeros(mr,mphi, mz);
r0 = 0.2;
x0 = r0/deltar + 0.5;

zrange = (mz/4):(3*mz/4);
for i = 1:mr
  binleft = (i-1)*deltar;
  binright = i*deltar;
###  r = (i-0.5)*deltar;
  r = (binright-binleft)/2.0;
  if (binright < r0)
    b(i,:,zrange) = -4*pi;
  elseif (binleft < r0)
    b(i,:,zrange) = -4*pi*(r0**2-binleft**2)/(binright**2-binleft**2);
  else
    b(i,:,zrange) = 0.0;
  endif
###  printf("i=%d, r=%f, x0=%f\n",i,r,x0);
endfor

blm = zeros(mr, mphi, mz);
tmp = zeros(mphi,mz);
for i=1:mr
  tmp(:,:) = b(i,:,:);
  blm(i,:,:) = fft2(tmp);
endfor

A = zeros(mr,mr);
rsoln = zeros(mr,1);
psilm = zeros(mr,mphi, mz);
psi = zeros(mr, mphi, mz);

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
    psilm(:,l+1,m+1) = A\blm(:,l+1,m+1);
  endfor
endfor

for i=1:mr
  tmp(:,:) = psilm(i,:,:);
  psi(i,:,:) = ifft2(tmp);
endfor

realpsi=real(psi);

endfunction
