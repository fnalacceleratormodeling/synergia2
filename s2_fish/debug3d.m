
mr = 12;
mphi = 64;
mz = 16;

b = zeros(mr,mphi, mz);

r0 = 0.4;
deltar = r0*2.0/(2.0*mr+1.0);
z0 = 2.0;

coords_r = zeros(1,mr);
coords_phi = zeros(1,mphi);
coords_z = zeros(1,mz);

for i = 1:mr
    r = (i-0.5)*deltar;
    coords_r(i) = r;
    for j = 1:mphi
        theta = (j-1.0)/(mphi-1.0)*2*pi;
        coords_phi(j) = theta;
        for k = 1:mz
            z = (k-1.0)/(mz-1.0)*2*z0 - z0;
            coords(k) = z;
            b(i,j,k) = (((18*r0**2-14*r**2)*sin(3*theta)**2+(18*r**2-18*r0**2)*cos(3* ...
                theta)**2)*cos(pi*z/ ...
                z0)**2*z0**2+(2*pi**2*r**4-2*pi**2*r**2*r0**2)*sin(3*theta)**2* ...
                sin(pi*z/ ...
                z0)**2+(2*pi**2*r**2*r0**2-2*pi**2*r**4)*sin(3*theta)**2* ...
                cos(pi*z/z0)**2)/(r**2*r0**2*z0**2);
        endfor
    endfor
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
      wavenumber_l = mod(l+mphi/2,mphi)-mphi/2;
      wavenumber_m = mod(l+mz/2,mz)-mz/2;
      A(i,i) += - wavenumber_l**2/r**2 - (2*pi*wavenumber_m/(2*z0))**2;
    endfor
    psilm(:,l+1,m+1) = A\blm(:,l+1,m+1);
  endfor
endfor

for i=1:mr
  tmp(:,:) = psilm(i,:,:);
  psi(i,:,:) = ifft2(tmp);
endfor

middle_mphi = mphi/2;
middle_mz = mz/2;

middle_mphi = mphi/2;
middle_mz = mz/2;
middle_mr = mr/2;
hold off
plot(coords_phi,psi(middle_mr,:,middle_mz),'r*');

analytic = zeros(1,mphi);
for j=1:mphi
    r = coords_r(middle_mr);
    theta = coords_phi(j);
    z = coords_z(middle_mz);
    analytic(j) = (1-r**2/r0**2)*sin(3*theta)**2*cos(pi*z/z0)**2;
endfor

hold on
plot(coords_phi,0.5-0.5*analytic);

