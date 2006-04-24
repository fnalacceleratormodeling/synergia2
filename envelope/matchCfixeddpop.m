function [C,r] = matchCfixeddpop(M,emitx,emity,deltapop,kenergy,particle,sigmaz)
  if strcmp(particle,"proton")
    m = 0.93827231;
  else if strcmp(particle,"electron")
      m = 0.00051099906;
    else
      error(sprintf("unknown particle %s",particle));
    endif
  endif
  gamma = kenergy/m + 1.0;
  beta = sqrt(1-1/gamma**2);


  sigmatp = deltapop*gamma*beta*m;
  if nargin > 5
    have_sigmaz = true;
    sigmat = sigmaz/beta;
  else
    have_deltapop = false;
  endif

  [B,lambdaM] = eig(M);
  for i = 1:6
    lambda(i) = lambdaM(i,i);
    E{i} = B(:,i)*B(:,i)';
  endfor
### Determine which eigenvectors go with which values. Eigenvectors are
### assumed to come in pairs (the pair could be 6,1). The x pair is
### chosen to be the pair with the largest x,x' determinant, likewise
### for y. The z pair is chosen to be the pair with the largest z,z
### component.
  nextind = [2,3,4,5,6,1];
  maxzcomp = 0.0;
  maxxdet = 0.0;
  maxydet = 0.0;
  for i = 1:6
    prod = lambda(i) * lambda(nextind(i));
    if (abs((real(prod) - 1.0) < 1.0e-10) && ...
	(abs(imag(prod)) < 1.0e-10))
      sumM = E{i} + E{nextind(i)};
      if (abs(det(sumM(1:2,1:2))) > maxxdet)
	maxxdet = abs(det(sumM(1:2,1:2)));
	Mx = sumM;
	ix1 = i;
	ix2 = nextind(i);
      endif
      if (abs(det(sumM(3:4,3:4))) > maxydet)
	maxydet = abs(det(sumM(3:4,3:4)));
	My = sumM;
	iy1 = i;
	iy2 = nextind(i);
      endif
      if (abs(sumM(6,6)) > maxzcomp)
	maxzcomp = abs(sumM(6,6));
	Mz = sumM;
	iz1 = i;
	iz2 = nextind(i);
      endif
    endif
  endfor
  if any(sort([ix1,ix2,iy1,iy2,iz1,iz2]) != [1,2,3,4,5,6])
    error("Failed to find proper order of eigenvalues.");
  endif

  sdetMx_x = sqrt(det(Mx(1:2,1:2)));
  sdetMy_x = sqrt(det(My(1:2,1:2)));
  sdetMz_x = sqrt(det(Mz(1:2,1:2)));
  sdetMx_y = sqrt(det(Mx(3:4,3:4)));
  sdetMy_y = sqrt(det(My(3:4,3:4)));
  sdetMz_y = sqrt(det(Mz(3:4,3:4)));
  Mx_zz = Mx(6,6);
  My_zz = My(6,6);
  Mz_zz = Mz(6,6);

### Solve for [ax,ay,az] using following Maxima code:
###
### eqx : ax*sdetMx_x + ay*sdetMy_x + az*sdetMz_x = emitx;
### eqy : ax*sdetMx_y + ay*sdetMy_y + az*sdetMz_y = emity;
### eqz : ax*Mx_zz + ay*My_zz + az*Mz_zz = sigmatp^2;
### ans : linsolve([eqx,eqy,eqz],[ax,ay,az]);
### load("f90");
### f90(ans[1]);
### f90(ans[2]);
### f90(ans[3]);
  ax = (sdetMy_x*(emity*Mz_zz-sdetMz_y*sigmatp**2)+sdetMz_x* ...
     (sdetMy_y*sigmatp**2-emity*My_zz)+emitx*(My_zz* ...
     sdetMz_y-Mz_zz*sdetMy_y))/ ...
     (sdetMx_x*(My_zz*sdetMz_y-Mz_zz*sdetMy_y)+sdetMy_x*(Mz_zz* ...
     sdetMx_y-Mx_zz*sdetMz_y)+(Mx_zz*sdetMy_y-My_zz*sdetMx_y)* ...
     sdetMz_x);
  ay = -(sdetMx_x*(emity*Mz_zz-sdetMz_y*sigmatp**2)+sdetMz_x* ...
     (sdetMx_y*sigmatp**2-emity*Mx_zz)+emitx*(Mx_zz* ...
     sdetMz_y-Mz_zz*sdetMx_y))/ ...
     (sdetMx_x*(My_zz*sdetMz_y-Mz_zz*sdetMy_y)+sdetMy_x*(Mz_zz* ...
     sdetMx_y-Mx_zz*sdetMz_y)+(Mx_zz*sdetMy_y-My_zz*sdetMx_y)* ...
     sdetMz_x);
  az = (sdetMx_x*(emity*My_zz-sdetMy_y*sigmatp**2)+sdetMy_x* ...
     (sdetMx_y*sigmatp**2-emity*Mx_zz)+emitx*(Mx_zz* ...
     sdetMy_y-My_zz*sdetMx_y))/ ...
     (sdetMx_x*(My_zz*sdetMz_y-Mz_zz*sdetMy_y)+sdetMy_x*(Mz_zz* ...
     sdetMx_y-Mx_zz*sdetMz_y)+(Mx_zz*sdetMy_y-My_zz*sdetMx_y)* ...
     sdetMz_x);

  C = ax*Mx + ay+My + az*Mz;

  if (abs(det(C(5:6,5:6))) > 1.0e-15)
    printf("nonsingular longitudinal components\n");
  else
    printf("***singular longitudinal components***\n");
    if !have_sigmaz
      error("sigmaz must be specified when longitudinal components \
	  are singular");
    endif
    C(5,5) = sigmat**2;
  endif
  r = zeros(6,6);
  for i = 1:6
    for j = 1:6
      r(i,j) = C(i,j)/sqrt(C(i,i)*C(j,j));
      if i == 6
	r(i,j) = -r(i,j);
      endif
      if j == 6
	r(i,j) = -r(i,j);
      endif
    endfor
  endfor

endfunction
