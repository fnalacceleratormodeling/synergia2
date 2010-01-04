function [C,r] = matchCfixedz(M,sigmax,sigmay,sigmaz,kenergy,particle,deltapop)
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

  sigmat = sigmaz/beta;
  if nargin > 5
    have_deltapop = true;
    sigmatp = deltapop*gamma*beta*m;
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
      if (abs(sumM(5,5)) > maxzcomp)
	maxzcomp = abs(sumM(5,5));
	Mz = sumM;
	iz1 = i;
	iz2 = nextind(i);
      endif
    endif
  endfor
  if any(sort([ix1,ix2,iy1,iy2,iz1,iz2]) != [1,2,3,4,5,6])
    error("Failed to find proper order of eigenvalues.");
  endif


### Solve for [ax,ay,az] using following Maxima code:
### eqx : ax*Mx_xx + ay*My_xx + az*Mz_xx = sigmax^2;
### eqy : ax*Mx_yy + ay*My_yy + az*Mz_yy = sigmay^2;
### eqz : ax*Mx_zz + ay*My_zz + az*Mz_zz = sigmat^2;
### ans : linsolve([eqx,eqy,eqz],[ax,ay,az]);
### load("f90");
### f90(ans[1]);
### f90(ans[2]);
### f90(ans[3]);
  ax = (My(1,1)*(Mz(5,5)*sigmay**2-Mz(3,3)*sigmat**2)+Mz(1,1)*(My(3,3)* ...
     sigmat**2-My(5,5)*sigmay**2)+(My(5,5)*Mz(3,3)-My(3,3)*Mz(5,5))...
     *sigmax**2)/ ...
     (Mx(1,1)*(My(5,5)*Mz(3,3)-My(3,3)*Mz(5,5))+...
     My(1,1)*(Mx(3,3)*Mz(5,5)-Mx(5,5)* ...
     Mz(3,3))+(Mx(5,5)*My(3,3)-Mx(3,3)*My(5,5))*Mz(1,1));
  ay = -(Mx(1,1)*(Mz(5,5)*sigmay**2-Mz(3,3)*sigmat**2)+Mz(1,1)*(Mx(3,3)* ...
     sigmat**2-Mx(5,5)*sigmay**2)+(Mx(5,5)*Mz(3,3)-Mx(3,3)*Mz(5,5))...
     *sigmax**2)/ ...
     (Mx(1,1)*(My(5,5)*Mz(3,3)-My(3,3)*Mz(5,5))+...
     My(1,1)*(Mx(3,3)*Mz(5,5)-Mx(5,5)* ...
     Mz(3,3))+(Mx(5,5)*My(3,3)-Mx(3,3)*My(5,5))*Mz(1,1));
  az = (Mx(1,1)*(My(5,5)*sigmay**2-My(3,3)*sigmat**2)+My(1,1)*(Mx(3,3)* ...
     sigmat**2-Mx(5,5)*sigmay**2)+(Mx(5,5)*My(3,3)-Mx(3,3)*My(5,5))...
     *sigmax**2)/ ...
     (Mx(1,1)*(My(5,5)*Mz(3,3)-My(3,3)*Mz(5,5))+...
     My(1,1)*(Mx(3,3)*Mz(5,5)-Mx(5,5)* ...
     Mz(3,3))+(Mx(5,5)*My(3,3)-Mx(3,3)*My(5,5))*Mz(1,1));

  C = ax*Mx + ay*My + az*Mz

  sigmax
  sigmay
  sigmaz
  sqrt(C(1,1))
  sqrt(C(3,3))
  sqrt(C(5,5))

  if (abs(det(C(5:6,5:6))) > 1.0e-15)
    printf("nonsingular longitudinal components\n");
  else
    printf("***singular longitudinal components***\n");
    if !have_deltapop
      error("deltapop must be specified when longitudinal components \
	  are singular");
    endif
    C(6,6) = sigmatp**2;
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
