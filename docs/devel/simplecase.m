L=17.0;
N=32;
for n=0:(N-1)
    zn(n+1) = n*L/N;
    f1n(n+1)=-(2*pi/L)^2*cos(2*pi/L* zn(n+1));
endfor

P=100;
for p = 0:(P-1)
    z(p+1) = p*L/P;
    phi1n(p+1) = cos(2*pi/L *z(p+1));
endfor

f1m = fft(f1n);
rhs_m(1) = 0.0;
for m=1:(N-1)
    k = mod(m+N/2,N)-N/2;
    rhs_m(m+1) = -(L/(2*pi*k))**2*f1m(m+1);
endfor
phin = ifft(rhs_m);
offset = phi1n(1)-phin(1);
phin += offset;


hold off
plot(z,phi1n,'r;expected;');
hold on
plot(zn,real(phin),'b*-;real(phi);');
plot(zn,imag(phin),'g*-;imag(phi);');
