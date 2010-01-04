function xp=envelope_rhs_gauss(x,s)

###From Alex Chao's "Physics of Collective Beam Instabilities, pg 37, equ 1.101
###we use the expression for the x,y rms beam evolution in the Gaussian beam
###case, so x(1)=sqrt(<x*x>) and x(2)=sqrt(<y*y>).  Note, this is an 
###approximation based on linearized force

global Q
global r0
global lambda
global A
global beta
global gamma
global eps_x
global eps_y

xi = 4.0*Q**2*r0*lambda/(A*(beta*gamma)**2);

xp(1) = x(3);
xp(2) = x(4);

xp(3) = xi/(2.*(x(1)+x(2))) - Kxy(s,1)*x(1) + eps_x**2/x(1)**3;
xp(4) = xi/(2.*(x(1)+x(2))) - Kxy(s,2)*x(2) + eps_y**2/x(2)**3;

endfunction


