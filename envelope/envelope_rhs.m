function xp=envelope_rhs(x,s)

### From Alex Chao's "Physics of Collective Beam Instabilities, pg 36, equ 1.98
### we use the expression for the x,y rms beam evolution in the KV case
### so, x(1)=sqrt(<x*x>) and x(2)=sqrt(<y*y>)

global Q
global r0
global lambda
global A
global beta
global gamma
global eps_x
global eps_y
global Kx
global Ky
global do_map

xi = 4.0*Q**2*r0*lambda/(A*beta**2*gamma**3);

xp(1) = x(2);
xp(3) = x(4);

xp(2) = xi/(4.*(x(1)+x(3))) - Kx*x(1) + eps_x**2/x(1)**3;
xp(4) = xi/(4.*(x(1)+x(3))) - Ky*x(3) + eps_y**2/x(3)**3;

if do_map
  sigma_x = x(1);
  sigma_y = x(3);
  K_s = xi/(4*(sigma_x+sigma_y)**2);

  xp(1+4) = x(5+4);
  xp(2+4) = x(6+4);
  xp(3+4) = x(7+4);
  xp(4+4) = x(8+4);

  c1 = 3*eps_x**2/sigma_x**4+K_s+Kx;
  xp(5+4) = -x(1+4)*c1-2*x(9+4)*K_s;
  xp(6+4) = -x(2+4)*c1-2*x(10+4)*K_s;
  xp(7+4) = -x(3+4)*c1-2*x(11+4)*K_s;
  xp(8+4) = -x(4+4)*c1-2*x(12+4)*K_s;

  xp(9+4) = x(13+4);
  xp(10+4) = x(14+4);
  xp(11+4) = x(15+4);
  xp(12+4) = x(16+4);

  c2 = 3*eps_y**2/sigma_y**4+K_s+Ky;
  xp(13+4) = -x(9+4)*c2-2*x(1+4)*K_s;
  xp(14+4) = -x(10+4)*c2-2*x(2+4)*K_s;
  xp(15+4) = -x(11+4)*c2-2*x(3+4)*K_s;
  xp(16+4) = -x(12+4)*c2-2*x(4+4)*K_s;
endif
  
endfunction


