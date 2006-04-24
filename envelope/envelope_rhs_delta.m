function xp=envelope_rhs_delta(x,s)

### From Alex Chao's "Physics of Collective Beam Instabilities, pg 36, equ 1.98
### we use the expression for the x,y rms beam evolution in the KV case
### so, x(1)=sqrt(<x*x>) and x(2)=sqrt(<y*y>)

### This file assumes that all Kx and Ky are delta functions (really K*L).
### Only terms proportional to delta functions are kept.

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

xp(1) = 0;
xp(3) = 0;

xp(2) = - Kx*x(1);
xp(4) = - Ky*x(3);

if do_map
  xp(1+4) = 0;
  xp(2+4) = 0;
  xp(3+4) = 0;
  xp(4+4) = 0;

  c1 = Kx;
  xp(5+4) = -x(1+4)*c1;
  xp(6+4) = -x(2+4)*c1;
  xp(7+4) = -x(3+4)*c1;
  xp(8+4) = -x(4+4)*c1;

  xp(9+4) = 0;
  xp(10+4) = 0;
  xp(11+4) = 0;
  xp(12+4) = 0;

  c2 = Ky;
  xp(13+4) = -x(9+4)*c2;
  xp(14+4) = -x(10+4)*c2;
  xp(15+4) = -x(11+4)*c2;
  xp(16+4) = -x(12+4)*c2;
endif
  
endfunction


