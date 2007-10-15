#!/usr/bin/env python

from matplotlib.pylab import *
from loadfile import loadfile

left = 1.54e4
x = loadfile("test6_out.dat")
right = max(x[:,0])

figure(1)
subplot(2,1,1)
ylabel("beta x")
plot(x[:,0],x[:,1],label="PT")
plot(x[:,0],x[:,7],label="CHEF")
newaxis = axis()
newaxis[0] = left
newaxis[1] = right
newaxis[3] = 20.0
axis(newaxis)
legend(loc=0)
subplot(2,1,2)
xlabel("s")
ylabel("PT-CHEF")
plot(x[:,0],x[:,1]-x[:,7],label="covariance")
plot(x[:,0],x[:,1]-x[:,5],label="twiss")
newaxis = axis()
newaxis[0] = left
newaxis[1] = right
axis(newaxis)
legend(loc=0)
savefig("betax_zoom.ps")
figure(2)
subplot(2,1,1)
ylabel("beta y")
plot(x[:,0],x[:,2],label="PT")
plot(x[:,0],x[:,8],label="CHEF")
newaxis = axis()
newaxis[0] = left
newaxis[1] = right
newaxis[3] = 20.0
axis(newaxis)
legend(loc=0)
subplot(2,1,2)
xlabel("s")
ylabel("PT-CHEF")
plot(x[:,0],x[:,2]-x[:,8],label="covariance")
plot(x[:,0],x[:,2]-x[:,6],label="twiss")
newaxis = axis()
newaxis[0] = left
newaxis[1] = right
axis(newaxis)
legend(loc=0)
savefig("betay_zoom.ps")

figure(3)
subplot(2,1,1)
ylabel("dispersion")
plot(x[:,0],x[:,3],label="PT")
plot(x[:,0],x[:,9],label="CHEF")
newaxis = axis()
newaxis[0] = left
newaxis[1] = right
axis(newaxis)
legend(loc=0)
subplot(2,1,2)
xlabel("s")
ylabel("PT-CHEF")
plot(x[:,0],x[:,3]-x[:,9],label="covariance")
newaxis = axis()
newaxis[0] = left
newaxis[1] = right
axis(newaxis)
legend(loc=0)
savefig("dispersion_zoom.ps")
figure(4)
subplot(2,1,1)
ylabel("dispersion'")
plot(x[:,0],x[:,4],label="PT")
plot(x[:,0],x[:,10],label="CHEF")
newaxis = axis()
newaxis[0] = left
newaxis[1] = right
axis(newaxis)
legend(loc=0)
subplot(2,1,2)
xlabel("s")
ylabel("PT-CHEF")
plot(x[:,0],x[:,4]-x[:,10],label="covariance")
newaxis = axis()
newaxis[0] = left
newaxis[1] = right
axis(newaxis)
legend(loc=0)
savefig("disp_prime_zoom.ps")


show()
