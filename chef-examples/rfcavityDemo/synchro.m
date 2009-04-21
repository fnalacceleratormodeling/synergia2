				# this script works with octave 3.0+

				# read synchro data file
sdata = load("lowv.dat");

plot(sdata(:,1),sdata(:,2),'.');
xlabel("cdt")
ylabel("ndp");
title("synchrotron motion")

				# z limit is .01436
				# dpop limit is 0.0001
zmax = .01436
dpopmax = .0001
betalong = zmax/dpopmax

#__gnuplot_set__("arrow from -0.005,0.0001 to +0.005,0.0001 nohead")
#__gnuplot_set__("label \"max dp/p = 0.0001\" at 0.006,0.0001")

text(-0.005,8e-5,"max dp/p = 0.0001")
text(0.008, 0.0, "max z = 0.01436")
text(-0.005,-8e-5, "beta_L = 143.6")
replot

figure
cdtt = abs(fft(sdata(2:end,1)));
n = rows(cdtt);


df = 1/n;
f = [0:n-1]*df;

f0 = 0.025;
f1 = 0.035;

i0 = floor(f0/df);
i1 = ceil(f1/df);

plot(f(i0:i1), cdtt(i0:i1))
xlabel("synchrotron frequency")
ylabel("power")
text(0.03, 150, "Peak power at 0.028741")
text(0.03, 130, "Synchrotron period = 34.794")
replot
				# the peak is at 0.028740
speak = .028740
sfreq = 1/speak
