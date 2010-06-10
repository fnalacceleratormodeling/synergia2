# n = load("noimp/booster_output.h5");
# i = load("imp/booster_output.h5");
no = load("noimpoffset2/booster_output.h5");
io = load("impoffset2/booster_output.h5");

# figure(0);
# hold off
# plot(io.s,io.mean(3,:));
# hold on
# plot(no.s,no.mean(3,:));

# nof = fft(no.mean(3,:));
# iof = fft(io.mean(3,:));

# plot(abs(nof(1:end/2)));
# hold on
# plot(abs(iof(1:end/2)));

# figure(1);
# hold off
tn1 = load("noimpoffset2/tracks/00/00/01/00.h5").points;
tn2 = load("noimpoffset2/tracks/00/00/02/00.h5").points;
# plot(tn1(:,3));
# hold on
# plot(tn2(:,3));

ti1 = load("impoffset3/tracks/00/00/01/00.h5").points;
ti2 = load("impoffset3/tracks/00/00/02/00.h5").points;
# plot(ti1(:,3));
# hold on
# plot(ti2(:,3));

figure(2);
hold off
fn1 = fft(tn1(:,3));
fn2 = fft(tn1(:,3));
plot(abs(fn1(1:end/2)));
hold on
plot(abs(fn2(1:end/2)));

fi1 = fft(ti1(:,3));
fi2 = fft(ti1(:,3));
plot(abs(fi1(1:end/2)));
plot(abs(fi2(1:end/2)));

