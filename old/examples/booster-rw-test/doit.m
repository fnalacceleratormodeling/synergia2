for i=0:10
printf("i=%02d\n",i);
x = load(sprintf("turn_%02d.h5",i)).particles;
plot(x(:,5),x(:,6),'*');
axis([-pi,pi,-0.002,0.002]);
yy = input("hit return");
endfor
