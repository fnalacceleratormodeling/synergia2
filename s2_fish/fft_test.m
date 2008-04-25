mr = 2;
mphi = 4;
mz = 4;

b = zeros(mr,mphi, mz);

for i = 1:mr
    for j = 1:mphi
        for k = 1:mz
            b(i,j,k) = i*100.0+j*10.0+k*.011;
        endfor
    endfor
endfor

b

blm = zeros(mr,mphi, mz);
tmp = zeros(mphi,mz);

for i=1:mr
    tmp(:,:) = b(i,:,:);
    blm(i,:,:) = fft2(tmp);
endfor

blm
newb = zeros(mr,mphi, mz);
for i=1:mr
    tmp(:,:) = blm(i,:,:);
    newb(i,:,:) = ifft2(tmp);
endfor

newb
max(max(max(newb - b)))