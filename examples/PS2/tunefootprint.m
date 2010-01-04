g = glob("tracks/*/*/*/*.h5");
total = length(g);
fulltracks = 0;
tunesxtmp = ones(1,total);
tunesytmp= ones(1,total);
for i=1:length(g)
    p = load(g{i}).points;
    if (size(p)(1) > 3)
        turns = p(end,7)/p(3,7);
    else
        turns = 0;
    endif
#    if (round(turns*3) == 3000)
        fulltracks += 1;
        fx = fft(p(:,1));
        fy = fft(p(:,3));
        [heightx,indexx]=max(fx(2:rows(fx)/2));
        [heighty,indexy]=max(fy(2:rows(fy)/2));
        tunesxtmp(fulltracks) = indexx*1.0/turns;
        tunesytmp(fulltracks) = indexy*1.0/turns;
#    endif
endfor
total
fulltracks
tunesx=tunesxtmp(1:fulltracks);
tunesy=tunesytmp(1:fulltracks);
