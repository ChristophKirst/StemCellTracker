%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Algorithms %%%
%%%%%%%%%%%%%%%%%%%%%%%


rl = [-1, -2]; ru = [1, 2];
rd = [0,0]; 

r = 0:0.02:3;

a = arrayfun(@(x) overlapDiskRectangle(rd, x, rl, ru), r);

figure(1); clf
plot(r, a)
hold on
plot(r, pi * r.^2 , 'r')
plot(r, 0.* r + prod(ru-rl), 'g')





%% Test DiskSegment


rl = [-1, -2]; ru = [1, 2];
rd = [0,0]; 

r = 0:0.02:3;

a = arrayfun(@(x) overlapDiskSegmentRectangle(rd, x, 3/4*pi, 5/4 * pi, rl, ru), r);

figure(1); clf
plot(r, a)
hold on
plot(r, pi * r.^2/4  , 'r')
plot(r, 1, 'g')