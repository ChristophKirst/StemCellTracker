%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Paritial Alignment of Black Images %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc
initialize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Overlap Statistics

img = loadTestImage();

img1 = img(1:230,:);
img2 = img(200:end, 20:end-10);

p = AlignmentPair(img1, img2);
p.align()

figure(1); clf
p.plot()

%%

[o1, o2] = p.overlap();

figure(2); clf
implottiling({o1; o2})

%%
clc
stats = overlapStatisticsImagePair(p, 'overlap.max', 100, 'aligned', true);

stats.from
stats.to

%% non-alinged statistics (default)
clc
stats = overlapStatisticsImagePair(p, 'overlap.max', 100);

stats.from
stats.to


%% one black image
clc
p = AlignmentPair(img1, 0.01 * rand(100,100));
p2 = AlignmentPair(img1, img2);

stats = overlapStatisticsImagePair(p);
stats.from
stats.to


stats2 = overlapStatisticsImagePair(p2);
stats2.from
stats2.to


%% overlap quality

q = overlapQuality(p, 'threshold.var', 1)
q = overlapQuality(p2, 'threshold.var', 1)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Connected Components

clear all
clear classes
close all
clc
initialize

%% create artifical disconnected alignment

a = Alignment({1,2; 3,4; 5,6}')

[a.pairs.from]
[a.pairs.to]
[a.pairs.orientation]

% disconnetct to sets visa 'black' overlaps
a.pairs(3).quality = 0;
a.pairs(4).quality = 0;
a.pairs(6).quality = 0;


a.nodes

%%

c = connectedAlignments(a, 'threshold.quality', 0.5, 'reduce', false)

%%
clc

for i = 1:length(c)
   [c(i).pairs.from]
   [c(i).pairs.to]
   c(1).nodes
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Real Images

clear all
clear classes
close all
clc
initialize


%% Load Images - 4 is 'black' as is overlap image 4 -> 5

tag = './Test/Images/hESCells_Tiling_Black/11_11_p<pos>t00000001z001c1.tif';

clear imgs
imax = -Inf;
imin = Inf;
ivar = 0;
for p = 1:5
   img = double(imread(tagexpr2string(tag, 'pos', p)))';
   imgs{p,1} = img;
   imax = max(imax, max(img(:)));
   imin = min(imin, min(img(:)));
   ivar = max(ivar, var(img(:)));
end

for i =1:numel(imgs)
   imgs{i} = imrescale(imgs{i}, 'class', 'double', 'data.min', imin, 'data.max', imax);
end


figure(1); clf
implottiling(imgs, 'link', false)

{imin, imax, ivar}


figure(2); clf
img = imgs{3};
{max(img(:)), min(img(:))}
subplot(1,2,1)
hist(img(:), 256)

img = imgs{4};
{max(img(:)), min(img(:))}
subplot(1,2,2)
hist(img(:), 256)


%% Identify Valid Pairs
a = Alignment(imgs)
[a.pairs.quality]

a.overlapQuality('threshold.max', 0.15, 'overlap.max', 80);
[a.pairs.quality]

%%
clc
sub = a.connectedComponents('threshold.quality', -eps, 'reduce', true);
for i = 1:length(sub)
   sub(i)
end

%% Align components

tic
sub(1).alignPairs('overlap.max', 100, 'overlap.min', 80, 'shift.max', 20)
sub(1).optimizePairwiseShifts();
toc

figure(1); clf
sub(1).plot

%%

sub(2).alignPairs('overlap.max', 100, 'overlap.min', 80, 'shift.max', 20)

figure(2); clf
sub(2).plot

%%

figure(3); clf
implot(sub(2).images{1})




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Large Tiling

clear all
clear classes
close all
clc
initialize


%% Load Images 

tag = './Test/Images/hESCells_Tiling_Large/11_150_p<pos,6>t00000001z001c01.tif';

clear imgs
imax = -Inf;
imin = Inf;
ivar = 0;
nt = 4;
for i = 1:nt
   for j = 1:nt
      p = (j-1)*27 + i;
      img = double(imread(tagexpr2string(tag, 'pos', p)))';
      imgs{i,j} = img;
      imax = max(imax, max(img(:)));
      imin = min(imin, min(img(:)));
      ivar = max(ivar, var(img(:)));
   end
end

for i =1:numel(imgs)
   imgs{i} = imrescale(imgs{i}, 'class', 'double', 'data.min', imin, 'data.max', imax);
end


figure(1); clf
implottiling(imgs, 'link', false)

{imin, imax, ivar}


figure(2); clf
img = imgs{2};
{max(img(:)), min(img(:))}
subplot(1,2,1)
hist(img(:), 256)

img = imgs{8};
{max(img(:)), min(img(:))}
subplot(1,2,2)
hist(img(:), 256)


%% Identify Valid Pairs
a = Alignment(imgs)
[a.pairs.quality]

a.overlapQuality('threshold.max', 0.05, 'overlap.max', 80);
[a.pairs.quality]


%%
clc
sub = a.connectedComponents('threshold.quality', -eps, 'reduce', false);
for i = 1:length(sub)
   sub(i)
end
length(sub)

var2char({sub.nodes})


%% Align components

for s = 1:length(sub)
   sub(s).alignPairs('overlap.max', 120, 'overlap.min', 50, 'shift.max', 20)
   sub(s).optimizePairwiseShifts();
   
   figure(1 + s); clf
   sub(s).plot
end







%%

tag = './Test/Images/hESCells_Tiling_Large/11_150_p<pos,6>t00000001z001c01.tif';

clear imgs
imax = -Inf;
imin = Inf;
ivar = 0;
nt = 2;
for i = 1:nt
   for j = 1:nt
      p = (j-1)*27 + i;
      %img = double(imread(tagexpr2string(tag, 'pos', p)))';
      %imgs{i,j} = img;
      %imax = max(imax, max(img(:)));
      %imin = min(imin, min(img(:)));
      %ivar = max(ivar, var(img(:)));
      
      
      info{i,j} = imreadBFInfo(tagexpr2string(tag, 'pos', p));
   end
end



for i = 1:length(info(:))
   info{i}
   info{i}.imetadata
end


%%

tag = './Test/Images/hESCells_Tiling_Large/11_150_p<pos,6>t00000001z001c01.tif';

clear imgs
imax = -Inf;
imin = Inf;
ivar = 0;
nt = 2;
for i = 1:nt
   for j = 1:nt
      p = (j-1)*27 + i;
      %img = double(imread(tagexpr2string(tag, 'pos', p)))';
      %imgs{i,j} = img;
      %imax = max(imax, max(img(:)));
      %imin = min(imin, min(img(:)));
      %ivar = max(ivar, var(img(:)));
      
      
      info{i,j} = imfinfo(tagexpr2string(tag, 'pos', p));
      info{i,j}
   end
end

