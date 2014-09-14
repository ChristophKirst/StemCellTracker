%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test ImageSourceTagged Class %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc

initialize
bfinitialize




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ImageSourceTiled 

%% prepare source

clc
texpr = tagexpr(dirr('./Test/Images/hESCells_Tiling/*.tif'), 'tagnames', {'tile'})
is = ImageSourceTagged(texpr);
is.print

is.itagranges

%% reudce to first 4 images

is.setTagRange('tile',  {37, 38, 33, 34})
is.print
is.cellsize

%% test id data access
clc
d = is.data('tile', 38);

figure(1);
implot(d);

%%

clc
c = is.celldata;
size(c)

figure(1);
implottiling(c)


%% ImageSourceTiled

clear all
clear classes
close all
clc

initialize
bfinitialize


%%
clc

texpr = tagexpr(dirr('./Test/Images/hESCells_Tiling/*.tif'), 'tagnames', {'tile'})
is = ImageSourceTagged(texpr);
is.tagrangesize

%%


is.setTagRange('tile',  {45, 46, 47, 48,  41, 42, 43, 44,  37, 38, 39, 40,  33, 34, 35, 36})
ist = ImageSourceTiled(is, [4,4])


%% get Tiles

tl = ist.tiles

figure(1);
implottiling(tl, 'tiling', [4,4])


%%
ts = ist.getTileSizes
ts{1}


%% align
clc
ist.align('overlap.max', 120,  'overlap.min', 80, 'shift.max', 30)

%% traditional alignment

imgs = reshape(tl, [4,4]);
imgs = cellfunc(@mat2gray, imgs);
shifts = alignImages(imgs, 'overlap.max', 120, 'overlap.min', 80, 'shift.max', 30);
figure(2); clf;
plotAlignedImages(imgs, shifts)

%% traditional pairwise alignment
clc
ip = ist.ialignment.ipairs

shifts = alignImages(imgs, 'pairs', ip, 'overlap.max', 120, 'overlap.min', 80, 'shift.max', 30);
var2char(shifts)
figure(2); clf;
plotAlignedImages(imgs, shifts)



%%

figure(3); clf;
ist.plotAlignedImages


%%

ist.imageShifts

%% stiching

clc
ist.align('overlap.max', 120,  'overlap.min', 80, 'shift.max', 30)

img = ist.stitch('method', 'Mean');

figure(4)
implot(img);


%% get the full data
ist.clearCache()

%%
ist.icache = 1;
img = ist.data;


figure(4)
implot(img);












