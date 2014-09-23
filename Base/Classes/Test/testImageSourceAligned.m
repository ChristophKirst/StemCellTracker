%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test ImageSourceAligned Class %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc

initialize
bfinitialize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setup source
clc
texpr = tagexpr(dirr('./Test/Images/hESCells_Tiling/*.tif'), 'tagnames', {'tile'})
is = ImageSourceTagged(texpr);
is.print


% setup tiling
 
is.setTagRange('tile',  {45, 46, 47, 48,  41, 42, 43, 44,  37, 38, 39, 40,  33, 34, 35, 36})
ist = ImageSourceTiled(is, 'tileshape', [4,4]);
ist.print


%% plot tiles

tl = ist.tiles
figure(1);
implottiling(tl, 'tiling', [4,4])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ImageSourceAligned

clc
ia = ImageSourceAligned(ist);
ia.print

%% align
clc
ia.align('overlap.max', 120,  'overlap.min', 80, 'shift.max', 30);
ia.print

%%
clc
figure(2); clf
ia.plotAlignedImages


%% stitch 
 
ia.stitch('method', 'Mean')

figure(3); clf
ia.plot


%% caching

ist.clearCache();
ia.clearCache()

%%
ist.icache = 1;
img = ia.data;

figure(4)
implot(img);


%% 
ia.clearCache()
img = ia.data;

figure(4)
implot(img);


%% manual stitching
sh= ia.imageShifts;
tiles= ist.tiles;

st = stitchImages(tiles(1:3), sh(1:3), 'method', 'Mean');
figure(1); clf
implot(st)






