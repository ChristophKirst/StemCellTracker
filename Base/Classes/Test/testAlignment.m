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
texpr = tagExpression(dirr('./Test/Images/hESCells_Tiling/*.tif'), 'tagnames', {'S'})

is = ImageSourceFiles(texpr);
is.printInfo


%% plot tiles to see how to align


tl = is.cell
figure(1);
implottiling(tl, 'tiling', [4,4])


%% reformat 
is.setReshape('S', 'UV', [4,4]);
is.setCellFormat('Uv');
is.printInfo

%%
% preview

figure(2); clf;
is.alignCellPreview


%% altenative tile

is.plottiling



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Alignement


clear all
clear classes
close all
clc

initialize
bfinitialize

obj = Alignment

texpr = tagExpression(dirr('./Test/Images/hESCells_Tiling/*.tif'), 'tagnames', {'S'})

is = ImageSourceFiles(texpr);
is.setReshape('S', 'UV', [4,4]);
is.setCellFormat('Uv');
is.printInfo

%%
clc
ia = Alignment(is);
ia.printInfo

%% align
clc
ia.align('overlap.max', 120,  'overlap.min', 80, 'shift.max', 30);
ia.printInfo

%%
clc
figure(2); clf
ia.plotAlignedImages


%% stitch 
 
ia.stitch('method', 'Mean');

figure(3); clf
ia.plot


%%
is.setCaching(true);
img = ia.data;

figure(4)
implot(img);


%% 

img = ia.data;

figure(4)
implot(img);


%% manual stitching
sh= ia.imageShifts;
tiles= is.cell;

st = stitchImages(tiles(1:3), sh(1:3), 'method', 'Mean');
figure(1); clf
implot(st)






