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
is.reformatCellFormat('Uv');
is.printInfo

%%

is.rawCellSize

%%
% preview

clc
figure(2); clf;
is.cellPreview


%% altenative tile

is.plottiling


%%

is.cellIndex


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Alignement


clear all
clear classes
close all
clc

initialize
bfinitialize

obj = Alignment;

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









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Alignment

clear all
close all
clear classes
clc

initialize
verbose = true;

%% ImageSource setup 

texp = tagExpression('./Test/Images/hESCells_GCamp_Vignetting/*.tif', 'tagnames', {'S', 'C'})

is = ImageSourceFiles(texp);
is.setReshape('S', 'UV', [3,3]);
is.setRawCellDataCaching(true);
is.setCellFormat('Uv');
is.printInfo


%% create 
clc
algn = Alignment(is);
algn.setCaching(false);
algn.printInfo


%% Background Intensity for Overlap Quality of Tiles

img1 = algn.sourceData(2);
nbins = 50;
th = thresholdFirstMin(img1, 'nbins', nbins, 'delta', 1/1000 * numel(img1))

if verbose
   figure(3); clf
   hist(img1(:), nbins)
end

%% Quality of Overlap between Neighbouring Tiles 

% parameter: see overlapQuality
algn.calculateOverlapQuality('threshold.max', th, 'overlap.max', 120);

hist(algn.overlapQuality, 256)


%% Align components
clc
clear subalgn
subalgn = algn.connectedComponents('threshold.quality', -eps);
nsubalgn = length(subalgn);
fprintf('Alignment: found %g connected components\n', nsubalgn);

if verbose  
   var2char({subalgn.anodes})
end


%%
for s = 1:nsubalgn
   fprintf('\n\nAligning component: %g / %g\n', s, nsubalgn)
   subalgn(s).align('alignment', 'Correlation', 'overlap.max', 100, 'overlap.min', 4, 'shift.max', 140);
   if verbose && s < 20 %%&& subalgn(s).nNodes < 75
      subalgn(s).printInfo 
      figure(100+s)
      
      subalgn(s).plotAlignedPreview('scale', 0.1)
   end
end








