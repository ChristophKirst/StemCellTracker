%%%%%%%%%%%%%%%%%%%%%%%
%%% Colony Detection %%
%%%%%%%%%%%%%%%%%%%%%%%

% no image caching is used
% usfull for laptop / computers with low ram

clear all
clear class
close all
clc

initialize
bfinitialize

verbose = true;

initializeParallelProcessing(4) % number of processors



%% Initialize Image Source

texp = './Test/Images/hESCells_Tiling_WildType/140902_RUES2_BMP4_DAPI_GFP_R_Cy5__9_p<P,6>t00000001z001c<C,2,s>.tif';

is = ImageSourceFiles(texp);
is.printInfo


%% Tiling / Format / Ranges
clc
is.setReshape('P', 'UV', [8, 19]);
is.printInfo

%% 
is.setCellDataFormat('XY', 'UvC');
is.printInfo

%%
is.setRange('U', 1:3, 'V', 1:3, 'C', 1);
is.cellSize

%% Preview

if verbose
   figure(1); clf
   is.plotPreviewStiched('overlap', 120, 'scale', 0.05, 'lines', false);
end

%% Plottiling 

% crashes matlab for large tilings!

if verbose && is.nCells < 20
   cd = is.cell;
   figure(2); clf
   implottiling(cd, 'link', false)
end


%% Alignment

algn = Alignment(is, 'UV');
algn.printInfo

%% Background Intensity for Overlap Quality of Tiles

clc
img1 = algn.sourceData(1);
nbins = 100;
th = thresholdFirstMin(img1, 'nbins', nbins, 'delta', 1/1000 * numel(img1))

%th = 240;

if verbose
   figure(3); clf
   imsubplot(2,1,1)
   implot(img1)
   imsubplot(2,1,2, 'marg', 0.1)
   hist(img1(:), nbins)
end

%% Overlap Quality between Neighbouring Tiles 

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
   
   as = cellfun(@length, {subalgn.anodes});   
   figure(5); clf
   hist(as, 256)
   {max(as), min(as)}
   slen = cellfun(@length, {subalgn.anodes});
   sum(slen)
end


%%
for s = 1:nsubalgn
   fprintf('\n\nAligning component: %g / %g\n', s, nsubalgn)
   subalgn(s).align('alignment', 'RMS', 'overlap.max', 170, 'overlap.min', 80, 'shift.max', 100);
   
   if verbose && s < 20 && subalgn(s).nNodes < 35
      subalgn(s).printInfo
      figure(100+s)
      
      %subalgn(s).plotAlignedPreview('scale', 0.5)
      subalgn(s).plotAlignedImages
   end
end

%% Merge Alignments

% align  connected components
subalgn.alignPositions;

% merge to single alignment
algnAll = subalgn.merge;
algnAll.printInfo

% mean primary overlap
subalgn.meanOverlapSizePrimary

%% Preview 

if verbose
   figure(5); clf
   algnAll.plotPreviewStiched
end


%% Colony Detection 

% detect by resampling and closing
scale = algnAll.source.previewScale;
roi = detectROIsByClosing(algnAll, 'scale', scale, 'threshold', th, 'strel', 1, 'radius', 100, 'dilate', 100)


%% Colony 

colonies = Colony(algnAll, roi);
ncolonies = length(colonies);


%%
figure(1); clf
colonies.plotPreview


%% Visualize

if verbose
   figure(10); clf
   for c = 1:min(ncolonies, 10)
      figure(10);
      img = colonies(c).data;
      imsubplot(10,4,c)
      if numel(img) > 0
         implot(img)
      end
   end
end


%%
save('./Test/Data/Colonies/colonies.mat', 'colonies')



