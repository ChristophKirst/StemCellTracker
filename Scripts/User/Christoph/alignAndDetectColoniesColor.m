%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Detection of Colonies 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% image data is read from disk once and cached
% useful for fast processing if sufficent memory is available

clear all
clear classes
close all
clc
initialize
bfinitialize

verbose = true;

% initializeParallelProcessing(4) % number of processors


%% Setup Image Source
clc

% infer the tag expression o he files from the file folder automatically
texp = tagExpression('/data/Science/Projects/StemCells/Experiment/H2B citrine tracked fate end point Snail/201410026.tif_Files/*.tif', 'tagnames', {'P', 'B'})
%texp='140902_RUES2_BMP4_DAPI_GFP_R_Cy5__9_p<pos,6>t00000001z001c<ch,2>.tif'

%%
is = ImageSourceFiles(texp)
is.printInfo

%%
texp = '/data/Science/Projects/StemCells/Experiment/H2B citrine tracked fate end point Snail/201410026.tif';

is = ImageSourceBF(texp)
is.printInfo

%%

is.setRange('C', 1);
is.printInfo

d = is.data(7);
figure(1)
implot(d)

%% Tiling

% set the tiling
is.setReshape('S', 'UV', [8, 19]);
is.setCellFormat('UvC');
is.setCaching(true);

is.printInfo

%% Color Keys

% set some keys for the colors
is.setKey('C', {'DAPI', 'GFP', 'R', 'Cy5'});

% restrict to DAPI first
clc
is.setRange('C', 'DAPI');
is.range
is.printInfo

%% Range

is.addRange('U', 1:8, 'V', 1:8);
is.printInfo


%% Preview

if verbose
   figure(1); clf
   is.plotPreviewStiched('overlap', 120, 'scale', 0.05, 'lines', false);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Alignment

% create 
clc
algn = Alignment(is, 'UV');
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

%% Quality of Overlap between neighbouring Tiles 

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
      figure(100+s); clf
      
      subalgn(s).plotPreviewStiched('scale', 0.05)
   end
end

%% Align Components
clc
subalgn.alignPositions;

% merge to single alignment

algnAll = subalgn.merge;
algnAll.printInfo

%%
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
   for c = 1:min(ncolonies, 20)
      figure(10);
      img = colonies(c).data;
      imsubplot(5,4,c)
      if numel(img) > 0
         implot(img)
      end
   end
end


%% Save

% make sure to clear cache before saving
colonies.clearCache();
save('./Test/Data/Colonies/colonies.mat', 'colonies')


%% Load

clear all %#ok<CLSCR>
close all
clc

verbose = true;

load('./Test/Data/Colonies/colonies.mat')


%% Change color cannel

as = colonies(1).source;
is = as.source;
is.printInfo

is.addRange('C', 'Cy5');
is.fileTagRange

%%
clc
figure(5); clf;
as.setPreviewScale(0.1);
as.plotPreviewStiched


%%
figure(6); clf
colonies.plotPreview

%%

if verbose
   figure(10); clf
   for c = 1:min(ncolonies, 10)
      figure(10);
      img = colonies(c).data;
      imsubplot(5,2,c)
      implot(img)
   end
end






