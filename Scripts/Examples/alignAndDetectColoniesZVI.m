%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Serial Detection of Colonies in a Zvi file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% image data is read from disk each time and not cached
% useful for large tilings that do not fit into memory

clear all
clear classes
close all
clc
initialize
bfinitialize

verbose = true;


%% Setup Image Source

fname = './Test/Images/hESCells_EpFgf52.zvi'

is = ImageSourceBF(fname);
is.printInfo

%% Tiling / Shape

%auto tiling detection
is.initializeFromTiling;
is.printInfo

%% 
clc
is.setCellDataFormat('XY', 'CUv');
is.printInfo


%% Range

is.setRange('C', {1});
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
   subalgn(s).align('alignment', 'RMS', 'overlap.max', 100, 'overlap.min', 4, 'shift.max', 140);
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
roi = detectROIsByClosing(algnAll, 'scale', scale, 'threshold', th, 'strel', 1, 'radius', 100, 'dilate', 100, 'check', true)


%% Colony 

colonies = Colony(algnAll, roi);
ncolonies = length(colonies);

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














%% Alignment

% create 
clc
isalgn = ImageSourceAligned(ist)
isalgn.setCache(false)
%var2char({[algn.pairs.from], [algn.pairs.to]})


%% Background Intensity for Overlap Quality of Tiles

img1 = isalgn.tile(3);
nbins = 50;
th = thresholdFirstMin(img1, 'nbins', nbins, 'delta', 1/1000 * numel(img1))

if verbose
   figure(3)
   hist(img1(:), nbins)
end

%% Quality of Overlap between Neighbouring Tiles 

% parameter: see overlapQuality
isalgn.overlapQuality('threshold.max', th, 'overlap.max', 120);
[isalgn.pairs.quality];


%% Connected Components based on overlap quality

clc
subalgn = isalgn.connectedComponents('threshold.quality', -eps);
nsubalgn = length(subalgn)

if verbose  
   var2char({subalgn.nodes})
end


%% Align components

clear subalgn
subalgn = isalgn.connectedComponents('threshold.quality', -eps);
nsubalgn = length(subalgn);

for s = 1:nsubalgn
   fprintf('\n\nAligning component: %g / %g\n', s, nsubalgn)
   subalgn(s).align('alignment', 'Correlation', 'overlap.max', 100, 'overlap.min', 40, 'shift.max', 140);
   if verbose && s < 20 % && subalgn(s).nnodes < 5
      subalgn(s).print%  
      figure(100+s)
      subalgn(s).plot
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Colony detection

% detect colonies in the aligned images
clc
colonies = [];

for s = 1:nsubalgn
   fprintf('\n\nDetecting colonies in component: %g / %g\n', s, nsubalgn)
   
   plt = verbose && s < 20 && subalgn(s).nnodes <= 32;
   if plt
      figure(200+s); clf
      implot(subalgn(s).data)
   end
   
   %rois = detectROIsByOpening(subalgn(s).data, 'threshold', th, 'output','ROIs', 'plot', true, 'strel', 50);
   [rois, pks] = detectROIsByPeakVolume(subalgn(s), 'radius', 200, 'dilate', 100, 'center', true, 'hmax', th, 'plot', plt);
   
   fprintf('Found %g regoins of interest\n', length(rois))
   
   for r = 1:length(rois)
      colonies = [colonies, Colony(subalgn(s), rois{r})]; %#ok<AGROW>
   end
end

ncolonies = length(colonies)
 
%% Visualize

if verbose
   figure(10); clf
   for c = 1:ncolonies
      figure(10);
      img = colonies(c).data;
      imsubplot(5,1,c)
      implot(img)
   end
end


if verbose
   figure(11); clf
   for c = 1:ncolonies
      figure(11);
      img = colonies(c).extract;
      imsubplot(5,1,c);
      implot(img)
   end
end

if verbose
   figure(12); clf
   for c = 1:ncolonies
      figure(12);
      img = colonies(c).mask;
      imsubplot(5,1,c);
      implot(img)
   end
end


%% Save Data

save('./Test/Data/Colonies/colonies.mat', colonies)














