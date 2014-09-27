%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Serial Detection of Colonies 
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

% infer the tag expression o he files from the file folder automatically
texp = tagexpr('./Test/Images/hESCells_Tiling_WildType/*.tif', 'tagnames', {'tile', 'ch'})
%texp='140902_RUES2_BMP4_DAPI_GFP_R_Cy5__9_p<pos,6>t00000001z001c<ch,2>.tif'

is = ImageSourceTagged(texp);
is.setDataFormat('py');
is.print

%
% (optional) restrict tiles to a certain subset
% - usefull to restrict alignment to a fixed z plane and channel etc
% - usefull to test on smaller subset
% is.setTagRange('tile', {9,10,17,18});
is.setTagRange('ch', {1}, 'tile', num2cell((1:(8*2)) + 8));
is.print

% 
% if verbose
%     figure(1); clf
%     a=[9,10,17,18]
%     is.plot('tile', a, 'tiling', [2,2])
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup Tiling

% parameter:
% tileshape   specifiess the shape of the tiles
% tileformat  specifies the ordering of the tiles, 
% in total the tiles are obtained via imuvwpermute(reshape(tiles,tileshape), tileformat)
% if the celldata in is returns the correct tiling already tileshape and tileformat are inferred from there

ist = ImageSourceTiled(is, 'tileshape', [8,2], 'tileformat', 'uv');
ist.print

% dont cache images if the number is large 
ist.setCache(true);

if verbose && ist.ntiles < 25 
   tiles = ist.tiles;
   figure(2); clf;
   implottiling(tiles)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Alignment

% create 
clc
isalgn = ImageSourceAligned(ist)
isalgn.setCache(false)
%var2char({[algn.pairs.from], [algn.pairs.to]})


%% Background Intensity for Overlap Quality of Tiles

img1 = isalgn.tile(2);
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
   subalgn(s).align('alignment', 'Correlation', 'overlap.max', 100, 'overlap.min', 4, 'shift.max', 140);
   if verbose && s < 20 && subalgn(s).nnodes < 5
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
   fprintf('\n\nDetecting colonies in compnent: %g / %g\n', s, nsubalgn)
   
   plt = verbose && s < 20 && subalgn(s).nnodes < 5;
   if plt
      figure(200+s); clf
      implot(subalgn(s).data)
   end
   
   %rois = findROIsByOpening(subalgn(s).data, 'threshold', th, 'output','ROIs', 'plot', true, 'strel', 50);
   rois = findROIsByPeakVolume(subalgn(s), 'radius', 100, 'dilate', 50, 'center', true, 'hmax', 500, 'plot', plt);
   
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
      imsubplot(10,5,c)
      implot(img)
   end
end


%% Save Data

save('./Test/Data/Colonies/colonies.mat', colonies)


















