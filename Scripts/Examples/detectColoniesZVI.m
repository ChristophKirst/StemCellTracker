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

% read info

info = imread_bf_info(fname)

%%

info = tileshapeFromZVI(info)
ts = info.icellsize

%%
is = ImageSourceTagged();
is.fromReadCommandAndFileExpression('filename', fname, 'readcommand', 'imread_bf(''<file>'', ''series'', <tile>, ''channel'', <ch>)',...
                                    'tagranges', setParameter('tile', num2cell(info.iseries), 'ch', num2cell(info.inimages)));
is.setDataFormat('py');
is.print

%
% (optional) restrict tiles to a certain subset
% - usefull to restrict alignment to a fixed z plane and channel etc
% - usefull to test on smaller subset
% is.setTagRange('tile', {9,10,17,18});
is.setTagRange('ch', {1});
is.print

% 
% if verbose
%     figure(1); clf
%     is.plot('tiling', ts)
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup Tiling

% parameter:
% tileshape   specifiess the shape of the tiles
% tileformat  specifies the ordering of the tiles, 
% in total the tiles are obtained via imuvwpermute(reshape(tiles,tileshape), tileformat)
% if the celldata in is returns the correct tiling already tileshape and tileformat are inferred from there

ist = ImageSourceTiled(is, 'tileshape', ts, 'tileformat', 'uv');
ist.print

% dont cache images if the number is large 
ist.setCache(true);

%%
if verbose && ist.ntiles <= 32
   tiles = ist.tiles;
   figure(2); clf;
   implottiling(tiles, 'link', false)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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














