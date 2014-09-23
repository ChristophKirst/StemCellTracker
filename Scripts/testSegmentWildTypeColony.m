%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Colony Analysis Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc
initialize
bfinitialize


verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup Image Source

% infer the tag expression o he files from the file folder automatically
texp = tagexpr('/home/ckirst/Desktop/test/*.tif', 'tagnames', 'tile')


%%
is = ImageSourceTagged(texp);
is.setDataFormat('xy');
is.print

%%
% (optional) restrict tiles to a certain subset
% - usefull to restrict alignment to a fixed z plane and channel etc
% - usefull to test on smaller subset
%is.setTagRange('tile', {1,2,3,4,5,6, 28,29,30,31,32,33, 55,56,57,58,59,60});
%is.print


if verbose
    figure(1); clf
    is.plot('tiling', [2,3])
end


%%

ii = imread_bf_info(tagexpr2string(texp, 'tile', 1))
ii.imetadata


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup Tiling

% parameter:
% tileshape   specifiess the shape of the tiles
% tileformat  specifies the ordering of the tiles, 
% in total the tiles are obtained via imuvwpermute(reshape(tiles,tileshape), tileformat)
% if the celldata in is returns the correct tiling already tileshape and tileformat are inferred from there

ist = ImageSourceTiled(is, 'tileshape', [2,3], 'tileformat', 'uv');
ist.print

% dont cache images if the number is large 
ist.icache = ist.ntiles < 100; 

if verbose && ist.ntiles < 20 
   tiles = ist.tiles;
   figure(2); clf;
   implottiling(tiles)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Alignment

% create 
clc
isalgn = ImageSourceAligned(ist)
%var2char({[algn.pairs.from], [algn.pairs.to]})


%% Background intensity for testing overlap quality of tiles

img1 = isalgn.tile(2);
nbins = 50;
th = thresholdFirstMin(img1, 'nbins', nbins, 'delta', 1/1000 * numel(img1))

if verbose
   figure(3)
   hist(img1(:), nbins)
end

%% Quality of Overlap between neighbouring tiles 

% parameter: see overlapQuality
isalgn.overlapQuality('threshold.max', th, 'overlap.max', 120);
[isalgn.pairs.quality];


%% Connected Components based on overlap quality

clc
subalgn = isalgn.connectedComponents('threshold.quality', -eps);
nsubalgn = length(subalgn)

if verbose 
%    for i = 1:nsubalgn
%       subalgn(i);
%    end   
   var2char({subalgn.nodes})
end


%% Align components


for s = 1:nsubalgn
   subalgn(s).align('alignment', 'Correlation', 'overlap.max', 110, 'overlap.min', 50, 'shift.max', 80);

   if verbose && s < 12
      subalgn(s).print
      figure(4); clf
      imsubplot(1, 1, s);
      subalgn(s).clearCache;
      subalgn(s).plot
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Colony detection

% detect colonies in the aligned images
clc
colonies = [];
for s = 1:nsubalgn
   figure(5)
   rois = findROIsByOpening(subalgn(s).data, 'threshold', th, 'output','ROIs', 'plot', true, 'strel', 50);
   
   for r = 1:length(rois)
       colonies = [colonies, Colony(subalgn(s), rois(r))]; %#ok<AGROW>
   end
end

ncolonies = length(colonies)

%%

figure(1); clf
for c = 1:ncolonies
   img = colonies(c).data;
   imsubplot(10,10,c)
   implot(img)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data analysis















