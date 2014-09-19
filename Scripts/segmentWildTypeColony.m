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
texp = tagexpr('./Test/Images/hESCells_Tiling_Large/*.tif', 'tagnames', 'tile')

%%
% setup a tagged image source
is = ImageSourceTagged(texp)
is.ntags


% (optional) restrict tiles to a certain subset
% - usefull to restrict alignment to a fixed z plane and channel etc
% - usefull to test on smaller subset
is.setTagRange('tile', {1,2,3,4, 28, 29, 30, 31, 45,46,47,48});
is.ntags

% if verbose
%     figure(1); clf
%     is.plot('tile', {1,2,3,4, 5,6,7,8}, 'tiling', [4,2])
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup Tiling

% parameter:
% tileshape   specifiess the shape of the tiles
% tileformat  specifies the ordering of the tiles, 
% in total the tiles are obtained via imuvwpermute(reshape(tiles,tileshape), tileformat)
% if the celldata in is returns the correct tiling already tileshape and tileformat are inferred from there

ist = ImageSourceTiled(is, 'tileshape', [4,3], 'tileformat', 'uy');
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
algn = Alignment(ist);
%var2char({[algn.pairs.from], [algn.pairs.to]})


%% Background intensity for testing overlap quality of tiles

img1 = ist.tile(1);
nbins = 50;
th = thresholdFirstMin(img1, 'nbins', nbins, 'delta', 1/1000 * numel(img1))

if verbose
   figure(3)
   hist(img1(:), nbins)
end

%% Quality of Overlap between neighbouring tiles 

% parameter: see overlapQuality
algn.overlapQuality('threshold.max', th, 'overlap.max', 120);
[algn.pairs.quality]


%% Connected Components based on overlap quality

clc
subalgn = algn.connectedComponents('threshold.quality', -eps);
nsubalgn = length(subalgn)

if verbose 
   %for i = 1:nsubalgn
   %   subalgn(i)
   %end   
   var2char({subalgn.nodes})
end


%% Align components

for s = 1:nsubalgn
   subalgn(s).alignPairs('overlap.max', 120, 'overlap.min', 100, 'shift.max', 10);
   subalgn(s).optimizePairwiseShifts();
   
   if verbose && s < 8
      figure(4)
      imsubplot(4, 2, s);
      subalgn(s).plotAlignedImages
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Colony detection

% detect colonies in the aligned images







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data analysis















