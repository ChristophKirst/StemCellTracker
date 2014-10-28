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
clc

% infer the tag expression o he files from the file folder automatically
texp = tagExpression('./Test/Images/hESCells_Tiling_WildType/*.tif', 'tagnames', {'S', 'C'});
%texp='140902_RUES2_BMP4_DAPI_GFP_R_Cy5__9_p<pos,6>t00000001z001c<ch,2>.tif'

is = ImageSourceFiles(texp);
is.printInfo

%%
% set the tiling
is.setReshape('S', 'UV', [8, 19]);
is.setCellFormat('UvC');
is.printInfo


%%

% set some keys for the colors
is.setKey('C', {'DAPI', 'GFP', 'R', 'Cy5'});


%%
% restrict to DAPI first
clc
is.setRange('C', 'DAPI');
is.range
is.printInfo

%% restric range to some sub set

is.addRange('U', 1:4, 'V', 1:3);
is.printInfo


%%

% % check formatting
% cd = is.cell(1:18)
% figure(1); clf
% implottiling(cd, 'tiling', [8,2])

 
%%

% make preview
preview = is.previewStitchedCells('overlap', 110, 'scale', 0.1, 'lines', true);
figure(2); clf
implot(preview)


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
      figure(100+s)
      
      subalgn(s).plotPreview('scale', 0.05)
   end
end

%% Align Components

subalgn.alignPositions;


%% merge to single alignment

algnAll = subalgn.merge;
algnAll.printInfo


%%
figure(5); clf
algnAll.plotPreview


%% Colony Detection from Preview

p = algnAll.preview;
s = algnAll.previewScale;

roi = detectROIsFromResampledImage(p, s, 'threshold', th, 'strel', 10, 'radius', 10, 'dilate', 8);

if verbose
   figure(1); clf
   implot(p)
   hold on
   roi.plotRescaled(s)
end

fprintf('Found %g regions of interest\n', length(roi))


%% Colony 

colonies = Colony(algnFull, roi);
colonies.plotPreview
ncolonies = length(colonies);



%% Visualize

if verbose
   figure(10); clf
   for c = 1:min(ncolonies, 10)
      figure(10);
      img = colonies(c).data;
      imsubplot(10,5,c)
      implot(img)
   end
end


%%

save('./Test/Data/Colonies/colonies.mat', 'colonies')


%% Change color cannel

a = colonies(1).source;
s = a.source;
s.printInfo

s.addRange('C', 4);
s.fileTagRange

%%

figure(5); clf;
a.clearPreview
a.plotPreview



%%

if verbose
   figure(10); clf
   for c = 1:min(ncolonies, 10)
      figure(10);
      img = colonies(c).data;
      imsubplot(10,5,c)
      implot(img)
   end
end














