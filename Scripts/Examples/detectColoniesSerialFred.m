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
texp = tagExpression('/var/run/media/ckirst/LACIE SHARE/H2B citrine tracked fate end point Snail/201410026.tif_Files/*.tif', 'tagnames', {'S', 'C'})


%%
texp = '/data/Science/Projects/StemCells/Experiment/H2B citrine tracked fate end point Snail/201410026.tif_Files/201410026_p<S,3>c2.tif';
fns = tagExpressionToFiles(texp, 'c', 2)

%%

imgs = cellfunc(@(x) imreadBF(x, 'S', 1, 'C', 2), fns);

%%

imgs = reshape(imgs, 13, 12);

%%

p = alignCellPreview(imgs, 'overlap', 100, 'scale', 0.1);
figure(1); clf
implot(p)


%%





%%

 fn = tagExpressionToString('/var/run/media/ckirst/LACIE SHARE/H2B citrine tracked fate end point Snail/201410026.tif_Files/201410026_p<S,3>c<C,1>.tif', 'S', 0, 'C', 0)
 
 
%info = imreadBFInfo(fn);

%info.printInfo
 
 
 
 %%
 
 data = imreadBF(fn);
 size(data)



 %%
 
%%
%/var/run/media/ckirst/LACIE SHARE/H2B citrine tracked fate end point Snail/201410026.tif_Files/

%texp='140902_RUES2_BMP4_DAPI_GFP_R_Cy5__9_p<pos,6>t00000001z001c<ch,2>.tif'

is = ImageSourceFiles('/var/run/media/ckirst/LACIE SHARE/H2B citrine tracked fate end point Snail/201410026.tif_Files/201410026_p<S,3>c<W,1>.tif')

is.printInfo

%%

is.irangekey = rmfield(is.irangekey, 'W');


%%
is.irangekey

%%

% set the data format
is.setDataFormat('XY');

% set some keys for the colors
%is.setRangeKey('C', {'DAPI', 'GFP', 'R', 'Cy5'});

% restrict to DAPI first
is.setRange('C', 2, 'W', 2);

% set the tiling and cell formats
is.setReshape('S', 'UV', [12, 13]);

is.setCellFormat('Uv');

%is.setCaching(false);
is.printInfo

%% restric range to some sub set

is.addRange('U', 1:4, 'V', 1:3)
is.printInfo

is.nCells


is.cleatCache
is.setCaching(false)


%%

is.cellIndex
is.rawCellIndex

%%

d = is.data(1);

%%


is.plottiling


%%

% % check formatting
% cd = is.cell(1:18)
% figure(1); clf
% implottiling(cd, 'tiling', [8,2])

 
%%

% make preview
preview = is.cellPreview('overlap', 110, 'scale', 0.1, 'lines', true);
figure(2); clf
implot(preview)


%%

figure(3)
%is.plottiling % this crashes for large tilings !



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Alignment

% create 
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
      
      subalgn(s).plotAlignedPreview('scale', 0.05)
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Colony detection

% detect colonies in the aligned images
clc
colonies = [];

for s = 1:nsubalgn
   fprintf('\n\nDetecting colonies in component: %g / %g\n', s, nsubalgn)
   
   plt = verbose && s < 20 && subalgn(s).nNodes < 15;
   if plt
      figure(200+s); clf
      implot(subalgn(s).data)
   end
   
   %rois = detectROIsByOpening(subalgn(s).data, 'threshold', th, 'output','ROIs', 'plot', true, 'strel', 50);
   rois = detectROIsByPeakVolume(subalgn(s), 'radius', 100, 'dilate', 50, 'center', true, 'hmax', 500, 'plot', plt);
   
   fprintf('Found %g regoins of interest\n', length(rois))
   
   for r = 1:length(rois)
      colonies = [colonies, Colony(subalgn(s), rois{r})]; %#ok<AGROW>
   end
end

ncolonies = length(colonies)
 
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

save('./Test/Data/Colonies/colonies.mat', colonies)


%%















