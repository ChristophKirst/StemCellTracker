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

initializeParallelProcessing(12)


tic

%% Overview Image

isOverview = ImageSourceBF('/home/ckirst/Data/Science/Projects/StemCells/Experiment/Cytoo_IF/131106_RUES2_Lam521_46hBMP4_Bra_Sox2_Cdx2/Cy_w1_s3009_t1.TIF')

figure(1)
isOverview.plot

%% Setup Image Source
clc

% infer the tag expression o he files from the file folder automatically
texp = '/home/ckirst/Data/Science/Projects/StemCells/Experiment/Cytoo_IF/131106_RUES2_Lam521_46hBMP4_Bra_Sox2_Cdx2/Cy_w4_s<S>_t1.TIF'

% fns = tagExpressionToFiles(texp);
% length(fns)

%%
clc
is = ImageSourceFiles();
is.fromFileExpression(texp, 'S', 1:3008);

% set the tiling and cell formats
is.setReshape('S', 'UV', [64, 47]);
is.setCellFormat('Vu');
is.printInfo


%%

% cd = is.cell('U', 1:4, 'V', 1:5);
% figure(2); clf
% implottiling(cd)


%% full preview
% 
%is.addRange('U', 1:30, 'V', 1:30);
% 
preview = is.preview('overlap', 110, 'scale', 0.05, 'lines', true);
% 
figure(2); clf
 implot(preview)

%% restric range to some sub set

%is.addRange('U', 1:4, 'V', 1:4)
%is.setRawCellDataCaching(true);
%figure(1); clf
%is.plottiling


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Alignment

% create 
clc
algn = Alignment(is, 'UV');
algn.setCaching(false);
algn.printInfo


%% Background Intensity for Overlap Quality of Tiles

clc
img1 = algn.sourceData(12);
nbins = 100;
th = thresholdFirstMin(img1, 'nbins', nbins, 'delta', 1/1000 * numel(img1))

%th = 240;

if verbose
   figure(3); clf
   imsubplot(2,1,1)
   implot(img1)
   imsubplot(2,1,2)
   hist(img1(:), nbins)
end


th = 240;

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
   
   as = cellfun(@length, {subalgn.anodes});   
   figure(5); clf
   hist(as, 256)
   {max(as), min(as)}
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

%% 

% merge to single alignment
algnAll = subalgn.merge;
algnAll.printInfo


%%

subalgn.meanOverlapSizePrimary

%%


figure(5)
hist([subalgn.nNodes], 256)
%%
length(subalgn)

sum([subalgn.anodes])

64*

%%
figure(5); clf
algnAll.clearPreview
algnAll.plotPreview




%% Colony Detection 

% detect by resampling and closing

scale = 0.025
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
   for c = 3:min(ncolonies, 40)
      figure(10);
      img = colonies(c).data;
      imsubplot(10,4,c)
      if numel(img) > 0
         implot(img)
      end
   end
end


%%

toc

%%
save('./Test/Data/Colonies/colonies.mat', 'colonies')


%%












