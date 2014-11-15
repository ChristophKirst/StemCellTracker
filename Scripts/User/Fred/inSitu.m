%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Analyze Wild Type Colonies - Multi channel %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initialize
bfinitialize

%%
clear all
clear classes
close all
clc
verbose = true;
%initializeParallelProcessing;



%% Setup Image Source

% infer the tag expression of the files from the file folder automatically

texp = '/var/run/media/ckirst/38cc9966-c6b8-4ff9-b338-90cd43814dda/for CK/1_p<P,6>t00000001z001c<C,2>.tif'

is = ImageSourceFiles(texp);
is.printInfo

%%

[s, f] = is.tileSizeAndFormat


%% 

is.setReshape('P','UV',[13,18]);
is.setCellFormat('UvC')
is.setRange('C',1, 'U', 1:3, 'V', 1:3)

is.printInfo

% %%
% 
% im=is.data(1);
% figure(4);
% implot(im)

%% Dangerous for large dataset, Restrict the range first!
% 
% figure(1);
% is.plottiling


%% Alignment

% create 
clc;
algn = Alignment(is, 'Uv');
algn.printInfo

%var2char({[algn.pairs.from], [algn.pairs.to]})

%% Background intensity for testing overlap quality of tiles

img1 = is.data(2);
nbins = 50;
th = thresholdFirstMin(img1, 'nbins', nbins, 'delta', 1/1000 * numel(img1))
%th=2000
if verbose
   figure(3)
   hist(img1(:), nbins)
end

%% Quality of Overlap between neighbouring tiles 

% parameter: see overlapQuality
algn.calculateOverlapQuality('threshold.max', th, 'overlap.max', 120);
[algn.pairs.quality];


%% Connected Components based on overlap quality

clc
clear subalgn
subalgn = algn.connectedComponents('threshold.quality', -eps);
nsubalgn = length(subalgn);
fprintf('Alignment: found %g connected components\n', nsubalgn);

if verbose  
   var2char({subalgn.anodes})
end

%% Align Components

for s = 1:nsubalgn
   fprintf('\n\nAligning component: %g / %g\n', s, nsubalgn)
   subalgn(s).align('alignment', 'RMS', 'overlap.max', 120, 'overlap.min', 40, 'shift.max', 140);
   if verbose && s < 20 %%&& subalgn(s).nNodes < 75
      subalgn(s).printInfo 
      figure(100+s)
      
      subalgn(s).plotPreview('scale', 0.05)
   end
end


%% Compose Final Alignment

% align sub cnnected components
subalgn.alignPositions;

% merge to single alignment
algnAll = subalgn.merge;
algnAll.printInfo


%%
figure(5); clf
algnAll.plotPreview



%% Colony Detection 

% detect by resampling and closing

scale = 0.025
roi = detectROIsByClosing(algnAll, 'scale', scale, 'threshold', th, 'strel', 1, 'radius', 100, 'dilate', 100)


%% Colony 

colonies = Colony(algnAll, roi);
ncolonies = length(colonies);

figure(1); clf
colonies.plotPreview


%% Visualize

if verbose
   figure(10); clf
   for c = 1:min(ncolonies, 10)
      figure(10);
      img = colonies(c).data;
      imsubplot(5,1,c)
      implot(img)
   end
end


%%

save('./Test/Data/Colonies/colonies.mat', 'colonies')


%%

clear all
close all

verbose = true;

%%

load('./Test/Data/Colonies/colonies.mat')

ncolonies = length(colonies);

%% Change color cannel

a = colonies(1).source;
s = a.source;
s.printInfo

s.addRange('C', 2);
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
      imsubplot(5,1,c)
      implot(img)
   end
end




