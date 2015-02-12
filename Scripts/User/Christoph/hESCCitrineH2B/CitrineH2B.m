%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Citrine Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clear classes
close all
clc
initialize
bfinitialize
verbose = true;
initializeParallelProcessing(12) % number of processors
%%
texpmip = '/data/Science/Projects/StemCells/Experiment/Unsorted/unsorted/Voyager/h2b_citrine/MIP/W1F<F,3>T<T,4>Z<Z,2>C<C,1>.tif';
%%
ismip = ImageSourceFiles(texpmip);
ismip.printInfo
%%
ismip.resetRange;
ismip.setReshape('F', 'UV', [10,15]);
ismip.setCellFormat('UvTZC');
ismip.printInfo
%%
ismip.setRange('T', 1)
ismip.printInfo
% %%
%
% ismip.setRange('T', 1, 'U', 1:4, 'V', 1:3)
% ismip.printInfo
%
% figure(1); clf;
% cd = ismip.cell
% implottiling(cd)
%%
figure(1); clf
ismip.plotPreviewStiched('overlap', 80, 'scale', 0.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Alignment
% create
clc
algn = Alignment(ismip, 'UV');
algn.setCaching(false);
algn.printInfo
%% Background Intensity for Overlap Quality of Tiles
d = ismip.data(2);
figure(2); clf;
imsubplot(2,1,1)
implot(d)
dd = d(:);
size(dd)
dd(dd < 2200) = [];
size(dd)
subplot(1,2,2)
hist(dd(:), 150)
th = thresholdFirstMin(dd, 'nbins', 150, 'delta', 100)
th = 2500;
%%
figure(30); clf
d = ismip.data(3);
implot(double(d > th))
%% Quality of Overlap between neighbouring Tiles
% parameter: see overlapQuality
algn.calculateOverlapQuality('threshold.max', th, 'overlap.max', 100);
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
subalgn(s).align('alignment', 'RMS', 'overlap.max', 95, 'overlap.min', 70, 'shift.max', 10);
if verbose && s < 20 %%&& subalgn(s).nNodes < 75
subalgn(s).printInfo
figure(100+s); clf
subalgn(s).plotPreviewStiched('scale', 0.5)
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
algnAll.plotPreviewStiched('scale', 0.5)
end
%% Colony Detection
% detect by resampling and closing
scale = algnAll.source.previewScale;
roi = detectROIsByClosing(algnAll, 'scale', scale, 'threshold', 1.0* th, 'strel', 1, 'radius', 90, 'dilate', 90, 'check', true)
%% Colony
colonies = Colony(algnAll, roi);
ncolonies = length(colonies);
figure(1); clf
colonies.plotPreview
%% bounding boxes
clc
coloniesC = colonies.copy;
for i = 1:length(colonies)
roi = colonies(i).roi;
%[c, r]= fitCircle(roi.p);
%roic = ROIDisk(c,r);
roic = roi.boundingBox;
coloniesC(i).roi = roic;
end
figure(2); clf
coloniesC.plotPreview
%% take full colonies only
figure(8); clf
hist(coloniesC.volume, 20)
coloniesD = coloniesC.copy;
ids = coloniesC.volume > 400000
coloniesD = coloniesD(ids)
% remove lower border
rois = [coloniesD.roi];
pp = [rois.p1];
ids = pp(2,:) > 0;
coloniesD = coloniesD(ids);
figure(9);
coloniesD.plotPreview
%% Save
% make sure to clear cache before saving
coloniesD.clearCache();
save('/home/ckirst/Data/Science/Projects/StemCells/Experiment/H2BCitrineTrackingAndFate/Analysis/coloniesMIP.mat', 'colonies', 'coloniesD')
ismip.clearCache();
save('/home/ckirst/Data/Science/Projects/StemCells/Experiment/H2BCitrineTrackingAndFate/Analysis/sourceMIP.mat', 'ismip')
%% Load
clear all %#ok<CLSCR>
close all
clc
verbose = true;
%%
load('/home/ckirst/Data/Science/Projects/StemCells/Experiment/H2BCitrineTrackingAndFate/Analysis/coloniesMIP.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract Data
cid = 30;
figure(9);
coloniesD(cid).plotPreview
%%
dd = coloniesD(cid).data;
figure(6); clf;
implot(dd)
%%
%% manual extraction as BF is slow for 1000000 files
texp = '/data/Science/Projects/StemCells/Experiment/unsorted/Voyager/h2b_citrine/Image/W1F<F,3>T<T,4>Z<Z,2>C<C,1>.tif';
texpout = '/home/ckirst/Data/Science/Projects/StemCells/Experiment/H2BCitrineTrackingAndFate/Analysis/ColonyData/ColonyW<Q,2>T<T,4>Z<Z,2>C<C,1>.tif';
cid = 10;
fnodes = coloniesD(cid).nodes;
fnodesraw = ismip.rawCellIndex(fnodes);
as = coloniesD(cid).source;
roi = coloniesD(cid).roi;
[ids, shids] = as.nodesFromROI(roi);
shifts = as.imageShifts;
shifts = shifts(shids);
cc = 1;
verbose = true;
for tt = 1:520
for zz = 1:12
csubids = ismip.cellIndexToSubIndex(fnodes);
csubids = csubids(:, [1,2]);
csubids = csubids - repmat(min(csubids), size(csubids,1), 1) + 1;
si = max(csubids)
csubids = num2cell(csubids);
n = length(fnodes);
imgs = cell(si);
for i = 1:n
fn = tagExpressionToString(texp, 'C', cc, 'T', tt, 'Z', zz, 'F', fnodesraw(i));
img = imread(fn);
imgs{csubids{i,:}} = imfrmtPermute(img, 'XY', 'Yx');
%imgs{csubids{i,:}} = imreadBF(fn);
end
%% stitch / extract
st = stitchImagesByInterpolate(imgs, shifts);
% figure(20);
% implot(st)
% size(st)
roiS = roi.copy;
roiS.shift(-shifts{1}+1)
dat = roiS.extractData(st);
if verbose
figure(21);
implot(dat)
drawnow
max(dat(:))
class(dat)
end
%% save
fnout = tagExpressionToString(texpout, 'Q', cid, 'C', cc, 'T', tt, 'Z', zz);
imwriteTIFF(int32(dat), fnout)
end
end
%%
%
% tdat = imread(fnout);
% size(tdat)
% max(tdat(:))
% min(tdat(:))
%
% figure(1); clf
% implot(tdat)