%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Citrine Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc

initialize
bfinitialize

initializeParallelProcessing(12) % number of processors

verbose = true;

%%

texpmip = '/data/Science/Projects/StemCells/Experiment/Unsorted/unsorted/Voyager/SMAD4_colonies/MIP/W1F<F,3>T<T,4>Z<Z,2>C<C,1>.tif';


%%
ismip = ImageSourceFiles(texpmip);

%%
ismip.printInfo

%%
ismip.resetRange;
ismip.setReshape('F', 'UVP', [6,6,2]);
ismip.setCellFormat('UvTZCP');
ismip.printInfo

%%

ismip.setRange('T', 1, 'P', 1, 'C', 2)
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

figure(1); clf;
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

%%
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
roi = detectROIsByClosing(algnAll, 'scale', scale, 'threshold', 1.25* th, 'strel', 1, 'radius', 20, 'dilate', 25, 'check', true)


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

coloniesD =  coloniesD(ids);


figure(9);
coloniesD.plotPreview






%% Estimate Flatfield / Background


fns = '/data/Science/Projects/StemCells/Experiment/Unsorted/unsorted/Voyager/SMAD4_colonies/MIP/W1F<F,3>T<T,4>Z<Z,2>C<C,1>.tif';

imgf = zeros([512 512]);
imgb = inf * ones([512 512]);

k = 1;

for t = 1:119 %119
   fprintf('t = %d\n', t);
   imgs = cell(72,1);
   
   parfor f = 1:72 % 72
      img = double(imread(tagExpressionToString(fns, 'C', 2, 'T', t, 'F', f, 'Z', 0)));
      imgs{f} = imfrmtPermute(img, 'XY', 'Yx');
   end

   imga = cat(3, imgs{:});
   imgmean = mean(imga, 3);

   imgb = min(imgb, min(imga,[], 3));
   imgf = ((k-1) * imgf + imgmean)/k;
      
   k = k + 1;

end


%%
figure(13); clf;

imgfs = filterGaussian(imgf, 100);

implottiling({mat2gray(imgf);mat2gray(imgb); mat2gray(filterGaussian(imgb, 15)); mat2gray(imgfs)});

figure(14); clf;

img = double(imread(tagExpressionToString(fns, 'C', 1, 'T', 16, 'F', 10, 'Z', 0)));
img = imfrmtPermute(img, 'XY', 'Yx');
implottiling({mat2gray(img); mat2gray((img - imgb)./ (imgfs - imgb))})



%% Estimate for each Z separately - channel 1 !

fns = '/data/Science/Projects/StemCells/Experiment/Unsorted/unsorted/Voyager/SMAD4_colonies/Image/W1F<F,3>T<T,4>Z<Z,2>C<C,1>.tif';


imgbZ = cell(1,12);
imgfZ = cell(1,12);

parfor z = 1:12

   imgf = zeros([512 512]);
   imgb = inf * ones([512 512]);

   k = 1;

   for t = 1:119 %119
      fprintf('t = %d\n', t);
      imgs = cell(72,1);
   
      for f = 1:72 % 72
         img = double(imread(tagExpressionToString(fns, 'C', 1, 'T', t, 'F', f, 'Z', z)));
         imgs{f} = imfrmtPermute(img, 'XY', 'Yx');
      end

      imga = cat(3, imgs{:});
      imgmean = mean(imga, 3);

      imgb = min(imgb, min(imga,[], 3));
      imgf = ((k-1) * imgf + imgmean)/k;
      
      k = k + 1;

   end
   
   imgbZ{z} = imgb;
   imgfZ{z} = imgf;
      
end

%%
for z = 1:12

   imgfZs{z} = filterGaussian(imgfZ{z}, 100);

   imgf = imgfZ{z};
   imgb = imgbZ{z};
   
   img = double(imread(tagExpressionToString(fns, 'C', 1, 'T', 16, 'F', 10, 'Z', z)));
   img = imfrmtPermute(img, 'XY', 'Yx');
   
   figure(13 + z); clf;
   implottiling({mat2gray(imgf), mat2gray(imgb); mat2gray(imgfZs{z}), 0* img; 
                 mat2gray(img), mat2gray((img - imgb)./ (imgfZs{z} - imgb))})

end

%% Save Correction Images:

save('/data/Science/Projects/StemCells/Experiment/Unsorted/SMAD4_Translocation/correctionC1.mat', 'imgfZ', 'imgbZ', 'imgfZs')


%% Estimate for each Z separately - channel 2 !


fns = '/data/Science/Projects/StemCells/Experiment/Unsorted/unsorted/Voyager/SMAD4_colonies/Image/W1F<F,3>T<T,4>Z<Z,2>C<C,1>.tif';


imgbZ = cell(1,12);
imgfZ = cell(1,12);

parfor z = 1:12

   imgf = zeros([512 512]);
   imgb = inf * ones([512 512]);

   k = 1;

   for t = 1:119 %119
      fprintf('t = %d\n', t);
      imgs = cell(72,1);
   
      for f = 1:72 % 72
         img = double(imread(tagExpressionToString(fns, 'C', 2, 'T', t, 'F', f, 'Z', z)));
         imgs{f} = imfrmtPermute(img, 'XY', 'Yx');
      end

      imga = cat(3, imgs{:});
      imgmean = mean(imga, 3);

      imgb = min(imgb, min(imga,[], 3));
      imgf = ((k-1) * imgf + imgmean)/k;
      
      k = k + 1;

   end
   
   imgbZ{z} = imgb;
   imgfZ{z} = imgf;
      
end

%%
for z = 1:12

   imgfZs{z} = filterGaussian(imgfZ{z}, 100);

   imgf = imgfZ{z};
   imgb = imgbZ{z};
   
   img = double(imread(tagExpressionToString(fns, 'C', 2, 'T', 16, 'F', 10, 'Z', z)));
   img = imfrmtPermute(img, 'XY', 'Yx');
   
   figure(13 + z); clf;
   implottiling({mat2gray(imgf), mat2gray(imgb); mat2gray(imgfZs{z}), 0* img; 
                 mat2gray(img), mat2gray((img - imgb)./ (imgfZs{z} - imgb))})

end

%% Save Correction Images:

save('/data/Science/Projects/StemCells/Experiment/Unsorted/SMAD4_Translocation/correctionC2.mat', 'imgfZ', 'imgbZ', 'imgfZs')



%% Plot Correction images


%%

load('/data/Science/Projects/StemCells/Experiment/Unsorted/SMAD4_Translocation/correctionC2.mat')

fns = '/data/Science/Projects/StemCells/Experiment/Unsorted/unsorted/Voyager/SMAD4_colonies/Image/W1F<F,3>T<T,4>Z<Z,2>C<C,1>.tif';

for z = 4:7

   imgf = imgfZ{z};
   imgb = imgbZ{z};
   imgbf = filterGaussian(imgb, 100);
   
   img = double(imread(tagExpressionToString(fns, 'C', 2, 'T', 16, 'F', 8, 'Z', z)));
   img = imfrmtPermute(img, 'XY', 'Yx');
   
   figure(13 + z); clf;
   implottiling({mat2gray(imgf), mat2gray(imgfZs{z}), mat2gray(imgb), mat2gray(imgbf);
                 mat2gray(img), mat2gray(img-imgbf), mat2gray((img)./ (imgfZs{z})), mat2gray((img - imgbf)./ (imgf - imgbf)), }')

end

%%
figure(18); clf;
imgbt = imgb;
imgbt(imgb < 1850) = 1850;
implot(mat2gray(imgbt))


%%




%% Save

% make sure to clear cache before saving
%coloniesD.clearCache();
%save('/home/ckirst/Data/Science/Projects/StemCells/Experiment/SMAD4_colonies/Analysis/coloniesMIP.mat', 'colonies', 'coloniesD')
%ismip.clearCache();
%save('/home/ckirst/Data/Science/Projects/StemCells/Experiment/SMAD4_colonies/Analysis/sourceMIP.mat', 'ismip')

%% Load
% 
% clear all %#ok<CLSCR>
% close all
% clc
% 
% verbose = true;

%load('/home/ckirst/Data/Science/Projects/StemCells/Experiment/Unsorted/H2BCitrineTrackingAndFate/Analysis/coloniesMIP.mat')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract Data 


cid = 2;

figure(9);
coloniesD(cid).plotPreview

colonies(1).source.asource


%%
clc
dd = coloniesD(cid).data;
figure(6); clf;
implot(dd)


%%
% 
% texpmip = '/data/Science/Projects/StemCells/Experiment/Unsorted/unsorted/Voyager/SMAD4_colonies/MIP/W1F<F,3>T<T,4>Z<Z,2>C<C,1>.tif';
% 
% %ismip = ImageSourceFiles(texpmip);
% %
% ismip.resetRange;
% ismip.setReshape('F', 'UVW', [6,6,2]);
% ismip.setCellFormat('UvTZCW');
% ismip.printInfo
% 
% %%
% ismip.setRange('T', 1, 'W', 1, 'C', 2)
% ismip.printInfo


%% manual extraction as BF is slow for 1000000 files 

texp =  '/data/Science/Projects/StemCells/Experiment/unsorted/Voyager/SMAD4_colonies/Image/W1F<F,3>T<T,4>Z<Z,2>C<C,1>.tif';
texpout = '/data/Science/Projects/StemCells/Experiment/SMAD4_Translocation/ColonyData/ColonyW<Q,2>T<T,4>Z<Z,2>C<C,1>.tif';

texp = '/data/Science/Projects/StemCells/Experiment/Unsorted/unsorted/Voyager/SMAD4_colonies/Image/W1F<F,3>T<T,4>Z<Z,2>C<C,1>.tif';
texpout = '/data/Science/Projects/StemCells/Experiment/Unsorted/SMAD4_Translocation/ColonyData/ColonyW<Q,2>T<T,4>Z<Z,2>C<C,1>.tif';

cid = 2;

as = coloniesD(cid).source;
roi = coloniesD(cid).roi;
[ids, shids] = as.nodesFromROI(roi);

shifts = as.imageShifts;
shifts = shifts(shids);

cc = 2;

load(tagExpressionToString('/data/Science/Projects/StemCells/Experiment/Unsorted/SMAD4_Translocation/correctionC<C,1>.mat', 'C', cc));

verbose = true;

%addpath(genpath('/home/ckirst/Science/Simulation/Matlab/StemCellTracker/Base/Utils/External/nu_corrector'));

%%
parfor tt = 1:119
   
   %%
   fprintf('Time: %d\n', tt)
   
   
   %%
   for zz = 1:12

      csubids = ismip.cellIndexToSubIndex(ids);
      csubids = csubids(:, [1,2]);
      csubids = csubids - repmat(min(csubids), size(csubids,1), 1) + 1;
      si = max(csubids);
      csubids = num2cell(csubids);
      
      n = length(ids);
      imgs = cell(si);

      
      %imgb = imgbZ{z};
      %imgfs = imgfZs{z};

      for i = 1:n
         fnoderaw =  ismip.rawCellIndex(ids(i));
         
         fn = tagExpressionToString(texp, ismip.rawCellIndexToRange(fnoderaw), 'C', cc, 'T', tt, 'Z', zz);
         img = imread(fn);
         img = imfrmtPermute(img, 'XY', 'Yx');

         imgs{csubids{i,:}} = (double(img) - double(imgbZ{zz})) ./ (imgfZs{zz} - imgbZ{zz} );
         %imgs{csubids{i,:}} = imreadBF(fn);

         % correct via de-vignetting
         %imgr = mat2gray(img);

         %[img,~]=vignCorrection_nonPara(round(255 * imgr));%you may need to input the second parameter and adapt it to the vignetting in the input image
         %imgs{csubids{i,:}} = img;
      end

      % stitch / extract
       st = stitchImages(imgs, shifts, 'method', 'Interpolate');
%       figure(20);      clf
%       %imsubplot(6,2,zz)
%       implot(st, 'color.scale', [0,11])
%       drawnow
      %size(st)
      
      roiS = roi.copy;
      roiS.shift(-shifts{1}+1)
      
      dat = roiS.extractData(st);
      
      %% save

      fnout = tagExpressionToString(texpout, 'Q', cid, 'C', cc, 'T', tt, 'Z', zz);
      imwriteTIFF(int32(dat / 11 * 2^16), fnout)
   end
end
  

%%
texp = '/data/Science/Projects/StemCells/Experiment/Unsorted/unsorted/Voyager/SMAD4_colonies/Image/W1F<F,3>T<T,4>Z<Z,2>C<C,1>.tif';
texpout = '~/Desktop/ColonyW<Q,2>T<T,4>Z<Z,2>C<C,1>.jpg';

cc = 2; tt = 1; zz = 5;
clear imgs;
for i = 1:n
   fnoderaw =  ismip.rawCellIndex(ids(i));
   
   fn = tagExpressionToString(texp, ismip.rawCellIndexToRange(fnoderaw), 'C', cc, 'T', tt, 'Z', zz);
   img = imread(fn);
   img = imfrmtPermute(img, 'XY', 'Yx');
   
   fnout = tagExpressionToString(texpout, 'Q', i, 'C', cc, 'T', tt, 'Z', zz);
   
   imwriteBF(double(mat2gray(img)), fnout)
   
   imgs{i} = double(mat2gray(img));
end


%%

sh = alignImagesOnGrid(imgs(1:2)')

st = stitchImagesByHugin(imgs(1:2), sh);
%st = stitchImages(imgs(1:2), sh);

figure(6); clf
implot(st)


%%
clc
st = align2ImagesOnGridByHugin(imgs(1:2))

      %%
%       
%       tdat = imread(fnout);
%       size(tdat)
%       max(tdat(:))
%       min(tdat(:))
% 
%       figure(1); clf
%       implot(tdat)
         



