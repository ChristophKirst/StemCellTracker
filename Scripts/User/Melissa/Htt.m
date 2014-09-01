%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Segmentation Htt     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
initialize

ilinitialize
bfinitialize

verbose = true;

addpath('./Scripts/User/Melissa')


%% Load Data

exp = Experiment('name', 'Melissa', 'description', 'some random segmentation for Melissa',...
                 'BaseDirectory', '/home/ckirst/Science/Projects/StemCells/Experiment/Data/Other/Melissa/', ...
                 'ImageDirectoryName', 'Images/',...
                 'ResultDirectoryName', 'Results/', ...
                 'ReadImageCommandFormat', 'imread_bf(''<file>'', ''series'', <tile>, ''channel'', <channel>)',...
                 'ReadImageFileFormat', '140513_RUES2<type,4,s>_500um_<sample>.zvi');

exp.info()

tags = exp.findImageTagRange


%%
%for sample = [8]
sample = 8
type = 'Ctrl';

%for sample = tags.sample

%% Plot Dapi Images

for t = 1:4
   imgs{t} = exp.readImage('type', type, 'sample', sample, 'tile', t, 'channel', 1);
end
imgs = [imgs(3), imgs(4); imgs(1), imgs(2)];


figure(1); clf; imcolormap('blue')
implottiling(imgs)

%% Stitch

sh = httalign(imgs, false);


% correct illumination etc

for t = 1:4
   imgsc{t} = imopen(imgs{t}, strel('disk', 20));
   imgsc{t}  = imgs{t} - imgsc{t};
end
imgsc = reshape(imgsc, 2,2);

img = stitchImagesByMean(imgsc, sh);

figure(3); imcolormap('blue');
implot(img)

imwrite(mat2gray(img), [exp.ImageDirectory('sample', sample) 'stitch.tif'])


%% Load from file

img = imread([exp.ImageDirectory('sample', sample) 'stitch.tif']);
img = img(:,:,1);
%figure(3); clf;
%implot(img)



%% Segmenting without classification

imgs = medianFilter(img, 3);
imgs = mat2gray(imgs);

%%
%imgss = diskFilter(imgs, 15, 0, 1, -1);
imgss = diskFilter(imgs, 10, 0, 1, -1);


%%
imgcl = imextendedmax(imgss, 0.05);

imgcl = imdilate(imgcl, strel('disk', 2));
imgcl = bwlabeln(imgcl);


if verbose 
   figure(31); clf
   implottiling({img, imgs; imgss, imoverlay(mat2gray(img), imgcl)})
end







%% Initialize the Classifier

illoadclassifier([exp.BaseDirectory 'classifier.h5'])


%% Run Classifier

% classify
imgcls = ilclassify(img);

%%
imgclsp = impqlpermute(imgcls, 'lpq', 'pql');
size(imgclsp)

%segment on max predictions
[~,imgseg] = max(imgclsp,[], 3);

%create label
imglab = bwlabeln(imgseg == 2);

figure(42); clf;
implottiling({imgclsp(:,:,1), imgclsp(:,:,2), imgclsp(:,:,3)})


%% Segmenta Cells




%% seeds
imgs = imgclsp(:,:,2)- 0.5 * imgclsp(:,:,1);
imgs(imgs < 0) = 0;
imgs = mat2gray(imgs);

%imgss = diskFilter(imgs, 15, 0, 1, -1);
imgss = sphereFilter(imgs, 15);

imgcl = imextendedmax(imgss, 0.2);


imgcl = imdilate(imgcl, strel('disk', 4));
imgcl = bwlabeln(imgcl);


if verbose 
   figure(31); clf
   implottiling({img, imgs; imgss, imoverlay(mat2gray(img), imgcl)})
end
















%%

%% seeds
imgs = imgclsp(:,:,2)- 0.5 * imgclsp(:,:,1);
imgs(imgs < 0) = 0;
imgs = mat2gray(imgs);

%imgss = diskFilter(imgs, 15, 0, 1, -1);
imgss = sphereFilter(imgs, 15);

imgcl = imextendedmax(imgss, 0.2);


imgcl = imdilate(imgcl, strel('disk', 4));
imgcl = bwlabeln(imgcl);


if verbose 
   figure(31); clf
   implottiling({img, imgs; imgss, imoverlay(mat2gray(img), imgcl)})
end


%% watershed

imgmask = ~(imgseg == 1);

imgf = mat2gray(medianFilter(img, 2)) - imgclsp(:,:,3);
imgf(imgf < 0 ) = 0;
imgf = mat2gray(imgf);

imgmin = imimposemin(max(imgf(:)) - imgf, imgcl);
imgws = watershed(imgmin);
imgws = immask(imgws, imgmask);
%imgws = imlabelseparate(imgws);

%%
imgws2 = imopen(imgws, strel('disk', 2));
imgws2 = bwlabeln(imgws2 > 0);
imgws2 = postProcessSegments(imgws2, 'volume.min', 250);


figure(43); clf;
implottiling({imgseg, imgws2, imoverlaylabel(img, imgws2)})


%%
stats = imstatistics(imgws2, 'PixelIdxList');

%% Read other images and stitch


for c = 1:3
   imgsr = exp.findImageFiles('channel', c, 'sample', sample);
   imgsr = cellfun(@(x) imread(x), imgsr, 'UniformOutput', false);
   imgsr = [imgsr(2), imgsr(4); imgsr(1), imgsr(3)];
   
   imgC{c} = stitchImages(imgsr, sh, 'method', 'Mean');
end


if verbose
   figure(10); clf;
   implottiling(imgC)
end

%% quantify



% median intensities 

for c = 1:3
   statsC{c} = imstatistics(imgws2, stats, 'MedianIntensity', imgC{c});
end


% color by median intensity of each

for c = 1:3
   statsI = statsC{c};
   imgR = zeros(size(img));
   for s = 1:length(statsI)
      imgR(statsI(s).PixelIdxList) = statsI(s).MedianIntensity;
   end

   imgRR{c} = mat2gray(imgR);
end



%%

cols = {'r', 'g', 'b'};
for c = 1:3
   imgRR{c} = imgray2color(imgRR{c}, cols{c});
   h = figure(42); clf; 
   implottiling(imgRR);
end
   

%% print

if ~isdir(exp.ResultDirectory)
   mkdir(exp.ResultDirectory);
end

print(h,'-dpdf', [exp.ResultDirectory, 'result_' num2str0(sample, 2) '.pdf']);


%% save

save([exp.ResultDirectory  'result_' num2str0(sample, 2) '.mat'], 'stats', 'statsC');


end











