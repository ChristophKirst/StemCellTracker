%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Segmentation Htt     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
initialize

verbose = true;



%% Load Data

exp = Experiment('name', 'Melissa', 'description', 'some random segmentation for Melissa',...
                 'BaseDirectory', '/home/ckirst/Science/Projects/StemCells/Experiment/Data/Other/Melissa/140805_Imaging/', ...
                 'ImageDirectoryName', '140513_RUES2_500um_<sample>.tif_Files/',...
                 'ResultDirectoryName', 'Results/', ...
                 'ReadImageCommandFormat', 'imread(''<file>'')',...
                 'ReadImageFileFormat', '140513_RUES2_500um_<sample>_p<tile,1>c<channel,1>.tif');

exp.info()

tags = exp.findImageTagRange


%%
%for sample = [8]
for sample = tags.sample

%% Plot Dapi 

imgs = exp.findImageFiles('channel', 0, 'sample', sample);
imgs = cellfun(@(x) imread(x), imgs, 'UniformOutput', false);
size(imgs)
imgs = reshape(imgs, 2,2);
imgs = [imgs(2), imgs(4); imgs(1), imgs(3)];


figure(1); clf; imcolormap('blue')
implottiling(imgs)


%% Stich

sh = alignImages(imgs, 'overlap.max', 150, 'overlap.min', 100);
figure(2); clf
plotAlignedImages(imgs, sh)

%%
img = stitchImagesByMean(imgs, sh);
figure(3); imcolormap('blue');
implot(img)

imwrite(mat2gray(img), [exp.ImageDirectory('sample', sample) 'stitch.tif'])


%% Load from file

img = imread([exp.ImageDirectory('sample', sample) 'stitch.tif']);
img = img(:,:,1);
%figure(3); clf;
%implot(img)



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
imgs = imgclsp(:,:,2)- 0.5 * imgclsp(:,:,3);
imgs(imgs < 0) = 0;
imgs = mat2gray(imgs);

%imgss = diskFilter(imgs, 15, 0, 1, -1);
imgss = sphereFilter(imgs, 15);

imgcl = imextendedmax(imgss, 0.05);


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











