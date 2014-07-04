%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Track Sox 17 Cells %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all

initialize()
ilinitialize()

%%
for t = 17:20

   
%%
close all


%% Load Data

exp = Experiment('name', 'Example', 'date', datestr(now), 'description', 'Sox 17 Tracking, Test',...
                 'BaseDirectory', '/home/ckirst/Science/Projects/StemCells/Experiment/Data/Sox17/20140530T121749/', ...
                 'ImageDirectoryName', 'Image/',...
                 'ResultDirectoryName', 'Results/', ...
                 'ReadImageCommandFormat', 'imread(''<file>'')',...
                 'ReadImageFileFormat', 'F<field,3>/W1F<field,3>T<time,4>Z<zpos,2>C1.tif');

exp.info()



%% Background


% img = imread_stack(exp.ImageFile('field', 25, 'time', 20, 'zpos', '1*'));
% size(img)
% 
% {min(img(:)), max(img(:))}
% 
% imgmean = mean(img(:,:,3:end), 3);
% 
% %imgop = imopen(imgback, strel('disk', 24));
% 
% imgback = medianFilter(imgmean, 10);
% 
% 
% {min(imgback(:)), max(imgback(:))}
% 
% figure(10); clf
% implottiling({img(:,:,1), imgmean, imgback})
% 
% zmax = exp.findImageTagRange('field', 10, 'time', 1);
% zmax = max([zmax.zpos])
% 
% imgback = repmat(imgback, [1 1 zmax]);

%
imgback = exp.loadData('background.mat');


%%

tmax = exp.findImageTagRange('field', 25, 'zpos', 1);
tmax = max([tmax.time])


%% Initialize the Classifier

illoadclassifier(fullfile(exp.ResultDirectory, 'classifier3D_2.h5'))


%% Read Image Stack

filenames = exp.ImageFile('field', 25, 'time', t, 'zpos', '**')

img = imread_stack(filenames);

%{max(img(:)), min(img(:))}

%img(img>3000) = 3000;

img(1,1) = 4100;
img = (double(img) - 1800) / 4100 * 255 ;
size(img)

{max(img(:)), min(img(:))}

figure(5); clf; imcolormap('igray')

implottiling(img)

%% Classify
imgprob  = ilclassify(img);


%% Select Class

[~,imgseg] = max(imgprob,[], 4);

 
% figure(1); clf
% implot(imgseg)
% 
% figure(2); clf;
% implottiling(imgseg)

%%
%figure(42); clf;
%implottiling(imoverlaylabel(img, imgseg))


%%
%imglab = bwlabeln(imgprob(:,:,:,2) > 0.75);
imglab = bwlabeln(imgseg==2);
imglab = postProcessSegments(imglab, 'volume.min', 10);
imglab = bwlabeln(imglab > 0);

%figure(45); clf; 
%implottiling(imoverlaylabel(img, imglab))


%% Postprocess and 3D seed generation

imglab2 = imglab;

% fill holes
for s = 1:size(imglab,3)
   imglab2(:,:,s) = imfill(imglab(:,:,s), 'holes');
end

% morphological transfrom
for s = 1:size(imglab,3)
   imglab2(:,:,s) = imerode(imglab2(:,:,s), strel('disk', 1));
end

imgpost = postProcessSegments(imglab2, 'volume.min', 5);

imglab3 = bwlabeln(imgpost > 0);

figure(42); clf;
set(gcf, 'Name', ['Seeds'])
implottiling(imoverlaylabel(img, imglab3))



%% Watershed

% masking, might want to use differnt method, either another backgroudnclassifier or thresholding...
imgmask = ~(imgseg == 1); 

imgf = medianFilter(img, 3);
imgf = mat2gray(imgf); 
% imgf = imgf - * imgcls(:,:,:,3); imgf(imgf < 0) = 0;

%%
imgmin = imimposemin(max(imgf(:)) - imgf, imglab);
imgws = watershed(imgmin);
imgws = immask(imgws, imgmask);
%imgws = imlabelseparate(imgws);
imgws = postProcessSegments(imgws, 'volume.min', 10);

figure(43); clf;
set(gcf, 'Name', ['Label'])
implottiling(imoverlaylabel(img, imgws))


%% 3D Visulaization

% figure(44); clf;
% implot(double(imgws))


%%

figure(45); clf;
implot3dsurface(imgws)


%% Create Objects and Frame

objs = label2DataObjects(imgws, img);
frame = Frame('objects', objs, 't', t, 'File', filenames);



%% Initialize matfile to store results

matout = matfile(exp.ResultFile('segmentation.mat'), 'Writable', true);

matout.(['T', num2str(t,4)]) = frame;

matout

%% 

%beep on
%beep





%% time loop
end



