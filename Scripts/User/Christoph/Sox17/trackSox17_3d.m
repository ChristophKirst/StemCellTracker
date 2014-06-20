%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Track Sox 17 Cells %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

initialize()
ilinitialize()


%% Load Data

exp = Experiment('name', 'Example', 'date', datestr(now), 'description', 'Sox 17 Tracking, Test',...
                 'BaseDirectory', '/home/ckirst/Science/Projects/StemCells/Experiment/Sox17/20140530T121749/', ...
                 'ImageDirectoryName', 'Image/F025',...
                 'ResultDirectoryName', 'Results/', ...
                 'ReadImageCommandFormat', 'imread(''<file>'')',...
                 'ReadImageFileFormat', 'W1F025T<time,4>Z<zpos,2>C1.tif');

exp.info()


%% Load Test Image and plot

%exp.ReadImageCommand('time', 1, 'zpos', 6)

img = exp.readImage('time', 4, 'zpos', 1);
size(img)

figure(1); clf; imcolormap('g')
implot(img)



%% Background


img = imread_stack(exp.fileName('time', 20, 'zpos', '1*'));
size(img)


{min(img(:)), max(img(:))}

imgmean = mean(img(:,:,3:end), 3);

%imgop = imopen(imgback, strel('disk', 24));

imgback = medianFilter(imgmean, 10);


{min(imgback(:)), max(imgback(:))}

figure(10); clf
implottiling({img(:,:,1), imgmean, imgback})

zmax = exp.fileTagRange('time', 1);
zmax = max([zmax.zpos])

imgback = repmat(imgback, [1 1 zmax]);


%% Plot Max-Projection for several times 

trange= (1:35) + 3;

img = zeros(512,512, length(trange));

ti = 1;
for t = trange
   t
   imgin = imread_stack(exp.fileName('time', t, 'zpos', '**')) - imgback;
   img(:,:,ti) = max(imgin, [], 3); 
   ti = ti + 1;
end
   
%img(img < 0) = 0;

figure(5); clf; imcolormap('igray')
implottiling(img)



%% Plot a Single Stack

t = 33;

img = imread_stack(exp.fileName('time', t, 'zpos', '**'));
size(img)

img = img-imgback;


figure(42 + t);clf
%imgs = img - min(img(:));
%imgs = mat2gray(imgs);
%imgs(imgs < 0.1) = 0;
imgs = img;
imgs(imgs < 0) = 0;
imgs = medianFilter(imgs,3);
implot(imgs)



%% Initialize the Classifier

illoadclassifier(fullfile(exp.ResultDirectory, 'classifier3D_2.h5'))


%% Initialize matfile to store results

matout = matfile(fullfile(exp.ResultDirectory, 'segmentation.mat'), 'Writable', true)


%% Loop over Time
tagrange = exp.fileTagRange;

%for t = tagrange.time
%for t = [1, 2]

t = 34;

%% Read Image Stack

filenames = exp.fileName('time', t, 'zpos', '**');

img = imread_stack(filenames);
img = mat2gray(img) * 255;
size(img)

%% Classify
imgprob  = ilclassify(img);


%%

figure(1); clf
[~,imgseg] = max(imgprob,[], 4);
implot(imgseg)

figure(2); clf;
implottiling(imgseg)

%%
figure(42); clf;
implottiling(imoverlaylabel(img, imgseg))


%%
%imglab = bwlabeln(imgprob(:,:,:,2) > 0.75);
imglab = bwlabeln(imgseg==2);
imglab = postProcessSegments(imglab, 'volume.min', 10);
imglab = bwlabeln(imglab > 0);

figure(45); clf; 
implottiling(imoverlaylabel(img, imglab))


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
implottiling(imoverlaylabel(img, imgws))


%% 3D Visulaization

figure(44); clf;
implot(double(imgws))


%%
figure(45); clf;
implot3dsurface(imgws)


%% Create Objects and Frame

objs = label2DataObjects(imgws, img);
frame = Frame('objects', objs, 't', t, 'File', filenames);


%%  Save Data To File

matout.(['T', num2str(t,4)]) = frame;


%% End Loop over times









%% Some Tests

test = imopen(imglab(:,:,5), strel('disk', 1));

figure(43); clf;
implottiling({imoverlaylabel(img(:,:,5), imglab(:,:,5)), imoverlaylabel(img(:,:,5), test)})


%%
figure(44); clf;
implottiling(imgprob(:,:,:,2))

%%

imglabn = bwlabeln(imgprob(:,:,:,2) > 0.75);
imglabn = postProcessSegments(imglabn, 'volume.min', 10);
imglabn = bwlabeln(imglabn > 0);

figure(45); clf; 
implottiling(imoverlaylabel(img, imglabn))





%% Time Loop
end


%%

% create a test file

hdf5write('/home/ckirst/Desktop/test.h5', '/results/test', [1 3 4 5])

%%

i = hdf5info('/home/ckirst/Desktop/test.h5')

%%

%%hdf5write('/home/ckirst/Desktop/test.h5', '/results/test2', [1 30], 'WriteMode', 'append')


hdf5write('/home/ckirst/Desktop/test.h5', '/doll/test7', [1 30], 'WriteMode', 'append')

%%



