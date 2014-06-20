%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Track Sox 17 Cells %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

initialize()
ilinitialize()


%% Load Data

exp = Experiment('name', 'Example', 'date', datestr(now), 'description', 'Sox 17 Tracking, Test',...
                 'BaseDirectory', '/home/ckirst/Science/Projects/StemCells/Experiment/Sox17/20140530T121749/', ...
                 'ImageDirectoryName', 'Image/',...
                 'ResultDirectoryName', 'Results/', ...
                 'ReadImageCommandFormat', 'imread(''<file>'')',...
                 'ReadImageFileFormat', 'W1F010T<time,4>Z<zpos,2>C1.tif');

exp.info()



%% Load Test Image an plot

%exp.ReadImageCommand('time', 1, 'zpos', 6)

img = exp.readImage('time', 1, 'zpos', 1);
size(img)

figure(1); clf; imcolormap('g')
implot(img)


%% Initialize the Classifier

illoadclassifier(fullfile(exp.ResultDirectory, 'classifier.h5'))



%% Loop over Time
tagrange = exp.fileTagRange;

%for t = tagrange.time
for t = [1, 2]

%% Read Image Stack

t = 1;

img = imread_stack(exp.fileName('time', t, 'zpos', '**'));
size(img)

   
   
figure(42 + t);clf

%imgs = img - min(img(:));

%imgs = mat2gray(imgs);



imgs(imgs < 0.05) = 0;

implot(imgs)


%% Classify

% classify
siz = size(img);
imgprob = zeros([siz, 3]);
imgcls = zeros(siz);
imglab = zeros(siz);

siz(3) = 7;

for z = 1:siz(3)
   res  = ilclassify(img(:,:,z));
   imgprob(:,:,z,:) = reshape(res, [siz(1), siz(2), 1, 3]);

   %segment on max predictions
   [~,imgseg] = max(res,[], 3);
   imgcls(:,:,z) = imgseg;

   %create label
   imglab(:,:,z) = bwlabeln(imgseg == 2);
   
end

figure(2); clf;
implottiling(imgprob)


figure(3); clf;
implottiling(imgcls)

figure(42); clf;
implottiling(imoverlaylabel(img, imglab))

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


%% Postprocess and 3D seed generation

imglab2 = imglabn;

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
imgmask = ~(imgcls == 1); 

imgf = medianFilter(img, 3);
imgf = mat2gray(imgf); 
% imgf = imgf - * imgcls(:,:,:,3); imgf(imgf < 0) = 0;

%%
imgmin = imimposemin(max(imgf(:)) - imgf, imglabn);
imgws = watershed(imgmin);
imgws = immask(imgws, imgmask);
%imgws = imlabelseparate(imgws);
imgws = postProcessSegments(imgws, 'volume.min', 10);

figure(43); clf;
implottiling(imoverlaylabel(img, imgws))


%% 3D Visulaization

figure(44); clf;
implot(double(imgws))

%% Save Data To File




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

figure(45); clf;
implot3dsurface(imgws)



