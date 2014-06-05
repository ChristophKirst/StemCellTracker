%%%%%%%%%%%%%%%%%%%%
%%% Test Ilastik %%%
%%%%%%%%%%%%%%%%%%%%

initialize()
ilinitialize()


%% Initialize the Classifier

illoadclassifier('./Interface/Ilastik/Test/classifier.h5')


%% Run Classifier

img = imread('./Test/Images/hESCells_DAPI.tif');

% classify
imgcls = ilclassify(img);

%segment on max predictions
[~,imgseg] = max(imgcls,[], 3);

%create label
imglab = bwlabeln(imgseg == 2);

figure(42); clf;
implottiling({imgseg, imoverlaylabel(img, imglab)})


%% Postprocess and Water shed

imgmask = ~(imgseg == 1);

imgpost = postProcessSegments(imglab, 'volume.min', 10);

figure(42); clf;
implottiling({imgseg, imoverlaylabel(img, imgpost)})

imgf = mat2gray(img); imgf = imgf - 2 * imgcls(:,:,3);
imgf(imgf < 0 ) = 0;
imgmin = imimposemin(max(img(:)) - img, imgpost);
imgws = watershed(imgmin);
imgws = immask(imgws, imgmask);
%imgws = imlabelseparate(imgws);
imgws = postProcessSegments(imgws, 'volume.min', 10);

figure(43); clf;
implottiling({imgseg, imgws, imoverlaylabel(img, imgws)})






%% Test Basic Routines 

initialize
ilinitialize

%% test run
py('eval', 'ilc = IlastikClassifier()')
py('eval', 'res = ilc.run()')

res = pyget('res');
size(res)

res = squeeze(res);
size(res)

%% segmentation of the resulting data

back = res(:,:,1);
nucl = res(:,:,2);
bord = res(:,:,3);

size(res)
[~,seg] = max(cat(3, back, nucl, 1.5*bord),[], 3);
size(seg)


img = imread('./Test/Images/hESCells_DAPI.tif');


imglab = bwlabeln(seg == 2);

figure(42); clf;
implot(seg)

figure(43); clf;
implot(imoverlaylabel(img', imglab))











