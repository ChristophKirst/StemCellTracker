%%%%%%%%%%%%%%%%%%%%
%%% Test Ilastik %%%
%%%%%%%%%%%%%%%%%%%%

initialize()
ilinitialize()


%% Initialize the Classifier

ilc = illoadclassifier('./Interface/Ilastik/Test/classifier.h5')


%% Run Classifier

img = imread('./Test/Images/hESCells_DAPI.tif');

figure(1); clf
implot(img);

%%

% classify
imgcls = ilclassify(ilc, img);


%%
%segment on max predictions
[~,imgseg] = max(imgcls,[], 3);

%create label
imglab = bwlabeln(imgseg == 2);

figure(42); clf;
implottiling({imgseg; imoverlaylabel(img, imglab)})


%% Postprocess and Watershed

imgmask = ~(imgseg == 1);

imgpost = postProcessSegments(imglab, 'volume.min', 10);

figure(42); clf;
implottiling({imgseg; imoverlaylabel(img, imgpost)})

imgf = mat2gray(img); imgf = imgf - 3 * imgcls(:,:,3);
imgf(imgf < 0 ) = 0;
imgmin = imimposemin(max(img(:)) - img, imgpost);
imgws = watershed(imgmin);
imgws = immask(imgws, imgmask);
%imgws = imlabelseparate(imgws);
imgws = postProcessSegments(imgws, 'volume.min', 10);

figure(43); clf;
implottiling({imgseg; imgws; imoverlaylabel(img, imgws)})



