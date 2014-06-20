%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Ilastik 3D %%%
%%%%%%%%%%%%%%%%%%%%%%%

initialize()
ilinitialize()


%% Initialize the Classifier

illoadclassifier('./Interface/Ilastik/Test/classifier2.h5')


%% Run Classifier

img = imread_stack('./Test/Images/hESCells_Stack/*.tif');

figure(1); clf;
implottiling(img)


%% Classify

% classify
siz = size(img);
imgprob = zeros([siz, 3]);
imgcls = zeros(siz);
imglab = zeros(siz);

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



%%

figure(45); clf;
implot3dsurface(imgws)



