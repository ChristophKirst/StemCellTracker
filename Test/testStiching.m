%%%%%%%%%%%%%%%%%%%%%
%%% Test Stiching %%%
%%%%%%%%%%%%%%%%%%%%%

initialize

%% 2D 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Test Image

img = imread('./Test/Images/hESCells_DAPI.tif');
img = mat2gray(img);

size(img)

img1 = img(1:303, 1:400);
img2 = img(260:end, 58:end);

img1 = img(1:303, 1:400);
img2 = img(260:end, 18:end);


figure(1)
clf
implot(img)

figure(2)
clf
imsubplot(1,2,1)
implot(img1);
imsubplot(1,2,2);
implot(img2);


%% 

%extract the border regoinsÖ

bb = 50;
img1b = img1(end-bb:end, :);
img2b = img2(1:bb, :);

figure(11)
implottiling({img1b, img2b})


[optimizer, metric] = imregconfig('monomodal')
%[optimizer, metric] = imregconfig('multimodal')


%%% for gradient descfent
%optimizer.InitialRadius = 0.009;
%optimizer.Epsilon = 1.5e-4;
%optimizer.GrowthFactor = 1.01;
%optimizer.GradientMagnitudeTolerance = 100;

optimizer.MaximumIterations = 3000;
optimizer.RelaxationFactor = 0.95;
optimizer.MaximumStepLength = 0.5;
optimizer.MinimumStepLength = 0.05;

%%% for evolutionary


%optimizer.InitialRadius = 15;
%optimizer.MaximumIterations = 3000;
%metric.NumberOfHistogramBins = 256;
%metric.NumberOfSpatialSamples = 1000;



%trsf0 = affine2d([1 0 0; 0 1 0; 60 size(img1, 1)-50 1]);
trsf0 = affine2d([1 0 0 ; 0 1 0; 0 0 1]);
%movingRegistered = imregister(img2b, img1b, 'translation', optimizer, metric, 'InitialTransform', trsf0,'DisplayOptimization',false);


tform = imregtform(img2b, img1b, 'translation', optimizer, metric, 'InitialTransform', trsf0,'DisplayOptimization',false)

%add shifts for full images

%T = tform.T;



     %moving = varargin{1};
     %fixed = varargin{2};
 %    dim =2;
%      if (dim == 2)
%         Rmoving = imref2d(size(img2));
%         Rfixed = imref2d(size(img1));
%      else
%          Rmoving = imref3d(size(img2));
%          Rfixed = imref3d(size(img1));
%      end
% 
%      
%      
%      
% 
% [movingRegistered,Rreg] = imwarp(img2,Rmoving,tform,'OutputView',Rfixed);



sh = tform.T;
sh = sh(3,[2 1])

%pixel shifts
sh = round(sh);

%absolute shift due to img1 size
sh = sh + size(img1) - [bb 0]

img2sh = zeros(size(img2) + sh);
img2sh(sh(1):end, sh(2):end) = img2;


% move the image


figure(3)
clf
imshowpair(img1b', img2sh','Scaling','joint');
set(gca, 'YDir', 'normal')

%%
size(movingRegistered)
size(img2)

clf
imsubplot(1,2,1)
implot(img2)
imsubplot(1,2,2)
implot(movingRegistered)



