%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Analze Sox 17 Cells %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

initialize()


%% Load Data


exp = Experiment('name', 'Example', 'description', 'Sox 17 Tracking, Test',...
                 'BaseDirectory', '/home/ckirst/Science/Projects/StemCells/Experiment/Data/Sox17/20140530T121749/', ...
                 'ImageDirectoryName', 'Image/',...
                 'ResultDirectoryName', 'Results/', ...
                 'ReadImageCommandFormat', 'imread(''<file>'')',...
                 'ReadImageFileFormat', 'F<field,3>/W1F<field,3>T<time,4>Z<zpos,2>C<channel,1>.tif');

exp.info()


exp.findImageTagRange


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Segment Sox 17 Cells %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Inspect Data

t = 1;
%%
img = imread_stack(exp.ImageFile('field', 44, 'time', t, 'zpos', '**', 'channel', 1));
img = impqlpermute(img, 'ypl', 'pql');
size(img)

figure(50); clf; imcolormap('g')
set(gcf, 'Name', ['3D'])
implottiling(mat2gray(img));


figure(51); clf; imcolormap('g')
set(gcf, 'Name', ['Max Proj'])
implot(max(img, [], 3));

t = t+1;

%%

%% Plot Max Propjections in Time
tt = 1:10;

img = zeros(512,512,length(tt));
ti = 1;
for t = tt
   fprintf('time" %d/%d\n', ti, length(tt))
   imgl = imread_stack(exp.ImageFile('field', 44, 'time', t, 'zpos', '5', 'channel', 1));
   imgl = impqlpermute(imgl, 'ypl', 'pql');
   
   img(:,:,ti) = max(imgl, [], 3);
   ti = ti + 1;
end


figure(51); clf; imcolormap('g')
set(gcf, 'Name', ['Time Evolution: Max Proj'])
implottiling(img);



%% Simple Segmentation of 2D Images


verbose = true;

t = 1;
ti = 1;

fields = {44; 48; 52};
gfppics = strfun(@(x) exp.ReadImageCommand('field', x, 'time', t, 'zpos', '5', 'channel', 1), fields);
dicpics = strfun(@(x) exp.ReadImageCommand('field', x, 'time', t, 'zpos', '8', 'channel', 3), fields);

%%
for p = 1:length(fields(:))
   fprintf('time %d/%d\n', ti, length(tt))
   
   % read tiles
   img{p} = eval(gfppics{p});
   %img = impqlpermute(img, 'ypl', 'pql');
   img{p} = max(img{p}, [], 3);

   imgdic{p} = eval(dicpics{p});
   %imgdic = impqlpermute(imgdic, 'yp', 'qp');

   % filter

   imgf{p} = medianFilter(img{p}, 3);
   
   if verbose
      figure(11 + p); clf;
      implot(imoverlayalpha(gray2rgb(mat2gray(imgdic{p})), imgray2color(mat2gray(imgf{p}), 'g')))
   end
end


%%






{min(img(:)), max(img(:))}

imgmean = mean(img(:,:,5:end), 3);

%imgop = imopen(imgback, strel('disk', 24));

imgback = medianFilter(imgmean, 10);


{min(imgback(:)), max(imgback(:))}

figure(10); clf
implottiling({img(:,:,1), imgmean, imgback})

zmax = exp.findImageTagRange('field', 10, 'time', 1);
zmax = max([zmax.zpos])

imgback = repmat(imgback, [1 1 zmax]);

%
imgback = exp.loadData('background.mat');


%%

tmax = exp.findImageTagRange('field', 25, 'zpos', 1);
tmax = max([tmax.time])


%% Initialize the Classifier

illoadclassifier(fullfile(exp.ResultDirectory, 'classifier3D_2.h5'))


%% Read Image Stack

filenames = exp.ImageFile('field', 25, 'time', t, 'zpos', '1')

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

































%%


% relabel joined segments


for t = 1:length(frame)
   t
   %objs = frame(t).objects();
   imglab = frame(t).labeledImage();
   imglab = bwlabeln(imglab > 0);
   imglab = postProcessSegments(imglab, 'volume.min', 50);
   
   objs = label2DataObjects(imglab, imglab > 0, 'time', t);
   framenew(i) = Frame('objects', objs, 't', t);
end
   
  


%% Tracking


param.optimize = true;
param.print.optimization = true;

[match, cost] = matchObjects(frame(1), frame(2), param);

figure(18)
clf
subplot(1,2,1)
plotMatchedObjects(match)
title('matches')
subplot(1,2,2)
plotMatchedCostMatrix(match, cost)
title('cost matrix')



%% match frames

%method
param.optimize = true;       % optimize match using an additional optimal coordinate transformation and second matching

% cost matrix
param.cost.creation  = [];   % costs for creating new objects ([])
param.cost.deletion = [];    % costs for deleting objects ([])

param.cutoff.cost = [];      % overall cost cutoff ([] = automatic)
param.cutoff.dist = [];      % cutoff for spatial distance ([] = automatic)
param.cutoff.time = -1;      % cutoff for temporal distances (-1 = ignore time differences)

param.weight.dist = 1.0;          % weight for distances in space (1.0)

param.weight.time.forward = 0.0;  % weight for positive time distances (0.0)
param.weight.time.backward = 0.0; % weight for negative time distances (0.0)
param.weight.volume = 100.0;        % weight for distances in volume (0.0)
param.weight.intensity = 100.0;     % weight for distances in intensity (0.0)
param.weight.type = Inf;          % weight for different types of objects (Inf)

% printing
param.print.load = 1;           % print info about loading the data files
param.print.save = 1;           % print info about saving the data files

param.print.match.frames  = 1;  % print info about matching the frames
param.print.match.objects = 1;  % print results on matching the objects

param.print.estimates = 1;      % print automatic determined estimates for cutoffs

% figures
param.figure.match = 0;          % generate figures for each match
param.figure.trajectories = 1;   % generate figure displaying the trajectories
param.figure.stats = 1;          % generate figure on statistics for the trajectories


[match, cost] = matchObjects(frame(1), frame(2), param);

figure(18)
clf
subplot(1,2,1)
plotMatchedObjects(match)
title('matches')
subplot(1,2,2)
plotMatchedCostMatrix(match, cost)
title('cost matrix')


%% full matching
[matches, traj] = trackObjects(frame, param);


%%

figure(42); clf;

plotMatchedTrajectories(frame, traj)


%% Plot object moving in time

imgtraj = traj(1).labeledImage();
for t = 2:10
   imgtraj = imgtraj + traj(t).labeledImage();
end
figure(42); clf
implot(imgtraj)













%% Plot Max-Projection for several times 

trange= (10:20);

img = zeros(512,512, length(trange));
imgback = exp.loadData('background.mat');

ti = 1;
for t = trange
   t
   imgin = imread_stack(exp.ImageFile('field', 25, 'time', t, 'zpos', '**')) - imgback;
   img(:,:,ti) = max(imgin, [], 3); 
   ti = ti + 1;
end
   
%img(img < 0) = 0;

figure(5); clf; imcolormap('igray')
implottiling(img)


%% Max projection of Trajectories

%imgtrajmax = gray2rgb(iminvert(mat2gray(img)));
imgtrajmax = iminvert(mat2gray(img));

%trajmax = 15;
trajmax = length(traj)

cols = imcolormap('jet',trajmax);
imgtrajobj = zeros(512,512, length(frame));


for tid = 1:trajmax
   objs = traj(tid).objects;
   fids = traj(tid).frameids;
   
   for i = 1:length(objs)
      imgobj = objs(i).labeledImage();
      imgobj = tid * (imgobj > 0);
      imgobj = max(imgobj, [], 3);
      
      t = fids(i);
      %imgtrajmax(:, :, fids(i), :) = squeeze(imgtrajmax(:, :, fids(i), :)) + imgray2color(imgobj, cols(tid,:));
      imgtrajobj(:,:,t) = max(imgtrajobj(:,:,t), imgobj);
   end
end

figure(42); clf
implottiling(imoverlaylabel(imgtrajmax, imgtrajobj))


%% Max-projection of Segmentation

%trajmax = 15;
tmax = length(frame);

imgobj = zeros(512,512, tmax);
%imgtrajmax = gray2rgb(iminvert(mat2gray(img)));
imgmax = iminvert(mat2gray(img));

for tid = 1:tmax

   imgframe = frame(i).labeledImage();
   imgframe = max(imgframe, [], 3);
      
   t = fids(i);
   %imgtrajmax(:, :, fids(i), :) = squeeze(imgtrajmax(:, :, fids(i), :)) + imgray2color(imgobj, cols(tid,:));
   imgmax(:,:,t) = max(imgmax(:,:,t), imgframe);
end

figure(45); clf
implottiling(imoverlaylabel(imgtrajmax, imgtrajobj))


%%

figure(43); clf;
imgseg = frame(1).labeledImage();
implottiling(imoverlaylabel(iminvert(mat2gray(img)), imgseg))



%% Generate Time Series

ts = TimeSeries('frames', frame, 'trajectories', traj);


%% Save 
exp.result = ts;
exp.saveExperiment('tracking.mat');



%% plot some statistics

figure(42); clf; 
subplot(1,2,1)
hist(double([frame(1).objects.intensity]), 15)
subplot(1,2,2)
hist(double([frame(1).objects.volume]), 15)

