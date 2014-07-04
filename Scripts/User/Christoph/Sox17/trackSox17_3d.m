%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Track Sox 17 Cells %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

initialize()


%% Load Data


exp = Experiment('name', 'Example', 'date', datestr(now), 'description', 'Sox 17 Tracking, Test',...
                 'BaseDirectory', '/home/ckirst/Science/Projects/StemCells/Experiment/Data/Sox17/20140530T121749/', ...
                 'ImageDirectoryName', 'Image/',...
                 'ResultDirectoryName', 'Results/', ...
                 'ReadImageCommandFormat', 'imread(''<file>'')',...
                 'ReadImageFileFormat', 'F<field,3>/W1F<field,3>T<time,4>Z<zpos,2>C1.tif');

exp.info()


%% Load Frames

matout = matfile(exp.ResultFile('segmentation.mat'), 'Writable', true);

timax = 11;
for ti = 1:timax;
   fprintf('loading frame %g / %g\n', ti, timax)
   frame(ti) = matout.(['T', num2str(ti+9,4)]);
end

frame


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

