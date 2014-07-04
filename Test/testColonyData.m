%% Gernerate Some TestData for Colony Class integration

initialize

%% Image sources

exp = Experiment('name', 'Example', 'description', 'Sox 17 Tracking, Test',...
                 'BaseDirectory', '/home/ckirst/Science/Projects/StemCells/Experiment/Data/Sox17/20140530T121749/', ...
                 'ImageDirectoryName', 'Image/',...
                 'ResultDirectoryName', 'ResultsColony/', ...
                 'ReadImageCommandFormat', 'imread(''<file>'')',...
                 'ReadImageFileFormat', 'F<field,3>/W1F<field,3>T<time,4>Z<zpos,2>C<channel,1>.tif');

exp.info()

exp.findImageTagRange



%% Determine Shifts once
%load grid of images -> in one of the next versions ImageSources will be integrated in the Experiment class directly...

t = 1;
fields = [44, 48, 52];
imgsdic = {}; imgsgfp = {};
for fi = 1:length(fields);
   f = fields(fi);
   imgsgfp{fi} = ImageSourceFile('filename', exp.ImageFile('field', f, 'time', t, 'zpos', 5, 'channel', 1), 'format', 'pq');
   imgsdic{fi} = ImageSourceFile('filename', exp.ImageFile('field', f, 'time', t, 'zpos', 8, 'channel', 3), 'format', 'pq');   
end

stichparam = parseParameter('method', 'Mean');
alignparam = parseParameter('method', 'SequentialShifts', 'overlap.max', 120, 'overlap.min', 90');

is = ImageSourceTiled('images', imgsgfp, 'alignparam', alignparam, 'stichparam', stichparam)

%%
is.align;
is

img = is.getData;
   
figure(1); clf; imcolormap('gray')
implot(img)

shifts = is.shifts;



%% Loop over Time 

verbose = true;
fields = [44, 48, 52];

tt = 1:10;
%tt = [1,2];

ti = 1;
for t = tt;

%% Simple Segmentation of 2D Images


fprintf('time %d/%d\n', ti, length(tt))

%% align 
imgsdic = {}; imgsgfp = {};
for fi = 1:length(fields);
   f = fields(fi);
   imgsgfp{fi} = ImageSourceFile('filename', exp.ImageFile('field', f, 'time', t, 'zpos', 5, 'channel', 1), 'format', 'pq');
   imgsdic{fi} = ImageSourceFile('filename', exp.ImageFile('field', f, 'time', t, 'zpos', 8, 'channel', 3), 'format', 'pq');   
end

isgfp = ImageSourceTiled('images', imgsgfp, 'alignparam', alignparam, 'stichparam', stichparam, 'shifts', shifts);
isgfp.initialize();

isdic = ImageSourceTiled('images', imgsdic, 'alignparam', alignparam, 'stichparam', stichparam, 'shifts', shifts);
isdic.initialize();

%% Filter

imgdic = isdic.getData();
imggfp = isgfp.getData();
imgf = medianFilter(imggfp, 3);

if verbose
   figure(11); clf;
   implot(imoverlayalpha(gray2rgb(mat2gray(imgdic)), imgray2color(mat2gray(imgf), 'g')))
end


%% Simple 2D segmentation

imgl = logFilter(imgf, [15, 15]);
imgl = iminvert(mat2gray(imgl));

imgmax = imextendedmax(imgl, 0.05);
imglab = bwlabeln(imgmax);
%imglab = postProcessSegments(imglab, 'volume.min', 10);

if verbose
   figure(12); clf; imcolormap('jet')
   implottiling({imgray2color(mat2gray(imgf), 'g'); imgray2color(imgl, 'g'); imoverlaylabel(imgf, imglab)})
end


%% Create Objects and Frame

objs = label2DataObjects(imglab, img);
frame(ti) = Frame('objects', objs, 't', t, 'source', isgfp);

if verbose
   imgt = frame(ti).getImage;
   figure(5); clf;
   implot(imgt)
end


ti = ti + 1;
end




%% Tracking

%method
param.optimize = false;       % optimize match using an additional optimal coordinate transformation and second matching

% cost matrix
param.cost.creation  = [];   % costs for creating new objects ([])
param.cost.deletion = [];    % costs for deleting objects ([])

param.cutoff.cost = [];      % overall cost cutoff ([] = automatic)
param.cutoff.dist = 20;      % cutoff for spatial distance ([] = automatic)
param.cutoff.time = -1;      % cutoff for temporal distances (-1 = ignore time differences)

param.weight.dist = 1.0;          % weight for distances in space (1.0)

param.weight.time.forward = 0.0;  % weight for positive time distances (0.0)
param.weight.time.backward = 0.0; % weight for negative time distances (0.0)
param.weight.volume = 0.0;        % weight for distances in volume (0.0)
param.weight.intensity = 10.0;     % weight for distances in intensity (0.0)
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




%% Generate Time Series

ts = TimeSeries('frames', frame, 'trajectories', traj);


%% Save 
exp.result = ts;
exp.saveExperiment('tracking.mat');



