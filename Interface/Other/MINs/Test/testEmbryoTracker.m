function testEmbryoTracker()

%% General Tests with Artificial Data

%
% test the Tracking 
%

path(path, './Classes')
path(path, './Interface')
path(path, './Tracker')


%% Test classes 

o1 = Object(1, 1, [0 1]', 100, [], [4 8 6]');
o2 = Object(2, 1, [0 2]', 90, [], [1 5 6]');
o3 = Object(3, 1, [0 3]', 10, [], [1 0 6]');

o4 = Object(4, 2, [1 1]', 90, [], [1 5 6]'); 
o5 = Object(5, 2, [1 3]', 10, [], [1 0 6]'); 


o6 = Object(6, 3, [2 1]', 100, [], [1 5 6]'); 
o7 = Object(7, 3, [2 2]', 90, [],  [1 0 6]'); 

f1 = Frame('test image 1', [o1, o2, o3]);
f2 = Frame('test image 2', [o4, o5]);
f3 = Frame('test image 3', [o6, o7]);


%% trajectories

t1 = Trajectory([o1, o4, o6]);
t2 = Trajectory([o3, o5]);
t3 = Trajectory([o2, o7]);

traj = [t1, t2, t3];

[~, ~, objs] = traj.timeSlice(3);

[objs.id]
f3.id

%% slicing

[trajpos, tpos, objs, preobjs, postobjs] = traj.timeSlicePrePost(2)

[preobjs.id]
[objs.id]
[postobjs.id]


%% Test tracking


data0 = f1;
data1 = f2;

[match, cost] = matchObjects(data0, data1);

figure(1)
clf
subplot(1,2,1)
plotMatchedObjects(match)
title('matches')
subplot(1,2,2)
plotMatchedCostMatrix(match, cost)
title('cost matrix')


%% Test optimal transformation

% some points
X0 = [1 2 3 4 8 3 3 3 11 -5; 3 4 5 6 1 3 9 10 -3 -2; 2 1 0 5 15 9 -1 -5 -8 10];

% some transformation
theta = 1.2;
Rtest = [1, 0, 0; 0, cos(theta), -sin(theta); 0, sin(theta), cos(theta)]
Ttest = [0.4; 1; 0]
Ctest = 1.0

n = size(X0,2);
X1 = zeros(3,n);
for i = 1:n
   X1(:,i) = (Ctest * Rtest * X0(:,i)) + Ttest;
end

[R, T, C] = optimalTransformation(X0, X1)

Xt = zeros(3,n);
for i = 1:size(X0,2)
   Xt(:,i) = (C * R * X0(:,i)) + T;
end

figure(2)
clf
hold on
h0 = scatter3(X0(1, :), X0(2, :), X0(3, :));
h1 = scatter3(X1(1, :), X1(2, :), X1(3, :), 100, [0 1 0]);
ht = scatter3(Xt(1, :), Xt(2, :), Xt(3, :), 50, [1 0 0]);
xlabel('x'); ylabel('y'); zlabel('z');
legend([h0, h1, ht], 'r0', 'r1', 'rotated r0');
title('Test optimalTransformation.m');
hold off








%% Tests with real Data

%
% Test code using typical sample Data
%


%% match embryo data - single step between two frames

data0 = loadEmbryoDataFile('./Test/Data/Timelapse_11152013_channel=0001_frame=0001_statistics.csv')
data1 = loadEmbryoDataFile('./Test/Data/Timelapse_11152013_channel=0001_frame=0002_statistics.csv')


[match, cost] = matchObjects(data0, data1);


%% plot result

figure(2)
clf
subplot(1,2,1)
plotMatchedObjects(match)
title('matches')
subplot(1,2,2)
plotMatchedCostMatrix(match, cost)
title('cost matrix')




%% extract data for single embryo for optimized tracking

data0 = loadEmbryoDataFile('./Test/Data/Timelapse_11152013_channel=0001_frame=0001_statistics.csv');
data1 = loadEmbryoDataFile('./Test/Data/Timelapse_11152013_channel=0001_frame=0002_statistics.csv');

dat = [data0.r];
indx = dat(2,:) > 340;
data0.objects = data0.objects(indx)

dat = [data1.r];
indx = dat(2,:) > 340;
data1.objects = data1.objects(indx)


%% match single embryo

[match, cost] = matchObjects(data0, data1);

%stats = matchedStatistics(data0, data1, match, cost);

figure(4)
clf
subplot(1,2,1)
plotMatchedObjects(match)
title('matches')
subplot(1,2,2)
plotMatchedCostMatrix(match, cost)
title('cost matrix')



%% find optimal rotation

[X0, X1] = match.toCoordinates();

disp 'optimal transformation:'
[R, T, C] = optimalTransformation(X0,X1)

n = size(X0,2);
Xt = zeros(3,n);
for i = 1:n
   Xt(:,i) = (C * R * X0(:,i)) + T;
end


figure(5)
clf
hold on
grid on
h0 = scatter3(X0(1, :), X0(2, :), X0(3, :));
h1 = scatter3(X1(1, :), X1(2, :), X1(3, :), 100, [0 1 0]);
ht = scatter3(Xt(1, :), Xt(2, :), Xt(3, :), 50, [1 0 0]);
xlabel('x'); ylabel('y'); zlabel('z');
legend([h0, h1, ht], 'r0', 'r1', 'rotated r0');
title('Test optimalTransformation.m');
hold off

%% second match after rotation

data0t = data0.copy; % deep copy
data0t  = data0t.transformCoordinates(R, T, C);

[matcht, costt] = matchObjects(data0t, data1);
matcht.objects0 = data0.objects;

figure(6)
clf
subplot(1,2,1)
plotMatchedObjects(matcht)
title('matches')
subplot(1,2,2)
plotMatchedCostMatrix(matcht, costt)
title('cost matrix')

%% compare results

stats = match.statistics(cost);
statst = matcht.statistics(costt);


figure(7)
clf
subplot(2,2,1)
hist(stats.dist.values);
title(sprintf('initial match spatial distances:\nn: %d mean: %g std: %g', length(stats.dist.values), stats.dist.mean, stats.dist.std));

subplot(2,2,2)
hist(statst.dist.values); 
title(sprintf('rotated match spatial distances:\nn: %d mean: %g std: %g', length(statst.dist.values), statst.dist.mean, statst.dist.std));

subplot(2,2,3)
hist(stats.cost.values);
title(sprintf('initial match costs:\nn: %d mean: %g std: %g', length(stats.cost.values), stats.cost.mean, stats.cost.std));

subplot(2,2,4)
hist(statst.dist.values); 
title(sprintf('rotated match costs:\nn: %d mean: %g std: %g', length(statst.dist.values), statst.cost.mean, statst.cost.std));





%% otimized Matching (= match, rotation, match)


param.optimize = true;
param.print.optimization = true;

[match, cost] = matchObjects(data0, data1, param);

figure(8)
clf
subplot(1,2,1)
plotMatchedObjects(match)
title('matches')
subplot(1,2,2)
plotMatchedCostMatrix(match, cost)
title('cost matrix')




%% track full list of time frames of single embryo

param.load.min = 2;       % at least one object in frame
param.load.change = 0.2;  % at most 20% change in number of objects

param.print.load = true;
param.print.match.objects = true;
param.print.match.optimization = false;

param.optimize = true;


frames = loadEmbryoData('./Test/Data', param);

for t = 1:length(frames)
   
   dat = frames(t).r;
   indx = dat(2,:) > 340;
   frames(t).objects = frames(t).objects(indx);
   
end

[matches, costs] = matchFrames(frames, param);


%% plot the result 

for t =1:length(matches)
   figure
   clf
   subplot(1,2,1)
   plotMatchedObjects(matches(t))
   title('matches')
   subplot(1,2,2)
   plotMatchedCostMatrix(matches(t), costs{t})
   title('cost matrix')
end


%% determine trajectories

trajs = findTrajectories(matches);

figure
clf
plotMatchedTrajectories(frames, trajs)

%% some statistics

stats = statisticsTrajectory(trajs);

figure
subplot(1,2,1)
hist(stats.length.values)
title(sprintf('trajectory time lengths:\nmean:%g std:%g', stats.length.mean, stats.length.std))
xlabel('time');

subplot(1,2,2)
hist(stats.dist.values)
title(sprintf('trajectory spatial lengths:\nmean:%g std:%g', stats.dist.mean, stats.dist.std))
xlabel('distance')


%% saving data

saveEmbryoData('./Test/Out', frames, trajs)



%% run full Tracker 

runTracker('./Test/Data', './Test/Out')


%% run full Tacker with standard parameter and a filter 

path(path, './Test')

param = setParameter();

param.filter =  @testFilter;

runTracker('./Test/Data', './Test/Out', param)


%% run full Tacker with test parameter 
path(path, './Test')

param = setParameterTest();

param.filter = @testFilter;

[frames, matches, trajs] = runTracker('./Test/Data.all', './Test/Out', param);




end