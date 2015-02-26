%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Mouse Nanog visualize MINs results %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc
initialize
bfinitialize
ijinitialize

verbose = true;

addpath('./Scripts/User/Christoph/mESCNanog');
addpath('./Interface/User/EmbryoData');


initializeParallelProcessing(12) % number of processors


%plotslices = 1:2:40;
tiling = [8,5];



%% Data

datadir = '/data/Science/Projects/StemCells/Experiment/Mouse/Nanog/12Aug14FGFonNanogH2BGFP-700/';
%dexp = 'T<F,3>/T<F,3>.tif';
%fulldexp = fullfile(datadir, dexp);
fulldexp = fullfile(datadir, '12Aug14FGFonNanogH2BGFP-700_movie.lsm');
dns = tagExpressionToFiles(fulldexp);

isd = ImageSourceBF(fulldexp);
isd.printInfo


%% Plot Data
f = 1;

img = isd.data('T', f, 'C', 1);

figure(1); clf
implottiling(img, 'tiling', tiling)

%% Create Filtered Data

resultdir = '/data/Science/Projects/StemCells/Experiment/Mouse/Nanog/12Aug14FGFonNanogH2BGFP-700_MINS_1/';

filterdir = fullfile(resultdir, 'filtered');
mkdir(filterdir)

outfile = fullfile(filterdir, 'T<T,3>.tif');

tmax = 10;


%%
parfor t = 1:tmax

   % load data
   img = isd.data('T', t, 'C', 1); 
   
   imax = 65;
   imax = 255;
   img = imclip(img, 0, imax);
   img = double(img) / imax;
   %img = imclose(img, strel('disk', 1));
   %img = filterMedian(img, [4,4,2]);
   
   figure(1); clf; 
   hist(img(:), 512);
   
   figure(2); clf;
   implottiling(img, 'tiling', tiling)
   
   
   %%
   %img = imclip(img, 0, 1);
   imgf = filterBM(img, 'sigma', 5, 'depth', 3);  % sigma is noise level out of 0..255, %img should be 0..1 though
   imgf = double(imgf);
   %imgf = img;
   
  
   figure(3); clf;
   imgp = zeros([size(img,1), size(img,2), 2* size(img,3)]);
%    k = 1;
%    for i = 1:size(img,3)
%       imgp(:,:,k) = img(:,:,i);
%       imgp(:,:,k+1) = imgf(:,:,i);
%       k = k + 2;
%    end
%    implottiling(imgp, 'tiling', [16,5]);
      

   %% gradient
   
%    imgg = imgradient3D(imgf);
%    figure(5); clf;
%    implottiling(imgg, 'tiling', tiling);

   %%
   %% 
   %imgfm = filterFunction(imgf, [3,3,3], 'median');
   
   %figure(5); clf;
   %implottiling(imgfm, 'tiling', tiling);

   %%
   imwriteStack(imgf, tagExpressionToString(outfile, 'T', t));
end

%% Illastik Segmentation

ildir = fullfile(resultdir, 'ilastik');
mkdir(ildir)

ilfile = fullfile(ildir, 'T<T,3>Z<Z,2>.tif');

t = 1;

imgf = imreadStack(tagExpressionToString(outfile, 'T', t));

zz = 1;
imgsh = zeros(size(imgf,1), size(imgf,2), 2* size(imgf,3));
for z = 1:size(imgf,3)
   for k = 1:2
      imgsh(:,:,zz) = imgf(:,:,z);
      imwriteBF(double(imgf(:,:,z))/255, tagExpressionToString(ilfile, 'T', t, 'Z', zz));
      zz = zz + 1;
   end
end

%% Illastik Segmentation

ilinitialize

cl = illoadclassifier(fullfile(resultdir, 'classifier.h5'))

%%
imgsg = ilclassify(cl, imgsh);


%%
figure(7);
implottiling(imgsg, 'tiling', tiling)

%%
figure(8)
implottiling(imgsg(:,:,:,3) > 0.5, 'tiling', tiling)

%%

%[~, imgseg] = max(imgsg(:, :, 1:2:end, :), [], 4);

imgseg = imgsg(:,:,1:2:end,2) > 0.5;
%imgseg = imgseg - (imgsg(:,:,1:2:end,3) > 0.5);
%imgseg(imgseg < 0) = 0;

imgseg = imerode(imgseg, strel('disk', 4));
imgseg = bwlabeln(imgseg);


figure(7);
colormap jet
implottiling(imoverlaylabel(imgf(:,:,:), imgseg(:,:,:)), 'tiling', tiling)


%%

imgm = filterMedian(img, [5,5,3]);

%%

% some pseudo psf:

psf = zeros(11,11,7);
h = (size(psf,3)+1)/2;
m = size();
psf(:,:,h);
for z = 1:h
   psf(:,:,h+z)


%%
psf = ones(7,7, 5);
[imgd, psfd] = deconvblind(imgf, psf);

%%
figure(18); clf;
implottiling(imgd, 'tiling', tiling)

figure(19); clf;
implottiling(psfd)


%% Segmentation 

for t = 1:tmax

   imgf = imreadStack(tagExpressionToString(outfile, 'T', t));
   imgf = imfrmtReformat(imgf, 'XYZ', 'yXZ');
   
%    imgs = isd.data('T', t, 'C', 1); 
%    imgs = double(imgs) / 255;
%    imgs = imclip(imgs, 0, 1);
%    
   
   size(imgf)
   
   imglog = log(imgf + eps + 1);
   %imglog(imglog < 0) = 0;
  
   if verbose
      figure(15); clf
      set(gcf, 'Name', ['Time: ', num2str(t), ' Mask']);
      implottiling(imgf, 'tiling', tiling);
      %implottiling(imglog, 'tiling', tiling);
   end

   
   %% Masking
      
   imgmask = 4*mat2gray(imgf) > 0.075;

   if verbose
      figure(15); clf
      set(gcf, 'Name', ['Time: ', num2str(t), ' Mask']);
      implottiling(imoverlay(12*mat2gray(imgf), imgmask, 'r', false), 'tiling', tiling);
   end

   
   %%
   imgmaske = imerode(imgmask, strel('disk', 4));
   if verbose
      figure(16); clf
      set(gcf, 'Name', ['Time: ', num2str(t), ' Mask']);
      implottiling(imoverlay(mat2gray(imgf), imgmaske, 'r', false), 'tiling', tiling);
   end
   
   
   %% gradient
   
   imgfm = immask(imgf, imgmask);
   
   imggrad = imgradient3D(imgfm);
   imggrad = mat2gray(imggrad);
   
   if verbose 
      figure(7); clf
      implottiling(imggrad, 'tiling', tiling)
   end
    
   %%
   
   imgfmg  = imgfm - 0. * imclip(imggrad, 0, 1);
   imgfmg(imgfmg < 0 ) = 0;

   if verbose
      figure(15); clf
      set(gcf, 'Name', ['Time: ', num2str(t), ' Mask']);
      %implottiling(imgfmg(:,:,:), 'tiling', tiling);
   end
    
   
   %% Seeding
   %imgf = filterMedian(imgs, [5, 5, 4]);
   %imgf = log(mat2gray(imgf)+eps) + 5;
   %imgf(imgf <0) = 0;
   
   imgfm = immask(imgfmg, imgmask);

   %imgf = filterLoG(max(imgfm(:)) - imgfm, [20,20,8]);
   imgf = filterDisk(imgfm, [10, 10, 4]);
 
   imgmax = imextendedmax(mat2gray(imgf), 0.05);
   %stackmax = imregionalmax(stackf);
   %imgmax = immask(imgmax, imgmask);
   %
   %stackmax = imerode(stackmax, strel('disk',3));
   %stackmax = imdilate(stackmax, strel('disk',3));
   
   if verbose
      %figure(4); clf;
      %implottiling(imgf, 'tiling', tiling);
      
      
      % figure(10)
      % clf
      % set(gcf, 'Name', ['PreSeeding: ' filenames]);
      % is = imgf - mean(imgf(:));
      % is(is < 0) = 0;
      % implot3d(is);
      %
      % figure(11)
      % clf
      % set(gcf, 'Name', ['Seeding: ' filenames]);
      % implot3d(mat2gray(imgmax));
      % overlay in imagej
      stackc = gray2rgb(imgmax);
      %size(stackc);
      stackc(:,:,:,2:3) = 0;
      stackc = stackc + gray2rgb(mat2gray(imgf));
      stackc(stackc > 1) = 1;
      
      figure(15); clf
      set(gcf, 'Name', ['Time: ', num2str(t), ' Seeding']);
      implottiling(stackc, 'tiling', tiling);
      
   end

   %% Watershed
   % 3d water shed on image + gradient
%    imgf = imgf;
%   
%    mi = max(imgf(:));
%    mg = max(imggrad(:));
%    imgws = imgf - 0.25 * mi / mg * imggrad;
%    imgws(imgws < 0) = 0;
%    
%    %imgmaxws = imdilate(imgmax, strel('disk', 4));
%    imgmaxws = imgseg  > 0;
   
   %imgmaxws = imdilate(imgmax, strel('disk', 4));
   imgmaxws =imgmax;

   %imgws = imgf;
   imgwsmin = imimposemin(iminvert(imgfmg), imgmaxws);
   imgseg = watershed(imgwsmin);
   imgseg = immask(imgseg, imgmask);
   
   if verbose
      %figure(13)
      %set(gcf, 'Name', ['Segmentation: ' filenames])
      %clf
      %implot3d(ws)
      %figure(14)
      %clf
      %set(gcf, 'Name', ['Segmentation: ' filenames])
      %imsurfaceplot3d(ws)
      %ijplot3d(imcolorize(ws) + gray2rgb(mat2gray(stackraw)))
      figure(14); clf
      set(gcf, 'Name', ['Time: ', num2str(t), ' Watershed']);
      imgovl = imoverlaylabel(4*imgf, impixelsurface(imgseg), false);
      implottiling(imgovl, 'tiling', tiling);
   end
   
    
   
   %% Create Frames
   objs  = label2DataObjects(imgseg, imgf);
   
   frame(t) = Frame('t', t, 'objects', objs);
    
end

%% Track result
clc

frames = frame;

param.load.min = 2;       % at least one object in frame
param.load.change = 0.2;  % at most 20% change in number of objects

param.print.load = true;
param.print.match.objects = true;
param.print.match.optimization = false;

param.optimize = true;

%vns = dir([resultdir, '*.csv']);
%frm = loadEmbryoDataFile([resultdir, vns(1).name])
 
%%
 
for t = 1:length(frames)
   
   dat = frames(t).r;
   %indx = dat(2,:) > 340;
   %frames(t).objects = frames(t).objects(indx);
   
end

[matches, costs] = matchFrames(frames, param);


%% plot the result 
% 
% for t =1:length(matches)
%    figure
%    clf
%    subplot(1,2,1)
%    plotMatchedObjects(matches(t))
%    title('matches')
%    subplot(1,2,2)
%    plotMatchedCostMatrix(matches(t), costs{t})
%    title('cost matrix')
% end


%% Determine trajectories

traj = findTrajectories(matches);

figure
clf
plotMatchedTrajectories(frames, traj)

%% Plot in 5d

moviedir = fullfile(resultdir, 'tracking');
mkdir(moviedir)

outfile = fullfile(moviedir, 'T<T,3>Z<Z,2>.tiff');

ids = {traj.objids};
ll = cellfun(@length, ids);

% only keep full length traj
%pos = ll == is.dataSize('F');

tmax = isd.cellSize('F');
ntraj = length(traj);

cols = colorcube(ntraj);
cols = cols(randperm(size(cols,1)),:);

imgmovie = zeros([isd.dataSize, 3, tmax]);

for t = 1:tmax
   
   % load data
   img = isd.data('T', t, 'C', 1);
   imgseg = isd.data('F', t);
   
   % determine color map
   [tpos, topos, tobjs] = traj.frameSlice(t);
   
   %maps is [tobjs.id] -> tpos
   colm = zeros(length(tpos), 3);
   colm([tobjs.id], :) = cols(tpos, :);

   imgmovie(:,:,:,:, t) = imoverlaylabel(mat2gray(img), impixelsurface(imgseg), false, 'color.map', colm,  'color.shuffle', 'noshuffle');
   
   imgwrt = imfrmtReformat(imgmovie(:,:,:,:,t), 'XYZC', 'XYCZ');
   
   for z = 1:size(imgwrt,4)
      imwrite(imgwrt(:,:,:,z), tagExpressionToString(outfile, 'T', t, 'Z', z));
   end
end


%%
ijplot5d(moviedir)


















'Z', z
%%
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














