%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Brain Slices %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

initialize()
bfinitialize
%initializeParallelProcessing(12)

addpath('./Scripts/User/Zeeshan/BrainSlices/2015_11');

datadir = '/home/ckirst/Data/Science/Projects/StemCells/Experiment/BrainSections/';
savedir   = '/home/ckirst/Data/Science/Projects/StemCells/Experiment/BrainSections/Analysis/2016_02/Data/'
resultdir = '/home/ckirst/Data/Science/Projects/StemCells/Experiment/BrainSections/Analysis/2016_02/';

dataname = 'DAPT D90 CTIP2 SATB2 TBR1 NURR1';

datafield = ' 1';

verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Source %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

is = ImageSourceBF([datadir  dataname '.lsm']);
clc
is.printInfo

%%

% is.setReshape('S', 'Uv', [1, 1]);
% is.setCellFormat('UV')
% is.setRange('C', 1);
% clc; is.printInfo

%%
% is.setRange('C', 1);
% is.plotPreviewStiched('overlap', 102, 'scale', 0.1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inspect Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%region = struct('U', [8:11], 'V', [9:12], 'C', 1);
%region = struct('U', 10, 'V', 11, 'C', 1, 'X', 1:300, 'Y', 1:300);
%region = struct('U', 10, 'V', 11, 'C', 1);
% region = struct('U', 7:10, 'V', 1:4, 'C', 1);
% 
% imgs = is.cell(region);
% size(imgs)
% 
% if verbose 
%    figure(10); clf
%    implottiling(imgs, 'link', false)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Align
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc
% sh = alignImages(imgs, 'alignment', 'RMS', 'overlap.max', 150);
% stmeth = 'Interpolate';
% %stmeth = 'Pyramid';


% %%
% if verbose 
%    img = stitchImages(imgs, sh, 'method', stmeth);
%    figure(1); clf; implot(img)
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stich Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nuclear marker
is.resetRange();

nch = 5;

imgsDAPI  = is.cell('C', 1);
imgsTBR1  = is.cell('C', 2);
imgsSATB2 = is.cell('C', 3);
imgsCTIP2 = is.cell('C', 4);
imgsNURR1 = is.cell('C', 5);

nch = 4;
for i = 1:nch
    imgsRaw{i}  = is.cell(region, 'C', i);
end

imgsAllRaw = {imgsRaw{3}, imgsRaw{4}, imgsRaw{2}, imgs{5}, imgs{1}};



channellabel = {'CTIP2', 'TBR1', 'SATB2', 'NURR1', 'DAPI'};
channelcolors  = {[1,0,0], [0,1,0], [1,1,0], [1, 0, 1], [0,0,1]};

%% Preprocess

imgs = cell(1,nch);
parfor i = 1:nch
   imgs{i} = imgsRaw{i};
   %thrs = 0 * [0, 720, 430, 0, 600];
   %imgs{i} = cellfunc(@(x) (x > thrs(i)) .* x,  imgsAll{i});
   %imgs{i} = cellfunc(@(x) x - imopen(x, strel('disk', 100)), imgs{i});   
   %imgs{i} = cellfunc(@(x) x - imopen(x, strel('disk', 10)), imgs{i});   
   %imgs{i} = cellfunc(@(x) filterAutoContrast(x/max(x(:))), imgs{i});
   imgs{i} = cellfunc(@(x) filterMedian(x, 3), imgs{i});
   imgs{i} = cellfunc(@(x) x - imopen(x, strel('disk', 75)), imgs{i});   
   imgs{i} = cellfunc(@(x) x - imopen(x, strel('disk', 10)), imgs{i});   
end

%%
img = cell(1,nch);
imgRaw = cell(1,nch);

for i = 1:nch
   img{i} = mat2gray(imgs{i}{1});
   imgRaw{i} = mat2gray(imgsRaw{i}{1});
   
   %img{i} = img{i} - imopen(img{i}, strel('disk', 20));
   %img{i} = filterAutoContrast(img{i}/max(img{i}(:)));
   
   if verbose
      figure(16);
      set(gcf, 'Name', 'Preprocessed Data')
      subplot(nch+1,1,i); 
      hist(img{i}(:), 256);
      title(channellabel{i});
   end
   
   img{i} = mat2gray(img{i});
   %img{i} = mat2gray(imclip(img{i}, 0, clip(i)));
end

if verbose
   figure(1); clf; colormap jet
   set(gcf, 'Name', 'Stitched and Preprocessed Data')
   implottiling({img{:}; imgRaw{:}}', 'titles', {channellabel{:}, channellabel{:}})
end



%% Color and Overall Intensity Image

nch = 5;

weights = [0.5, 0.5, 0.5, 0, 1]; % exclude NURR1 as it is fuzzy
weightsfull = 2*[0.5, 0.5, 0.5, 0.3, 0.5]; %exclude DAPI

imgC = zeros([size(img{1}), 3]);
imgCfull = zeros([size(img{1}), 3]);

for i = 1:nch
   for c = 1:3
      imgC(:,:,c) = imgC(:,:,c) + channelcolors{i}(c) * weights(i) * img{i}; % - imopen(imgsWs{i}, strel('disk', 10)));
      imgCfull(:,:,c) = imgCfull(:,:,c) + channelcolors{i}(c) * weightsfull(i) * img{i}; 
   end
end
    
%imgC = imclip(imgC, 0, 2500);
imgC = imgC / max(imgC(:));
imgC = filterAutoContrast(imgC);

%imgCfull = imgCfull / max(imgCfull(:));
imgCfull = filterAutoContrast(imgCfull);

imgI = sum(imgC,3);

if verbose
   figure(2); clf;
   set(gcf, 'Name', 'Image Data')
   implottiling({imgC, imgCfull, img{5}, imgI}, 'titles' , {'color', 'color full', 'DAPI','intensity'});
end

%% Save preporcessing state

close all
clear('imgs', 'imgsRaw', 'i', 'c')

%%
savefile = [savedir, dataname, datafield, '_Raw.mat'];
save(savefile)
%load(savefile)

%% Restrict to sub regoin
sub = false;
%subregion = {1:2500, 1:2500};
%subregion = {1:1000, 1:1000};

if sub
   imgCsub = imgC(subregion{1}, subregion{2}, :);
   
   if verbose
      figure(4); clf
      implot(imgCsub)
   end
    
   imgC = imgCsub;
   for i = 1: length(img)
      img{i} = img{i}(subregion{:});
      imgRaw{i} = imgRaw{i}(subregion{:});
   end
end
clear('imgCsub')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if verbose
  for i = 1:nch
    figure(16);
    set(gcf, 'Name', 'Preprocessed Data')
    subplot(nch+1,1,i); 
    hist(img{i}(:), 256);
    title(channellabel{i});
  end
end

%%

%figure(1); clf; colormap jet
%implottiling(imgsSt, 'titles', chlabel)

%imgm =  imgsSt{5} - imopen(imgsSt{5} , strel('disk', 10));
imgM = img{1}+ img{2};
imgMRaw = imgRaw{1}+ imgRaw{2};

if verbose
    figure(12); clf;
    set(gcf, 'Name', 'Histograms')
    subplot(1,4,1);
    hist(imgm(:), 256)
    title('imgM');
    
    subplot(1,4,2);
    hist(imgRaw{5}(:), 256)
    title('DAPI Raw');
    
    subplot(1,4,3);
    hist(img{5}(:), 256)
    title('DAPI');
    
    subplot(1,4,4);
    hist(imgMRaw(:), 256)
    title('imgMRaw');
end
    
%%
if verbose 
    figure(15); clf
    implottiling({img{nch}, imopen(img{nch}, strel('disk', 5)); mat2gray(img{1}+img{2}), mat2gray(imgRaw{1} + imgRaw{2})}, 'titles', {'DAPI', 'open'; 'imgM', 'imgMRaw'})
end

%%

%imgmaskLo = imgMRaw > 0.05;
imgmaskLo = imgM > 0.05;
imgmaskLo = imopen(imgmaskLo, strel('disk', 2));
%imgmaskLo = postProcessSegments(bwlabeln(imgmaskLo), 'volume.min', 15, 'fillholes', false) > 0;

imgmask = imopen(imgmaskLo, strel('disk', 2));
%imgmask = postProcessSegments(bwlabeln(imgmask), 'volume.min', 15, 'fillholes', false) > 0;

if verbose
   %max(img(:))   
   figure(21); clf;
   set(gcf, 'Name', 'Masking')
   implottiling({mat2gray(imgM),  imgmask;
                 img{nch}, imgRaw{nch}
                 }, 'titles', {'imgm', 'Mask'})
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cell Centers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%imgBM = filterBM(mat2gray(imgC), 'profile', 'np', 'sigma', 9);
%imgCf = imgBM;
%imgI  = sum(imgBM,3);
%imgI = imgC(:,:,3);
%imgI  = imgC(:,:,1) + imgC(:,:,2);
imgI = img{1}+ img{2};

if verbose
   figure(3); clf;
   set(gcf, 'Name', 'BM')
   implottiling({imgC; img{nch}; mat2gray(sum(imgC,3)); mat2gray(imgI)})
end


%% DoG filter + Cell Detection

%imgBMI = sum(imgBM, 3);
%imgI = imgBM(:,:,3);
%imgI = imgCf(:,:,3) + imgCf(:,:,1) + imgCf(:,:,2);

disk = strel('disk',3);
disk = disk.getnhood;
%imgV = 1 ./ filterStd(imgI, disk) .* imgI ;
imgV = imgI - 1.5 * filterStd(imgI, disk);
%imgV(imgV > 100) = 100;

if verbose
   figure(31); clf 
   set(gcf, 'Name', 'Variance')
   implottiling({imgC, mat2gray(imgV), mat2gray(filterDoG(imgV, 5)),  mat2gray(imgI)}');
end

%%
%imgf = imgV;
imgf = filterDoG(imgV, 8); % .* imgI;
%imgf = filterSphere(imgV, 4);
%imgf = filterDisk(imgI,12,1);
%imgF = imgBMI;


imgf2 = imgf; % - imopen(imgf, strel('disk', 1));

imgmax = imextendedmax(mat2gray(imgf2), 0.01);
%imgmax  = immask(imgmax, imgmask);
%imglab = bwlabeln(imgmax);

if verbose 
   figure(32); clf   
   set(gcf, 'Name', 'Maxima')
   %implottiling({imgC, imgsSt{5}; imgI,  imgf;
   %              imoverlay(mat2gray(imgf2), imgmax), imoverlay(imgC, imgmax); 10 * imgI, 100 * imgsSt{5}}')
   implottiling({imoverlay(mat2gray(imgf2), imgmax), imoverlay(imgC, imgmax)}')
end


%figure(67); clf;
%implottiling(bwlabeln(imgmax))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cell Shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Watershed

%imgWs = imimposemin(iminvert(imgf), imgmax);
imgWs = imimposemin(iminvert(imgMRaw), imgmax);
imgWs = watershed(imgWs);
imgWs = immask(imgWs, imgmask);
imgWs = imlabelseparate(imgWs);
imgWs = postProcessSegments(bwlabeln(imgWs), 'volume.min', 15);

if verbose
   %figure(7); clf;
   %implot(imgWs)

   %imglab = bwlabel(imgWs);

   figure(8); clf;
   set(gcf, 'Name', 'Cells from Watershed')
   imgsurf = impixelsurface(imgWs);
   implottiling({ imoverlay(imgC, imgmax); imoverlaylabel(imgM, imgsurf, false);imoverlaylabel(imgC, imgsurf, false)})
   
end
   
imgS = imgWs;      


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filter Cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mode = 'MeanIntensity';

stats = imstatistics(imgS, {mode, 'Volume'}, immask(img{5}, imgmask));

if verbose
   figure(7); clf;
   subplot(1,2,1)
   hist([stats.(mode)], 128)
   title(mode)
   set(gcf, 'Name', 'WS Stats')

   subplot(1,2,2);
   hist([stats.('Volume')], 128)
   title('Volume')
   
   %hist(imgI(:), 256)
end

%% Remove weak cells

[imgSP, stats] = postProcessSegments(imgS, imgRaw{5}, 'intensity.mean.min', 0.1, 'volume.min', 30, 'volume.max', 220, 'fillholes', false);

if verbose
   imgsurf = impixelsurface(imgSP);
   figure(5); clf;
   implottiling({imoverlaylabel(imgC, imgsurf, false), imoverlaylabel(img{5}, imgsurf, false)}');
end

fprintf('watershed    : %d cells\n', max(imgS(:)))
fprintf('postprocessed: %d cells\n', max(imgSP(:)))

imglab = imrelabel(imgSP);


%% Cell Center Label 

%imglabc = imlabelapplybw(imglab, @bwulterode);

stats = imstatistics(imglab, stats, {'Centroid', 'PixelIdxList'});

centers = fix([stats.Centroid]);

imglabc = zeros(size(imglab));
ind = imsub2ind(size(imglab), centers');
imglabc(ind) = 1:length(ind);
imglabc = imdilate(imglabc, strel('disk', 4));

if verbose
    h = figure(78); clf;
    implottiling({imoverlaylabel(img{5}, imglabc > 0, false, 'color.map', [[0,0,0]; [1,0,0]]), 
                  imoverlaylabel(imgC, imglabc > 0, false, 'color.map', [[0,0,0]; [1,0,0]])}');
              
    %saveas(h, [datadir, dataname, datafield, '_CellDetection.pdf'])
end

%% Cell Outline

imglabcSurf = impixelsurface(imglabc);

if verbose
    h = figure(78); clf;
    set(gcf, 'Name', 'Cell Centers')
    implottiling({imoverlaylabel(img{5}, imglabcSurf > 0, false, 'color.map', [[0,0,0]; [1,0,0]]), 
                  imoverlaylabel(imgC, imglabcSurf > 0, false, 'color.map', [[0,0,0]; [1,0,0]])});
              
    %saveas(h, [datadir, dataname, datafield, '_CellDetection.pdf'])
end

statsSurf = imstatistics(imglabcSurf, {'PixelIdxList'});


%% Save Image Preporcessing Result


close all
clear('imgV', 'imgm', 'imgsurf', 'imglabcSurf', 'imgf', 'imgf2', 'imgS', 'imgSP', 'imgWs', 'disk', 'i', 'imgmax', 'centers')

%%

savefile = [savedir, dataname, datafield, '_ImageProcessing.mat'];
save(savefile)
%load(savefile)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Statistics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stats = imstatistics(imglabc, {'PixelIdxList', 'Centroid'});

%% Intensities

%mode = 'MedianIntensity';
mode = 'MeanIntensity';

clear statsCh
parfor c = 1:nch
   statsChRaw{c} = imstatistics(imglabc, mode,  imgRaw{c});
   statsCh{c} = imstatistics(imglabc, stats, mode,  img{c});
end

clear statsChN
parfor c = 1:nch
   statsChN{c} = statsCh{c};
   statsChRawN{c} = statsChRaw{c};
   if c < 5
      for i = 1:length(statsCh{c})
         statsChN{c}(i).(mode) = statsCh{c}(i).(mode) / statsCh{5}(i).(mode);
         statsChRawN{c}(i).(mode) = statsChRaw{c}(i).(mode) / statsChRaw{5}(i).(mode);
      end
   end
end


%% Measurements

% image sets to measure on -> some might have bad background -> used background reomed image for those

imgsStMeasure = {imgRaw{1}, imgRaw{2}, imgRaw{3}, img{4}, imgRaw{5}};
statsMeasure = {statsChRaw{1}, statsChRaw{2}, statsChRaw{3}, statsCh{4}, statsChRaw{5}};


%% Scatter Plot Colored Expression

if verbose 
    clip = 0.85;

    cols  = {[0,1,0], [0,0,1], [1,0,1], [1,0,0], [0,0,1]};
    %cols  = {0 * [0,1,0], [-1,-1,0],0* [0.75,0,0.25], [0,-1,-1] + 0* [0.25, 0.25, 0.25], [0,0,1]};
 
    xy = [stats.Centroid]';
    ncells = size(xy,1);

    %cdat = ones(ncells, 3);
    
    figure(79); clf;
    set(gcf, 'Name', 'Separate Expressions')
    for c = 1:nch-1
        cdat = zeros(ncells,3);
       fi = [statsMeasure{c}.(mode)]';
       fis = sort(fi);
       ncs = length(fis);
       fith = int64(clip * ncs);
       fi = imclip(fi, 0, fis(fith)) / fis(fith);
       cdat = cdat + fi * cols{c};
       
       subplot(2,2,c)
       %scatter(xy(:,1), xy(:,2), 30, cdat, 'filled');
       colormap jet
       scatter(xy(:,1), xy(:,2), 20, fi, 'filled');
       title(channellabel{c});
    end
end

%%

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 20  20/ 1.5656 * 4];
fig.PaperPositionMode = 'manual';
set(gcf, 'papersize', [20, 20/ 1.5656 * 4]);

saveas(fig, [datadir dataname datafield '_Quantification.png']);


%% Expression in Space - Separate Channels

xy = [stats.Centroid]';

ha = figure(21); clf;
set(gcf, 'Name', 'Scatter Expressions')
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 30 20];
fig.PaperPositionMode = 'manual';
set(gcf, 'papersize', [30,20]);

for c = 1:nch
   fi = [statsMeasure{c}.(mode)]';
   
   fis = sort(fi);
   ncs = length(fis);
   p90 = int64(0.95 * ncs);
   fi = imclip(fi, 0, fis(p90));
   
   figure(21);
   subplot(3,nch,c+2*nch);
   [nb,xb] = hist(fi, 256);
   bh = bar(xb,nb);
   set(bh,'facecolor','k');
   
   subplot(3,nch,c+nch);
   %imcolormap(cm{c});
   colormap jet
   scatter(xy(:,1), xy(:,2), 10, fi, 'filled');
   xlim([min(xy(:,1)), max(xy(:,1))]); ylim([min(xy(:,2)), max(xy(:,2))]);
   title(['M ' channellabel{c}]);
   %freezecolormap(gca)
   
   subplot(3,nch,c)
   %implot(imgsStRaw{c})
   
   fi = [statsChRaw{c}.(mode)]';
   fis = sort(fi);
   ncs = length(fis);
   p90 = int64(0.95 * ncs);
   fi = imclip(fi, 0, fis(p90));
   colormap jet
   scatter(xy(:,1), xy(:,2), 10, fi, 'filled');
   xlim([min(xy(:,1)), max(xy(:,1))]); ylim([min(xy(:,2)), max(xy(:,2))]);
   title(['R ' channellabel{c}]);
   
   %h = figure(22); clf
   %colormap jet
   %scatter(xy(:,1), xy(:,2), 10, fi, 'filled');
   %xlim([0, size(imglab,1)]); ylim([0, size(imglab,2)]);
   %title(chlabel{c}); 
   %colorbar('Ticks', [])
   %saveas(h, [datadir dataname datafield '_Quantification_' lab{c} '.pdf']);
end

%saveas(ha, [datadir dataname datafield '_Quantification_Normalized.pdf']);


%% Flourescence Expression

fi = cell(1,3);
for c = 1:nch
   fi{c} = [statsMeasure{c}.(mode)]';
   fis = sort(fi{c});
   ncs = length(fis);
   p95 = int64(0.97 * ncs);
   fi{c} = imclip(fi{c}, 0, fis(p95));
   %fi{c} = mat2gray(fi{c});
end
 
pairs = {[1,2], [1,3], [1,4], [2, 3], [2,4], [3,4]};

np = length(pairs);

figure(22); clf;
for n = 1:np
   subplot(2, np/2,n)
   %fi = imclip(fi, 0, cl(c));
   %scatter(fi{pairs{n}(1)}, fi{pairs{n}(2)}, 10, 'b', 'filled');
   scattercloud(fi{pairs{n}(1)}, fi{pairs{n}(2)}, 100, .5, 'r+', [jet(256); repmat([0.5,0,0],1000,1)] );
   %scattercloud(fi{pairs{n}(1)}, fi{pairs{n}(2)}, 75, 1, 'r+', [flipud(gray(256)); repmat([0,0,0], 400,1)] );
   %xlim([0,1]); ylim([0,1]);
   xlabel(channellabel{pairs{n}(1)}); ylabel(channellabel{pairs{n}(2)});
   %freezecolormap(gca)
end
set(gcf, 'Name', 'Normalized Expression')

if verbose > 1
   for n = 1:np
      h = figure(50+n); clf;
      %fi = imclip(fi, 0, cl(c));
      %scatter(fi{pairs{n}(1)}, fi{pairs{n}(2)}, 10, 'b', 'filled');
      scattercloud(fi{pairs{n}(1)}, fi{pairs{n}(2)});
      xlim([0,1]); ylim([0,1]);
      xlabel(channellabel{pairs{n}(1)}); ylabel(channellabel{pairs{n}(2)});

      %saveas(h, [datadir dataname datafield '_Quantification_Scatter_' chlabel{pairs{n}(1)} ,'_' chlabel{pairs{n}(2)} '.pdf']);
      %freezecolormap(gca)
   end
end

%%

savefile = [datadir, dataname, datafield, '_Results.mat'];
save(savefile)
%load(savefile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Classify Cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Mask for ROIs

%h  = figure(300); clf;
%implot(imgsStRaw{4});
%
%p = impoly
%impoly(gca, xy)
%%

nx = size(imgC,1);
ny = size(imgC,2);

%xy  = [699.92088607595 1369.45632911392;732.346202531646 678.391772151899;716.133544303798 341.979113924051;509.422151898735 202.144936708861;343.242405063291 262.942405063292;233.806962025317 495.999367088608;162.876582278481 718.923417721519;162.876582278481 881.050000000001;136.531012658228 1067.49556962025;211.514556962025 1120.18670886076;185.168987341772 1207.32974683544;126.398101265823 1389.72215189873;294.604430379747 1489.0246835443;545.900632911392 1493.07784810127;665.468987341772 1472.81202531646];
%imggood = poly2mask(xy(:,2), xy(:,1), size(imgI,1), size(imgI,2));
imggood = ~logical(zeros(nx,ny));

if verbose > 1
   figure(6); clf;
   set(gcf, 'Name', 'Bad Image Region');
   implot(imggood)
end

xy = fix([stats.Centroid]);
xy(1,xy(1,:) > nx) = nx; xy(1,xy(1,:) < 1) = 1;
xy(2,xy(2,:) > ny) = ny; xy(2,xy(2,:) < 1) = 1;
   
isgood = zeros(1, size(xy,2));
for i = 1:size(xy, 2)
   isgood(i) = imggood(xy(1,i), xy(2,i));
end
isgood = isgood > 0;

fprintf('%d / %d cells in ROI\n', sum(isgood), size(xy,2))
 

%% histograms

figure(10); clf
set(gcf, 'Name', 'Histograms');
for c= 1:nch
   subplot(4,nch,c); 
   %implot(imgsStRaw{c})
   dat = [statsChRaw{c}.(mode)];
   hist(dat(isgood), 256)
   title(['Raw ' channellabel{c}])
   
   subplot(4, nch, c+nch);
   dat = [statsChRawN{c}.(mode)];
   hist(dat(isgood), 256)
   title(['Raw Norm. ' channellabel{c}])
   
   subplot(4, nch, c+2*nch); 
   dat = [statsCh{c}.(mode)];
   hist(dat(isgood), 256)
   title(['Corr ' channellabel{c}])
   
   subplot(4, nch, c+3*nch); 
   dat = [statsChN{c}.(mode)];
   hist(dat(isgood), 256)
   title(['Corr Norm.' channellabel{c}])
   
end

%% Classify

%cth = {0.16, 0.35, 0.18, 0.18, 0}; % measured
cth = {0.15, 0.1, 0.175, 0.1, 0}; % measured
neuroBaseColor = {[1,0,0], [0,1,0], [1,1,0], [1,0,1],[0,0,1]};

neuroClass = zeros(nch, length(isgood));
for c = 1:nch
   neuroClass(c,:) = double(and([statsMeasure{c}.(mode)] > cth{c}, isgood));
end
neuroClassTotal = int32(fix(neuroClass(1,:) + 2 * neuroClass(2,:) + 4 * neuroClass(3,:) + 8 * neuroClass(4,:)));

if verbose
    %imgClass1 = cell(1,nch);
    imgClass1 = cell(1,2);
    for c = 1:nch
       neuroColor1 = {[0,0,0], neuroBaseColor{c}};
       neuroColor1 = neuroColor1(neuroClass(c,:)+1);

       R = imgsStMeasure{c}; G = R; B = R;
       for i = 1:length(stats);
          if neuroClass(c,i) > 0
             R(statsSurf(i).PixelIdxList) =  neuroColor1{i}(1);
             G(statsSurf(i).PixelIdxList) =  neuroColor1{i}(2);
             B(statsSurf(i).PixelIdxList) =  neuroColor1{i}(3);
          else
             R(statsSurf(i).PixelIdxList) =  neuroBaseColor{5}(1);
             G(statsSurf(i).PixelIdxList) =  neuroBaseColor{5}(2);
             B(statsSurf(i).PixelIdxList) =  neuroBaseColor{5}(3);
          end
       end
       imgClass1{c} = cat(3, R, G, B);

       R = imgRaw{c}; G = R; B = R;
       for i = 1:length(stats);
          if neuroClass(c,i) > 0
             R(statsSurf(i).PixelIdxList) =  neuroColor1{i}(1);
             G(statsSurf(i).PixelIdxList) =  neuroColor1{i}(2);
             B(statsSurf(i).PixelIdxList) =  neuroColor1{i}(3);
          else
             R(statsSurf(i).PixelIdxList) =  neuroBaseColor{5}(1);
             G(statsSurf(i).PixelIdxList) =  neuroBaseColor{5}(2);
             B(statsSurf(i).PixelIdxList) =  neuroBaseColor{5}(3);
          end
       end
       imgClassRaw{c} = cat(3, R, G, B);
    end


    h = figure(7); clf;
    implottiling({imgClassRaw{:}; imgClass1{:}}', 'titles', [channellabel(1:nch)', channellabel(1:nch)']);
    set(gcf, 'Name', 'Classifications');

    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 30];
    fig.PaperPositionMode = 'manual';
    set(gcf, 'papersize', [8,30]);

    %saveas(h, [datadir dataname datafield '_Classification' '.pdf']);
end



%%

cc = num2cell(sum(neuroClass, 2));
tb = table(cc{:}, 'VariableNames', channellabel)

writetable(tb, [datadir dataname datafield '_Class_Count.txt'])



%% Image + Classes + Combined Classes

if verbose
    % neuroClassColor = num2cell([0.285714, 0.285714, 0.285714;
    %    0.929412, 0.290196, 0.596078;
    %    0.462745, 0.392157, 0.67451;
    %    0.709804, 0.827451, 0.2;
    %    0.984314, 0.690196, 0.25098;
    %    0.976471, 0.929412, 0.196078;
    %    0.584314, 0.752941, 0.705882;
    %    0, 0.682353, 0.803922
    %    ], 2);

    %neuroClassColor = num2cell(colorcube(max(neuroClassTotal(:))), 2);

    ncls = 2^(nch-1);
    neuroClassColor = cell(1, ncls);
    for i = 1:ncls
       cb = dec2bin(i-1);
       cb = cb(end:-1:1);
       %cb = [cb(end:-1,1), '0000'];
       %cb = cb(1:nch);

       col = [0,0,0]; k = 0;
       for l = 1:length(cb)
          if cb(l) == '1'
             col = col + neuroBaseColor{l};
             k = k + 1;
          end
       end
       if k > 1
          col = col / k;
       end

       neuroClassColor{i} = col;
    end

    R = imgRaw{c}; G = R; B = R;
    for p = 1:length(stats);
       nc = neuroClassColor{neuroClassTotal(p)+1};

       R(stats(p).PixelIdxList) =  nc(1);
       G(stats(p).PixelIdxList) =  nc(2);
       B(stats(p).PixelIdxList) =  nc(3);
    end

    imgClass = cat(3, R, G, B);

    % 
    % figure(7); clf;
    % implottiling({imgsSt{1:nch}, imgCf, imgClass}', 'tiling', [3,2]);


    h = figure(7); clf;
    implottiling({imgRaw{1:nch}, imgC; imgClass1{:}, imgClass}, 'titles', [[channellabel(1:nch)'; {'merge'}], [channellabel(1:nch)'; {'merge'}]]');

    %saveas(h, [datadir dataname datafield '_Classification_Channels' '.pdf']);
    %saveas(h, [datadir dataname datafield '_Classification_Image' '.pdf']);
end


%% Co-localization - NURR1 cells

if verbose
   neuroNurrColor = num2cell([0.285714, 0.285714, 0.285714;
       0.584314, 0.752941, 0.705882;
       0.984314, 0.690196, 0.25098;
       0.709804, 0.827451, 0.2;
       0.462745, 0.392157, 0.67451;
       0.976471, 0.929412, 0.196078;
       0.929412, 0.290196, 0.596078;
       0, 0.682353, 0.803922
       ], 2);

    nurr1ids = find(neuroClass(4,:));
    %nurr1xy = xy(nurr1ids,:);
    nurr1xy = xy(:, nurr1ids);
    nurr1NeuroClassTotal = neuroClassTotal(nurr1ids);
    
    ord = [9,10,11,12,13,14,15,16];
    h = figure(66); clf; 
    hold on
    k = 1;
    for i = ord
        ids = nurr1NeuroClassTotal + 1 == i;
        %scatter(nurr1xy(ids,1), nurr1xy(ids,2), 20, neuroNurrColor{k}, 'filled');
        scatter(nurr1xy(2,ids), size(imgI,1) - nurr1xy(1,ids), 20, neuroNurrColor{k}, 'filled');
        k = k + 1;
    end
end

%%

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition =  [0 0  20/ 1.5656  20 ];
fig.PaperPositionMode = 'manual';
set(gcf, 'papersize',  [20/ 1.5656 , 20]);

saveas(fig, [datadir dataname datafield '_Colocalization_Nurr1+.png']);

%% Class in Space

if verbose
    
    isgood = logical(isgood);

    xy1 = [stats.Centroid]';
    xy1 = xy1(isgood, :)

    h = figure(27)

    cdat = (cell2mat({[0,0,0], neuroClassColor{:}}'))

    scatter(xy1(:,1), xy1(:,2), 20, cdat(neuroClassTotal(isgood) + 2, :), 'filled');
    xlim([0, size(imglab,1)]); ylim([0, size(imglab,2)]);
    %title(chlabel{c});
    %freezecolormap(gca)

    pbaspect([1,1,1])

    %saveas(h, [datadir dataname datafield '_Quantification_Space.pdf']);
end


%% Histogram of Cell Classes
 
ncls = length(neuroClassColor);

nstat = zeros(1,ncls);
for i = 1:ncls
   id = num2cell(dec2bin(i-1));
   id = id(end:-1:1);
   for k = length(id)+1:nch
      id{k} = '0';
   end

   clslab{i} = '';
   for c = 1:nch
      if id{c} == '1'
         clslab{i} = [clslab{i}, '_', channellabel{c}];
      end
   end
   clslab{i} = clslab{i}(2:end);
   if isempty(clslab{i})
      clslab{i} = 'None';
   end
   
   nstat(i) = sum(neuroClassTotal + 1 == i);
end

nc = num2cell(nstat);
tb = table(nc{:}, 'VariableNames', clslab)

writetable(tb, [datadir dataname datafield '_Counts_All.txt'])

% neuroClassColor = num2cell([0.285714, 0.285714, 0.285714;
%    0.929412, 0.290196, 0.596078;
%    0.462745, 0.392157, 0.67451;
%    0.709804, 0.827451, 0.2;
%    0.984314, 0.690196, 0.25098;
%    0.976471, 0.929412, 0.196078;
%    0.584314, 0.752941, 0.705882;
%    0, 0.682353, 0.803922
%    ], 2);

h = figure(10); clf; hold on

k = 1;
%ord = [1, 5, 3, 6, 7, 8, 4, 2];
ord = 1:ncls
for i = ord
   bar(k, nstat(i), 'FaceColor', neuroClassColor{i});
   
   k = k + 1;
end

set(gca, 'XTickLabel', '')  
xlabetxt = strrep(clslab(ord), '_', '+');
n = length(clslab);
ypos = -max(ylim)/50;
text(1:n,repmat(ypos,n,1), xlabetxt','horizontalalignment','right','Rotation',35,'FontSize',15)

%saveas(h, [datadir dataname datafield '_Classification_Statistics' '.pdf']);


%% Histogram Classes NURR1 positive 

h = figure(11); clf; hold on

neuroNurrColor = num2cell([0.285714, 0.285714, 0.285714;
   0.584314, 0.752941, 0.705882;
   0.984314, 0.690196, 0.25098;
   0.709804, 0.827451, 0.2;
   0.462745, 0.392157, 0.67451;
   0.976471, 0.929412, 0.196078;
   0.929412, 0.290196, 0.596078;
   0, 0.682353, 0.803922
   ], 2);

k = 1;
%ord = [1, 5, 3, 6, 7, 8, 4, 2];
%ord = 1:ncls
ord = [9,10,11,12,13,14,15,16];
for i = ord
   bar(k, nstat(i), 'FaceColor', neuroNurrColor{k});
   k = k + 1;
end

set(gca, 'XTickLabel', '')  
xlabetxt = strrep(clslab(ord), '_', '+');
n = length(xlabetxt);
ypos = -max(ylim)/50;
text(1:n,repmat(ypos,n,1), xlabetxt','horizontalalignment','right','Rotation',35,'FontSize',15)

%saveas(h, [datadir dataname datafield '_Classification_Statistics' '.pdf']);



%% Save the Statistics to File

statsfile = [datadir dataname datafield '_statistics' '.mat'];
save(statsfile, 'chlabel', 'stats', 'statsCh', 'statsChN', 'statsMeasure', 'neuroClass', 'lineshires', 'dinner', 'douter');
%load(statsfile)

%% Save For Mathematica

centers = [stats.Centroid];

intensities = zeros(length(statsMeasure), size(centers,2));
intensities_normalized = zeros(length(statsMeasure), size(centers,2));
for c = 1:length(statsMeasure)
   intensities(c, :) = [statsMeasure{c}.(mode)];
   %intensities_normalized(c, :) = [statsChN{c}.MedianIntensity];
end

statsfilemathematica = [datadir dataname datafield '_statistics_mathematica' '.mat'];
save(statsfilemathematica, '-v6', 'chlabel', 'centers', 'intensities', 'intensities_normalized', 'neuroClass');

%%
