%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Brain Slices %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

initialize()
bfinitialize
initializeParallelProcessing(12)

addpath('./Scripts/User/Zeeshan/BrainSlices/2015_11');

%datadir = '/data/Science/Projects/StemCells/Experiment/BrainSections_2015_11/';
datadir = '/home/ckirst/Data/Science/Projects/StemCells/Experiment/BrainSections_2015_11/Analysis/';
%datadir = '/home/ckirst/Data/Science/Projects/StemCells/Experiment/BrainSections_2015_11/';
dataname = 'PCW15 Sample K Slide 73 CTIP2 SATB2 TLE4 NURR1 - new';

datafield = ' 1';

verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Source %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

is = ImageSourceBF([datadir  dataname '.lsm']);
clc
is.printInfo

%%

is.setReshape('S', 'Uv', [15, 19]);
is.setCellFormat('UV')
is.setRange('C', 1);
clc; is.printInfo

%%
is.setRange('C', 1);
is.plotPreviewStiched('overlap', 102, 'scale', 0.1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inspect Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%region = struct('U', [8:11], 'V', [9:12], 'C', 1);
%region = struct('U', 10, 'V', 11, 'C', 1, 'X', 1:300, 'Y', 1:300);
%region = struct('U', 10, 'V', 11, 'C', 1);
region = struct('U', 1:6, 'V', 8:10, 'C', 1);

imgs = is.cell(region);
size(imgs)

if verbose 
   figure(10); clf
   implottiling(imgs, 'link', false)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Align
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
sh = alignImages(imgs, 'alignment', 'RMS', 'overlap.max', 150);
stmeth = 'Interpolate';
%stmeth = 'Pyramid';


%%
if verbose 
   img = stitchImages(imgs, sh, 'method', stmeth);
   figure(1); clf; implot(img)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stich Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nuclear marker
is.resetRange();

nch = 5;

imgsDAPI  = is.cell(region, 'C', 1);
imgsSATB2  = is.cell(region, 'C', 2);
imgsCTIP2 = is.cell(region, 'C', 3);
imgsTLE4 = is.cell(region, 'C', 4);
imgsNURR1 = is.cell(region, 'C', 5);

imgsAllRaw = {imgsCTIP2, imgsTLE4, imgsSATB2, imgsNURR1, imgsDAPI};

chlabel = {'CTIP2', 'TLE4', 'SATB2', 'NURR1', 'DAPI'};

chcols  = {[1,0,0], [0,1,0], [1,1,0], [1, 0, 1], [0,0,1]};
cm      = {'r',      'g',    'y',     'm',       'b'};

%% Test Preprocessing
% 
% for i  = [4]
%    imgsAll{i} = imgsAllRaw{i}(1:2);
%    imgsAll{i} = cellfunc(@(x) filterMedian(x, 3), imgsAll{i});
%    imgsAll{i} = cellfunc(@(x) x - imopen(x, strel('disk', 75)), imgsAll{i});   
%    imgsAll{i} = cellfunc(@(x) x - imopen(x, strel('disk', 10)), imgsAll{i});   
% end
%    
% figure(12); clf;
% implottiling({imgsAllRaw{i}{1:2}; imgsAll{i}{:}})

%% Preprocess

parfor i = 1:nch
   imgsAll{i} = imgsAllRaw{i};
   %thrs = 0 * [0, 720, 430, 0, 600];
   %imgsAll{i} = cellfunc(@(x) (x > thrs(i)) .* x,  imgsAll{i});
   %imgsAll{i} = cellfunc(@(x) x - imopen(x, strel('disk', 100)), imgsAll{i});   
   %imgsAll{i} = cellfunc(@(x) x - imopen(x, strel('disk', 10)), imgsAll{i});   
   %imgsAll{i} = cellfunc(@(x) filterAutoContrast(x/max(x(:))), imgsAll{i});
   imgsAll{i} = cellfunc(@(x) filterMedian(x, 3), imgsAll{i});
   imgsAll{i} = cellfunc(@(x) x - imopen(x, strel('disk', 75)), imgsAll{i});   
   imgsAll{i} = cellfunc(@(x) x - imopen(x, strel('disk', 10)), imgsAll{i});   
end

%%

imgsSt = cell(1,nch);
imgsStRaw = cell(1,nch);

for i = 1:nch
   imgsStRaw{i} = stitchImages(imgsAllRaw{i}, sh, 'method', stmeth);
   imgsStRaw{i} = mat2gray(imgsStRaw{i});
   
   imgSt = stitchImages(imgsAll{i}, sh, 'method', stmeth);
   
   %imgSt = imgSt - imopen(imgSt, strel('disk', 20));
   %imgsS = filterAutoContrast(imgsS/max(imgsS(:)));
   
   if verbose
      figure(16);
      subplot(nch+1,1,i); 
      hist(imgSt(:), 256);
   end
   
   imgsSt{i} = mat2gray(imgSt);
   %imgsSt{i} = mat2gray(imclip(imgSt, 0, cl(i)));
end

if verbose
   figure(1); clf; colormap jet
   set(gcf, 'Name', 'Stitched and Preprocessed Data')
   implottiling({imgsSt{:}; imgsStRaw{:}}', 'titles', {chlabel{:}, chlabel{:}})
end

%% Save preporcessing state

savefile = [datadir, dataname, datafield, '.mat'];
%save(savefile)
load(savefile)

%%

nch = 5;

weights = [0.5, 0.5, 0.5, 0, 1]; % exclude NURR1 as it is fuzzy
chcols  = {[1,0,0], [0,1,0], [1,1,0], [1, 0, 1], [0,0,1]};

imgsWs = cellfunc(@(x,y) x * y, imgsSt(1:nch), num2cell(weights));

imgC = zeros([size(imgsWs{1}), 3]);
for i = 1:nch
   for c = 1:3
      imgC(:,:,c) = imgC(:,:,c) + chcols{i}(c) * (imgsWs{i}); % - imopen(imgsWs{i}, strel('disk', 10)));
   end
end
    
%imgC = imclip(imgC, 0, 2500);
imgC = imgC / max(imgC(:));
imgC = filterAutoContrast(imgC);

imgI = sum(imgC,3);

if verbose
   figure(2); clf;
   set(gcf, 'Name', 'Data')
   implottiling({imgC, imgsSt{5}, imgI}', 'titles' , {'color','DAPI','intensity'});
end

%% Restrict to sub regoin
sub = true;
subregion = {700:5200, 1:500};
%subregion = {700:5200, ':'};

if sub
   imgCsub = imgC(subregion{1}, subregion{2}, :);
   
   if verbose
      figure(4); clf
      implot(imgCsub)
   end
    
   imgC = imgCsub;
   for i = 1: length(imgsSt)
      imgsSt{i} = imgsSt{i}(subregion{:});
      imgsStRaw{i} = imgsStRaw{i}(subregion{:});
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(7); clf;
%implottiling({imgsSt{5}, imgsStRaw{5}, imgsSt{5} - imgsStRaw{5}})
hist(imgsSt{5}(:), 256)

%%
%figure(1); clf; colormap jet
%implottiling(imgsSt, 'titles', chlabel)

%imgm =  imgsSt{5} - imopen(imgsSt{5} , strel('disk', 10));
imgm = imgsSt{5};


imgmaskLo = imgm > 0.08;
imgmaskLo = imopen(imgmaskLo, strel('disk', 2));
imgmaskLo = postProcessSegments(bwlabeln(imgmaskLo), 'volume.min', 15, 'fillholes', false) > 0;


%figure(66); clf;
%implottiling({255*(postProcessSegments(bwlabeln(imgmaskLo), 'volume.min', 10, 'fillholes', false)>0); 2*(bwlabeln(imgmaskLo)); 255*imgmaskLo}')

%get rid of blod vessels
imgmaskHi1 = imgsStRaw{4} > 0.5;
imgmaskHi1 = imdilate(imgmaskHi1, strel('disk', 3));
imgmaskHi1 = postProcessSegments(bwlabeln(imgmaskHi1), 'volume.min', 500, 'fillholes', false) > 0;
imgmaskHi1 = not(imgmaskHi1);

imgmaskHi2 = imgsStRaw{1} > 1.9;
imgmaskHi2 = imdilate(imgmaskHi2, strel('disk', 3));
%imgmaskHi2 = postProcessSegments(bwlabeln(imgmaskHi2), 'volume.min', 100, 'fillholes', false) > 0;
imgmaskHi2 = not(imgmaskHi2);

imgmaskHi = and(imgmaskHi1, imgmaskHi2);

imgmask = and(imgmaskLo, imgmaskHi);
imgmask = imopen(imgmask, strel('disk', 2));
imgmask = postProcessSegments(bwlabeln(imgmask), 'volume.min', 15, 'fillholes', false) > 0;

if verbose
   %max(img(:))   
   figure(21); clf;
   set(gcf, 'Name', 'Masking')
   implottiling({imgsStRaw{5},imgmaskLo; 
                 imgsStRaw{4},imgmaskHi1;
                 imgsStRaw{1},imgmaskHi2;
                 mat2gray(imgm),  imgmask
                 }, 'titles', {'I', 'Lo', 'I','Hi','I', 'Mask'})
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cell Centers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%imgBM = filterBM(mat2gray(imgC), 'profile', 'np', 'sigma', 9);
%imgCf = imgBM;
%imgI  = sum(imgBM,3);
imgI = imgC(:,:,3);

if verbose > 1
   figure(3); clf;
   set(gcf, 'Name', 'BM')
   implottiling({imgC; mat2gray(sum(imgBM,3)); mat2gray(imgI)}')
end


%% DoG filter + Cell Detection

%imgBMI = sum(imgBM, 3);
%imgI = imgBM(:,:,3);
%imgI = imgCf(:,:,3) + imgCf(:,:,1) + imgCf(:,:,2);

disk = strel('disk',3);
disk = disk.getnhood;
%imgV = 1 ./ filterStd(imgI, disk) .* imgI ;
imgV = imgI - 1 * filterStd(imgI, disk);
%imgV(imgV > 100) = 100;

if verbose
   figure(31); clf 
   set(gcf, 'Name', 'Variance')
   implottiling({imgC, mat2gray(imgV), mat2gray(filterDoG(imgV, 5)),  mat2gray(imgI)});
end

%%
%imgf = imgV;
imgf = filterDoG(imgV, 10); % .* imgI;
%imgf = filterSphere(imgV, 4);
%imgf = filterDisk(imgI,12,1);
%imgF = imgBMI;


imgf2 = imgf; % - imopen(imgf, strel('disk', 1));

imgmax = imextendedmax(mat2gray(imgf2), 0.01);

if verbose 
   figure(32); clf   
   set(gcf, 'Name', 'Maxima')
   %implottiling({imgC, imgsSt{5}; imgI,  imgf;
   %              imoverlay(mat2gray(imgf2), imgmax), imoverlay(imgC, imgmax); 10 * imgI, 100 * imgsSt{5}}')
   implottiling({imoverlay(mat2gray(imgf2), imgmax), imoverlay(imgC, imgmax)})
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cell Shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Watershed

%imgWs = imimposemin(iminvert(imgf), imgmax);
imgWs = imimposemin(iminvert(immask(imgsStRaw{5}, imgmask)), imgmax);
imgWs = watershed(imgWs);
imgWs = immask(imgWs, imgmask);
imgWs = imlabelseparate(imgWs);

%[imgWs, stats] = postProcessSegments(bwlabeln(imgWs), 'volume.min', 0);

if verbose
   %figure(7); clf;
   %hist([stats.('Volume')], 256)

   figure(8); clf;
   set(gcf, 'Name', 'Cells from Watershed')
   imgsurf = impixelsurface(imgWs);
   implottiling({imoverlaylabel(imgsStRaw{5}, imgsurf, false); imoverlaylabel(imgC, imgsurf, false)}')
   
end
   
imgS = imgWs;      


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filter Cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mode = 'MeanIntensity';

stats = imstatistics(imgS, {mode, 'Volume'}, immask(imgsSt{5}, imgmask));

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

[imgSP, stats] = postProcessSegments(imgS, imgsStRaw{5}, 'intensity.mean.min', 0.1, 'volume.min', 10, 'volume.max', 50000, 'fillholes', false);

if verbose
   imgsurf = impixelsurface(imgSP);
   figure(5); clf;
   implottiling({imoverlaylabel(imgC, imgsurf, false), imoverlaylabel(imgsSt{5}, imgsurf, false)});
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
    implottiling({imoverlaylabel(imgsSt{5}, imglabc > 0, false, 'color.map', [[0,0,0]; [1,0,0]]), 
                  imoverlaylabel(imgC, imglabc > 0, false, 'color.map', [[0,0,0]; [1,0,0]])}');
              
    %saveas(h, [datadir, dataname, datafield, '_CellDetection.pdf'])
end

%% Cell Outline

imglabcSurf = impixelsurface(imglabc);

if verbose
    h = figure(78); clf;
    set(gcf, 'Name', 'Cell Centers')
    implottiling({imoverlaylabel(imgsSt{5}, imglabcSurf > 0, false, 'color.map', [[0,0,0]; [1,0,0]]), 
                  imoverlaylabel(imgC, imglabcSurf > 0, false, 'color.map', [[0,0,0]; [1,0,0]])}');
              
    %saveas(h, [datadir, dataname, datafield, '_CellDetection.pdf'])
end

statsSurf = imstatistics(imglabcSurf, {'PixelIdxList'});


%% Save Image Preporcessing Result

savefile = [datadir, dataname, datafield, '_ImageProcessing.mat'];
%save(savefile)
load(savefile)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Statistics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stats = imstatistics(imglabc, {'PixelIdxList', 'Centroid'});

%% Intensities

%mode = 'MedianIntensity';
mode = 'MeanIntensity';

clear statsCh
for c = 1:nch
   statsChRaw{c} = imstatistics(imglabc, mode,  imgsStRaw{c});
   statsCh{c} = imstatistics(imglabc, stats, mode,  imgsSt{c});
end

clear statsChN
for c = 1:nch
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

imgsStMeasure = {imgsSt{1}, imgsStRaw{2}, imgsStRaw{3}, imgsSt{4}, imgsStRaw{5}};
statsMeasure = {statsCh{1}, statsChRaw{2}, statsChRaw{3}, statsCh{4}, statsChRaw{5}};


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
       
       imsubplot(1,4,c)
       %scatter(xy(:,1), xy(:,2), 30, cdat, 'filled');
       colormap jet
       scatter(xy(:,1), xy(:,2), 30, fi, 'filled');
       title(chlabel{c});
    end
end


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
   title(['M ' chlabel{c}]);
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
   title(['R ' chlabel{c}]);
   
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
   p95 = int64(1 * ncs);
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
   scattercloud(fi{pairs{n}(1)}, fi{pairs{n}(2)}, 100, .5, 'r+', [jet(256); repmat([0.5,0,0],600,1)] );
   %scattercloud(fi{pairs{n}(1)}, fi{pairs{n}(2)}, 75, 1, 'r+', [flipud(gray(256)); repmat([0,0,0], 400,1)] );
   %xlim([0,1]); ylim([0,1]);
   xlabel(chlabel{pairs{n}(1)}); ylabel(chlabel{pairs{n}(2)});
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
      xlabel(chlabel{pairs{n}(1)}); ylabel(chlabel{pairs{n}(2)});

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
   title(['Raw ' chlabel{c}])
   
   subplot(4, nch, c+nch);
   dat = [statsChRawN{c}.(mode)];
   hist(dat(isgood), 256)
   title(['Raw Norm. ' chlabel{c}])
   
   subplot(4, nch, c+2*nch); 
   dat = [statsCh{c}.(mode)];
   hist(dat(isgood), 256)
   title(['Corr ' chlabel{c}])
   
   subplot(4, nch, c+3*nch); 
   dat = [statsChN{c}.(mode)];
   hist(dat(isgood), 256)
   title(['Corr Norm.' chlabel{c}])
   
end

%% Classify

%cth = {0.16, 0.35, 0.18, 0.18, 0}; % measured
cth = {0.16, 0.25, 0.18, 0.18, 0}; % measured
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

       R = imgsStRaw{c}; G = R; B = R;
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
    implottiling({imgClassRaw{:}; imgClass1{:}}, 'titles', [chlabel(1:nch)', chlabel(1:nch)']');
    set(gcf, 'Name', 'Classifications');

    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 30];
    fig.PaperPositionMode = 'manual';
    set(gcf, 'papersize', [8,30]);

    %saveas(h, [datadir dataname datafield '_Classification' '.pdf']);
end



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

    R = imgsStRaw{c}; G = R; B = R;
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
    implottiling({imgsStRaw{1:nch}, imgC; imgClass1{:}, imgClass}, 'titles', [[chlabel(1:nch)'; {'merge'}], [chlabel(1:nch)'; {'merge'}]]');

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
    nurr1xy = xy(nurr1ids,:);
    nurr1NeuroClassTotal = neuroClassTotal(nurr1ids);
    
    ord = [9,10,11,12,13,14,15,16];
    h = figure(66); clf; 
    hold on
    k = 1;
    for i = ord
        ids = nurr1NeuroClassTotal + 1 == i;
        scatter(nurr1xy(ids,1), nurr1xy(ids,2), 10, neuroNurrColor{k});
        k = k + 1;
    end
end


%% Class in Space

if verbose
    
    isgood = logical(isgood);

    xy = [stats.Centroid]';
    xy = xy(isgood, :)

    h = figure(27)

    cdat = (cell2mat({[0,0,0], neuroClassColor{:}}'))

    scatter(xy(:,1), xy(:,2), 5, cdat(neuroClassTotal(isgood) + 2, :), 'filled');
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
         clslab{i} = [clslab{i}, '_', chlabel{c}];
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

% writetable(tb, [datadir dataname datafield '_Counts.txt'])

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



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zones and distances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
figure(101); clf
implot(imgC)
%
%figure(100); clf
%implot(permute(imgC,[2,1,3]));

%%
p  = impoly

%%
outer = [200.287974683544 2873.68037974683;167.534810126582 2749.94620253164;167.534810126582 2688.07911392405;145.699367088607 2604.37658227848;145.699367088607 2527.95253164557;142.060126582278 2469.7246835443;127.503164556962 2411.49683544304;91.1107594936707 2324.15506329114;83.8322784810125 2211.33860759494;83.8322784810125 2083.96518987342;69.2753164556962 1938.39556962025;76.5537974683543 1800.10443037975;47.4398734177214 1567.19303797468;32.882911392405 1403.42721518987;25.6044303797466 1185.07278481013;14.6867088607594 1035.86392405063;11.0474683544303 875.737341772152;-7.14873417721537 457.224683544304;-7.14873417721537 -19.5158227848101];
inner = [4498.23101265823 2851.84493670886;4450.92088607595 2702.63607594937;4407.25 2560.70569620253;4389.05379746835 2455.16772151899;4356.30063291139 2324.15506329114;4301.71202531645 2069.4082278481;4283.51582278481 1934.75632911392;4250.76265822785 1760.07278481013;4221.64873417721 1643.61708860759;4207.0917721519 1519.8829113924;4170.69936708861 1385.23101265823;4163.42088607595 1268.7753164557;4141.58544303797 1086.81329113924;4119.75 1014.02848101266;4134.30696202532 904.851265822785;4119.75 788.395569620253;4116.11075949367 668.300632911392;4119.75 489.977848101266;4097.91455696202 391.71835443038;4105.19303797468 289.819620253164;4116.11075949367 180.642405063291;4105.19303797468 64.1867088607596;4112.47151898734 2.31962025316489];
ISVZ_OSVZ = [4061.52215189873 2877.31962025316;4043.32594936709 2717.19303797468;4017.85126582278 2578.90189873418;3992.37658227848 2426.05379746835;3966.90189873418 2287.76265822785;3930.50949367089 2164.02848101266;3890.47784810126 2000.26265822785;3872.28164556962 1905.64240506329;3846.80696202531 1781.9082278481;3843.16772151899 1632.69936708861;3792.21835443038 1501.68670886076;3766.74367088607 1370.67405063291;3744.9082278481 1225.10443037975;3741.26898734177 1061.33860759494;3737.62974683544 923.04746835443;3748.54746835443 752.003164556962;3763.10443037975 653.743670886076;3766.74367088607 544.566455696202;3737.62974683544 377.161392405063;3723.07278481013 217.034810126582;3704.87658227848 136.971518987342;3737.62974683544 -19.5158227848101];
SP_CP = [1055.50949367089 2866.40189873418;1015.47784810127 2724.47151898734;971.806962025316 2535.23101265823;939.053797468354 2305.95886075949;920.857594936709 2123.99683544304;920.857594936709 1898.36392405063;866.268987341772 1578.11075949367;858.990506329113 1327.00316455696;840.794303797468 1003.11075949367;844.433544303797 806.591772151899;858.990506329113 624.629746835443;851.712025316455 439.028481012658;877.186708860759 202.477848101266;858.990506329113 56.9082278481014;858.990506329113 -4.95886075949375];
CP_M = [371.332278481013 2877.31962025316;338.57911392405 2757.2246835443;320.382911392405 2640.76898734177;280.351265822785 2473.36392405063;273.072784810126 2364.18670886076;251.237341772152 2196.78164556962;229.401898734177 2065.76898734177;218.48417721519 1956.5917721519;233.041139240506 1847.41455696203;283.990506329114 1727.31962025316;287.629746835443 1636.33860759494;276.712025316456 1530.80063291139;262.155063291139 1385.23101265823;251.237341772152 1246.93987341772;214.844936708861 1050.42088607595;207.566455696202 930.325949367088;196.648734177215 788.395569620253;160.256329113924 646.465189873418;131.142405063291 570.041139240506;152.977848101266 471.78164556962;145.699367088607 340.768987341772;160.256329113924 220.674050632912;167.534810126582 93.3006329113923;163.895569620253 -8.59810126582261];

VZ_ISVZ = [4327.437890625 2869.17734375;4281.724609375 2725.0046875;4243.044140625 2598.4140625;4193.814453125 2405.01171875;4155.133984375 2208.09296875;4105.904296875 1990.07578125;4077.773046875 1817.771875;4021.510546875 1641.9515625;4003.928515625 1511.84453125;3993.379296875 1413.38515625;3979.313671875 1283.278125;3961.731640625 1065.2609375;3954.698828125 840.2109375;3951.182421875 678.45625;3951.182421875 460.4390625;3940.633203125 326.815625;3947.666015625 133.41328125;3926.567578125 -10.7593750000005];
IZ_SP = [2017.158984375 2883.24296875;1996.060546875 2788.3;1887.051953125 2510.50390625;1823.756640625 2236.22421875;1774.526953125 1859.96875;1742.879296875 1532.94296875;1721.780859375 1170.753125;1707.715234375 833.178125;1711.231640625 502.6359375;1728.813671875 263.5203125;1721.780859375 52.5359374999994;1714.748046875 -17.7921875000006];
OSVZ_IZ = [2407.480078125 2879.7265625;2386.381640625 2732.0375;2326.602734375 2415.5609375;2273.856640625 2116.66640625;2242.208984375 1845.903125;2196.495703125 1504.81171875;2200.012109375 1272.72890625;2200.012109375 910.5390625;2210.561328125 622.19375;2185.946484375 383.078125;2178.913671875 150.9953125;2175.397265625 17.3718749999994;2150.782421875 -45.9234375000005];

lines = {inner, VZ_ISVZ, ISVZ_OSVZ, OSVZ_IZ, IZ_SP, SP_CP, CP_M, outer};

nlines = length(lines);
lineshires = cell(1, nlines);
for l = 1:nlines
   lineshires{l} = [interp1(lines{l}(:,2), lines{l}(:,1), 1:size(imgC, 2)); (1:size(imgC,2))];
end
   
if verbose
   h = figure(101); clf
   implot(imgC)
   hold on
   for l = 1:nlines
      plot(lineshires{l}(1,:), lineshires{l}(2,:), 'r', 'LineWidth',1)
   end
   hold off
   
   %saveas(h, [datadir dataname datafield '_Annotation' '.pdf']);
end


%% Calculate Relative Cell Position

xy = [stats.Centroid];

ncells = size(xy,2);

dinner = zeros(1,ncells);
douter = zeros(1,ncells);
for i = 1:ncells
   if mod(i, 500) == 0
      fprintf('%d / %d\n', i, ncells)
   end
   dinner(i) = minimalDistance(xy(:,i), lineshires{1});
   douter(i) = minimalDistance(xy(:,i), lineshires{end});
end


%% Plot Expressions as relative distances

drel = dinner ./ (dinner + douter);

% non normalized
if verbose
   h  = figure(200); clf;
   for c = 1:nch-1
      dat = [statsCh{c}.(mode)];
      subplot(1,nch-1,c)
      plot(drel, dat, ['.' cm{c}])
      title(chlabel{c})
      xlabel('rel. position')
      ylabel('intensity [AU]')
   end
   
   fig = gcf;
   fig.PaperUnits = 'inches';
   fig.PaperPosition = [0 0 20 8];
   fig.PaperPositionMode = 'manual';
   set(gcf, 'papersize', [20,8]);
   
   saveas(h, [datadir dataname datafield '_RelativePosition_vs_Expression' '.pdf']);
end

%%
if verbose
   h  = figure(200); clf;
   for c = 1:nch-1
      dat = [statsChN{c}.(mode)];
      subplot(1,nch-1,c)
      plot(drel, dat, ['.' cm{c}])
      title(chlabel{c})
      xlabel('rel. position')
      ylabel('intensity [AU]')
   end
   
   fig = gcf;
   fig.PaperUnits = 'inches';
   fig.PaperPosition = [0 0 20 8];
   fig.PaperPositionMode = 'manual';
   set(gcf, 'papersize', [20,8]);
   
   saveas(h, [datadir dataname datafield '_RelativePosition_vs_Expression_DAPINormalized' '.pdf']);
end




%% Save the Statistics to File

statsfile = [datadir dataname datafield '_statistics' '.mat'];
save(statsfile, 'chlabel', 'stats', 'statsCh', 'statsChN',  'neuroClass', 'lineshires', 'dinner', 'douter');
%load(statsfile)

%% Save For Mathematica

centers = [stats.Centroid];

intensities = zeros(length(statsCh), size(centers,2));
intensities_normalized = zeros(length(statsCh), size(centers,2));
for c = 1:length(statsCh)
   intensities(c, :) = [statsCh{c}.MedianIntensity];
   intensities_normalized(c, :) = [statsChN{c}.MedianIntensity];
end


statsfilemathematica = [datadir dataname datafield '_statistics_mathematica' '.mat'];
save(statsfilemathematica, '-v6', 'chlabel', 'centers', 'dinner', 'douter', 'lineshires', 'intensities', 'intensities_normalized', 'neuroClass');