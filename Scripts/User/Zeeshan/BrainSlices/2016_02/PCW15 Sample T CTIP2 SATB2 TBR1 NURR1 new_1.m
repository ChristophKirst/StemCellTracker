%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Brain Slices %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

initialize()
bfinitialize
initializeParallelProcessing(12)

addpath('./Scripts/User/Zeeshan/BrainSlices/2015_11');

datadir = '/home/ckirst/Data/Science/Projects/StemCells/Experiment/BrainSections/';

savedir   = '/home/ckirst/Data/Science/Projects/StemCells/Experiment/BrainSections/Analysis/2016_02/Data/'
resultdir = '/home/ckirst/Data/Science/Projects/StemCells/Experiment/BrainSections/Analysis/2016_02/';
dataname = 'PCW15 Sample T CTIP2 TBR1 SATB2 NURR1 new'

datafield = ' 1';

verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Source %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

is = ImageSourceBF([datadir  dataname '.czi']);
clc
is.printInfo

%%

is.setReshape('S', 'Uv', [13, 19]);
is.setCellFormat('UV')
is.setRange('C', 1);
clc; is.printInfo

%%

if verbose
    is.setRange('C', 1);
    is.plotPreviewStiched('overlap', 102, 'scale', 0.1);
end

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
shift = alignImages(imgs, 'alignment', 'RMS', 'overlap.max', 150);
stitchmethod = 'Interpolate';
%stitchmethod = 'Pyramid';


%%
if verbose 
   img = stitchImages(imgs, shift, 'method', stitchmethod);
   figure(1); clf; implot(img)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stich Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

is.resetRange();

nch = 5;
for i = 1:5
    imgsRaw{i}  = is.cell(region, 'C', i);
end

imgsRaw = {imgsRaw{2}, imgsRaw{1}, imgsRaw{4}, imgsRaw{3}, imgsRaw{5}};
imgsRaw = {imgsRaw{2}, imgsRaw{1}, imgsRaw{4}, imgsRaw{3}, imgsRaw{5}};

channellabel  = {'CTIP2', 'TBR1', 'SATB2', 'NURR1', 'DAPI'};
channelcolors = {[1,0,0], [0,1,0], [1,1,0], [1, 0, 1], [0,0,1]};

%% Preprocess

imgs = cell(1,nch);
parfor i = 1:nch
   imgs{i} = imgsRaw{i};
   %thrs = 0 * [0, 720, 430, 0, 600];
   %imgs{i} = cellfunc(@(x) (x > thrs(i)) .* x,  imgs{i});
   %imgs{i} = cellfunc(@(x) x - imopen(x, strel('disk', 100)), imgs{i});   
   %imgs{i} = cellfunc(@(x) x - imopen(x, strel('disk', 10)), imgs{i});   
   %imgs{i} = cellfunc(@(x) filterAutoContrast(x/max(x(:))), imgs{i});
   imgs{i} = cellfunc(@(x) filterMedian(x, 3), imgs{i});
   imgs{i} = cellfunc(@(x) x - imopen(x, strel('disk', 75)), imgs{i});   
   imgs{i} = cellfunc(@(x) x - imopen(x, strel('disk', 10)), imgs{i});   
end

% figure(12); clf;
% implottiling({imgsRaw{i}{1:2}; imgs{i}{1:2}})

%% Stitch

img = cell(1,nch);
imgRaw = cell(1,nch);

for i = 1:nch
   imgRaw{i} = stitchImages(imgsRaw{i}, shift, 'method', stitchmethod);
   imgRaw{i} = mat2gray(imgRaw{i});
   
   img{i} = stitchImages(imgs{i}, shift, 'method', stitchmethod);
   
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

%% Restrict to sub Region

figure(99); clf;
implot(imgC)

%%

sub = true;
%subregion = {50:4900, 1:500};
subregion = {50:4900, ':'};


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

%imgm =  imgsSt{5} - imopen(imgsSt{5} , strel('disk', 10));
imgm = img{5};

imgmaskLo = imgm > 0.08;
imgmaskLo = imopen(imgmaskLo, strel('disk', 2));
imgmaskLo = postProcessSegments(bwlabeln(imgmaskLo), 'volume.min', 15, 'fillholes', false) > 0;

%figure(66); clf;
%implottiling({255*(postProcessSegments(bwlabeln(imgmaskLo), 'volume.min', 10, 'fillholes', false)>0); 2*(bwlabeln(imgmaskLo)); 255*imgmaskLo}')

%get rid of blod vessels
imgmaskHi1 = imgRaw{4} > 0.5;
imgmaskHi1 = imdilate(imgmaskHi1, strel('disk', 3));
imgmaskHi1 = postProcessSegments(bwlabeln(imgmaskHi1), 'volume.min', 500, 'fillholes', false) > 0;
imgmaskHi1 = not(imgmaskHi1);

imgmaskHi2 = imgRaw{1} > 1.9;
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
   implottiling({imgRaw{5},imgmaskLo; 
                 imgRaw{4},imgmaskHi1;
                 imgRaw{1},imgmaskHi2;
                 mat2gray(imgm),  imgmask
                 }, 'titles', {'I', 'Lo', 'I','Hi','I', 'Mask'})
end

clear('imgmaskHi1', 'imgmaskHi2', 'imgmaskHi', 'imgmaskLo');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cell Centers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%imgBM = filterBM(mat2gray(imgC), 'profile', 'np', 'sigma', 9);
% if verbose > 1
%    figure(3); clf;
%    set(gcf, 'Name', 'BM')
%    implottiling({imgC; mat2gray(sum(imgBM,3)); mat2gray(imgI)}')
% end

%imgCf = imgBM;
%imgI  = sum(imgBM,3);
imgI = imgC(:,:,3);


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
imgWs = imimposemin(iminvert(immask(imgRaw{5}, imgmask)), imgmax);
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
   implottiling({imoverlaylabel(imgRaw{5}, imgsurf, false); imoverlaylabel(imgC, imgsurf, false)}')
   
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

[imgSP, stats] = postProcessSegments(imgS, imgRaw{5}, 'intensity.mean.min', 0.1, 'volume.min', 10, 'volume.max', 50000, 'fillholes', false);

if verbose
   imgsurf = impixelsurface(imgSP);
   figure(5); clf;
   implottiling({imoverlaylabel(imgC, imgsurf, false), imoverlaylabel(img{5}, imgsurf, false)});
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
    fig = figure(78); clf;
    implottiling({imoverlaylabel(img{5}, imglabc > 0, false, 'color.map', [[0,0,0]; [1,0,0]]), 
                  imoverlaylabel(imgC, imglabc > 0, false, 'color.map', [[0,0,0]; [1,0,0]])}');
              
    %saveas(h, [datadir, dataname, datafield, '_CellDetection.pdf'])
end

%% Cell Outline

imglabcSurf = impixelsurface(imglabc);

if verbose
    fig = figure(78); clf;
    set(gcf, 'Name', 'Cell Centers')
    implottiling({imoverlaylabel(img{5}, imglabcSurf > 0, false, 'color.map', [[0,0,0]; [1,0,0]]), 
                  imoverlaylabel(imgC, imglabcSurf > 0, false, 'color.map', [[0,0,0]; [1,0,0]])}');
              
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

% image sets to measure on -> some might have bad background -> used background removed image for those

%imgMeasure = {img{1}, img{2}, img{3}, img{4}, imgRaw{5}};
%statsMeasure = {statsCh{1}, statsCh{2}, statsCh{3}, statsCh{4}, statsChRaw{5}};

imgMeasure = {imgRaw{1}, img{2}, img{3}, img{4}, imgRaw{5}};
statsMeasure = {statsChRaw{1}, statsCh{2}, statsCh{3}, statsCh{4}, statsChRaw{5}};


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
       scatter(xy(:,1), xy(:,2), 2, fi, 'filled');
       title(channellabel{c});
       axis off
    end 
    
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 20  20/ 1.5656 * 4];
    fig.PaperPositionMode = 'manual';
    set(gcf, 'papersize', [20, 20/ 1.5656 * 4]);

    saveas(fig, [resultdir dataname datafield '_Intensities_Space.png']);  
end


%% Expression in Space - Separate Channels

xy = [stats.Centroid]';

if verbose > 1
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
end

%% Flourescence Expression

fi = cell(1,3);
for c = 1:nch
   fi{c} = [statsMeasure{c}.(mode)]';
   fis = sort(fi{c});
   ncs = length(fis);
   p95 = int64(0.95 * ncs);
   fi{c} = imclip(fi{c}, 0, fis(p95));
   %fi{c} = mat2gray(fi{c});
end


if verbose
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
   set(gcf, 'Name', 'Normalized Expression');
   
   saveas(gcf, [resultdir dataname datafield '_Intensities_Cloud.png']);
  
%    for n = 1:np
%       h = figure(50+n); clf;
%       %fi = imclip(fi, 0, cl(c));
%       %scatter(fi{pairs{n}(1)}, fi{pairs{n}(2)}, 10, 'b', 'filled');
%       scattercloud(fi{pairs{n}(1)}, fi{pairs{n}(2)});
%       xlim([0,1]); ylim([0,1]);
%       xlabel(channellabel{pairs{n}(1)}); ylabel(channellabel{pairs{n}(2)});
% 
%       saveas(h, [resultdir dataname datafield '_Intensities_Cloud' chlabel{pairs{n}(1)} ,'_' chlabel{pairs{n}(2)} '.png']);
%       %freezecolormap(gca)
%    end
end

%% Save Measurement Results

close all
clear('cdat', 'fis', 'ind', 'fi', 'nb', 'cdat', 'bh', 'fig', 'p90', 'p95', 'np', 'ncs', 'ncells', 'n', 'fith', 'clip', 'c', 'ha', 'h', 'cols', 'pairs', 'xb', 'xy')

%%
savefile = [savedir, dataname, datafield, '_Statistics.mat'];
save(savefile)
%load(savefile)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Zones and Distances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
figure(101); clf
implot(imgC)
%
%figure(100); clf
%implot(permute(imgC,[2,1,3]));

%%
%p  = impoly

%%
line_outer = [362.1921875 3104.953125;237.6359375 2664.215625;156.1953125 2276.175;117.8703125 1916.878125;89.1265625 1428.234375;65.1734375 1040.19375;74.7546875 700.059375;84.3359375 345.553125;117.8703125 0.628124999999727;132.2421875 -32.9062500000005]
line_inner = [4898.9140625 3100.1625;4740.8234375 2611.51875;4597.1046875 2003.109375;4558.7796875 1672.55625;4515.6640625 1212.65625;4482.1296875 757.546875;4472.5484375 484.48125;4530.0359375 29.3718749999998;4530.0359375 0.628124999999727]
line_ISVZ_OSVZ = [4482.1296875 3109.74375;4347.9921875 2760.028125;4237.8078125 2314.5;4170.7390625 1959.99375;4118.0421875 1562.371875;4098.8796875 1164.75;4094.0890625 656.94375;4141.9953125 158.71875;4204.2734375 -23.3250000000003];
line_SP_CP = [1143.0640625 3085.790625;1032.8796875 2616.309375;932.2765625 2156.409375;893.9515625 1682.1375;850.8359375 1370.746875;812.5109375 987.496875;788.5578125 618.61875;759.8140625 268.903125;774.1859375 -28.1156250000004]
line_CP_M = [525.0734375 3100.1625;453.2140625 2875.003125;424.4703125 2721.703125;410.0984375 2544.45;304.7046875 2242.640625;261.5890625 1936.040625;247.2171875 1653.39375;228.0546875 1274.934375;194.5203125 973.125;175.3578125 738.384375;180.1484375 565.921875;194.5203125 403.040625;218.4734375 197.04375;242.4265625 -37.6968750000001]
line_VZ_ISVZ = [4630.6390625 3085.790625;4482.1296875 2740.865625;4400.6890625 2343.24375;4309.6671875 1801.903125;4285.7140625 1562.371875;4252.1796875 1260.5625;4242.5984375 838.9875;4242.5984375 494.0625;4280.9234375 82.0687499999999;4304.8765625 -23.3250000000003]
line_IZ_SP = [2086.8171875 3119.325;1991.0046875 2621.1;1981.4234375 2405.521875;1952.6796875 2041.434375;1914.3546875 1739.625;1876.0296875 1365.95625;1871.2390625 968.334375;1895.1921875 609.0375;1890.4015625 393.459375;1890.4015625 177.88125;1909.5640625 10.2093749999999]
line_OSVZ_IZ = [2700.0171875 3133.696875;2604.2046875 2846.259375;2537.1359375 2539.659375;2503.6015625 2237.85;2455.6953125 1912.0875;2441.3234375 1567.1625;2436.5328125 1351.584375;2436.5328125 997.078125;2455.6953125 733.59375;2431.7421875 446.15625;2426.9515625 149.1375;2474.8578125 -28.1156250000004];

zonelines = {line_inner, line_VZ_ISVZ, line_ISVZ_OSVZ, line_OSVZ_IZ, line_IZ_SP, line_SP_CP, line_CP_M, line_outer};

if sub
  for i = 1:2
    if subregion{i} == ':'
        zoneshift(i) = 1;
    else
        zoneshift(i) = subregion{i}(1);
    end
  end
  zonelines = cellfunc(@(x) x - repmat(zoneshift, size(x,1),1), zonelines);  
end

zonenames = {'VZ', 'ISVZ', 'OSVZ', 'IZ', 'SP', 'CP', 'MZ'};
                   
zonecolors = num2cell([0.709804, 0.827451, 0.2;
                       0.929412, 0.290196, 0.596078;
                       0.584314, 0.752941, 0.705882;
                       0.984314, 0.690196, 0.25098;
                       0, 174/256, 205 / 256;
                       0.462745, 0.392157, 0.67451;
                       0.976471, 0.929412, 0.196078], 2);
                  
pixelsize = 0.83 * 0.83; % um^2;
pixelsize =  pixelsize / (1000 * 1000); % mm^2

nzonelines = length(zonelines);
zonelineshires = cell(1, nzonelines);

nzones = length(zonelines) - 1;

for l = 1:nzonelines
   zonelineshires{l} = [interp1(zonelines{l}(:,2), zonelines{l}(:,1), 1:size(imgC, 2)); (1:size(imgC,2))];
   
   % remove all nans and infs
   for i = 1:2
       ids = zonelineshires{l}(i,:) < inf;
       ids = and(ids, ~isnan(zonelineshires{l}(i,:)));
       %ids = and(ids, lineshires{l}(i,:) > 0);
       zonelineshires{l} =  zonelineshires{l}(:, ids);
   end
end
   
if verbose
   fig = figure(101); clf
   implot(imgC)
   hold on
   for l = 1:nzonelines
      plot(zonelineshires{l}(1,:), zonelineshires{l}(2,:), 'r', 'LineWidth',1)
   end
   hold off
   
   saveas(fig, [resultdir dataname datafield '_Annotation' '.png']);
end


%% Calculate Relative Cell Position

xy = [stats.Centroid];

ncells = size(xy,2);

distance_inner = zeros(1,ncells);
distance_outer = zeros(1,ncells);
parfor i = 1:ncells
   if mod(i, 500) == 0
      fprintf('%d / %d\n', i, ncells)
   end
   distance_inner(i) = minimalDistance(xy(:,i), zonelineshires{1});
   distance_outer(i) = minimalDistance(xy(:,i), zonelineshires{end});
end

%% Measure Densities in Zones

% Area of zones
area = zeros(nzones, 1);
for i = [1:length(zonelineshires)-1]
    pp = ROIPolygon([zonelineshires{i}, zonelineshires{i+1}(:, end:-1:1)]);
    area(i) = pp.volume;
end

if verbose
    figure(16); clf; hold on
    set(gcf, 'Name', 'Area of Zones');
    
    for i = 1:nzones
       bar(i, area(i), 0.75, 'FaceColor', zonecolors{i}); 
    end
    set(gca, 'XTickLabel', '')  
    title('Area of Zones')
    text(1:nzones,repmat(-max(ylim)/50, nzones, 1), zonenames','horizontalalignment','right','Rotation',35,'FontSize',15)
    
    saveas(gcf, [resultdir dataname datafield '_Zones_Areas.png']);
end

%% Cells in Zones

zonesid = cell(nzones,1);
zonesdensity = zeros(nzones,1);
%imgROIAll = zeros(size(imgI));

if verbose 
    figure(123); clf; hold on
    set(gcf, 'Name', 'Cells in Zones');
    for a = 1:nzones
        roi = [zonelineshires{a}, zonelineshires{a+1}(:, end:-1:1)];
        imgROI = poly2mask(roi(2,:), roi(1,:), size(imgI,1), size(imgI,2));

        roiind = imsub2ind(size(imgI), fix(xy)');
        roiind = imgROI(roiind);

        xyROI = fix(xy(:, roiind));
        scatter(xyROI(1,:), xyROI(2,:), 20, zonecolors{a}, 'filled');

        zonesid{a} = roiind;
        %zonesdensity(a) = sum(roiind) / area(a);
        zonesdensity(a) = sum(roiind) / area(a) / pixelsize;

        %imgROIAll = imgROIAll + a * imgROI;
    end
    saveas(gcf, [resultdir dataname datafield '_Zones_Space.png']);


    figure(16); clf; hold on
    set(gcf, 'Name', 'Cell Density in Zones [1/mm^2]');
    
    for i = 1:nzones
       bar(i, zonesdensity(i), 0.75, 'FaceColor', zonecolors{i}); 
    end
    
    set(gca, 'XTickLabel', '')  
    title('Cell Density in Zones [#/mm^2]')
    text(1:nzones,repmat(-max(ylim)/50, nzones, 1), zonenames','horizontalalignment','right','Rotation',35,'FontSize',15)
    
    saveas(fig, [resultdir dataname datafield '_Zones_Density.png']);
end

nc = num2cell(zonesdensity);
tb = table(nc{:}, 'VariableNames', zonenames)

writetable(tb, [resultdir dataname datafield '_Zones_Density.txt'])




%% Plot Expressions Vs Relative Distances

distance_rel = distance_inner ./ (distance_inner + distance_outer);

class_thresholds = {0.12, 0.17, 0.11, 0.25, 0};

% non normalized
if verbose
   fig  = figure(200); clf;
   for c = 1:nch-1
      dat = [statsMeasure{c}.(mode)];
      subplot(1,nch-1,c)
      plot(dat, distance_rel, '.', 'color', channelcolors{c}); hold on
      plot( [class_thresholds{c}, class_thresholds{c}],[0, 1], 'k')
      title(channellabel{c})
      ylabel('rel. position')
      xlabel('intensity [AU]')
   end
   
   fig = gcf;
   fig.PaperUnits = 'inches';
   fig.PaperPosition = [0 0 30 8];
   fig.PaperPositionMode = 'manual';
   set(gcf, 'papersize', [30,8]);
   
   saveas(fig, [resultdir dataname datafield '_Intensities_RelativePosition.png']);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Classify Cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Mask for ROIs

%h  = figure(300); clf;
%implot(img{4});
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

fig = figure(10); clf
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

fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 20  15];
fig.PaperPositionMode = 'manual';
set(gcf, 'papersize', [20, 15]);

saveas(gcf, [resultdir dataname datafield '_Intensities_Histograms.png']);

%% Classify

%class_thresholds = {0.16, 0.35, 0.18, 0.18, 0}; % measured
%class_thresholds = {0.16, 0.25, 0.16, 0.18, 0}; % measured

%classcolors = {[1,0,0], [0,1,0], [1,1,0], [1,0,1],[0,0,1]};
classcolors = channelcolors;

cellclass = zeros(nch, length(isgood));
for c = 1:nch
   cellclass(c,:) = double(and([statsMeasure{c}.(mode)] > class_thresholds{c}, isgood));
end
cellclassfull = int32(fix(cellclass(1,:) + 2 * cellclass(2,:) + 4 * cellclass(3,:) + 8 * cellclass(4,:)))

%nclassesfull = num2cell(colorcube(max(neuroClassTotal(:))), 2);
nclassesfull = 2^(nch-1);

for i = 1:nclassesfull
    id = num2cell(dec2bin(i-1));
    id = id(end:-1:1);
    for k = length(id)+1:nch
      id{k} = '0';
    end

    cellclassfulllabel{i} = '';
    for c = 1:nch
      if id{c} == '1'
         cellclassfulllabel{i} = [cellclassfulllabel{i}, '_', channellabel{c}];
      end
    end
    cellclassfulllabel{i} = cellclassfulllabel{i}(2:end);
    if isempty(cellclassfulllabel{i})
      cellclassfulllabel{i} = 'None';
    end
end



%%
classfullcolor = cell(1, nclassesfull);
for i = 1:nclassesfull
   cb = dec2bin(i-1);
   cb = cb(end:-1:1);
   %cb = [cb(end:-1,1), '0000'];
   %cb = cb(1:nch);
   col = [0,0,0]; k = 0;
   for l = 1:length(cb)
      if cb(l) == '1'
         col = col + classcolors{l};
         k = k + 1;
      end
   end
   if k > 1
      col = col / k;
   end
   classfullcolor{i} = col;
end

if verbose
    %imgClass1 = cell(1,nch);
    imgClass1 = cell(1,2);
    for c = 1:nch
       ccol = {[0,0,0], classcolors{c}};
       ccol = ccol(cellclass(c,:)+1);

       R = imgMeasure{c}; G = R; B = R;
       for i = 1:length(stats);
          if cellclass(c,i) > 0
             R(statsSurf(i).PixelIdxList) =  ccol{i}(1);
             G(statsSurf(i).PixelIdxList) =  ccol{i}(2);
             B(statsSurf(i).PixelIdxList) =  ccol{i}(3);
          else
             R(statsSurf(i).PixelIdxList) =  classcolors{5}(1);
             G(statsSurf(i).PixelIdxList) =  classcolors{5}(2);
             B(statsSurf(i).PixelIdxList) =  classcolors{5}(3);
          end
       end
       imgClass1{c} = cat(3, R, G, B);

       R = imgRaw{c}; G = R; B = R;
       for i = 1:length(stats);
          if cellclass(c,i) > 0
             R(statsSurf(i).PixelIdxList) =  ccol{i}(1);
             G(statsSurf(i).PixelIdxList) =  ccol{i}(2);
             B(statsSurf(i).PixelIdxList) =  ccol{i}(3);
          else
             R(statsSurf(i).PixelIdxList) =  classcolors{5}(1);
             G(statsSurf(i).PixelIdxList) =  classcolors{5}(2);
             B(statsSurf(i).PixelIdxList) =  classcolors{5}(3);
          end
       end
       imgClassRaw{c} = cat(3, R, G, B);
    end


    fig = figure(7); clf;
    xlabel('')
    implottiling({imgClassRaw{:}; imgClass1{:}}, 'titles', [channellabel(1:nch)', channellabel(1:nch)']');
    set(gcf, 'Name', 'Classification');

    figfac = nch * ny / (2 * nx) ; 
    
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 10 figfac*10 + 0.5 * nch];
    fig.PaperPositionMode = 'manual';
    set(gcf, 'papersize', [10, figfac*10 + 0.5 * nch]);

    saveas(fig, [resultdir dataname datafield '_Classification.png']);
end


%% Image + Classes + Combined Classes

if verbose > 1
    R = imgRaw{c}; G = R; B = R;
    for p = 1:length(stats);
       nc = classfullcolor{cellclassfull(p)+1};

       R(stats(p).PixelIdxList) =  nc(1);
       G(stats(p).PixelIdxList) =  nc(2);
       B(stats(p).PixelIdxList) =  nc(3);
    end

    imgClass = cat(3, R, G, B);

    fig = figure(7); clf;
    implottiling({imgRaw{1:nch}, imgC; imgClass1{:}, imgClass}, 'titles', [[channellabel(1:nch)'; {'merge'}], [channellabel(1:nch)'; {'merge'}]]');

    %saveas(h, [resultdir dataname datafield '_Classification_Channels' '.pdf']);
end


%% Class in Space

if verbose > 1
    isgood = logical(isgood);

    xy1 = [stats.Centroid]';
    xy1 = xy1(isgood, :)

    fig = figure(27)

    cdat = (cell2mat({[0,0,0], classfullcolor{:}}'))

    scatter(xy1(:,1), xy1(:,2), 20, cdat(cellclassfull(isgood) + 2, :), 'filled');
    xlim([0, size(imglab,1)]); ylim([0, size(imglab,2)]);
    %title(chlabel{c});
    %freezecolormap(gca)

    pbaspect([1,1,1])

    %saveas(h, [resultdir dataname datafield '_Quantification_Space.pdf']);
end


%% Histogram of Cell Classes
 
nclassesfull = length(classfullcolor);

nstat = zeros(1,nclassesfull);
for i = 1:nclassesfull   
   nstat(i) = sum(cellclassfull + 1 == i);
end

nc = num2cell(nstat);
tb = table(nc{:}, 'VariableNames', cellclassfulllabel)

writetable(tb, [resultdir dataname datafield '_Classification_Counts_All.txt'])


fig = figure(10); clf; hold on
k = 1;
%ord = [1, 5, 3, 6, 7, 8, 4, 2];
ord = 1:nclassesfull
for i = ord
   bar(k, nstat(i), 'FaceColor', classfullcolor{i});
   
   k = k + 1;
end

set(gca, 'XTickLabel', '')  
xlabetxt = strrep(cellclassfulllabel(ord), '_', '+');
n = length(cellclassfulllabel);
ypos = -max(ylim)/50;
text(1:n,repmat(ypos,n,1), xlabetxt','horizontalalignment','right','Rotation',35,'FontSize',15)

saveas(fig, [resultdir dataname datafield '_Classification_Counts_All.png']);


%% Measure Classes in Zones

zonecounts = zeros(nzones, nclassesfull);

fig = figure(10); clf;
for z = 1:nzones;
    for i = 1:nclassesfull
       zonecounts(z, i) = sum(cellclassfull(zonesid{z}) + 1 == i);
    end

    subplot(1,nzones,z); hold on

    k = 1;
    %ord = [1, 5, 3, 6, 7, 8, 4, 2];
    ord = 1:nclassesfull;
    for i = ord
       bar(k, zonecounts(z, i), 'FaceColor', classfullcolor{i});
       k = k + 1;
    end

    set(gca, 'XTickLabel', '')  
    xlabetxt = strrep(cellclassfulllabel(ord), '_', '+');
    n = length(cellclassfulllabel);
    ypos = -max(ylim)/50;
    text(1:n,repmat(ypos,n,1), xlabetxt','horizontalalignment','right','Rotation', 35, 'FontSize', 6)
    title(zonenames{z}) 
end

nc = num2cell(zonecounts,1);
tb = table(nc{:}, 'VariableNames', cellclassfulllabel, 'RowNames', zonenames)

writetable(tb, [resultdir dataname datafield '_Classification_Counts_Zone.csv'], 'WriteRowNames', true)

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition =  [0 0  35  7 ];
fig.PaperPositionMode = 'manual';
set(gcf, 'papersize',  [35 , 7]);
saveas(fig, [resultdir dataname datafield '_Classification_Counts_Zones.png']);

%% Measure Classes in Zones with Error Bars

nregions = 3;
regionbounds = (0:nregions) / nregions * ny;
for i = 1:nregions
    regionid{i} = and(regionbounds(i) <= xy(2,:), xy(2,:) < regionbounds(i+1))';
end

regioncounts = zeros(nzones, nregions, nclassesfull);

regionnames = cell(nregions, nzones);
for r = 1:nregions
  for z = 1:nzones
    for i = 1:nclassesfull
       regioncounts(z, r, i) = sum(cellclassfull(and(regionid{r}, zonesid{z})) + 1 == i);
    end
    regionnames{r,z} = [zonenames{z} '_' num2str(r)]
  end
end

rg = permute(regioncounts, [3,2,1]);
rg = rg(:,:); rg = rg';
rg = num2cell(rg, 1);
tb = table(rg{:}, 'VariableNames', cellclassfulllabel, 'RowNames', {regionnames{:}})

writetable(tb, [resultdir dataname datafield '_Classification_Counts_Regions.csv'], 'WriteRowNames', true)


% mean and std
regioncountsmean = squeeze(mean(regioncounts, 2))
regioncountsstd = squeeze(std(regioncounts,[], 2))

rg = num2cell(regioncountsmean, 1);
tb = table(rg{:}, 'VariableNames', cellclassfulllabel, 'RowNames', zonenames)
writetable(tb, [resultdir dataname datafield '_Classification_Counts_Mean.csv'], 'WriteRowNames', true)

rg = num2cell(regioncountsstd, 1);
tb = table(rg{:}, 'VariableNames', cellclassfulllabel, 'RowNames', zonenames)
writetable(tb, [resultdir dataname datafield '_Classification_Counts_Std.csv'], 'WriteRowNames', true)


%% mean and std / normalized

for r = 1:nregions
    regioncountsN(:,r,:) = regioncounts(:,r,:) / total(regioncounts(:,r,:));
end

regioncountsNmean = squeeze(mean(regioncountsN, 2))
regioncountsNstd = squeeze(std(regioncountsN,[], 2))

rg = num2cell(regioncountsNmean, 1);
tb = table(rg{:}, 'VariableNames', cellclassfulllabel, 'RowNames', zonenames)
writetable(tb, [resultdir dataname datafield '_Classification_Counts_Mean_Normalized.csv'], 'WriteRowNames', true)

rg = num2cell(regioncountsNstd, 1);
tb = table(rg{:}, 'VariableNames', cellclassfulllabel, 'RowNames', zonenames)
writetable(tb, [resultdir dataname datafield '_Classification_Counts_Std_Normalized.csv'], 'WriteRowNames', true)


%%
fig = figure(17); clf; 
for z = 1:nzones
    subplot(1,nzones,z); hold on
    
    k = 1;
    %ord = [1, 5, 3, 6, 7, 8, 4, 2];
    ord = 1:nclassesfull;
    for i = ord
       hb = bar(k, regioncountsmean(z,i), 'FaceColor', classfullcolor{i});
       xData = hb.XData+hb.XOffset;
       elo =regioncountsmean(z,i) - regioncountsstd(z,i);
       if elo < 0; elo = regioncountsmean(z,i); else elo = regioncountsstd(z,i); end
       errorbar(xData,regioncountsmean(z,i),elo, regioncountsstd(z,i), 'k.') 
       k = k + 1;
    end

    set(gca, 'XTickLabel', '')  
    xlabetxt = strrep(cellclassfulllabel(ord), '_', '+');
    n = length(cellclassfulllabel);
    ypos = -max(ylim)/50;
    text(1:n,repmat(ypos,n,1), xlabetxt','horizontalalignment','right','Rotation', 35, 'FontSize', 6)
    title(zonenames{z})
end

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition =  [0 0  25  8 ];
fig.PaperPositionMode = 'manual';
set(gcf, 'papersize',  [25 , 8]);
saveas(fig, [resultdir dataname datafield '_Classification_Counts_Zones_Errors_Normalized.png']);


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

    nurr1ids = find(cellclass(4,:));
    %nurr1xy = xy(nurr1ids,:);
    nurr1xy = xy(:, nurr1ids);
    nurr1NeuroClassTotal = cellclassfull(nurr1ids);
    
    ord = [9,10,11,12,13,14,15,16];
    fig = figure(66); clf; 
    hold on
    k = 1;
    for i = ord
        ids = nurr1NeuroClassTotal + 1 == i;
        %scatter(nurr1xy(ids,1), nurr1xy(ids,2), 20, neuroNurrColor{k}, 'filled');
        scatter(nurr1xy(2,ids), size(imgI,1) - nurr1xy(1,ids), 20, neuroNurrColor{k}, 'filled');
        k = k + 1;
    end

    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition =  0.5 * [0 0  30/ 1.5656  20 ];
    fig.PaperPositionMode = 'manual';
    set(gcf, 'papersize',  0.5 * [30/ 1.5656 , 20]);
    title('NURR1 positive cells in space')
    saveas(fig, [resultdir dataname datafield '_Classification_NURR1.png']);
end


%% Save the Statistics to File

statsfile = [resultdir dataname datafield '_Results.mat'];
save(statsfile, 'channellabel', 'stats', 'statsCh', 'statsChN', 'statsMeasure', 'cellclass', 'zonelineshires', 'distance_inner', 'distance_outer');
%load(statsfile)

%% Save For Mathematica

centers = [stats.Centroid];

intensities = zeros(length(statsMeasure), size(centers,2));
intensities_normalized = zeros(length(statsMeasure), size(centers,2));
for c = 1:length(statsMeasure)
   intensities(c, :) = [statsMeasure{c}.(mode)];
   %intensities_normalized(c, :) = [statsChN{c}.MedianIntensity];
end

statsfilemathematica = [resultdir dataname datafield '_Results_Mathematica' '.mat'];

save(statsfilemathematica, '-v6', 'channellabel', 'centers', 'distance_inner', 'distance_outer', 'zonelineshires', 'intensities', 'cellclass');

