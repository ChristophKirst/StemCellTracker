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
dataname = 'PCW15 Sample F CTIP2 SATB2 LHX2 MEF2C'

datafield = ' 1';

verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Source %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

is = ImageSourceBF([datadir  dataname '.lsm']);
clc
is.printInfo

%%

is.setReshape('S', 'Uv', [13, 17]);
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
region = struct('U', 7:11, 'V', 9:11, 'C', 1);

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

imgsRaw = {imgsRaw{3}, imgsRaw{4}, imgsRaw{2}, imgsRaw{5}, imgsRaw{1}};

channellabel  = {'CTIP2', 'MEF2C', 'SATB2', 'NURR1', 'DAPI'};
channelcolors = {[1,0,0], [0,1,0], [1,1,0], [1, 0, 1], [0,0,1]};

%% Preprocess

imgs = imgsRaw;

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

%% Save preprocessing state

close all
clear('imgs', 'imgsRaw', 'i', 'c')

%%
savefile = [savedir, dataname, datafield, '_Raw.mat'];
save(savefile)
%load(savefile)

%% Restrict to sub Region

figure(99); clf;
implot(imgCfull)

%%

sub = true;
%subregion = {1000:4100, 1200:1700};
subregion = {80:4300, ':'};

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
imgmaskHi1 = imgRaw{4} > 0.65;
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
                 }, 'titles', {'I5', 'I4', 'I1','igm','Lo', 'Hi1', 'Hi2', 'Mask'})
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

if verbose > 1
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

imgMeasure = {imgRaw{1}, imgRaw{2}, imgRaw{3}, img{4}, imgRaw{5}};
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
implot(imgCfull)
%
%figure(100); clf
%implot(permute(imgC,[2,1,3]));

%%
%p  = impoly

%%
line_inner = [434.598828125 2918.32265625;654.395703125 2812.0875;778.947265625 2731.4953125;852.212890625 2581.30078125;914.488671875 2478.72890625;958.448046875 2365.1671875;1028.050390625 2222.29921875;1061.019921875 2075.76796875;1086.662890625 2013.4921875;1156.265234375 1987.84921875;1170.918359375 1914.58359375;1141.612109375 1852.3078125;1123.295703125 1768.05234375;1097.652734375 1680.13359375;1104.979296875 1522.6125;1082.999609375 1306.47890625;1013.397265625 1130.64140625;885.182421875 903.517968750001;804.590234375 756.98671875;698.355078125 610.45546875;500.537890625 423.628125000001;350.343359375 255.1171875;130.546484375 -15.9656249999994];
line_VZ_ISVZ = [896.172265625 2932.97578125;991.417578125 2793.77109375;1053.693359375 2687.5359375;1152.601953125 2497.0453125;1207.551171875 2328.534375;1280.816796875 2156.36015625;1335.766015625 1965.86953125;1324.776171875 1819.33828125;1295.469921875 1603.2046875;1251.510546875 1387.07109375;1207.551171875 1229.55;1141.612109375 1035.39609375;1082.999609375 910.844531250001;936.468359375 665.4046875;863.202734375 592.1390625;727.661328125 438.281250000001;599.446484375 317.392968750001;445.588671875 174.525000000001;302.720703125 -15.9656249999994];
line_ISVZ_OSVZ = [1196.561328125 2914.659375;1258.837109375 2753.475;1310.123046875 2632.58671875;1368.735546875 2519.025;1434.674609375 2379.8203125;1463.980859375 2200.31953125;1500.613671875 2039.13515625;1526.256640625 1844.98125;1496.950390625 1724.09296875;1467.644140625 1573.8984375;1431.011328125 1427.3671875;1416.358203125 1310.1421875;1405.368359375 1207.5703125;1365.072265625 1090.3453125;1269.826953125 841.242187500001;1148.938671875 643.425000000001;929.141796875 438.281250000001;786.273828125 313.729687500001;614.099609375 189.178125000001;482.221484375 -8.63906249999946];
line_OSVZ_IZ = [2555.638671875 2914.659375;2584.944921875 2599.6171875;2617.914453125 2409.1265625;2614.251171875 2189.3296875;2577.618359375 2057.4515625;2566.628515625 1771.715625;2529.995703125 1504.29609375;2519.005859375 1302.815625;2467.719921875 1119.6515625;2456.730078125 973.120312500001;2409.107421875 789.956250000001;2310.198828125 599.465625000001;2214.953515625 474.914062500001;2061.095703125 423.628125000001;1918.227734375 321.056250000001;1797.339453125 262.443750000001;1705.757421875 196.504687500001;1639.818359375 123.239062500001;1570.216015625 -1.31249999999955];
line_IZ_SP = [2976.916015625 2900.00625;3002.558984375 2577.6375;2991.569140625 2200.31953125;2932.956640625 1987.84921875;2943.946484375 1830.328125;2921.966796875 1537.265625;2914.640234375 1361.428125;2881.670703125 1134.3046875;2779.098828125 866.88515625;2749.792578125 669.067968750001;2669.200390625 540.853125000001;2555.638671875 390.65859375;2394.454296875 291.750000000001;2203.963671875 244.127343750001;1998.819921875 174.525000000001;1881.594921875 101.259375000001;1746.053515625 -1.31249999999955]
line_SP_CP = [3596.010546875 2936.6390625;3544.724609375 2533.678125;3515.418359375 2214.97265625;3508.091796875 1965.86953125;3493.438671875 1683.796875;3464.132421875 1478.653125;3445.816015625 1288.1625;3409.183203125 1017.0796875;3339.580859375 848.56875;3258.988671875 621.445312500001;3193.049609375 489.567187500001;3119.783984375 343.035937500001;3009.885546875 222.147656250001;2804.741796875 97.5960937500008;2573.955078125 71.9531250000005;2372.474609375 57.3000000000007;2181.983984375 6.01406250000082];
line_CP_M = [4141.839453125 2907.3328125;4149.166015625 2771.79140625;4134.512890625 2617.93359375;4119.859765625 2511.6984375;4116.196484375 2343.1875;4053.920703125 2068.44140625;4031.941015625 1885.27734375;4013.624609375 1713.103125;3998.971484375 1529.9390625;3991.644921875 1335.78515625;3969.665234375 1167.27421875;3940.358984375 1035.39609375;3885.409765625 852.232031250001;3826.797265625 669.067968750001;3764.521484375 522.536718750001;3680.266015625 365.015625000001;3552.051171875 104.922656250001;3471.458984375 -1.31249999999955];
line_outer = [4237.084765625 2903.66953125;4229.758203125 2782.78125;4218.768359375 2698.52578125;4215.105078125 2573.97421875;4207.778515625 2471.40234375;4189.462109375 2207.64609375;4163.819140625 2105.07421875;4138.176171875 1987.84921875;4112.533203125 1870.62421875;4097.880078125 1731.41953125;4079.563671875 1570.23515625;4064.910546875 1365.09140625;4061.247265625 1196.58046875;3969.665234375 888.864843750001;3896.399609375 654.414843750001;3775.511328125 408.975;3676.602734375 178.188281250001;3555.714453125 -8.63906249999946];


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
   %zonelineshires{l} = [interp1(zonelines{l}(:,2), zonelines{l}(:,1), 1:size(imgC, 2)); (1:size(imgC,2))];
   par = parametrizePolygon(zonelines{l}(:,2), zonelines{l}(:,1));
   ny = size(imgC, 2);
   samplepoints = 1:ny / ny * max(par);
   zonelineshires{l} = [interp1(par, zonelines{l}(:,1), samplepoints); interp1(par, zonelines{l}(:,2), samplepoints)];

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

class_thresholds = {0.26, 0.23, 0.17, 0.22, 0};

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
cellclassfull = int32(fix(cellclass(1,:) + 2 * cellclass(2,:) + 4 * cellclass(3,:) + 8 * cellclass(4,:)));

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

nregions = 4;
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
saveas(fig, [resultdir dataname datafield '_Classification_Counts_Zones_Errors.png']);




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
fig = figure(18); clf; 
for z = 1:nzones
    subplot(1,nzones,z); hold on
    
    k = 1;
    %ord = [1, 5, 3, 6, 7, 8, 4, 2];
    ord = 1:nclassesfull;
    for i = ord
       hb = bar(k, regioncountsNmean(z,i), 'FaceColor', classfullcolor{i});
       xData = hb.XData+hb.XOffset;
       elo =regioncountsNmean(z,i) - regioncountsNstd(z,i);
       if elo < 0; elo = regioncountsNmean(z,i); else elo = regioncountsNstd(z,i); end
       errorbar(xData,regioncountsNmean(z,i),elo, regioncountsNstd(z,i), 'k.') 
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

