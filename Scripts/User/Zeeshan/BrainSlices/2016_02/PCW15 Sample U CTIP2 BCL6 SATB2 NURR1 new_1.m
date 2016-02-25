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
dataname = 'PCW15 Sample U CTIP2 BCL6 SATB2 NURR1 new'

datafield = ' 1';

verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Source %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

is = ImageSourceBF([datadir  dataname '.lsm']);
clc
is.printInfo

%%

is.setReshape('S', 'Uv', [13, 18]);
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
region = struct('U', 1:7, 'V', 10:13, 'C', 1);

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
imgsRaw = {imgsRaw{1}, imgsRaw{2}, imgsRaw{3}, imgsRaw{4}, imgsRaw{5}};

channellabel  = {'CTIP2', 'BCL6', 'SATB2', 'NURR1', 'DAPI'};
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
implot(imgCfull)

%%

sub = true;
%subregion = {430:6500, 2100:3300};
%subregion = {700:5200, ':'};
subregion = {':', ':'};

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
imgmaskHi1 = imgRaw{4} > 0.75;
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
implot(imgCfull)
%
%figure(100); clf
%implot(permute(imgC,[2,1,3]));

%%
%p  = impoly(gca, line_outer)

%%
line_outer = [1964.780859375 4128.1640625;1864.191796875 3966.10390625;1713.308203125 3714.63125;1573.601171875 3457.5703125;1500.953515625 3273.15703125;1568.012890625 3066.390625;1640.660546875 2954.625;1785.955859375 2809.3296875;1746.837890625 2731.09375;1551.248046875 2770.21171875;1316.540234375 2814.91796875;1160.068359375 2596.975;1025.949609375 2311.97265625;914.183984375 2032.55859375;841.536328125 1753.14453125;835.948046875 1501.671875;735.358984375 1116.08046875;757.712109375 914.902343750001;724.182421875 702.54765625;618.005078125 590.78203125;545.357421875 249.896875000001;472.709765625 -7.16406249999953];line_inner = [7150.705859375 4133.75234375;7150.705859375 3675.51328125;6999.822265625 3658.7484375;6765.114453125 3563.74765625;6496.876953125 3334.628125;6155.991796875 2982.56640625;5781.576953125 2552.26875;5541.280859375 2177.85390625;5401.573828125 1920.79296875;5289.808203125 1697.26171875;5110.983203125 1194.31640625;4993.629296875 870.196093750001;4809.216015625 484.6046875;4747.744921875 210.778906250001;4708.626953125 -29.5171874999996]
line_ISVZ_OSVZ = [7066.881640625 4128.1640625;6731.584765625 4016.3984375;6463.347265625 3776.10234375;6150.403515625 3479.9234375;5798.341796875 3122.2734375;5613.928515625 2921.0953125;5390.397265625 2580.21015625;5340.102734375 2451.6796875;5178.042578125 2345.50234375;5088.630078125 2177.85390625;4988.041015625 1982.2640625;4831.569140625 1574.31953125;4697.450390625 1289.3171875;4591.273046875 1021.0796875;4529.801953125 825.489843750001;4451.566015625 568.428906250001;4401.271484375 372.8390625;4345.388671875 104.601562500001;4323.035546875 15.1890625000005;4334.212109375 -46.2820312499998];
line_SP_CP = [2713.610546875 4116.9875;2573.903515625 3865.51484375;2411.843359375 3546.9828125;2311.254296875 3351.39296875;2272.136328125 3094.33203125;2328.019140625 2893.15390625;2406.255078125 2759.03515625;2417.431640625 2619.328125;2389.490234375 2496.3859375;2294.489453125 2356.67890625;2171.547265625 2239.325;1953.604296875 2166.67734375;1836.250390625 2105.20625;1707.719921875 1870.4984375;1573.601171875 1490.4953125;1473.012109375 1266.9640625;1422.717578125 1004.31484375;1344.481640625 663.4296875;1327.716796875 417.5453125;1288.598828125 126.954687500001;1266.245703125 -1.57578124999964]
line_CP_M = [2149.194140625 4122.57578125;2059.781640625 4066.69296875;2015.075390625 3977.28046875;1948.016015625 3820.80859375;1841.838671875 3753.74921875;1774.779296875 3546.9828125;1679.778515625 3312.275;1696.543359375 3150.21484375;1746.837890625 2954.625;1931.251171875 2809.3296875;1920.074609375 2680.79921875;1746.837890625 2636.09296875;1467.423828125 2664.034375;1249.480859375 2529.915625;1120.950390625 2267.26640625;1087.420703125 1993.440625;1037.126171875 1920.79296875;1025.949609375 1775.49765625;970.066796875 1619.02578125;891.830859375 1468.1421875;858.301171875 1272.55234375;841.536328125 1093.72734375;813.594921875 847.842968750001;808.006640625 713.724218750001;713.005859375 540.487500000001;662.711328125 423.133593750001;668.299609375 233.13203125;573.298828125 -18.3406249999998];
line_VZ_ISVZ = [7150.705859375 4122.57578125;7150.705859375 3826.396875;6955.116015625 3781.690625;6653.348828125 3641.98359375;6429.817578125 3496.68828125;6223.051171875 3267.56875;5949.225390625 2982.56640625;5653.046484375 2624.91640625;5507.751171875 2418.15;5328.926171875 2127.559375;5183.630859375 1876.08671875;5133.336328125 1674.90859375;5043.923828125 1468.1421875;5021.570703125 1339.61171875;4887.451953125 1021.0796875;4758.921484375 741.665625000001;4641.567578125 400.78046875;4563.331640625 149.3078125;4535.390234375 43.1304687500005];
line_IZ_SP = [3630.088671875 4122.57578125;3479.205078125 3871.103125;3305.968359375 3569.3359375;3210.967578125 3356.98125;3127.143359375 3183.74453125;3037.730859375 2932.271875;3043.319140625 2786.9765625;3071.260546875 2680.79921875;3071.260546875 2529.915625;3004.201171875 2272.8546875;2881.258984375 2054.91171875;2758.316796875 1926.38125;2674.492578125 1836.96875;2573.903515625 1641.37890625;2462.137890625 1440.20078125;2333.607421875 1155.1984375;2216.253515625 780.78359375;2160.370703125 467.839843750001;2115.664453125 -57.4585937499996];
line_OSVZ_IZ = [4518.625390625 4128.1640625;4200.093359375 3910.22109375;3915.091015625 3491.1;3702.736328125 3066.390625;3551.852734375 2580.21015625;3317.144921875 2105.20625;3216.555859375 1915.2046875;3138.319921875 1635.790625;3054.495703125 1373.14140625;3004.201171875 1144.021875;2903.612109375 881.372656250001;2864.494140625 579.605468750002;2875.670703125 361.662500000001;2786.258203125 -51.8703124999984];

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

class_thresholds = {0.15, 0.17, 0.11, 0.29, 0};

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

