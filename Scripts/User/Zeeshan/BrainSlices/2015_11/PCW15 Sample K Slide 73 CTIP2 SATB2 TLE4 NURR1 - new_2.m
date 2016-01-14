%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Brain Slices %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

initialize()
bfinitialize
%initializeParallelProcessing(12)

addpath('./Scripts/User/Zeeshan/BrainSlices/2015_11');

%datadir = '/data/Science/Projects/StemCells/Experiment/BrainSections_2015_11/';
datadir = '/home/ckirst/Science/Projects/StemCells/Experiment/BrainSections_2015_11/';

dataname = 'PCW15 Sample K Slide 73 CTIP2 SATB2 TLE4 NURR1 - new';

datafield = ' 2';

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
region = struct('U', 7:10, 'V', 1:4, 'C', 1);

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
imgsTLE4  = is.cell(region, 'C', 2);
imgsSATB2 = is.cell(region, 'C', 3);
imgsCTIP2 = is.cell(region, 'C', 4);
imgsNURR1 = is.cell(region, 'C', 5);

imgsAllRaw = {imgsSATB2, imgsCTIP2, imgsTLE4, imgsNURR1, imgsDAPI};

chlabel = {'SATB2', 'CTIP2', 'TLE4', 'NURR1', 'DAPI'};

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
save(savefile)
%load(savefile)

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
subregion = {600:3780, ':'};

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

%figure(1); clf; colormap jet
%implottiling(imgsSt, 'titles', chlabel)

%imgm =  imgsSt{5} - imopen(imgsSt{5} , strel('disk', 10));
imgm = imgsSt{5};

imgmaskLo = imgm > 0.1;
imgmaskLo = imopen(imgmaskLo, strel('disk', 2));
imgmaskLo = postProcessSegments(bwlabeln(imgmaskLo), 'volume.min', 15, 'fillholes', false) > 0;


%figure(66); clf;
%implottiling({255*(postProcessSegments(bwlabeln(imgmaskLo), 'volume.min', 10, 'fillholes', false)>0); 2*(bwlabeln(imgmaskLo)); 255*imgmaskLo}')

%get rid of blod vessels
imgmaskHi1 = imgsStRaw{4} > 0.5;
imgmaskHi1 = imdilate(imgmaskHi1, strel('disk', 3));
imgmaskHi1 = postProcessSegments(bwlabeln(imgmaskHi1), 'volume.min', 500, 'fillholes', false) > 0;
imgmaskHi1 = not(imgmaskHi1);

imgmaskHi2 = imgsStRaw{1} > 0.6;
imgmaskHi2 = imdilate(imgmaskHi2, strel('disk', 3));
imgmaskHi2 = postProcessSegments(bwlabeln(imgmaskHi2), 'volume.min', 500, 'fillholes', false) > 0;
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
   implottiling({imgC, mat2gray(imgV), mat2gray(filterDoG(imgV, 5)),  mat2gray(imgI)}');
end

%%
%imgf = imgV;
imgf = filterDoG(imgV, 8); % .* imgI;
%imgf = filterSphere(imgV, 4);
%imgf = filterDisk(imgI,12,1);
%imgF = imgBMI;


imgf2 = imgf; % - imopen(imgf, strel('disk', 1));

imgmax = imextendedmax(mat2gray(imgf2), 0.025);
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
imgWs = imimposemin(iminvert(imgsSt{5}), imgmax);
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
   implottiling({imoverlaylabel(imgsSt{5}, imgsurf, false);imoverlaylabel(imgC, imgsurf, false)})
   
end
   
imgS = imgWs;      


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
max(imgS(:))

imgSP = imgS;

mode = 'MedianIntensity';

stats = imstatistics(imgSP, {mode}, imgI);

if verbose
   figure(7); clf; 
   hist([stats.(mode)], 256)
   set(gcf, 'Name', mode)
   %hist(imgI(:), 256)
end

%% Remove weak cells

imgSP = postProcessSegments(imgS, imgI, 'intensity.median.min', 0.1, 'volume.min', 15, 'fillholes', false);

if verbose > 1
   imgSp1 = impixelsurface(imgSP);
   imgSp1 = imoverlaylabel(imgC, imgSp1, false);
   figure(5); clf;
   implot(imgSp1);
end

fprintf('watershed    : %d cells\n', max(imgS(:)))
fprintf('postprocessed: %d cells\n', max(imgSP(:)))

%% Remove weak DAPI Cells


%statsDAPI = imstatistics(imgSP, {mode},  imgsSt{5});

if verbose > 1
   figure(6); clf;
   hist([statsDAPI.(mode)], 256)
   title('DAPI intensity')
end

%%
imgSP2 = postProcessSegments(imgSP, imgsSt{5}, 'intensity.median.min', 0.08);

if verbose > 1
   imgSp2 = impixelsurface(imgSP2);
   imgSp2 = imoverlaylabel(imgC, imgSp2, false);
   figure(5); clf;
   implot(imgSp2);
end

fprintf('postprocessed : %d cells\n', max(imgSP(:)))
fprintf('dapi corrected: %d cells\n', max(imgSP2(:)))

%statsDAPI = imstatistics(imgSP2, {mode},  imgsSt{5});
%figure(6); clf;
%hist([statsDAPI.(mode)], 256)
%title('DAPI intensity')


%% Plot Results

imglab = imgSP2;
fprintf('final: %d cells\n', max(imglab(:)))

%%
if verbose 
   
   
   statsF = imstatistics(imglab, {'PixelIdxList', 'Centroid'});
   % mode = 'MedianIntensity';

   figure(7); clf; 
   for c = 1:nch
      statsCh{c} = imstatistics(imglab, statsF, mode, imgsSt{c});
      
      imsubplot(1, nch, c)
      imgch = zeros(size(imglab));
      for i = 1:length(statsF)
         imgch(statsF(i).PixelIdxList) = statsCh{c}(i).(mode);
      end
      implot(imgch);
      %imcolormap(cm(c));
      title(chlabel{c});
      drawnow
   end
end

%% Save Image Preporcessing Result

savefile = [datadir, dataname, datafield, '_ImageProcessing.mat'];
save(savefile)
%load(savefile)


%%

%savefile = [datadir, dataname, datafield, '_Results.mat'];
%save(savefile)
%load(savefile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Statistics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stats = imstatistics(imglab, {'PixelIdxList', 'Centroid'});

mode = 'MedianIntensity';

clear statsCh
for c = 1:nch
   %statsCh{c} = imstatistics(imglab, stats, mode,  imgsStRaw{c});
   statsCh{c} = imstatistics(imglab, stats, mode,  imgsSt{c});
end
  

%% Cell Center Image 

%imglabc = imlabelapplybw(imglab, @bwulterode);
cc = fix([stats.Centroid]);
imglabc = zeros(size(imglab));
ind = imsub2ind(size(imglab), cc');
imglabc(ind) = 1:length(ind);
imglabc = imdilate(imglabc, strel('disk', 2));

if verbose > 1
    figure(78); clf;
    implottiling({imglabc})
end


%% Plot Detected Cells

if verbose > 1 
   h = figure(7); clf;
   %subreg = {[600:900], [550:950], ':'};
   %subreg = {[200:300], [150:250], ':'};
   subreg = {':', ':', ':'};
   imgCfsub = imgCf(subreg{:});
   imgCsub = imgC(subreg{:});
   imglabsub = imglabc(subreg{:});

   implot(imoverlaylabel(imgCsub, imglabsub > 0, false, 'color.map', [[0,0,0]; 0.5*[1,0,0]]));
   axis off
   xlabel([]); ylabel([])

   %saveas(h, [datadir, dataname, datafield, '_Segmentation_Seeds.pdf'])
end

%% Raw Image Data Color Plot

if verbose
   h = figure(11); clf;

   %implot(filterAutoContrast(imgCf));
   implot(imgC);
   axis off
   xlabel([]); ylabel([])

   saveas(h, [datadir, dataname, datafield, '_Raw.pdf'])
end


%% 

if verbose
   h = figure(7); clf;
   implottiling({imgC, imoverlaylabel(imgC, imglabc > 0,false, 'color.map', [[0,0,0]; [1,0,0]])}');
   
   saveas(h, [datadir, dataname, datafield, '_CellDetection.pdf'])
end


%% Flourescence in Space

xy = [stats.Centroid]';
cl = 100 * [0.1, 0.6, 0.5, 0.5, 0.6];
%ct = [0.110, 0.0, 0.100]

figure(21); clf;
for c = 1:nch
   fi = [statsCh{c}.(mode)]';
   fi = imclip(fi, 0, cl(c));
   
   figure(21);
   subplot(2,nch,c+nch);
   hist(fi, 256)
   h = findobj(gca,'Type','patch');
   h.FaceColor = 'k';
   %h.EdgeColor = 'w';
   
   subplot(2,nch,c);
   %imcolormap(cm{c});
   colormap jet
   scatter(xy(:,1), xy(:,2), 10, fi, 'filled');
   %xlim([0, size(imglab,1)]); ylim([0, size(imglab,2)]);
   title(chlabel{c});
   %freezecolormap(gca)
   
   %h = figure(22); clf
   %colormap jet
   %scatter(xy(:,1), xy(:,2), 10, fi, 'filled');
   %xlim([0, size(imglab,1)]); ylim([0, size(imglab,2)]);
   %title(chlabel{c}); 
   %colorbar('Ticks', [])
   
   %saveas(h, [datadir dataname datafield '_Quantification_' lab{c} '.pdf']);
end


%% Normalized Intensities - Flourescence in Space

clear statsChN
for c = 1:nch
   statsChN{c} = statsCh{c};
   if c < 5
      for i = 1:length(statsCh{c})
         statsChN{c}(i).(mode) = statsCh{c}(i).(mode) / statsCh{5}(i).(mode);
      end
   end
end


xy = [stats.Centroid]';

ha = figure(21); clf;
set(gcf, 'Name', 'Normalized Expression')

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 30 20];
fig.PaperPositionMode = 'manual';
set(gcf, 'papersize', [30,20]);


for c = 1:nch
   fi = [statsChN{c}.(mode)]';
   
   fis = sort(fi);
   ncs = length(fis);
   p90 = int64(0.9 * ncs);
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
   title(chlabel{c});
   %freezecolormap(gca)
   
   subplot(3,nch,c)
   implot(imgsStRaw{c})
      
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
   fi{c} = [statsChN{c}.(mode)]';
   
   fis = sort(fi{c});
   ncs = length(fis);
   p95 = int64(0.99 * ncs);
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
   scattercloud(fi{pairs{n}(1)}, fi{pairs{n}(2)}, 50);
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
   subplot(3,nch,c); 
   d = imgsSt{c};
   %hist(d(:), 256);
   implot(imgsSt{c})
   title(['Raw ', chlabel{c}])
   
   subplot(3, nch, c+nch);
   dat = [statsCh{c}.(mode)];
   hist(dat(isgood), 256)
   title(chlabel{c})
   
   subplot(3, nch, c+2*nch); 
   dat = [statsChN{c}.(mode)];
   hist(dat(isgood), 256)
   title(['Normalized ' chlabel{c}])
   
end

%% Classify

cth = {0.8, 1, 0.6, 0.8, 0};

neuroClass = zeros(nch, length(isgood));
for c = 1:nch
   neuroClass(c,:) = double(and([statsChN{c}.(mode)] > cth{c}, isgood));
end
neuroClassTotal = int32(fix(neuroClass(1,:) + 2 * neuroClass(2,:) + 4 * neuroClass(3,:) + 8 * neuroClass(4,:)));

neuroBaseColor = {[1,0,0], [0,1,0], [1,1,0], [1,0,1],[0,0,1]};

%%
imgClass1 = cell(1,nch);
for c = 1:nch
   neuroColor1 = {[0,0,0], neuroBaseColor{c}};
   neuroColor1 = neuroColor1(neuroClass(c,:)+1);
   
   R = imgsSt{c}; G = R; B = R;
   for i = 1:length(stats);
      if neuroClass(c,i) > 0
         R(stats(i).PixelIdxList) =  neuroColor1{i}(1);
         G(stats(i).PixelIdxList) =  neuroColor1{i}(2);
         B(stats(i).PixelIdxList) =  neuroColor1{i}(3);
      end
   end
   imgClass1{c} = cat(3, R, G, B);
end

clear R
clear B
clear G

%%
h = figure(7); clf;
implottiling({imgsStRaw{1:nch}; imgClass1{:}}, 'titles', [chlabel(1:nch)', chlabel(1:nch)']');
set(gcf, 'Name', 'Classifications');

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 30];
fig.PaperPositionMode = 'manual';
set(gcf, 'papersize', [8,30]);

saveas(h, [datadir dataname datafield '_Classification' '.pdf']);


%%



%%
%neuroClassColor = {[0,0,0]; [0.6,0,0]; [0,0.6,0]; [0.33,0.33,0]; [0,0.33,0]; [0.33,0,0.33]; [0,0.33,0.33]; [0.5,0.5,0.5]};
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
   

%%
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


%% Class in Space
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
   
   nstat(i) = sum(neuroClassTotal + 1== i);
end

nc = num2cell(nstat);
tb = table(nc{:}, 'VariableNames', clslab)


%% save numbers
writetable(tb, [datadir dataname datafield '_Counts.txt'])


%% hist

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

saveas(h, [datadir dataname datafield '_Classification_Statistics' '.pdf']);




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
outer = [-3.73734177215101 424.963924050633;125.76012658228 400.982911392405;327.200632911393 362.613291139241;562.214556962026 338.632278481013;778.043670886077 285.874050632911;989.076582278482 242.708227848101;1190.5170886076 237.912025316456;1406.34620253165 266.689240506329;1569.4170886076 290.670253164557;1785.24620253165 271.485443037975;2005.87151898734 261.893037974683;2240.88544303798 285.874050632911;2461.51075949367 338.632278481013;2614.98924050633 338.632278481013;2797.24493670886 377.001898734177;2950.72341772152 453.741139240507;3080.22088607595 492.110759493671;3190.5335443038 554.461392405063];
inner = [3180.94113924051 2918.98924050633;3041.85126582279 2899.80443037975;2845.20696202532 2871.02721518987;2725.30189873418 2847.04620253165;2682.13607594937 2832.65759493671;2528.65759493671 2837.45379746835;2384.77151898734 2861.43481012658;2293.64367088608 2847.04620253165;2111.38797468355 2880.61962025316;1900.35506329114 2909.39683544304;1732.48797468355 2928.58164556962;1636.56392405063 2957.35886075949;1521.45506329114 2986.13607594937;1363.18037974684 3019.70949367089;1286.44113924051 3038.89430379747;1214.49810126582 3029.30189873418;1075.4082278481 3110.83734177215;873.967721518988 3273.9082278481;758.858860759495 3365.03607594937;495.067721518988 3504.12594936709;303.219620253165 3633.62341772152;125.76012658228 3724.75126582279;-8.53354430379659 3806.28670886076];
vz =  [3195.32974683544 2760.71455696203;2989.09303797469 2736.7335443038;2821.22594936709 2703.16012658228;2658.15506329114 2674.38291139241;2403.95632911393 2674.38291139241;2125.77658227848 2712.75253164557;1943.52088607595 2760.71455696203;1622.1753164557 2770.30696202532;1526.25126582279 2818.26898734177;1343.99556962025 2851.84240506329;1166.53607594937 2875.82341772152;965.095569620255 2966.95126582278;835.598101265824 3058.07911392405;619.768987341773 3211.55759493671;480.679113924052 3278.70443037975;260.053797468355 3451.36772151899;120.963924050634 3561.68037974684;20.2436708860769 3648.01202531646];
cortexL = [3180.94113924051 1417.77784810127;2907.55759493671 1398.59303797468;2591.0082278481 1317.05759493671;2284.05126582279 1254.70696202532;1962.70569620253 1211.54113924051;1698.91455696203 1211.54113924051;1248.07151898734 1235.52215189873;998.668987341774 1302.66898734177;768.451265822786 1355.42721518987;490.271518987343 1451.35126582279;188.110759493672 1547.2753164557;-8.53354430379659 1657.58797468354];
cortexU = [3190.5335443038 655.18164556962;3070.62848101266 592.831012658228;2897.96518987342 540.072784810126;2706.1170886076 492.110759493671;2499.88037974684 453.741139240507;2312.82848101266 405.779113924051;2149.75759493671 372.205696202532;1857.18924050633 333.836075949367;1718.09936708861 362.613291139241;1555.02848101266 367.409493670886;1367.97658227848 357.817088607595;1180.9246835443 343.428481012658;1017.85379746836 333.836075949367;897.948734177216 348.224683544304;730.081645569621 381.798101265823;600.584177215191 410.575316455696;466.290506329115 453.741139240507;312.812025316457 463.333544303797;144.944936708862 516.091772151899;15.4474683544313 535.276582278481];

lines = {inner, vz, cortexL, cortexU, outer};

nlines = length(lines);
lineshires = cell(1, nlines);
for l = 1:nlines
   lineshires{l} = [(1:size(imgC,1)); interp1(lines{l}(:,1), lines{l}(:,2), 1:size(imgC, 1))];
end
   
if verbose
   h = figure(101); clf
   implot(imgC)
   hold on
   for l = 1:nlines
      plot(lineshires{l}(1,:), lineshires{l}(2,:), 'r', 'LineWidth',1)
   end
   hold off
   
   saveas(h, [datadir dataname datafield '_Annotation' '.pdf']);
end


%% Calculate Relative Cell 

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


%% Plot Expressions as Relative Distances

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