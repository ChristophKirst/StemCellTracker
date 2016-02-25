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

dataname = 'PCW10 Sample H CTIP2 TLE4 SATB2 NURR1';

datafield = ' 1';

verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Source %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

is = ImageSourceBF([datadir  dataname '.lsm']);
clc
is.printInfo

%%

is.setReshape('S', 'Uv', [7, 16]);
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
region = struct('U', 5:7, 'V', 10:15, 'C', 1);

imgs = is.cell(region);
size(imgs)

if verbose 
   figure(1); clf
   implottiling(imgs, 'link', false)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Align
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
sh = alignImages(imgs, 'alignment', 'RMS', 'overlap.max', 150);
stmeth = 'Interpolate';
%stmeth = 'Pyramid';

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
subregion = {750:2750, 550:5350};

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
imgmaskHi1 = imgsStRaw{4} > 0.6;
imgmaskHi1 = imdilate(imgmaskHi1, strel('disk', 3));
imgmaskHi1 = postProcessSegments(bwlabeln(imgmaskHi1), 'volume.min', 100, 'fillholes', false) > 0;
imgmaskHi1 = not(imgmaskHi1);

imgmaskHi2 = imgsStRaw{1} > 0.5;
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


statsDAPI = imstatistics(imgSP, {mode},  imgsSt{5});

if verbose 
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
%save(savefile)
load(savefile)


%%

savefile = [datadir, dataname, datafield, '_Results.mat'];
%save(savefile)
load(savefile)

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

if verbose 
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

cth = {0.2, 0.4, 0.5, 0.5, 0};

neuroClass = zeros(nch, length(isgood));
for c = 1:nch
   neuroClass(c,:) = double(and([statsChN{c}.(mode)] > cth{c}, isgood));
end
neuroClassTotal = int32(fix(neuroClass(1,:) + 2 * neuroClass(2,:) + 4 * neuroClass(3,:) + 8 * neuroClass(4,:)));

neuroBaseColor = {[1,0,0], [0,1,0], [1,1,0], [1,0,1],[0,0,1]};


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


h = figure(7); clf;
implottiling({imgsStRaw{1:nch}; imgClass1{:}}, 'titles', [chlabel(1:nch)', chlabel(1:nch)']');
set(gcf, 'Name', 'Classifications');

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 30];
fig.PaperPositionMode = 'manual';
set(gcf, 'papersize', [8,30]);

%saveas(h, [datadir dataname datafield '_Classification' '.pdf']);


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
%implot(imgsStRaw{5})
%
%figure(100); clf
implot(permute(imgC,[2,1,3]));

%%

p  = impoly

%%
inner = [-1.76890359168283 113.858695652174;116.214083175803 136.547731568998;216.04584120983 177.387996219281;302.264177693761 181.925803402646;402.095935727788 191.001417769376;501.927693761814 195.539224952741;606.297258979206 222.76606805293;751.507088846881 249.992911153119;892.179111531191 290.833175803402;973.859640831758 340.749054820415;1019.23771266541 358.900283553875;1105.45604914934 367.975897920604;1227.97684310019 390.664933837429;1341.42202268431 408.816162570888;1445.7915879017 436.043005671077;1550.16115311909 458.732041587901;1690.8331758034 472.345463137996;1786.12712665406 490.496691871455;1858.7320415879 522.261342155009;1935.8747637051 531.336956521739;2013.01748582231 544.950378071833;2090.16020793951 540.412570888469;2208.143194707 531.336956521739;2285.2859168242 526.799149338374;2425.95793950851 517.723534971644;2575.70557655955 531.336956521739;2739.06663516068 535.874763705104;2829.82277882798 549.488185255198;2911.50330812855 554.025992438563;3038.56190926276 563.101606805293;3165.62051039698 576.715028355387;3251.83884688091 590.328449905482;3328.98156899811 613.017485822306;3401.58648393195 640.244328922495;3528.64508506616 690.160207939508;3651.16587901701 717.387051039697;3714.69517958412 740.076086956522;3796.37570888469 762.765122873346;3900.74527410208 821.756616257089;4005.11483931947 848.983459357278;4086.79536862004 871.672495274102;4182.0893194707 907.974952741021;4263.76984877127 948.815217391304;4331.83695652174 980.579867674858;4386.29064272212 1025.95793950851;4431.66871455577 1071.33601134215;4490.66020793951 1112.17627599244;4567.80293005671 1116.7140831758;4708.47495274102 1143.94092627599;4812.84451795841 1175.70557655955];
outer = [-1.76890359168283 1611.33506616257;138.903119092627 1638.56190926276;220.583648393194 1643.09971644613;284.112948960302 1647.63752362949;379.406899810964 1661.25094517958;452.011814744801 1652.17533081285;538.230151228733 1656.71313799622;633.524102079395 1656.71313799622;719.742438563327 1647.63752362949;805.960775047259 1652.17533081285;869.490075614367 1638.56190926276;923.943761814745 1620.4106805293;960.246219281664 1634.0241020794;1055.54017013233 1634.0241020794;1096.38043478261 1638.56190926276;1178.06096408318 1611.33506616257;1291.5061436673 1602.25945179584;1355.0354442344 1602.25945179584;1373.18667296786 1593.18383742911;1427.64035916824 1597.72164461248;1532.00992438563 1584.10822306238;1609.15264650284 1565.95699432892;1672.68194706994 1538.73015122873;1772.51370510397 1516.04111531191;1872.345463138 1506.96550094518;1963.10160680529 1511.50330812854;2022.09310018904 1506.96550094518;2099.23582230624 1488.81427221172;2171.84073724008 1484.27646502836;2244.44565217391 1470.66304347826;2321.58837429112 1457.04962192817;2376.04206049149 1457.04962192817;2457.72258979206 1443.43620037807;2525.78969754253 1434.36058601134;2607.4702268431 1429.82277882798;2711.83979206049 1425.28497164461;2798.05812854442 1416.20935727788;2902.42769376182 1407.13374291115;2997.72164461248 1402.59593572779;3052.17533081286 1420.74716446125;3111.1668241966 1438.89839319471;3188.3095463138 1447.97400756144;3306.29253308129 1443.43620037807;3378.89744801512 1457.04962192817;3442.42674858223 1466.1252362949;3537.72069943289 1479.73865784499;3619.40122873346 1502.42769376182;3669.31710775047 1493.35207939509;3732.84640831758 1506.96550094518;3787.30009451796 1547.80576559546;3864.44281663516 1588.64603024575;3918.89650283554 1606.79725897921;3986.96361058601 1615.87287334594;4055.03071833649 1638.56190926276;4114.02221172023 1652.17533081285;4159.40028355388 1674.86436672968;4213.85396975425 1706.62901701323;4277.38327032136 1733.85586011342;4340.91257088847 1756.54489603025;4404.44187145558 1779.23393194707;4445.28213610586 1788.3095463138;4508.81143667297 1815.53638941399;4563.26512287335 1820.07419659735;4608.643194707 1847.30103969754;4631.33223062382 1856.37665406427;4694.86153119093 1856.37665406427;4767.46644612476 1860.91446124764;4808.30671077505 1865.452268431];
cortexU = outer;
vz = [-1.76890359168283 363.43809073724;116.214083175803 390.664933837429;225.121455576559 404.278355387523;275.037334593572 399.740548204159;352.180056710775 408.816162570888;447.474007561436 426.967391304348;529.154536862004 440.580812854442;583.608223062382 431.505198487713;656.213137996219 440.580812854442;710.666824196597 449.656427221172;783.271739130435 490.496691871455;864.952268431002 490.496691871455;969.321833648393 522.261342155009;1032.8511342155 554.025992438563;1146.29631379962 576.715028355387;1282.43052930057 599.404064272212;1364.11105860113 617.555293005671;1473.01843100189 635.70652173913;1545.62334593573 653.85775047259;1631.84168241966 662.933364839319;1722.59782608696 694.698015122873;1831.50519848771 712.849243856333;1890.49669187146 731.000472589792;1990.32844990548 740.076086956522;2081.08459357278 749.151701323251;2194.5297731569 721.924858223062;2280.74810964083 708.311436672968;2389.65548204159 708.311436672968;2457.72258979206 721.924858223062;2530.3275047259 740.076086956522;2612.00803402647 721.924858223062;2689.15075614367 721.924858223062;2793.52032136106 726.462665406427;2884.27646502836 735.538279773157;2961.41918714556 740.076086956522;3038.56190926276 749.151701323251;3188.3095463138 771.840737240075;3269.99007561437 803.605387523629;3347.13279773157 821.756616257089;3415.19990548204 835.370037807183;3496.88043478261 858.059073724007;3560.40973534972 871.672495274102;3633.01465028356 898.899338374291;3732.84640831758 935.20179584121;3850.82939508507 966.966446124764;3937.047731569 989.655482041588;4014.1904536862 1035.03355387524;4073.18194706994 1062.26039697543;4186.62712665407 1116.7140831758;4263.76984877127 1157.55434782609;4327.29914933838 1193.85680529301;4395.36625708885 1248.31049149338;4454.35775047259 1289.15075614367;4522.42485822306 1311.83979206049;4595.0297731569 1325.45321361059;4685.7859168242 1339.06663516068;4772.00425330813 1352.68005671078;4808.30671077505 1361.7556710775];
cortexL = [-1.76890359168283 1180.24338374291;75.3738185255195 1189.31899810964;193.356805293005 1193.85680529301;275.037334593572 1202.93241965974;470.163043478261 1221.08364839319;556.381379962193 1225.62145557656;683.439981096408 1234.69706994329;842.263232514178 1239.23487712665;910.33034026465 1243.77268431002;1023.77551984877 1248.31049149338;1105.45604914934 1248.31049149338;1182.59877126654 1252.84829867675;1264.27930056711 1261.92391304348;1382.26228733459 1266.46172022684;1463.94281663516 1270.99952741021;1581.92580340265 1252.84829867675;1704.4465973535 1239.23487712665;1858.7320415879 1202.93241965974;2044.78213610586 1175.70557655955;2144.61389413989 1157.55434782609;2285.2859168242 1134.86531190926;2403.26890359168 1121.25189035917;2503.10066162571 1121.25189035917;2621.0836483932 1107.63846880907;2698.2263705104 1107.63846880907;2811.67155009452 1107.63846880907;2911.50330812855 1107.63846880907;2988.64603024575 1125.78969754253;3065.78875236295 1134.86531190926;3138.39366729679 1143.94092627599;3224.61200378072 1148.47873345936;3306.29253308129 1166.62996219282;3433.3511342155 1184.78119092628;3560.40973534972 1202.93241965974;3637.55245746692 1221.08364839319;3737.38421550095 1252.84829867675;3814.52693761815 1275.53733459357;3932.50992438563 1339.06663516068;4009.65264650284 1361.7556710775;4109.48440453686 1393.52032136106;4173.01370510397 1425.28497164461;4227.46739130435 1475.20085066163;4286.45888468809 1506.96550094518;4363.60160680529 1534.19234404537;4440.7443289225 1565.95699432892;4526.96266540643 1588.64603024575;4595.0297731569 1638.56190926276;4649.48345935728 1652.17533081285;4735.70179584121 1674.86436672968;4781.07986767486 1693.01559546314;4808.30671077505 1706.62901701323];

lines = {inner, vz, cortexL, cortexU, outer};

for l = 1:length(lines)
   lines{l} = lines{l}(:,[2,1]);
end

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