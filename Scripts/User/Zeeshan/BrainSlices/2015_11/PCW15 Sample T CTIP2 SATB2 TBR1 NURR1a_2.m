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

dataname = 'PCW15 Sample T CTIP2 SATB2 TBR1 NURR1a';

datafield = ' 2';

verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Source %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

is = ImageSourceBF([datadir  dataname '.lsm']);
clc
is.printInfo

%%

is.setReshape('S', 'Uv', [14, 20]);
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
region = struct('U', 1:7, 'V', 12:13, 'C', 1);

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
imgsTBR1  = is.cell(region, 'C', 2);
imgsSATB2 = is.cell(region, 'C', 3);
imgsCTIP2 = is.cell(region, 'C', 4);
imgsNURR1 = is.cell(region, 'C', 5);

imgsAllRaw = {imgsSATB2, imgsCTIP2, imgsTBR1, imgsNURR1, imgsDAPI};

chlabel = {'SATB2', 'CTIP2', 'TBR1', 'NURR1', 'DAPI'};

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
subregion = {650:6500, ':'};

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
imgmaskHi1 = imgsStRaw{4} > 0.7;
imgmaskHi1 = imdilate(imgmaskHi1, strel('disk', 3));
imgmaskHi1 = postProcessSegments(bwlabeln(imgmaskHi1), 'volume.min', 100, 'fillholes', false) > 0;
imgmaskHi1 = not(imgmaskHi1);

imgmaskHi2 = imgsStRaw{1} > 0.9;
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


%statsDAPI = imstatistics(imgSP, {mode},  imgsSt{5});

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

cth = {0.12, 0.2, 0.25, 0.7, 0};

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






%% No Zones etc ...

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
outer = [934.37936827957 2052.95396505376;902.922379032258 1903.53326612903;875.39751344086 1777.70530913978;867.533266129032 1702.99495967742;816.41565860215 1588.96337365591;718.112567204301 1435.61055107527;643.402217741935 1333.37533602151;549.03125 1203.61525537634;450.728158602151 1058.12668010753;379.949932795699 936.230846774193;309.171706989247 786.810147849462;250.189852150538 653.117943548387;179.411626344086 448.64751344086;112.565524193548 259.905577956989;61.4479166666666 102.62063172043;22.1266801075268 -7.47883064516145];
inner = [5822.00907258065 2068.68245967742;5751.23084677419 1962.51512096774;5652.92775537634 1840.61928763441;5597.87802419355 1754.1125672043;5546.76041666667 1620.42036290323;5472.0500672043 1470.99966397849;5401.27184139785 1317.64684139785;5358.0184811828 1199.68313172043;5283.30813172043 1054.19455645161;5204.66565860215 900.841733870968;5129.95530913979 719.964045698925;5063.10920698925 578.407594086021;5019.85584677419 448.64751344086;4976.60248655914 287.430443548387;4945.14549731183 193.059475806451;4909.7563844086 82.9600134408597;4878.29939516129 4.3175403225805];
vz = [5291.17237903226 2072.61458333333;5243.98689516129 1879.94052419355;5133.8874327957 1612.5561155914;5023.78797043011 1329.44321236559;4968.73823924731 1191.8188844086;4905.82426075269 1010.94119623656;4838.97815860215 861.520497311828;4752.47143817204 657.050067204301;4677.76108870968 480.104502688172;4614.84711021505 354.276545698925;4500.81552419355 177.330981182795;4433.96942204301 27.9102822580644];
cortexL = [1866.29267473118 2060.81821236559;1819.10719086022 1891.73689516129;1783.71807795699 1722.65557795699;1767.98958333333 1569.30275537634;1736.53259408602 1404.15356182796;1701.1434811828 1274.3934811828;1689.34711021505 1215.41162634409;1532.06216397849 1148.56552419355;1370.84509408602 999.144825268817;1229.28864247312 892.97748655914;1130.98555107527 763.217405913978;1009.08971774194 535.154233870967;997.293346774193 428.98689516129;965.836357526882 303.158938172043;895.05813172043 165.534610215053;855.73689516129 0.385416666666515];
cortexU = [1103.46068548387 2064.75033602151;1107.39280913978 1946.78662634409;1115.25705645161 1844.55141129032;1119.18918010753 1754.1125672043;1138.8497983871 1679.40221774194;1130.98555107527 1577.16700268817;1009.08971774194 1506.38877688172;859.669018817204 1376.62869623656;710.248319892473 1215.41162634409;623.741599462366 1117.10853494624;533.302755376344 955.891465053763;458.592405913978 751.421034946236;403.542674731183 605.932459677419;368.153561827957 499.765120967741;313.103830645161 428.98689516129;297.375336021505 358.208669354839;246.257728494624 232.380712365591;230.529233870968 126.213373655914;187.275873655914 -7.47883064516145];

lines = {inner, vz, cortexL, cortexU, outer};

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