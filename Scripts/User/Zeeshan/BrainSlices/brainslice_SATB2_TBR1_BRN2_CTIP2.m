%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Brain Slices %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

initialize()
bfinitialize
initializeParallelProcessing(12)

addpath('./Scripts/User/Zeeshan/BrainSlices/');

datadir = '/data/Science/Projects/StemCells/Experiment/BrainSections/';
dataname = 'Sample A Slide 33 SATB2 TBR1 BRN2 CTIP2';

verbose = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Source %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%is = ImageSourceBF('/data/Server/smb/upload/Brain sections/Sample B Slide 33 SATB2 CUX2 NURR1 CTIP2.lsm');
is = ImageSourceBF([datadir  dataname '.lsm']);
clc
is.printInfo

%%

is.setReshape('S', 'Uv', [15, 21]);
is.setCellFormat('UV')
is.setRange('C', 1);
clc; is.printInfo

%%
%is.setRange('C', 1);
%is.plotPreviewStiched('overlap', 102, 'scale', 0.1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inspect Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%region = struct('U', [8:11], 'V', [9:12], 'C', 1);
%region = struct('U', 10, 'V', 11, 'C', 1, 'X', 1:300, 'Y', 1:300);
%region = struct('U', 10, 'V', 11, 'C', 1);
region = struct('U', [9:10], 'V', [10:11], 'C', 1);

imgs = is.cell(region);
size(imgs);

figure(1); clf
implottiling(imgs, 'link', false)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Align
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
sh = alignImages(imgs, 'alignment', 'RMS', 'overlap.max', 150);
stmeth = 'Interpolate';
%stmeth = 'Pyramid';

img = stitchImages(imgs, sh, 'method', stmeth);
figure(1); clf; implot(img)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stich Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nuclear marker
is.resetRange();

%region2 = {':',':'};
region2 = {500:1900, 450:1850};

imgsDAPI  = is.cell(region, 'C', 1);
imgsBRN2  = is.cell(region, 'C', 2);
imgsCTIP2 = is.cell(region, 'C', 3);
imgsSATB2 = is.cell(region, 'C', 4);
imgsTBR1  = is.cell(region, 'C', 5);

nch = 4;

% corect illumination and background

imgsAll = {imgsSATB2, imgsCTIP2, imgsBRN2, imgsTBR1, imgsDAPI};

chlabel = {'SATB2', 'CTIP2', 'BRN1', 'TBR1', 'DAPI'};
chcols  = {[1,0,0], [0,1,0], [0,0,1], [1, 0, 1]};

%%
% figure(6); clf;
% k = 4; l = 2;
% img1= mat2gray(imopen(imgsAll{l}{k}, strel('disk', 30)));
% img2= mat2gray(imopen(imgsAll{l}{k}, strel('disk', 50)));
% img3= mat2gray(imopen(imgsAll{l}{k}, strel('disk', 150)));
% img4= mat2gray(imclose(imgsAll{l}{k}, strel('disk', 150)));
% implottiling({img1, img2; img3, img4})

parfor i = 1:3
   imgsAll{i} = cellfunc(@(x) x - imopen(x, strel('disk', 30)), imgsAll{i});
   %imgsAll{i} = cellfunc(@(x) filterAutoContrast(x/max(x(:))), imgsAll{i});
   imgsAll{i} = cellfunc(@(x) x - imopen(x, strel('disk', 150)), imgsAll{i});   
   %imgsAll{i} = cellfunc(@(x) x - imclose(x, strel('disk', 150)), imgsAll{i});
end


%%
imgsSt =cell(1,nch);
%th = [200, 50, 100, 100];
th = [0,0,0, 0,0];
cl = [1500, 2000, 2000, 3000, 4000];

for i = 1:nch+1
   imgSt = stitchImages(imgsAll{i}, sh, 'method', stmeth);
   imgSt = imgSt(region2{:});
   
   %imgSt = imgSt - imopen(imgSt, strel('disk', 20));
   
   %imgsS = filterAutoContrast(imgsS/max(imgsS(:)));
   
   figure(16);
   subplot(nch+1,1,i); 
   hist(imgSt(:), 256);
   
   %imgSt(imgSt < th(i)) = 0;
   imgsSt{i} = mat2gray(imgSt);
   %imgsSt{i} = mat2gray(imclip(imgSt, 0, cl(i)));
   
end

figure(1); clf; colormap jet
implottiling(imgsSt', 'titles', chlabel)

%%

weights = [2, 2, 1, 1];
imgsWs = cellfunc(@(x,y) x * y, imgsSt(1:nch), num2cell(weights));

imgC = zeros([size(imgsWs{1}), 3]);
for i = 1:nch
   for c = 1:3
      imgC(:,:,c) = imgC(:,:,c) + chcols{i}(c) * imgsWs{i};
   end
end
    
%imgC = imclip(imgC, 0, 2500);
imgC = imgC / max(imgC(:));
imgC = filterAutoContrast(imgC);

figure(2)
implot(imgC)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filter Image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%imgBM = filterBM(imgC, 'profile', 'np', 'sigma', 30);

imgBM = filterBM(imgC, 'profile', 'np', 'sigma', 15);
figure(3); clf;
implottiling({imgC; imgBM})

%%
% 
% imgsM = {imgBM(:,:,1), imgBM(:,:,2), imgBM(:,:,3)};
% imgCm = cellfunc(@(x) filterMedian(x, 5), imgsM);
% imgCm = cat(3, imgCm{:});
% figure(3);
% implottiling({imgC; imgCm})

%imgBM = imgBM(300:500, 400:600, :);

clc
imgCf = imgBM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imgI = imgCf(:,:,1) + imgCf(:,:,2) + imgCf(:,:,3);
imgI = imgI - imopen(imgI, strel('disk', 30));
imgI = mat2gray(imgI);

%imgmask = img > 0.125;
%imgmask = img > 0.25;
imgmask = imgI > 0.025;
%imgmask = imopen(imgmask, strel('disk', 3));
%imgmask = postProcessSegments(bwlabeln(imgmask), 'volume.min', 50) > 0;

%get rid of blod vessels
% figure(7); clf; 
% imgsA = imgsSt{1};
% for i = 2:length(imgsSt)
%    imgsA= imgsA + imgsSt{i};
% end
% implot(imclose(imgsSt{5}> 0.9, strel('disk', 5)))

imgmask = and(imgmask, ~imclose(imgsSt{5} > 0.9, strel('disk', 5)));

if verbose
   %max(img(:))
   
   figure(21); clf;
   set(gcf, 'Name', ['Masking'])
   implottiling({imgI;  imgmask})
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SLIC Segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of pixels: typical cell 9x9
npxl = fix(1.1 * numel(imgI) / (8*8))

%segment
imgS = segmentBySLIC(imgCf, 'superpixel', npxl, 'compactness', 10);

if verbose
   imgSp = impixelsurface(imgS);
   figure(5); clf;
   implot(imoverlaylabel(imgCf, imgSp, false));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
max(imgS(:))

imgSP = imgS;

stats = imstatistics(imgSP, {'MedianIntensity'}, imgI);
figure(7); clf; 
hist([stats.MedianIntensity], 256)
%hist(imgI(:), 256)

%%
imgSP = postProcessSegments(imgS, imgI, 'intensity.median.min', 0.025, 'volume.min', 15, 'fillholes', false);

if verbose
   imgSp1 = impixelsurface(imgSP);
   imgSp1 = imoverlaylabel(imgCf, imgSp1, false);
   figure(5); clf;
   implot(imgSp1);
end

%%

imgSP2 = immask(imgSP, imgmask);
imgSP2 = imlabelapplybw(imgSP2, @(x) imopen(x, strel('disk', 2)));
imgSP2 = imlabelseparate(imgSP2);

% stats = imstatistics(imgSP, {'MinIntensity'}, imgI);
% figure(7); clf; 
% hist([stats.MinIntensity], 56)
%hist(imgI(:), 256)

%%

imgSP2 = postProcessSegments(imgSP2, imgI, 'intensity.median.min', 0.0, 'volume.min', 15, 'fillholes', false);
%imgSP = postProcessSegments(imgSP, 'volume.min', 7, 'fillholes', false);
imgSP2 = imrelabel(imgSP2);
fprintf('detected: %d cells', max(imgSP2(:)))

if verbose
   imgSp = impixelsurface(imgSP2);
   figure(6); clf;
   implottiling({imgSp1; imoverlaylabel(imgCf, imgSp, false);  imoverlaylabel(imgCf, imgmask, false)});
end

%%

imgSPP = imgSP2;
%stats = imstatistics(imgSPP, {'Volume', 'PixelIdxList', 'MedianIntensity', 'Perimeter', 'Extent', 'FilledArea'}, imgI);

%
% figure(78); clf; hist([stats.MinIntensity], 256);
% 
% if verbose 
%    figure(15);clf;
%    %scatter(bbarea, [stats.Area])
%    scatter([stats.Perimeter]/2/pi, sqrt([stats.Volume]/pi))
%    
%    figure(16); clf
%    hist([stats.MedianIntensity], 256);
%    
%    figure(17); clf;
%    scatter([stats.Volume], [stats.FilledArea]);
% end
% 
% %bbarea = [stats.BoundingBox]; bbarea = reshape(bbarea, 4, []);
% %bbarea = bbarea(3:4, :); bbarea = prod(bbarea, 1);
% 
% ids = zeros(1,length(stats));
% %ids = [stats.Extent] < 0.5; %* pi /4;
% ids = or(ids, [stats.FilledArea] > ([stats.Volume] + 10));
% ids = or(ids, 0.65 * ([stats.Perimeter]) / (2 * pi) >= sqrt([stats.Volume] / pi));
% %ids = or(ids, [stats.Volume] <= 3);
% %ids = or(ids, [stats.MedianIntensity] < 0.045);
% %ids = ~ids;
% 
% imgSPP = imgSP;
% for i = find(ids)
%    imgSPP(stats(i).PixelIdxList) = 0;
% end
% imgSPP = imrelabel(imgSPP);
% max(imgSPP(:))

if verbose 
   %figure(6); clf;
   %implottiling({imoverlaylabel(imgCf, impixelsurface(imgSPP), false); imoverlaylabel(imgCf, imgSP, false)});
   
   statsF = imstatistics(imgSPP, {'PixelIdxList', 'Centroid'});
   
%    mode = 'MedianIntensity';
%    statsR = imstatistics(imgSP, stats, mode,  imgCf(:,:,1));
%    statsG = imstatistics(imgSP, stats, mode,  imgCf(:,:,2));
%    statsB = imstatistics(imgSP, stats, mode,  imgCf(:,:,3));


   mode = 'MedianIntensity';

   figure(7); clf; 
   for c = 1:nch
      statsCh{c} = imstatistics(imgSPP, statsF, mode,  imgsSt{c});
      
      imsubplot(nch, 1, c)
      imgch = zeros(size(imgSPP));
      for i = 1:length(statsF)
         imgch(statsF(i).PixelIdxList) = statsCh{c}(i).(mode);
      end
      
      implot(imgch);
      title(chlabel{c});
      drawnow
   end
end

imglab = imgSPP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Statistics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stats = imstatistics(imglab, {'PixelIdxList', 'Centroid'});

mode = 'MedianIntensity';

clear statsCh
for c = 1:nch
   statsCh{c} = imstatistics(imglab, stats, mode,  imgsSt{c});
end
  

%%
%imglabc = imlabelapplybw(imglab, @bwulterode);
cc = fix([stats.Centroid]);
imglabc = zeros(size(imglab));
ind = imsub2ind(size(imglab), cc');
imglabc(ind) = 1;
imglabc = imdilate(imglabc, strel('disk', 2));


%%
h = figure(7); clf;
subreg = {[600:900], [550:950], ':'};
imgCfsub = imgCf(subreg{:});
imgCsub = imgC(subreg{:});
imglabsub = imglabc(subreg{:});

implot(imoverlaylabel(imgCsub, imglabsub > 0,false, 'color.map', [[0,0,0]; 0.2*[1,1,1]]));
axis off
xlabel([]); ylabel([])

%saveas(h, [datadir, dataname, '_Segmentation_Seeds.pdf'])

h = figure(8); clf;
imgCsub = imgC(subreg{:});

implottiling({imgCsub; imgCsub});
axis off
xlabel([]); ylabel([])

%saveas(h, [datadir, dataname, '_Segmentation_Segments.pdf'])

%% 
h = figure(9); clf;
imgCsub = imgCf(subreg{:});

implot(imgCsub);
axis off
xlabel([]); ylabel([])

%saveas(h, [datadir, dataname, '_Segmentation_Filtered.pdf'])


%%
h = figure(10); clf;
imgCsub = imgCf(subreg{:});

implot(imgCsub);
axis off
xlabel([]); ylabel([])

%saveas(h, [datadir, dataname, '_Segmentation_Raw.pdf'])

%%
h = figure(11); clf;

%implot(filterAutoContrast(imgCf));
implot(imgCf);
axis off
xlabel([]); ylabel([])

%saveas(h, [datadir, dataname, '_Raw.pdf'])


%% Individual channels

for c = 1:3

   h = figure(11+c); clf;

   %implot(filterAutoContrast(imgCf));
   imgCfc = imgCf;
   imgCfc(:,:,setdiff([1,2,3], c)) = 0;
   implot(imgCfc);
   axis off
   xlabel([]); ylabel([])

   %saveas(h, [datadir, dataname, '_Raw_' lab{c} '.pdf'])
end

%% 

figure(7); clf;
implottiling({imgCf, imgC, imoverlaylabel(imgCf, imglabc> 0,false, 'color.map', [[0,0,0]; [1,0,0]])}');


%% Flourescence in Space

xy = [stats.Centroid]';
cm = {'r', 'g', 'b', 'm'};
cl = [0.5, 0.6, 0.4, 0.4];
%ct = [0.110, 0.0, 0.100]

figure(21); clf;
for c = 1:nch
   fi = [statsCh{c}.(mode)]';
   fi = imclip(fi, 0, cl(c));
   
   figure(21);
   subplot(2,nch,c+nch);
   hist(fi, 256)
 
   subplot(2,nch,c);
   %imcolormap(cm{c});
   colormap jet
   scatter(xy(:,1), xy(:,2), 10, fi, 'filled');
   xlim([0, size(imglab,1)]); ylim([0, size(imglab,2)]);
   title(chlabel{c});
   %freezecolormap(gca)
   
   h = figure(22); clf
   colormap jet
   scatter(xy(:,1), xy(:,2), 10, fi, 'filled');
   xlim([0, size(imglab,1)]); ylim([0, size(imglab,2)]);
   title(chlabel{c}); 
   colorbar('Ticks', [])
   
   %saveas(h, [datadir dataname '_Quantification_' lab{c} '.pdf']);
end



%% Flourescence Expression

figure(22); clf;
fi = cell(1,3);
for c = 1:nch
   fi{c} = [statsCh{c}.(mode)]';
   fi{c} = imclip(fi{c},0, cl(c));
   fi{c} = mat2gray(fi{c});
end
 
pairs = {[1,2], [1,3], [1,4], [2, 3], [2,4], [3,4]};

np = length(pairs);

for n = 1:np
   subplot(1, np,n)
   %fi = imclip(fi, 0, cl(c));
   %scatter(fi{pairs{n}(1)}, fi{pairs{n}(2)}, 10, 'b', 'filled');
   scattercloud(fi{pairs{n}(1)}, fi{pairs{n}(2)}, 50);
   xlim([0,1]); ylim([0,1]);
   xlabel(chlabel{pairs{n}(1)}); ylabel(chlabel{pairs{n}(2)});
   %freezecolormap(gca)
end


for n = 1:np
   h = figure(50+n); clf;
   %fi = imclip(fi, 0, cl(c));
   %scatter(fi{pairs{n}(1)}, fi{pairs{n}(2)}, 10, 'b', 'filled');
   scattercloud(fi{pairs{n}(1)}, fi{pairs{n}(2)});
   xlim([0,1]); ylim([0,1]);
   xlabel(chlabel{pairs{n}(1)}); ylabel(lab{pairs{n}(2)});
   
   %saveas(h, [datadir dataname '_Quantification_Scatter_' chlabel{pairs{n}(1)} ,'_' chlabel{pairs{n}(2)} '.pdf']);
   %freezecolormap(gca)
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Classify Cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Mask for non valid regoins

xp = [0.5, 44.8458, 67.84, 90.8341, 177.8834, 241.9385, 332.2726, 471.8798, 616.4144, 790.5129, 910.4109, 984.3206, 1028.6665, 1081.2245, 1166.6313, 1212.6196, 1243.8259, 1230.6864, 1252.0381, 1281.602, 1298.0264, 1312.8083, 1319.3781, 1311.1659, 1334.16, 1362.0815, 1403.1424, 1401.5, 1334.16, 1279.9596, 1253.6805, 1207.6923, 1156.7767, 1073.0123, 943.2597, 885.7743, 825.0041, 747.8095, 718.2456, 670.6149, 621.3417, 578.6383, 526.0803, 499.8013, 374.976, -1.1424, 0.5];
yp = [613.9508, 528.544, 502.2649, 426.7128, 361.0152, 305.1723, 239.4748, 193.4865, 173.7773, 185.2743, 203.3411, 219.7655, 247.687, 278.8933, 323.2392, 375.7972, 398.7913, 413.5733, 439.8523, 482.5557, 526.9015, 559.7503, 592.5991, 622.163, 646.7995, 663.2239, 697.7151, 76.8734, 91.6553, 103.1524, 96.5826, 101.51, 99.8675, 101.51, 96.5826, 88.3705, 67.0188, 47.3095, 52.2368, 47.3095, 17.7456, 21.0305, 16.1032, 4.6061, -5.2485, -3.6061, 613.9508];

imgbad = ~poly2mask(xp',yp', size(imgI,1), size(imgI,2))';

figure(6); clf;
implot(imgbad)

%% histograms
figure(10); clf
for c= 1:nch
   subplot(2,nch,c); 
   d = imgsSt{c};
   hist(d(:), 256);
   
   subplot(2, nch, c+nch)
   hist([statsCh{c}.(mode)], 256)
   
   title(chlabel{c})
end


%% classify

cth = {0.06, 0.1, 0.075, 0.2};

clear neuroClass
for c = 1:nch
   xy = fix([statsCh{c}.Centroid]);
   xy(1,xy(1,:) > size(imgI,1)) = size(imgI,1); xy(1,xy(1,:) < 1) = 1;
   xy(2,xy(2,:) > size(imgI,2)) = size(imgI,2); xy(2,xy(2,:) < 1) = 1;
   
   isgood = zeros(1, size(xy,2));
   for i = 1:size(xy, 2)
      isgood(i) = imgbad(xy(1,i), xy(2,i));
   end

   neuroClass(c,:) = double(and([statsCh{c}.(mode)] > cth{c}, isgood));
end
neuroClassTotal = fix(neuroClass(1,:) + 2 * neuroClass(2,:) + 4 * neuroClass(3,:) + 8 * neuroClass(4,:))+1;

neuroBaseColor = {[1,0,0], [0,1,0], [0,0,1], [1, 1,1]};


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


figure(7); clf;

implottiling({imgsSt{1:nch}; imgClass1{:}}', 'titles', [chlabel(1:nch), chlabel(1:nch)]);

saveas(h, [datadir dataname '_Classification_Channels' '.pdf']);


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

ncls = 2^nch;
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
   
R = zeros(size(imgsSt{1})); G = R; B = R;
for p = 1:length(stats);
   nc = neuroClassColor{neuroClassTotal(p)};
   
   R(stats(p).PixelIdxList) =  nc(1);
   G(stats(p).PixelIdxList) =  nc(2);
   B(stats(p).PixelIdxList) =  nc(3);
end

imgClass = cat(3, R, G, B);

% 
% figure(7); clf;
% implottiling({imgsSt{1:nch}, imgCf, imgClass}', 'tiling', [3,2]);


h = figure(7); clf;
implottiling({imgsSt{1:nch}, imgC; imgClass1{:}, imgClass}', 'titles', [chlabel(1:nch), 'merge', chlabel(1:nch), 'merge']);

saveas(h, [datadir dataname '_Classification_Channels' '.pdf']);
%saveas(h, [datadir dataname '_Classification_Image' '.pdf']);


%% Class in Space
isgood = logical(isgood);

xy = [stats.Centroid]';
xy = xy(isgood, :)

h = figure(27)

cdat = (cell2mat({[0,0,0], neuroClassColor{:}}'))
   
scatter(xy(:,1), xy(:,2), 5, cdat(neuroClassTotal(isgood) + 1, :), 'filled');
xlim([0, size(imglab,1)]); ylim([0, size(imglab,2)]);
%title(chlabel{c});
%freezecolormap(gca)

pbaspect([1,1,1])

saveas(h, [datadir dataname '_Quantification_Space.pdf']);



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
   
   nstat(i) = sum(neuroClassTotal == i);
end

nc = num2cell(nstat);
tb = table(nc{:}, 'VariableNames', clslab)


%% save numbers
writetable(tb, [datadir dataname '_Counts.txt'])


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
saveas(h, [datadir dataname '_Classification_Statistics' '.pdf']);



%%

