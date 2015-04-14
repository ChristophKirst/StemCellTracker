%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Brain Slices %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

initialize()
%%

bfinitialize
initializeParallelProcessing(12)

addpath('./Scripts/User/Zeeshan/BrainSlices/');

datadir = '/data/Science/Projects/StemCells/Experiment/BrainSections/';
dataname = 'PCW10 Sample A Slide 53 SATB2 TBR1 BRN2 CTIP2';

verbose = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Source %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%is = ImageSourceBF('/data/Server/smb/upload/Brain sections/Sample B Slide 33 SATB2 CUX2 NURR1 CTIP2.lsm');
is = ImageSourceBF([datadir  dataname '.czi']);
clc
is.printInfo

%%

figure(1); clf;
implottiling(is.data)


%% NURR1 

imgsCTIP2 = {is.data('C', 3)};
imgsSATB2 = {is.data('C', 4)};
imgsNURR1 = {is.data('C', 5)};
imgsTBR1  = {is.data('C', 2)};
imgsDAPI  = {is.data('C', 1)};

imgsRaw = {imgsSATB2, imgsCTIP2, imgsNURR1, imgsTBR1, imgsDAPI};

nch = length(imgsRaw);

% correct illumination and background
parfor i = 1:nch
   %imgsAll{i} = cellfunc(@(x) filterAutoContrast(x/max(x(:))), imgsAll{i});
   imgsAll{i} = cellfunc(@(x) x - imopen(x, strel('disk', 10)), imgsRaw{i});
   imgsAll{i} = cellfunc(@(x) x - imopen(x, strel('disk', 150)), imgsRaw{i});
end

%NURR1 is noisy
imgsAll{3} =  cellfunc(@(x) imopen(x, strel('disk', 2)), imgsAll{i});

chlabel = {'SATB2', 'CTIP2', 'NURR1', 'TBR1', 'DAPI'};
chcols  = {[1,0,0], [0,1,0], [0,0,1], [0.5, 0.5, 0], [1,1,1]};

% select relevant channels only
nch = 3;

%%

for i = 1:nch+1
   figure(17);
   subplot(nch+1,1,i); 
   hist(imgsAll{i}{1}(:), 256);
   
   figure(18); clf;
   implottiling(cellfunc(@(x) x{1}, imgsAll)')
end

%%

subregion = {':', ':', ':'};

imgsSt = cell(1,nch);
imgsStRaw = cell(1,nch);

%th = [200, 50, 100];
th = [285,560,285, 900,0];
cl = [700, 5000, 1000, 3000, 5000];

for i = 1:nch +1 
   fprintf('processing image channel: %d / %d\n', i, nch+1)
   
   imgsS = imgsAll{i}{1};
   
   %max(imgsS(:))
   
   imgsS = imgsS(subregion{:});
   %imgsS = mat2gray(imgsS);
   %imgsS = filterAutoContrast(imgsS/max(imgsS(:)));
   
   figure(16);
   subplot(nch+1,1,i); 
   hist(imgsS(:), 256);
   title(chlabel{i})
   
   imgsS(imgsS < th(i)) = th(i);
   imgsSt{i} = imclip(imgsS, 0, cl(i));  
   imgsSt{i} = mat2gray(imgsSt{i});
  
   imgsS = imgsRaw{i}{1};
   imgsS = imgsS(subregion{:});
   %imgsS = mat2gray(imgsS);
   %imgsS = filterAutoContrast(imgsS/max(imgsS(:)));
   
   imgsS(imgsS < th(i)) = th(i);
   imgsStRaw{i} = imclip(imgsS, 0, cl(i)); 
   imgsStRaw{i} = mat2gray(imgsStRaw{i});
   
end

if verbose
   figure(1); clf; colormap jet
   implottiling([imgsSt; imgsStRaw]')
end


%
weights = [1, 0.5, 0.5];
weights= [1,1,1];
imgsWs = cellfunc(@(x,y) x * y, imgsSt(1:nch), num2cell(weights));
%imgsWs = cellfunc(@(x) filterMedian(x, 2), imgsSt(1:nch));
imgsWs = cellfunc(@(x) imopen(x, strel('disk', 2)), imgsWs(1:nch));

imgC = zeros([size(imgsWs{1}), 3]);
imgCRaw = zeros(size(imgC));
for i = 1:nch
   for c = 1:3
      imgC(:,:,c) = imgC(:,:,c) + chcols{i}(c) * imgsWs{i};
      imgCRaw(:,:,c) = imgCRaw(:,:,c) + chcols{i}(c) * imgsSt{i};
   end
end
    
%imgC = imclip(imgC, 0, 2500);
imgC = imgC / max(imgC(:));
%imgC = filterAutoContrast(imgC);


if verbose
   figure(2)
   implottiling({imgC; imgCRaw})
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filter Image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%imgBM = filterBM(imgC, 'profile', 'np', 'sigma', 30);

% imgBM = filterBM(imgC, 'profile', 'np', 'sigma', 15);
% figure(3); clf;
% implottiling({imgCRaw; imgC; imgBM})

%%
% 
% imgsM = {imgBM(:,:,1), imgBM(:,:,2), imgBM(:,:,3)};
% imgCm = cellfunc(@(x) filterMedian(x, 5), imgsM);
% imgCm = cat(3, imgCm{:});
% figure(3);
% implottiling({imgC; imgCm})

%imgBM = imgBM(300:500, 400:600, :);

clc
imgCf = imgC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imgHSV = rgb2hsv(imgCf);

imgI = imgCf(:,:,1) + imgCf(:,:,2) + imgCf(:,:,3) + 0* imgsAll{end}{1};
%imgI = imgI - imopen(imgI, strel('disk', 30));
imgI = mat2gray(imgI);

%imgmask = img > 0.125;
imgmask = imgHSV(:,:,3) < 0.9;
%imgmask = and(imgmask, mat2gray(imgsSt{4}) < 215/255);
%imgmask = and(and(and(imgI > 0.05, imgI < 240/255), mat2gray(imgsSt{4}) < 215/255), imgHSV(:,:,3) < 0.9);
%imgmask = ~imclose(~imgmask, strel('disk', 2));
%imgmask = imopen(imgmask, strel('disk', 3));
imgmask = ~(postProcessSegments(bwlabeln(~imgmask), 'volume.min', 50) > 0);
imgmask = ~imdilate(~imgmask, strel('disk', 5));

imgmask = and(imgmask, imgI > 0.005);
%get rid of blod vessels
% figure(7); clf; 
% imgsA = imgsSt{1};
% for i = 2:length(imgsSt)
%    imgsA= imgsA + imgsSt{i};
% end
% implot(imclose(imgsSt{5}> 0.9, strel('disk', 5)))

%imgmask = and(imgmask, ~imclose(mat2gray(imgsSt{4}) > 0.975, strel('disk', 5)));
%imgmask = and(imgmask, ~imclose(mat2gray(imgI) > 0.975, strel('disk', 5)));



if verbose
   %max(img(:))
   
   figure(21); clf;
   colormap gray
   set(gcf, 'Name', ['Masking'])
   implottiling({imgI; mat2gray(imgsSt{1}); imgmask})
    
   imgSeg = immask(imgCf, imgmask); 
   figure(22); clf;
   implottiling({imgCf; imgSeg})
end


  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % number of pixels: typical cell 9x9
% npxl = fix(%%
% 
% imgSeg = immask(imgCf, imgmask);
% 
% %segment
% imgS = segmentBySLIC(imgCf, 'superpixel', npxl, 'compactness', 40);
% 
% if verbose
%    imgSp = impixelsurface(imgS);
%    figure(5); clf;
%    implot(imoverlaylabel(imgCf, imgSp, false));
% end

%% HSV

imgHSV = rgb2hsv(imgC);

if verbose
   figure(7);
   colormap jet
   implottiling({imgHSV(:,:,1), imgHSV(:,:,2); imgHSV(:,:,3), imgC}')
end

%% Seeds

imgBMI = sum(imgCf,3);
imgHSV = rgb2hsv(imgCf);

imgIf = imgHSV(:,:,3);
%imgIf = filterSphere(imgIf, 4);
%imgIf = filterDisk(imgIf, 8, 1, 1, -1);
imgIf = filterDoG(imgIf, 8);
imgIf = mat2gray(imgIf);
imgImax = imextendedmax(imgIf, 0.01);   
%imgImax = imregionalmax(imgIf);

stats = imstatistics(bwlabel(imgImax), {'PixelIdxList', 'MaxIntensity'}, imgIf);
statsF = stats([stats.MaxIntensity] < 0.35);
pxl = {statsF.PixelIdxList};
pxl = cat(1, pxl{:});
imgImax(pxl) = 0;

%figure(10); clf
%implottiling({imgC; imoverlay(imgIf, imgImax); imoverlay(imgC, imgImax, [1,1,1])})

imgmax = cell(1,nch);
for i = 1:nch
   imgmax{i} = imgCf(:,:,i);
   %imgmax{i} = filterSphere(imgmax{i}, 3); 
   %imgmax{i} = filterDisk(imgmax{i}, 8, 1, 1, -1);
   imgmax{i} = filterDoG(imgmax{i}, 8);
   imgmax{i} = mat2gray(imgmax{i});
   imgmax{i} = imextendedmax(imgmax{i}, 0.01);   
   
   stats = imstatistics(bwlabel(imgmax{i}), {'PixelIdxList', 'MinIntensity'}, imgCf(:,:,i));

   statsF = stats([stats.MinIntensity] < 0.1);
   pxl = {statsF.PixelIdxList};
   pxl = cat(1, pxl{:});
   imgmax{i}(pxl) = 0;
end

imgmaxAll = imgImax + imgmax{1} + imgmax{2} +0 * imgmax{3} > 0;

if verbose
   figure(6); clf
   implottiling({imoverlay(imgIf, imgImax), imoverlay(imgC, imgImax, [1,1,1]); ...
              imoverlay(imgCf(:,:,1), imgmax{1}),  imoverlay(imgCf(:,:,2), imgmax{2}); ...
              imoverlay(imgCf(:,:,3), imgmax{3}),  imoverlay(imgC, imgmaxAll, [1,1,1])}')
end

           
%% Watershed

imgWs = imimposemin(iminvert(imgIf), imgmaxAll);
imgWs = watershed(imgWs);
imgWs = immask(imgWs, imgmask);
imgWs = imlabelseparate(imgWs);


if verbose
   figure(7); clf;
   implot(imgWs)

   %imglab = bwlabel(imgWs);

   figure(8); clf;
   implottiling(imoverlay(imoverlaylabel(imgC, impixelsurface(imgWs), false), imgmaxAll, [1,1,1]))
end
   
imgS = imgWs;      


%% Postprocess

imgSP = postProcessSegments(imgS, imgIf, 'intensity.median.min', 0.025, 'volume.min', 15, 'fillholes', false);
imgSP = imrelabel(imgSP);

if verbose
   imgSp1 = impixelsurface(imgSP);
   imgSp1 = imoverlaylabel(imgCf, imgSp1, false);
   figure(5); clf;
   implot(imgSp1);
end

imglab = imgSP;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Statistics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stats = imstatistics(imglab, {'PixelIdxList', 'Centroid', 'SurfacePixelIdxList'});

mode = 'MedianIntensity';

clear statsCh
for c = 1:nch
   statsCh{c} = imstatistics(imglab, stats, mode,  imgsSt{c});
end

%% Center Label

%imglabc = imlabelapplybw(imglab, @bwulterode);
cc = fix([stats.Centroid]);
imglabc = zeros(size(imglab));
ind = imsub2ind(size(imglab), cc');
imglabc(ind) = 1;
imglabc = imdilate(imglabc, strel('disk', 2));


%% Plot Result

figure(7); clf;
implottiling({imgCRaw, imgCf, imoverlaylabel(imgCRaw, imglabc > 0,false, 'color.map', [[0,0,0]; [1,1,1]])}');


%% Plot Seeding 

%subreg = {[600:900], [550:950], ':'};
subreg = {':', ':', ':'};
imgCfsub = imgCf(subreg{:});
imgCsub = imgC(subreg{:});
imglabsub = imglabc(subreg{:});

h = figure(7); clf;
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

%% Plot Filtered

h = figure(9); clf;
imgCsub = imgC(subreg{:});

implot(imgCsub);
axis off
xlabel([]); ylabel([])

%saveas(h, [datadir, dataname, '_Segmentation_Filtered.pdf'])


%% Plot Raw
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


%% Histogram of Intensities 


xy = [stats.Centroid]';
cm = {'r', 'g', 'b', 'm'};
cl = [0.5, 0.4, 0.4, 0.4];
%ct = [0.110, 0.0, 0.100]

figure(21); clf;
for c = 1:nch
   fi = [statsCh{c}.(mode)]';
   %fi = imclip(fi, 0, cl(c));
   
   figure(21);
   subplot(1,nch,c);
   hist(fi, 256)
   title(chlabel{c}); 
end



%% Flourescence in Space

xy = [stats.Centroid]';
cm = {'r', 'g', 'b', 'm'};
cl = [0.08, 0.08, 0.005, 0.4];
%ct = [0.110, 0.0, 0.100]

figure(21); clf;
for c = 1:nch
   fi = [statsCh{c}.(mode)]';
   %fi = imclip(fi, 0, cl(c));
   
   figure(21);
   subplot(nch,4,(c-1)*4+1);
   hist(fi, 256)
   title(chlabel{c});
 
   subplot(nch,4,(c-1)*4+2);
   %imcolormap(cm{c});
   colormap jet
   scatter(xy(:,1), xy(:,2), 10, fi, 'filled');
   xlim([0, size(imglab,1)]); ylim([0, size(imglab,2)]);
   axis equal

   %freezecolormap(gca)
   
   subplot(nch,4,(c-1)*4+3);
   %imcolormap(cm{c});
   scatter(xy(:,1), xy(:,2), 10, fi > cl(c), 'filled');
   xlim([0, size(imglab,1)]); ylim([0, size(imglab,2)]);
   %freezecolormap(gca)
   axis equal

   subplot(nch,4,(c-1)*4+4);
   implot(imgCf(:,:,c))
   
%    h = figure(22); clf
%    colormap jet
%    scatter(xy(:,1), xy(:,2), 10, fi, 'filled');
%    xlim([0, size(imglab,1)]); ylim([0, size(imglab,2)]);
%    title(chlabel{c}); 
%    colorbar('Ticks', [])
%    
   %saveas(h, [datadir dataname '_Quantification_' lab{c} '.pdf']);
end



%% Flourescence Expression

figure(22); clf;
fi = cell(1,3);
for c = 1:nch
   fi{c} = [statsCh{c}.(mode)]';
   %fi{c} = imclip(fi{c},0, cl(c));
   fi{c} = mat2gray(fi{c});
end
 
%pairs = {[1,2], [1,3], [1,4], [2, 3], [2,4], [3,4]};
pairs = {[1,2], [1,3], [2, 3]};

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


%% Color Settings

colorscheme = 'Z';

if strcmp(colorscheme, 'rgb')

   neuroBaseColor = {[1,0,0], [0,1,0], [0,0,1], [1, 1,1]};
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

else % zeeshan color scheme
   %neuroClassColor = {[0,0,0]; [0.6,0,0]; [0,0.6,0]; [0.33,0.33,0]; [0,0.33,0]; [0.33,0,0.33]; [0,0.33,0.33]; [0.5,0.5,0.5]};
   neuroClassColor = num2cell([0.285714, 0.285714, 0.285714;
      0.929412, 0.290196, 0.596078;
      0.462745, 0.392157, 0.67451;
      0.709804, 0.827451, 0.2;
      0.984314, 0.690196, 0.25098;
      0.976471, 0.929412, 0.196078;
      0.584314, 0.752941, 0.705882;
      0, 0.682353, 0.803922
      ], 2);
   
   neuroBaseColor = neuroClassColor([1,2,4]+1);
end


%neuroClassColor = num2cell(colorcube(max(neuroClassTotal(:))), 2);


%% Mask for non valid regoins

%%
% h = figure(5); clf;
% implot(imgC);
% p = impoly

%%

xp = [1028.536, 893.6171, 737.9414, 568.4279, 417.3649, 317.0405, 201.7252, 96.7883, -3.536, -3.536, 152.1396, 207.491, 337.7973, 455.4189, 581.1126, 712.5721, 809.4369, 940.8964, 1014.6982, 1017.0045, 1028.536]
yp = [707.9595, 687.2027, 683.7432, 680.2838, 686.0495, 702.1937, 720.6441, 737.9414, 742.5541, 1015.8514, 1010.0856, 974.3378, 980.1036, 981.2568, 976.6441, 975.491, 990.482, 991.6351, 1004.3198, 899.3829, 902.8423]
%p = impoly(gca, [xp', yp'])

%%
imgroi1 = poly2mask(yp',xp', size(imgI,1), size(imgI,2));

imgroi = {imgroi1, imgroi1, imgroi1};

size(imgroi)
size(imgI)

%%
%imgbad = ones(size(imgI));

figure(6); clf;
implottiling(imgroi')

%% Histograms

figure(10); clf
for c= 1:nch
   subplot(2,nch,c); 
   d = imgsSt{c};
   hist(d(:), 256);
   
   subplot(2, nch, c+nch)
   hist([statsCh{c}.(mode)], 256)
   
   title(chlabel{c})
end


%% Classify

%cth = {0.06, 0.1, 0.25, 0.2};
%cl = [0.08, 0.08, 0.08, 0.4];

cl = [0.08, 0.08, 0.02, 0.4];
cth = num2cell(cl);


clear neuroClass
for c = 1:nch
   xy = fix([statsCh{c}.Centroid]);
   xy(1,xy(1,:) > size(imgI,1)) = size(imgI,1); xy(1,xy(1,:) < 1) = 1;
   xy(2,xy(2,:) > size(imgI,2)) = size(imgI,2); xy(2,xy(2,:) < 1) = 1;
   
   isgood = zeros(1, size(xy,2));
   for i = 1:size(xy, 2)
      isgood(i) = imgroi{c}(xy(1,i), xy(2,i));
   end
  
   neuroClass(c,:) = double(and([statsCh{c}.(mode)] > cth{c}, isgood));
end
neuroClassTotal = fix(neuroClass(1,:) + 2 * neuroClass(2,:) + 4 * neuroClass(3,:))+1; % + 8 * neuroClass(4,:))+1;


imgClass1 = cell(1,nch);
for c = 1:nch
   neuroColor1 = {[0,0,0], neuroBaseColor{c}};
   neuroColor1 = neuroColor1(neuroClass(c,:)+1);
   

   R = imgsSt{c}; G = R; B = R;
   for i = 1:length(stats);
      if neuroClass(c,i) > 0
         R(stats(i).SurfacePixelIdxList) =  neuroColor1{i}(1);
         G(stats(i).SurfacePixelIdxList) =  neuroColor1{i}(2);
         B(stats(i).SurfacePixelIdxList) =  neuroColor1{i}(3);
      end
   end
   imgClass1{c} = cat(3, R, G, B);
end


figure(7); clf;

implottiling({imgsSt{1:nch}; imgClass1{:}}', 'titles', [chlabel(1:nch), chlabel(1:nch)]);

saveas(h, [datadir dataname '_Classification_Channels' '.pdf']);



%%

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

%saveas(h, [datadir dataname '_Classification_Channels' '.pdf']);
%saveas(h, [datadir dataname '_Classification_Image' '.pdf']);

for i = 1:8
   imsubplot(4,2,i)
   axis off
   xlabel('');
   ylabel('');
end
   

saveas(h, [datadir dataname '_Classification_Channels' '.pdf']);


%% Class in Space

xy = fix([statsCh{1}.Centroid]);
xy(1,xy(1,:) > size(imgI,1)) = size(imgI,1); xy(1,xy(1,:) < 1) = 1;
xy(2,xy(2,:) > size(imgI,2)) = size(imgI,2); xy(2,xy(2,:) < 1) = 1;

isgood = zeros(1, size(xy,2));
for i = 1:size(xy, 2)
   isgood(i) = imgroi{1}(xy(1,i), xy(2,i));
end

isgood = logical(isgood);

xy = [stats.Centroid]';
xy = xy(isgood, :)

h = figure(27)

cdat = (cell2mat({[0,0,0], neuroClassColor{:}}'))
   
scatter(xy(:,1), xy(:,2), 4, cdat(neuroClassTotal(isgood) + 1, :), 'filled');
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
   
   nstat(i) = sum(neuroClassTotal(isgood) == i);
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



%% Distance from Center

%figure(7); clf;
%implottiling({imoverlaylabel(imgCRaw, imglabc > 0,false, 'color.map', [[0,0,0]; [1,1,1]])}');
%p = impoly

%%
xb = [1021.6556, 790.3074, 674.6333, 540.9444, 409.1519, 311.4926, 251.7593, 187.2852, 106.6926, 55.4926, 1.4481];
yb = [692.1741, 675.1074, 681.7444, 677.9519, 684.5889, 699.7593, 705.4481, 714.9296, 733.8926, 742.4259, 743.3741];

%% use higher sampling

k = 20; n = length(xb);
xbb = interp1(1:n, xb, 1:1/k:n);
ybb = interp1(1:n, yb, 1:1/k:n);

figure(17); clf;
imgA = imgCf;
for i = 1:3
   imgAA = imgA(:,:,i);
   imgAA(~logical(imgroi{1})) = 0;
   imgA(:,:,i) = imgAA;
end
implot(imgA)

%% 

% get cells in roi

xy = [stats.Centroid];

clear neuroClassA
for c = 1:nch
   xy = fix([statsCh{c}.Centroid]);
   xy(1,xy(1,:) > size(imgI,1)) = size(imgI,1); xy(1,xy(1,:) < 1) = 1;
   xy(2,xy(2,:) > size(imgI,2)) = size(imgI,2); xy(2,xy(2,:) < 1) = 1;
   
   isgoodA = zeros(1, size(xy,2));
   for i = 1:size(xy, 2)
      isgoodA(i) = imgroi{c}(xy(1,i), xy(2,i));
   end
  
   neuroClassA(c,:) = double(and(neuroClass(c,:), isgoodA));
end

isgoodA2 = any(neuroClassA, 1);

xy = xy(:, isgoodA2);

dm = distanceMatrix(xy, [xbb; ybb]);
d = min(dm,[],2);

figure(18); clf
hist(d,256)

ncl = logical(neuroClassA(:,isgoodA2));

figure(8); clf
for i = 1:3
   subplot(1,3,i)
   hist(d(ncl(i,:)), 1:10:1300)
   title(chlabel{i})
   ylim([0, 200])
   xlim([0, 500])
end

%save numbers

dmat = [ncl; d'];

save([datadir dataname '_Quantification_Distance.mat'], 'dmat')



%% Plot 3D + dist


xyz = zeros(length(statsCh{c}),3);
for c = 1:nch
   xyz(:,c) = [statsCh{c}.(mode)]';
   %fi = imclip(fi, 0, cl(c));
end
   
 
figure(20); clf
plot3(xyz(isgoodA2,1),xyz(isgoodA2,2), xyz(isgoodA2,3), '.')
 
figure(21); clf
scatter3(xyz(isgoodA2,1),xyz(isgoodA2,2), xyz(isgoodA2,3), 5, d)

figure(22); clf
scatter(xyz(isgoodA2,1),xyz(isgoodA2,2), 5, xyz(isgoodA2,3))






%% ROI

h = figure(7)
implottiling({imgC; immask(imgC, imgroi{1})})

%saveas(h, [datadir dataname '_ROI' '.pdf']);


%% Save 

save( [datadir dataname '_workspace.mat'])
