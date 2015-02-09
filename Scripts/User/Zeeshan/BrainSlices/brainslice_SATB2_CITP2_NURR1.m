%%%%%%%%%%%%%%%%%%%%
%%% Brain Slices %%%
%%%%%%%%%%%%%%%%%%%%

clear all
close all

initialize()
bfinitialize
initializeParallelProcessing(12)

addpath('./Scripts/User/Zeeshan/BrainSlices/');

datadir = '/data/Science/Projects/StemCells/Experiment/BrainSections/';
dataname = 'Sample B Slide 33 SATB2 CUX2 NURR1 CTIP2';

lab = {'SATB2', 'CTIP2', 'NURR1'}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%is = ImageSourceBF('/data/Server/smb/upload/Brain sections/Sample B Slide 33 SATB2 CUX2 NURR1 CTIP2.lsm');
is = ImageSourceBF([datadir  dataname '.lsm']);
clc
is.printInfo

%%

is.setReshape('S', 'Uv', [16, 22]);
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
region = struct('U', 1+[9:10], 'V', 1+[10:11], 'C', 1);

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
%% Segment Cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% nuclear marker
is.resetRange();

%region2 = {':',':'};
region2 = {500:1900, 450:1850};

imgsCTIP2 = is.cell(region, 'C', 3);
imgsSATB2 = is.cell(region, 'C', 4);
imgsNURR1  = is.cell(region, 'C', 5);

% corect illumination and background

imgsAll = {imgsSATB2, imgsCTIP2, imgsNURR1};

parfor i = 1:3
   %imgsAll{i} = cellfunc(@(x) filterAutoContrast(x/max(x(:))), imgsAll{i});
   imgsAll{i} = cellfunc(@(x) x - imopen(x, strel('disk', 10)), imgsAll{i});
   imgsAll{i} = cellfunc(@(x) x - imopen(x, strel('disk', 150)), imgsAll{i});
end


%%
imgSt =cell(1,3);
th = [200, 50, 100];
for i = 1:3
   imgsS = stitchImages(imgsAll{i}, sh, 'method', stmeth);
   imgsS = imgsS(region2{:});
   %imgsS = filterAutoContrast(imgsS/max(imgsS(:)));
   
   figure(16);
   subplot(1,3,i); 
   hist(imgsS(:), 256);
   
   imgsS(imgsS < th(i)) = 0;
   imgSt{i} = imclip(imgsS, 0, 2000);
   
end

figure(1); clf; colormap jet
implottiling(imgSt')

imgC = cat(3, imgSt{:});
%imgC = imclip(imgC, 0, 2500);
imgC = imgC / max(imgC(:));

figure(2)
implot(imgC)

%%

for c = 1:3
   figure(5); 
   subplot(3,1,c)
   dd = imgC(:,:,c);
   hist(dd(:), 256);
end

%% filter image

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

%% segment
clc
imgCf = imgBM;
[imglab, imgI, imgmask] = brainslice_segment(imgCf, false);

%%
stats = imstatistics(imglab, {'PixelIdxList', 'Centroid'});

mode = 'MedianIntensity';
statsR = imstatistics(imglab, stats, mode,  imgCf(:,:,1));
statsG = imstatistics(imglab, stats, mode,  imgCf(:,:,2));
statsB = imstatistics(imglab, stats, mode,  imgCf(:,:,3));
  
si = size(imgCf(:,:,1));
R = zeros(si); G = R; B = R;
R2 = imgCf(:,:,1); G2 = imgCf(:,:,2); B2 = imgCf(:,:,3);

for i = 1:length(stats);
   R(stats(i).PixelIdxList) =  statsR(i).(mode);
   G(stats(i).PixelIdxList) =  statsG(i).(mode);
   B(stats(i).PixelIdxList) =  statsB(i).(mode);
   
   R2(stats(i).PixelIdxList) =  0*statsR(i).(mode);
   G2(stats(i).PixelIdxList) =  0*statsG(i).(mode);
   B2(stats(i).PixelIdxList) =  0*statsB(i).(mode);
end
imgCC = cat(3, R, G, B);
imgCC2 = cat(3, R2, G2, B2);
statsCC = {statsR, statsG, statsB};

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
imglabsub = imglabc(subreg{:});

implot(imoverlaylabel(imgCfsub, imglabsub > 0,false, 'color.map', [[0,0,0]; 0.2*[1,1,1]]));
axis off
xlabel([]); ylabel([])

saveas(h, [datadir, dataname, '_Segmentation_Seeds.pdf'])


h = figure(8); clf;
imgCCsub = imgCC(subreg{:});

implot(imgCCsub);
axis off
xlabel([]); ylabel([])

saveas(h, [datadir, dataname, '_Segmentation_Segments.pdf'])

%%
h = figure(9); clf;
imgCCsub = imgCf(subreg{:});

implot(imgCCsub);
axis off
xlabel([]); ylabel([])

saveas(h, [datadir, dataname, '_Segmentation_Filtered.pdf'])


%%
h = figure(10); clf;
imgCCsub = imgCf(subreg{:});

implot(imgCCsub);
axis off
xlabel([]); ylabel([])

saveas(h, [datadir, dataname, '_Segmentation_Raw.pdf'])

%%
h = figure(11); clf;

%implot(filterAutoContrast(imgCf));
implot(imgCf);
axis off
xlabel([]); ylabel([])

saveas(h, [datadir, dataname, '_Raw.pdf'])


%% individual channels

for c = 1:3

   h = figure(11+c); clf;

   %implot(filterAutoContrast(imgCf));
   imgCfc = imgCf;
   imgCfc(:,:,setdiff([1,2,3], c)) = 0;
   implot(imgCfc);
   axis off
   xlabel([]); ylabel([])

   saveas(h, [datadir, dataname, '_Raw_' lab{c} '.pdf'])
end

%% 

figure(7); clf;
implottiling({imgCf, imgCC; imgCC2, imoverlaylabel(imgCf, imglabc> 0,false, 'color.map', [[0,0,0]; [1,0,0]])});


%% Flourescence in Space

xy = [stats.Centroid]';
cm = {'r', 'g', 'b'};
cl = [0.5, 0.6, 0.4];
%ct = [0.110, 0.0, 0.100]

figure(21); clf;
for c = 1:3
   fi = [statsCC{c}.(mode)]';
   fi = imclip(fi, 0, cl(c));
   
   figure(21);
   subplot(2,3,c+3);
   hist(fi, 256)
 
   subplot(2,3,c);
   %imcolormap(cm{c});
   colormap jet
   scatter(xy(:,1), xy(:,2), 10, fi, 'filled');
   xlim([0, size(imglab,1)]); ylim([0, size(imglab,2)]);
   title(lab{c});
   %freezecolormap(gca)
   
   h = figure(22); clf
   colormap jet
   scatter(xy(:,1), xy(:,2), 10, fi, 'filled');
   xlim([0, size(imglab,1)]); ylim([0, size(imglab,2)]);
   title(lab{c}); 
   colorbar('Ticks', [])
   
   saveas(h, [datadir dataname '_Quantification_' lab{c} '.pdf']);
end



%% Flourescence Expression

figure(22); clf;
fi = cell(1,3);
for c = 1:3
   fi{c} = [statsCC{c}.(mode)]';
   fi{c} = imclip(fi{c},0, cl(c));
   fi{c} = mat2gray(fi{c});
end
 
pairs = {[1,2], [1,3], [2,3]};

np = length(pairs);

for n = 1:np
   subplot(1, np,n)
   %fi = imclip(fi, 0, cl(c));
   %scatter(fi{pairs{n}(1)}, fi{pairs{n}(2)}, 10, 'b', 'filled');
   scattercloud(fi{pairs{n}(1)}, fi{pairs{n}(2)});
   xlim([0,1]); ylim([0,1]);
   xlabel(lab{pairs{n}(1)}); ylabel(lab{pairs{n}(2)});
   %freezecolormap(gca)
end


for n = 1:np
   h = figure(50+n); clf;
   %fi = imclip(fi, 0, cl(c));
   %scatter(fi{pairs{n}(1)}, fi{pairs{n}(2)}, 10, 'b', 'filled');
   scattercloud(fi{pairs{n}(1)}, fi{pairs{n}(2)});
   xlim([0,1]); ylim([0,1]);
   xlabel(lab{pairs{n}(1)}); ylabel(lab{pairs{n}(2)});
   
   saveas(h, [datadir dataname '_Quantification_Scatter_' lab{pairs{n}(1)} ,'_' lab{pairs{n}(2)} '.pdf']);
   %freezecolormap(gca)
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Classify Cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% histograms
figure(10); clf
for c= 1:3
   subplot(2,3,c); 
   d = imgCf(:,:,c);
   hist(d(:), 256);
end

for c= 1:3
   subplot(2,3,c+3); 
   d = imgCC(:,:,c);
   hist(d(:), 128);
end


%% classify

cth = {0.1, 0.15, 0.1};

clear neuroClass
for c = 1:3
   neuroClass(c,:) = double([statsCC{c}.(mode)] > cth{c});
end
neuroClass = fix(neuroClass(1,:) + 2 * neuroClass(2,:) + 4 * neuroClass(3,:))+1;

neuroClassColor = {[0,0,0]; [0.6,0,0]; [0,0.6,0]; [0.33,0.33,0]; [0,0.33,0]; [0.33,0,0.33]; [0,0.33,0.33]; [0.5,0.5,0.5]};
neuroColor = reshape([neuroClassColor{neuroClass}], 3,[])';
size(neuroColor)

R = zeros(si); G = R; B = R;
for i = 1:length(stats);
   R(stats(i).PixelIdxList) =  neuroColor(i,1);
   G(stats(i).PixelIdxList) =  neuroColor(i,2);
   B(stats(i).PixelIdxList) =  neuroColor(i,3);
end
imgClass = cat(3, R, G, B);


figure(7); clf;
implottiling({imgCf; imgClass});

saveas(h, [datadir dataname '_Classification_Image' '.pdf']);



%% histogram of cell classes
 
ncls = length(neuroClassColor);

nstat = zeros(1,ncls);

for i = 1:ncls
   id = num2cell(dec2bin(i-1));
   id = id(end:-1:1);
   for k = length(id)+1:3
      id{k} = '0';
   end

   clslab{i} = '';
   for c = 1:3
      if id{c} == '1'
         clslab{i} = [clslab{i}, '_', lab{c}];
      end
   end
   clslab{i} = clslab{i}(2:end);
   if isempty(clslab{i})
      clslab{i} = 'None';
   end
   
   nstat(i) = sum(neuroClass == i);
end

nc = num2cell(nstat);
tb = table(nc{:}, 'VariableNames', clslab)

%% save numbers
writetable(tb, [datadir dataname '_Counts.txt'])







