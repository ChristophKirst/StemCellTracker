function [imglab, imgI, imgmask] = brainslice_segment(imgCf, varargin)

if nargin > 1
   verbose = varargin{1};
else
   verbose = false;
   %verbose = true;
end

%imgCfO = imgCf;
%imgCf = imgCfO;

%%

%imgCff = cellfunc(@(x) filterMedian(imgC(:,:,x),3), {1,2,3});
%imgCff = cat(3, imgCff{:});

%% mask

imgI = imgCf(:,:,1) + imgCf(:,:,2) + imgCf(:,:,3);
imgI = mat2gray(imgI);

%imgmask = img > 0.125;
%imgmask = img > 0.25;
imgmask = imgI > 0.075;
%imgmask = imopen(imgmask, strel('disk', 3));
%imgmask = postProcessSegments(bwlabeln(imgmask), 'volume.min', 50) > 0;

if verbose
   %max(img(:))
   
   figure(21); clf;
   set(gcf, 'Name', ['Masking'])
   implottiling({imoverlaylabel(imgI, double(imgmask));  imgmask})
end
   
%% Prepare Seeding Image

% %imgCff = cellfunc(@(x) filterMedian(imgCf(:,:,x),3), {1,2,3});
% imgCff = cellfunc(@(x) imgCf(:,:,x), {1,2,3});;
% imgCff = cat(3, imgCff{:});
% 
% imgG = mat2gray(imgradient(imgCff(:,:,1) + imgradient(imgCff(:,:,2)) + imgradient(imgCff(:,:,3))));
% imgG = imclip(imgG, 0, 0.5);
% imgT = imgI - 0.5 * imgG;
% imgT(imgT < 0) = 0;
% 
% figure(2); clf;
% implottiling({imgT; imgG});


%% Seeds

% imglab = zeros(size(imgI));
% for c = 1:3
%    imgf = imgCff(:,:,c) -0.2 * imgG;
%    imgf(imgf < 0) = 0;
%    %imgf = imgT;
%    %imgf = filterLoG(max(imgf(:)) - imgf, [8,8], []);
%    %imgf = filterDisk(imgf, 6, 1, 1, 0);
%    imgf = filterSphere(imgf, 7);
%    imgf = mat2gray(imgf);
% 
% 
%    % h-max detection (only local maxima with height > hmax are considered as maxima
%    %imgmax = imextendedmax(imgf,  0.00005);
%    imgmax = imextendedmax(imgf,  0.005);
% 
%    % constrain to maxima within mask
%    imgmax = immask(imgmax, imgmask);
% 
%    % Combine nearby points
%    %imgmax = imdilate(imgmax, strel('disk', 4));
%    imgmax = imdilate(imgmax, strel('disk', 1));
% 
%    % fill holes - combination of nearby points can lead to holes
%    imgmax = imfill(imgmax,'holes');
%    % shrink to single points - extended maxima usually give better segmentation results
%    % imgmax = bwmorph(imgmax,'shrink',inf);           
% 
%    % plot the results.
%    if verbose  
%       figure(22); clf
%       set(gcf, 'Name', ['Seeding'])
%       implottiling({imoverlaylabel(imgCf, imgmax, false); imoverlay(imgf, imgmax)});
%    end
% 
%    imglab = imgmax + imglab;
% end
% 
% imglab = imglab > 0;
% imgmax = imfill(imgmax,'holes');
% figure(19)
% implottiling({imoverlaylabel(imgCf, 1* imglab +0 * imgmask , false); imoverlaylabel(imgI, imglab, false)})
% 
% 
% imgl = bwlabeln(imglab);
% max(imgl(:))

%% SLIC Segmentation

% number of pixels: typical cell 9x9
npxl = fix(1.1 * numel(imgI) / (8*7))

%segment
imgS = segmentBySLIC(imgCf, 'superpixel', npxl, 'compactness', 10);

if verbose
   imgSp = impixelsurface(imgS);
   figure(5); clf;
   implot(imoverlaylabel(imgCf, imgSp, false));
end


%% Postprocess
clc
max(imgS(:))

imgSP = imgS;

stats = imstatistics(imgSP, {'MinIntensity'}, imgI);
figure(7); clf; 
hist([stats.MinIntensity], 256)
%hist(imgI(:), 256)


%%
imgSP = postProcessSegments(imgS, imgI, 'intensity.min', 0.085, 'volume.min', 15, 'fillholes', false);

if verbose
   imgSp = impixelsurface(imgSP);
   figure(5); clf;
   implot(imoverlaylabel(imgCf, imgSp, false));
end


%%


imgSP2 = immask(imgSP, imgmask);
imgSP2 = imlabelapplybw(imgSP2, @(x) imopen(x, strel('disk', 2)));
imgSP2 = imlabelseparate(imgSP2);

% stats = imstatistics(imgSP, {'MinIntensity'}, imgI);
% figure(7); clf; 
% hist([stats.MinIntensity], 56)
%hist(imgI(:), 256)

imgSP2 = postProcessSegments(imgSP2, imgI, 'intensity.min', 0.085, 'volume.min', 15, 'fillholes', false);
%imgSP = postProcessSegments(imgSP, 'volume.min', 7, 'fillholes', false);
imgSP2 = imrelabel(imgSP2);
max(imgSP2(:))

if verbose
   imgSp = impixelsurface(imgSP2);
   figure(5); clf;
   implot(imoverlaylabel(imgCf, imgSp, false));
end

%%

imgSPP = imgSP2;

stats = imstatistics(imgSPP, {'Volume', 'PixelIdxList', 'MedianIntensity', 'Perimeter', 'Extent', 'FilledArea'}, imgI);
% 
% %figure(78); clf; hist([stats.MinIntensity], 256);
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
   figure(6); clf;
   %implottiling({imoverlaylabel(imgCf, impixelsurface(imgSPP), false); imoverlaylabel(imgCf, imgSP, false)});
   
   statsF = imstatistics(imgSPP, {'PixelIdxList', 'Centroid'});
   
%    mode = 'MedianIntensity';
%    statsR = imstatistics(imgSP, stats, mode,  imgCf(:,:,1));
%    statsG = imstatistics(imgSP, stats, mode,  imgCf(:,:,2));
%    statsB = imstatistics(imgSP, stats, mode,  imgCf(:,:,3));
   
   mode = 'MedianIntensity';
   statsR = imstatistics(imgSPP, statsF, mode,  imgCf(:,:,1));
   statsG = imstatistics(imgSPP, statsF, mode,  imgCf(:,:,2));
   statsB = imstatistics(imgSPP, statsF, mode,  imgCf(:,:,3));
   
   si = size(imgI);
   R = zeros(si); G = R; B = R;
   for i = 1:length(stats);
      R(statsF(i).PixelIdxList) =  statsR(i).(mode);
      G(statsF(i).PixelIdxList) =  statsG(i).(mode);
      B(statsF(i).PixelIdxList) =  statsB(i).(mode);
   end
   imgCC = cat(3, R, G, B);
   
   figure(7); clf;
   implottiling({imgCf; imgCC});
   R = zeros(si); G = R; B = R;
   for i = 1:length(stats);
      R(statsF(i).PixelIdxList) =  statsR(i).(mode);
      G(statsF(i).PixelIdxList) =  statsG(i).(mode);
      B(statsF(i).PixelIdxList) =  statsB(i).(mode);
   end
   imgCC = cat(3, R, G, B);
   
   figure(7); clf;
   implottiling({imgCf; imgCC});
end

imglab = imgSPP;

%% Label

% imglab = bwlabeln(img);
% 
% if verbose
%    figure(23); clf;
%    set(gcf, 'Name', ['Labeling'])
%    implot(imoverlaylabel(imgCf, imglab))
% end
