function [stats, varargout] = fucci_classify(img, imglab, imgseg, varargin)

if nargin > 1
   verbose = varargin{1};
else
   verbose = false;
end

%% Classify Labels

stats = imstatistics(imglab, 'PixelIdxList');

mode = 'MeanIntensity';

se = 50;
imgR = img(:,:,1);
imgR = imgR - imopen(imgR, strel('disk', se));
imgG = img(:,:,2);
imgG = imgG - imopen(imgG, strel('disk', se));
imgB = img(:,:,3);
imgB = imgB - imopen(imgB, strel('disk', se));

statsR = imstatistics(imglab, stats, mode, imgR);
statsG = imstatistics(imglab, stats, mode, imgG);
statsB = imstatistics(imglab, stats, mode, imgB);

%[~, class] = max(cat(2, 1.5 * [statsR.(mode)]', [statsG.(mode)]'), [], 2);
[~, class] = max(cat(2, 1.5 * [statsR.(mode)]', [statsG.(mode)]'), [], 2);

%{max([statsB.(mode)]), min([statsB.(mode)])}
class([statsR.(mode)] + [statsG.(mode)] < 0.1) = 3;

cls = num2cell(class);
[stats.class] = cls{:};

if verbose || nargout > 1
   % plot classification
   pixl = {stats.PixelIdxList};
   imgcls = imglab;
   for i = 1:length(pixl)
      imgcls(pixl{i}) = class(i);
   end
   imgp = imoverlay(cat(3, imgR, imgG, imgB), imgcls == 1, 'r');
   imgp = imoverlay(imgp, imgcls == 2, 'g');
   imgp = imoverlay(imgp, imgcls == 3, 'b');
end

if verbose 
   figure(99); clf
   subplot(1,4,1)
   hist([statsR.(mode)])
   xlabel('R')
   
   subplot(1,4,2)
   hist(1.5 * [statsR.(mode)])
   xlabel('1.5 R')
   xlim([0,1]);
   
   subplot(1,4,3)
   hist([statsG.(mode)])
   xlabel('G')
   subplot(1,4,4)
   hist(class)
   xlabel('B')
   
   figure(100); clf;   
   implot(imgp);
end

statsR = num2cell([statsR.(mode)]);
[stats.RIntensity] = statsR{:};

statsG = num2cell([statsG.(mode)]);
[stats.GIntensity] = statsG{:};

statsB = num2cell([statsB.(mode)]);
[stats.BIntensity] = statsB{:};

if nargout > 1
   varargout{1} = imgp;
end

if nargout > 2
   varargout{2} = cat(3, imgR, imgG, imgB);
end

