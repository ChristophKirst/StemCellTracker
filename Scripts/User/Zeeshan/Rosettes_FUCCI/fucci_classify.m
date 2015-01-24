function [stats, varargout] = fucci_classify(img, imglab, imgseg, varargin)

if nargin > 1
   verbose = varargin{1};
else
   verbose = false;
end

%% Classify Labels

stats = imstatistics(imglab, 'PixelIdxList');

mode = 'MeanIntensity';

statsR = imstatistics(imglab, stats, mode, img(:,:,1));
statsG = imstatistics(imglab, stats, mode, img(:,:,2));
statsB = imstatistics(imglab, stats, mode, imgseg);

[~, class] = max(cat(2, 1.5 * [statsR.(mode)]', [statsG.(mode)]'), [], 2);


{max([statsB.(mode)]), min([statsB.(mode)])}
class([statsB.(mode)] < 0.45) = 3;

cls = num2cell(class);
[stats.class] = cls{:};

if verbose 
   figure(99); clf
   subplot(1,3,1)
   hist([statsR.(mode)])
   subplot(1,3,2)
   hist([statsG.(mode)])
   subplot(1,3,3)
   hist(class)

   % plot classification 
   pixl = {stats.PixelIdxList};

   imgcls = imglab;
   for i = 1:length(pixl)
      imgcls(pixl{i}) = class(i);
   end

   figure(100); clf;
   imgp = imoverlay(img(:,:,:), imgcls == 1, 'r');
   imgp = imoverlay(imgp, imgcls == 2, 'g');
   imgp = imoverlay(imgp, imgcls == 3, 'b');
   
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

