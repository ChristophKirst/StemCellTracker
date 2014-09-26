function imglab = findROIsByOpening(img, varargin)
%
% imglab = findROIsByOpening(img, param)
%
% descritption:
%     finds shapes in an image img using morphological opening and thresholding
%
% input:
%    img    image
%    param  paramteer struct
%           .threshold    threshold ([] = thresholdFirstMin(img))
%           .resize       resize factor before detection of objects ([]=false)
%           .filtersize   size of Gaussian filter ([] = no filter) 
%           .strel        structure element for morphological closing ([] = strel('disk', 20))
%           .output       'ROIs' or 'label' ('label')
%           .plot         plot result
%           
% output:
%    imglab labeled image or ROIMask objects

param = parseParameter(varargin);

rs = getParameter(param, 'resize', []);
if ~isempty(rs)
   imglab = imresize(img, rs);
else
   imglab= img;
end

fs = getParameter(param, 'filtersize', []);
if ~isempty(fs)
   imglab = gaussianFilter(imglab, fs);
end

se= getParameter(param, 'strel', strel('disk', 20));
if isnumeric(se)
   se = strel('disk', se);
end
imglab = imclose(imglab, se);


th = getParameter(param, 'threshold', []);
if isempty(th)
   th = thresholdFirstMin(mat2gray(img));
end

imglab = imglab > th;

if ~isempty(rs)
   imglab = imresize(imglab, 1/rs);
end

imglab = bwlabeln(imglab);

if getParameter(param, 'plot', false)
   implot(imoverlaylabel(img, imglab));
end


if isequal(getParameter(param, 'output', 'label'), 'ROIs')
   ids = regionprops(imglab, 'PixelIdxList');
   
   if ~isempty(ids)
      rois(length(ids)) = ROIMask;
   
      for i = 1:length(ids);
         imgm = zeros(size(imglab));
         imgm(ids(i).PixelIdxList) = 1;
         rois(i).mask = imgm;
      end
      imglab = rois;
   else
      imglab = [];
   end
end