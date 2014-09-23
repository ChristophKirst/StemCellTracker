function imglab = findROIsByOpening(img, varargin)
%
% imglab = findShapes(img, param)
%
% descritption:
%     finds shapes in aan image img
%
% input:
%    img    image
%    param  paramteer struct
%           .threshold    threshold ([] = thresholdFirstMin(img))
%           .resize       resize factor before detection of objects ([]=false)
%           .filtersize   size of Gaussian filter ([] = no filter) 
%           .strel        structure element for morphological closing ([] = strel('disk', 20))
%           

param = parseParameter(varargin);

rs = getParameter(param, 'resize', []);
if ~isempty(rs)
   img = imresize(img, rs);
end

fs = getParameter(param, 'filtersize', []);
if ~isempty(fs)
   img = gaussianFilter(img, fs);
end

se= getParameter(param, 'strel', strel('disk', 20));
img = imclose(img, se);


th = getParameter(param, 'threshold', []);
if isempty(th)
   th = thresholdFirstMin(mat2gray(imgro));
end

img = img > th;

if ~isempty(rs)
   img = imresize(img, 1/rs);
end

imglab = bwlabeln(img);


end