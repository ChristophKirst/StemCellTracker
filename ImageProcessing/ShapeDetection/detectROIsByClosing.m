function roi = detectROIsByClosing(imgs, varargin)
%
% roi = detectROIsByClosing(img, param)
% roi = detectROIsByClosing(imgs, shifts, param)
% roi = detectROIsByClosing(is, param)
%
% descritption:
%     finds shapes in an image img using morphological opening and thresholding
%
% input:
%    img    image as numeric array
%    imgs   images as cell arrays
%    is     ImageSource or Alignment class
%    param  paramteer struct
%           .threshold    threshold ([] = thresholdFirstMin(img))
%           .scale        rescale factor before detection of objects ([] = 1)
%           .filtersize   size of Gaussian filter ([] = no filter) 
%           .strel        structure element for morphological closing ([] = strel('disk', 20))
%           .radius       probe radius passed to detectAlphaVolume (100)
%           .preview      stitch scaled images and save in preview of the ImageSource (true if input is ImageSource)
%           
% output:
%    roi    region of interest
%
% note:
%    for Alignment classes rois coords are not shifted by the position of the Alignment

%prepare
if isnumeric(imgs)
   shifts = {zeros(1, ndims(imgs))};
   imgs = {imgs};
   
elseif iscell(imgs)
   if nargin < 2 || ~iscell(varargin{1})
      error('detectROIsByClosing: expects image shifts as second argument!')
   end
   
   shifts = varargin{1};
   varargin = varargin(2:end);
end

param = parseParameter(varargin);

rs = getParameter(param, 'scale', []);
fs = getParameter(param, 'filtersize', []);
se = getParameter(param, 'strel', strel('disk', 20));
if isnumeric(se)
   se = strel('disk', se);
end

th = getParameter(param, 'threshold', []);

% find rois
if isa(imgs, 'Alignment')
   
   pre = getParameter(param, 'preview', true);

   n = imgs.nNodes;
   pks = cell(1, n);
   source = imgs.source;
   is = source.dataSize;
   nds = imgs.nodes;
   rsf = [];

   if pre
      if isempty(rs)
         rs = imgs.asource.previewScale;
      else
         imgs.asource.setPreviewScale(rs);
      end
   end

   parfor i = 1:n
      if pre
         img = imgs.asource.preview(nds(i)); %#ok<PFBNS>
         img = img{1};
         rsf = (size(img)-1) ./ (is-1); %#ok<PFTIN>
      elseif ~isempty(rs)
         img = source.dataResample(rs, nds(i)); %#ok<PFBNS>
         rsf = (size(img)-1) ./ (is-1); 
      else
         img = source.data(nds(i));
      end
      
      pks{i} = detectPeaksByClosing(img, 'filtersize', fs, 'strel', se, 'threshold', th);
      %pks{i}
      
      if ~isempty(rs) 
%          size(pks{i})
%          size(repmat(rsf(:), 1, size(pks{i},2)))
         pks{i} = round( (pks{i}-1)./repmat(rsf(:), 1, size(pks{i},2)) + 1);
      end
   end

   pks = stitchPoints(pks, imgs.imageShifts, imgs.imageSizes);

else
   n = numel(imgs);
   pks = cell(1, n);
   
   parfor i = 1:n
      if ~isempty(rs)
         img = imresize(imgs{i}, rs);
      else
         img = imgs{i};
      end
      pks{i} = detectPeaksByClosing(img, 'filtersize', fs, 'strel', se, 'threshold', th);
      
      if ~isempty(rs)
         pks{i} = 1/rs * (pks{i}-1) + 1;
      end
   end
   
   pks = stitchPoints(pks, shifts, cellfunc(@size, imgs), param);
end

if size(pks,2) < 3
   roi = [];
   return
end
   
% get the polygons from the points
shps = points2shapes(pks, param);
nr = length(shps);

if nr == 0
   roi = [];
   return
end

%convert polygons to ROIs 
roi(nr) = ROIPolygon;
for i = 1:length(shps)
   roi(i) = ROIPolygon(shps{i});
   %roi{i} = roi{i}.shift(origin -1);
end

end
