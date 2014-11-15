function [img, scale] = stitchPreview(celldata, varargin)
%
% [img, scale] = stitchPreview(celldata, varargin)
%
% description:
%   stitches cells in a resmaples way using approximate overlaps
%
% input:
%   celldata    cell array of cell data
%   param       parameter struct with entires
%               .scale    scale for resmapling ([] = 0.1)
%               .shifts   shifts of the images in celldata ([] = infer from overlap and cell grid)
%               .size     size of the final image, everything outsie this size will be cut ([] = automatic full size)
%               .preview  if ture caches preview images in ImageSource (true)
%
% output:
%   img         stitched image
%   scale       (optional) scale used for down sampling
%
% notes:
%   if input is an image source with preview images of the same scale these are used without resampling
%   if the scale is diffrent or no chached images exists they are chached if 'preview' option is true


param = parseParameter(varargin);

scale = getParameter(param, 'scale',   []);
if isempty(scale)
   if isa(celldata, 'ImageSource')
      scale = celldata.previewScale;
   else
      scale = 0.1;
   end
end
%scale

sh    = getParameter(param, 'shifts',  []);
siz   = getParameter(param, 'size',    []);
pre   = getParameter(param, 'preview', true);

if isempty(celldata)
   if isempty(siz)
      img = [];
   else
      img = zeros(ceil(siz .* scale));
   end
   return
end

overlap = getParameter(param, 'overlap', []);
if isempty(overlap)
   overlap = 0;
end
overlap = padright(overlap, 2, overlap);

% rescale
if isa(celldata, 'Alignment')

   if isempty(siz)
      siz = ceil(celldata.dataSize .* scale);
   end
   if isempty(sh)
      sh = celldata.imageShifts;
      sh = cellfunc(@(x) ceil(x .* scale), sh);
   end
   
   if pre
      %scale
      celldata.setPreviewScale(scale);
      cdat = celldata.asource.preview(celldata.nodes);
      %size(cdat)
      %size(cdat{1})
      
   else
      % read sequentially
      nodes = celldata.nodes;
      cdat = cell(1,length(nodes));
      parfor i = 1:numel(cdat)
         cdat{i} = celldata.asource.dataResample(scale, nodes(i)); %#ok<PFBNS>
      end
   end
   
   celldata = cdat;
   
elseif isa(celldata, 'ImageSource')

   if pre
      celldata.setPreviewScale(scale);
      cdat = celldata.preview;
   else
      % read sequentially 
      cdat = cell(imfrmtAllocateSize(celldata.cellSize));
   
      parfor i = 1:numel(cdat)
         cdat{i} = celldata.dataResample(scale, i); %#ok<PFBNS>
      end
   end
   
   cs = celldata.cellSize;
   celldata = cdat;
   celldata = squeeze(reshape(celldata, cs));
   
else
   celldata = cellfunc(@(x) imresize(x, scale), celldata);
end

overlap  = ceil(scale .* overlap);

% shifts
if isempty(sh)
   isize = size(celldata{1}); % assume all images are same size
   sh = imgrid(size(celldata));
 
   for i = 1:length(sh)
      sh{i} = sh{i}(:);
   end
   sh = [sh{:}] - 1;
   sh = sh .* repmat(isize - overlap, size(sh,1), 1) + 1;

   sh = num2cell(sh,2);
end

if isempty(siz)
   [sh, siz] = absoluteShiftsAndSize(sh, cellfunc(@size, celldata));
end


%var2char(pos)

img = stitchImages(celldata, sh, 'method', 'Overwrite', param, 'size', siz);


if getParameter(param, 'lines', false)
   sh = cell2mat(sh);
   posx = unique(sh(:,1)); 
   
   mval = max(img(:));
   for p = posx(:)'
      img = impixelline(img, [p,1], [p, siz(2)], mval);
   end
   
   posy = unique(sh(:,2)); 
   for p = posy(:)'
      img = impixelline(img, [1, p], [siz(1), p], mval);
   end
end


end