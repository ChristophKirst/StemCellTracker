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
%
% output:
%   img         stitched image
%   scale       (optional) scale used for down sampling


param = parseParameter(varargin);

scale = getParameter(param, 'scale',    0.1);
sh    = getParameter(param, 'shifts', []);
siz   = getParameter(param, 'size',     []);

if isempty(celldata)
   if isempty(siz)
      img = [];
   else
      img = zeros(ceil(siz .* scale));
   end
   return
end

overlap = getParameter(param, 'overlap', []);
overlap = padright(overlap, 2, 0);

% rescale
if isa(celldata, 'Alignment')
   if isempty(siz)
      siz = ceil(celldata.dataSize .* scale);
   end
   if isempty(sh)
      sh = celldata.imageShifts;
      sh = cellfunc(@(x) ceil(x .* scale), sh);
   end
   
   % read sequentially 
   nodes = celldata.nodes;
   cdat = cell(1, length(nodes));
   parfor i = 1:numel(cdat)
      cdat{i} = celldata.asource.dataResample(scale, nodes(i)); %#ok<PFBNS>
   end
   celldata = cdat;
   
elseif isa(celldata, 'ImageSource')
   % read sequentially 
   cdat = cell(celldata.cellSize);
   
   parfor i = 1:numel(cdat)
      cdat{i} = celldata.dataResample(scale, i); %#ok<PFBNS>
   end
   celldata = cdat;
   
else
   celldata = cellfunc(@(x) imresize(x, scale), celldata);
end

overlap  = ceil(scale .* overlap);

% 
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