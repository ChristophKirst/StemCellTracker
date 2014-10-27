function [img, scale] = alignCellPreview(celldata, varargin)

param = parseParameter(varargin);

scale = getParameter(param, 'scale',    0.1);
pos   = getParameter(param, 'position', []);
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
   if isempty(pos)
      pos = celldata.imageShifts;
      pos = cellfunc(@(x) ceil(x .* scale), pos);
   end
   
   % read sequentially 
   nodes = celldata.nodes;
   cdat = cell(1, length(nodes));
   for i = 1:numel(cdat)
      cdat{i} = celldata.asource.dataResample(scale, nodes(i));
   end
   celldata = cdat;
   
elseif isa(celldata, 'ImageSource')
   % read sequentially 
   cdat = cell(celldata.cellSize);
   for i = 1:numel(cdat)
      cdat{i} = celldata.dataResample(scale, i);
   end
   celldata = cdat;
   
else
   celldata = cellfunc(@(x) imresize(x, scale), celldata);
end

overlap  = ceil(scale .* overlap);

% 
if isempty(pos)
   isize = size(celldata{1}); % assume all images are same size
   pos = imgrid(size(celldata));
 
   for i = 1:length(pos)
      pos{i} = pos{i}(:);
   end
   pos = [pos{:}] - 1;
   pos = pos .* repmat(isize - overlap, size(pos,1), 1) + 1;

   pos = num2cell(pos,2);
end

if isempty(siz)
   [pos, siz] = absoluteShiftsAndSize(pos, cellfunc(@size, celldata));
end


%var2char(pos)

img = stitchImages(celldata, pos, param, 'size', siz);


if getParameter(param, 'lines', false)
   pos = cell2mat(pos);
   posx = unique(pos(:,1)); 
   
   mval = max(img(:));
   for p = posx(:)'
      img = impixelline(img, [p,1], [p, siz(2)], mval);
   end
   
   posy = unique(pos(:,2)); 
   for p = posy(:)'
      img = impixelline(img, [1, p], [siz(1), p], mval);
   end

end


end