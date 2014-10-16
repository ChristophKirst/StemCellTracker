function img = alignCellPreview(celldata, varargin)

param = parseParameter(varargin);

scale = getParameter(param, 'scale',    0.5);
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
if isa(celldata, 'ImageSource')
   % read sequentially 
   
   cdat = cell(celldata.cellSize);
   for i = 1:numel(cdat)
      cdat{i} = celldata.dataResample(scale, i);
   end
   celldata = cdat;
end



celldata = cellfunc(@(x) imresize(x, scale), celldata);
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

end