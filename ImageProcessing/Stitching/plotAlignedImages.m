function varargout = plotAlignedImages(imgs, shifts, varargin)
%
% img = plotAlignedImages(imgs, shifts, param)
%
% description:
%    plots aligned images using the shifts and colors conveniently for 
%    inspection of the result
%
% input: 
%    imgs    images
%    shifts  shifts
%

%    param   (optional) paramter struct with entries
%            .tileformat     tile format which specifies different coloring of tiles
% 

%param = parseParameter(varargin{:});
%size(imgs)
%size(shifts)

if ~iscell(imgs) || ~iscell(shifts) || numel(imgs) ~= numel(shifts)
   error('plotAlignedImages: inconsistent input');
end

% determine image size and absolute shifts w.r.t composed image

imgsizes = cellfunc(@size, imgs);

[ashifts, asize] = absoluteShiftsAndSize(shifts, imgsizes);


% calculate connectivity structure for coloring
pairs = alignPairsFromOverlap(ashifts, imgsizes);

if isempty(pairs)
   cids = 1;
else
   edges = [[pairs.from]', [pairs.to]'];
   cids = graphVertexColoring(1:numel(imgsizes), edges);
end

ncols = min(max(cids), 8);
cids = mod(cids, 8);

if ncols <= 1
   cols = {[0, 1, 0], [1, 0, 1]};
elseif ncols <= 4
   cols = {[0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0], [0, 0, 0.5]};
else
   cols =  {...
         [0.25, 0,     0   ],... 
         [0,    0.25,  0   ],...
         [0,    0,     0.25],...
         [0.25  0.25,  0   ],...
         [0,    0.25,  0.25],... 
         [0.25  0,     0.25],... 
         [0.125,0.25,  0   ],...
         [0.125,0,     0.25]};
end

img = imgray2color(zeros(asize), 'white');

imax = cellfun(@(x) max(x(:)), imgs);
nrm = double(max(imax(:)));

for i = 1:numel(imgs)
   imga = imgray2color(double(imgs{i})/nrm, cols{cids(i)});
   imgb = imextract(img, [1 + ashifts{i}, 1],  [ashifts{i} + imgsizes{i}, 3]);
   img  = imreplace(img, imga + imgb, [1 + ashifts{i},1]);
end

implot(img);

if nargout > 0
   varargout{1} = img;
end

end