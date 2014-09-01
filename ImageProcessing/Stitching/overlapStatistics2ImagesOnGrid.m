function stats = overlapStatistics2ImagesOnGrid(imgs, varargin)
%
% stats = overlapStatistics2ImagesOnGrid(imgs, varargin)
%
% description:
%    for overlaps with no signal accurate alignment is not possible
%    this routine determines statistics in the overlap regoin 
%
% input:
%    imgs   images as prealigned cell array
%    param  parameter struct with entries
%           .overlap.max   maximal overlap
%
% output:
%    stats  struct with entries .from and .to
%           each having entries .var, .max, .min of potential overlap regoins
%


if numel(imgs) ~= 2
   error('overlapStatistics2ImagesOnGrid: expect precisle cell array with two images!');
end

param = parseParameter(varargin{:});
mo = getParameter(param, 'overlap.max', []);

img1 = imgs{1}; img2 = imgs{2};

% orientation

si = size(imgs);
dim = ndims(imgs);
pos = find(si == 2, 1);
idx = repmat({':'}, 1, dim);

% stats of img1 overlap region
si1 = size(img1);
sip  = si1(pos);

if isempty(mo)
   moi = sip;
else
   moi = mo;
end

idx{pos} = sip-moi+1:sip;

img1 = img1(idx{:});
fr.var = var(double(img1(:)));
fr.max = double(max(img1(:)));
fr.min = double(min(img1(:)));


% stats of img2 overlap region

si2 = size(img2);
sip  = si2(pos);

if isempty(mo)
   moi = sip;
else
   moi = mo;
end

idx{pos} = 1:moi;

img2 = img2(idx{:});
to.var = var(double(img2(:)));
to.max = double(max(img2(:)));
to.min = double(min(img2(:)));



stats.from = fr;
stats.to   = to;

end



