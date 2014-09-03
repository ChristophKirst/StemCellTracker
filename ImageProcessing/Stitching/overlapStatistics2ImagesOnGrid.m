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

[img1, img2] = overlap2ImagesOnGrid(imgs, varargin{:});

fr.var = var(double(img1(:)));
fr.max = double(max(img1(:)));
fr.min = double(min(img1(:)));

to.var = var(double(img2(:)));
to.max = double(max(img2(:)));
to.min = double(min(img2(:)));

stats.from = fr;
stats.to   = to;

end



