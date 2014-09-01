function stats = overlapStatisticsImagePair(pair, varargin)
%
% stats = overlapStatisticsImagePair(pair, varargin)
%
% description:
%    for overlaps with no signal accurate alignment is not possible
%    this routine determines statistics in the overlap regoin 
%
% input:
%    imgs   images as prealigned cell array
%    param  parameter struct with entries
%           .overlap.max   maximal overlap
%           .aligned       images are aligned if true and pair.shift is used to detemine overlap (false)
%
% output:
%    stats  struct with entries .from and .to
%           each having entries .var, .max, .min of potential overlap regoins
%

param = parseParameter(varargin{:});
al = getParameter(param, 'aligned', false);

if al
   stats = overlapStatistics2AlignedImages(pair.from, pair.to, pair.shift);
else
   c = pair.toCell();
   stats = overlapStatistics2ImagesOnGrid(c, param);
end


end



