function [info, pos] = tileshapeFromBF(fname, varargin)
%
% ts = tileshapeFromBF(fname)
%
% description:
%     tries to infer the tile shape form a bio-format file
%
% input:
%     fname    filename or bfreader
%
% output:
%     ts       the inferred tile shape
%     pos      the positions read form the meta data
%

pos = imreadBFPositions(fname. varargin);
pos = cell2mat([pos{:}]);

% use meanshift clustering to detect number of images in each tiling direction

is = [ireader.getSizeX, ireader.getSizeY];
vs = imreadBFVoxelsize('series', 1);
vs = vs{1};
is = is .* vs(1:2); 

ni = ireader.getSeriesCount;

cl = clusterMeanShift(xpos, is(1)/2);
ntx = length(cl);

if mod(ni, ntx) == 0
   info.icellformat = 'uv';
   info.icellsize   = [ntx, ni/ntx];
   return
end

cl = clusterMeanShift(ypos, is(2)/2);
nty = length(cl);

if mod(ni, nty) == 0
   info.icellformat = 'uv';
   info.icellsize   = [ni/nty, nty];
   return
end

warning('tileshapeFromZVI: cannot infer tile format, found u = %g, v = %g', ntx, nty)

end









