function [ts, pos] = imreadBFTileSize(fname, varargin)
%
% ts = imreadBFTileSize(fname, varargin)
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

ireader= imreadBFReader(fname);

pos = imreadBFStagePositions(ireader, varargin);
ppos = [pos{:}];

xpos = ppos(1,:);
ypos = ppos(2,:);

% use meanshift clustering to detect number of images in each tiling direction

is = [ireader.getSizeX; ireader.getSizeY];
vs = imreadBFVoxelSize(ireader, 'S', 1, varargin);
vs = vs{1};

is = is .* vs(1:2); 

ni = ireader.getSeriesCount;

cl = clusterMeanShift(xpos, is(1)/2);
ntx = length(cl);

if mod(ni, ntx) == 0
   ts = [ntx, ni/ntx];
   return
end

cl = clusterMeanShift(ypos, is(2)/2);
nty = length(cl);

if mod(ni, nty) == 0
   ts   = [ni/nty, nty];
   return
end

warning('imreadBFTileshape: cannot infer tile format, found U = %g, V = %g', ntx, nty)

end









