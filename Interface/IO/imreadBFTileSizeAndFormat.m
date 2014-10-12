function [ts, tf] = imreadBFTileSizeAndFormat(fname, varargin)
%
% [ts, tf] = imreadBFTileSizeAndFormat(fname, varargin)
%
% description:
%     tries to infer the tile format and shape for the meta data
%
% input:
%     fname    filename or bfreader
%
% output:
%     ts       the inferred tile shape
%     tf       the inffered tile format
%


% single color, time an z are usually good
% pos ordered according to S
[ts, pos] = imreadBFTileSize(fname, 'C', 1, 'Z', 1, 'T', 1);

tf = 'UV';

% from positions try to infer orientation
pos = reshape(pos,ts);

if ts(1) >1
   dxy = pos{2} - pos{1};
   
   if dxy(1) < 0
      tf(1) = 'u';
   end
end

if ts(2) >1
   dxy = pos{ts(1)+1} - pos{1};
   
   if dxy(2) > 0
      tf(2) = 'v';
   end
end

end









