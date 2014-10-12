function [info, pos] = tileshapeFromImageInfo(fname)
%
% ts = tileshapeFromImageInfo(fname)
%
% description:
%     tries to infer the tile shape form the meta data in ImageInfo
%
% input:
%     fname    filename, info struct / ImageInfo class
%
% output:
%     ts       the inferred tile shape
%     pos      the positions read form the meta data
%


if ischar(fname)
   info = imreadBFInfo(fname);
else
   info = fname;
end

mdata = info.imetadata;

% detect grid from info

parnames = mdata.Parameter;
parvals  = mdata.Value;

ids = ~cellfun(@isempty, strfind(parnames, 'X position for position #'));

s = length('X position for position #');

% x positions
tidsx = cellfun(@(x) str2double(x(s+1:end)), parnames(ids))';
xpos = cellfun(@str2double, parvals(ids))';

% y positions
ids =  ~cellfun(@isempty, strfind(parnames, 'Y position for position #'));
tidsy = cellfun(@(x) str2double(x(s+1:end)), parnames(ids))';
ypos = cellfun(@str2double, parvals(ids))';

% positions

if nargout > 1
   pos = [xpos(tidsx), ypos(tidsy)]';
end


% using meanshift clustering to detect number of images in each tiling direction

is = info.idatasize(1:2);
is = is .* info.iscale; 

ni = prod(info.icellsize);

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

warning('tileshapeFromImageInfo: cannot infer tile format, found u = %g, v = %g', ntx, nty)

end









