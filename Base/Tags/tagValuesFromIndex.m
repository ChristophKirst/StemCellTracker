function [tvs, ids] = tagValuesFromIndex(tagRange, i)
%
% [tvs, ids] = tagValuesFromIndex(tagRange, i)
%
% description:
%    returns tag values for index i as a cell, tvs{i} corresponds to the i-th name in tagRange
%
% input:
%    tagRange   the tag range
%    i          the index assumed to be scalar
%
% output:
%   tvs         tag values as cell
%   ids         the position in each dimension

si = tagRangeSize(tagRange);

if i > prod(si) || i < 1
   error('tagValuesFromIndex: index %g out of bounds %g:%g !', i, 1, prod(si));
end

ids = imind2sub(si, i);

names = fieldnames(tagRange);
nnames = length(names);

tvs = cell(1, nnames);

for i = 1:nnames
   tr = tagRange.(names{i});
   if iscell(tr)
      tvs{i} = tr{ids(i)};
   else
      tvs{i} = tr(ids(i));
   end
end

end