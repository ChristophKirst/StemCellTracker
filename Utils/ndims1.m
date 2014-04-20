function d = ndims1(val)
%
% d = ndims1(val)
%
% description:
%    returns dimension of val as ndims
%    returns 1 if ndims(val) == 2 and size(val,2) == 1
%  
% See also: ndims

si = size(val);
d = length(si);

if d == 2 && si(2) == 1
   d = 1;
end

% data dim
% si = size(val);
% si(si == 1) = [];
% d = length(si);

end
