function d = ndimsd(val)
%
% d = ndimsd(val)
%
% description:
%    returns the number of dimensions of size > 1
%  
% See also: ndims, ndims1

si = size(val);
si(si <= 1) = [];
d = length(si);

end
