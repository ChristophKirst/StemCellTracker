function a = extract(a, ind, dim)
%
% a = extract(a, ind, dim)
%
% description:
%    extract all data for indices ind in dimension dim of array a
%
% input:
%    a   array
%    ind indices to extract
%    dim dimensions
%

d = ndims(a);
inds = repmat({':'},1,d);
inds{dim} = ind;
a = a(inds{:});

end