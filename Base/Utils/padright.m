function res = padright(vec, n, val)
%
% res = padright(res, n, val)
% res = padright(res, n, type)
%
% description: 
%    padds the vector vec on the right to total size n, with values val, or
%    type as in padarray
%
% See also: padarray

res = vec(:)';

if nargin < 3
   val = 0;
end

l =length(res);
if n <= l
   res = res(1:n);
   return
end

if iscell(res)
   res = padcell(res, n, val);
else
   res = padarray(res, [0, n - l], val, 'post');
end

end
