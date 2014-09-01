function isize = hipto2sizes(pto)
%
% isize = hipto2sizes(pto)
%
% description:
%   convert pto parameter to realtive shifts of images in pixel
%
% input:
%   pto     pto file or struct as obtained by hiparsepto
%
% output:
%   isize   image sizes
%
% See also: hiparsepto

if ischar(pto)
   pto = hiparsepto(pto);
end

n = length(pto);
isize = cell(1,n);

if n == 0
   return
end

for i = 1:n
   isize{i} = [pto(i).h, pto(i).w];   
end

end






