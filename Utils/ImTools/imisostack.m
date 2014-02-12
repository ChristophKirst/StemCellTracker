function isostack = imisostack(stack, zslices)
%
% stack = imisostack(image, xyspacing, zspacing)
%
% description:
%     returns isotropic stack by copying z-slices to match xyspacing
%
% input:
%     stack    image stack
%     zslices  number of zslices that make stack isotropic
%
% output:
%     isostack  new isotropic stack


dim = size(stack);

if length(dim) ~=3
   error('imisostack: expect 3d stack of gray scale images.')
end

isostack = zeros(dim(1), dim(2), zslices * dim(3));

for i = 1:dim(3);
   isostack(:,:,(i-1)*zslices+1:i*zslices) = repmat(stack(:,:,i), 1,1,zslices);
end

end

