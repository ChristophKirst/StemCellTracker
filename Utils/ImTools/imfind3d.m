function [p,q,l] = imfind3d(image)
%
% [p,q,l] = imfind3d(image)
%
% descrition:
%     finds the coordinates of nonzero values in a 3d image
%
% input:
%     image      3d grayscale image
%
% output:
%     p,q,l      pixel coordinates

if ndims(image) ~= 3
   error('imfind3d: input is not a 3d bw or gray scale image!')
end

idx = find(image);
[p,q,l] = ind2sub(size(image), idx);

% if nargin > 1
%    if ischar(coords) && strcmp(coords, 'matlab') % swap h and w
%       pp = p; p = q; q = pp;
%    end
% end

end