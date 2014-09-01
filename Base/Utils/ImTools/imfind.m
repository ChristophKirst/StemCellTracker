function [p,varargout] = imfind(image)
%
% [p,q,l] = imfind(image);
% [p,q]   = imfind(image);
% pql     = imfind(image);
%
% descrition:
%     finds the coordinates of nonzero values in a 3d image
%
% input:
%     image      2d or 3d grayscale image
%
% output:
%     p,q,l      pixel coordinates
%     pql        matrix of coordinates

dim = ndims(image);

if dim < 2 || dim > 3
   error('imfind: input is not a 2d or 3d bw or gray scale image!')
end

idx = find(image);
if dim ==1
      [p,q] = ind2sub(size(image), idx);
else
      [p,q,l] = ind2sub(size(image), idx);
end

if nargout <= 1
   if dim == 2 
      p = [p, q];
   else
      p = [p, q, l];
   end
elseif dim == 2 && nargout == 2
   varargout{1} = q;
elseif dim ==3 && nargout == 3
   varargout{1} = q; 
   varargout{2} = l;
else
   error('imfind: inconsistent number of output arguments');
end

% if nargin > 1
%    if ischar(coords) && strcmp(coords, 'matlab') % swap h and w
%       pp = p; p = q; q = pp;
%    end
% end

end