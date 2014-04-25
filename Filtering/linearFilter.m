function imgf = linearFilter(img, ker, varargin)
%
% imgf = linearFilter(img, ker, varargin)
% imgf = linearFilter(img, type, ksize, varargin)
% imgf = linearFilter(img, ..., padding)
%
% description:
%    apply the filter kernel ker to img
%
% input:
%    img          image to be filtered
%    ker          filter kernel matrix or type
%    type         special type
%    ksize        h x w (x l) size of the fitler
%    padding      padding of array at borders as used in imfilter
%
% output:
%    imgf         filtered image
%
% See also: fspecial2, fspecial3, imfilter

dim = ndims(img);

if dim < 2 || dim > 4
   error('linearFilter: img must be a 2d or 3d gray scale img')
end

if ischar(ker)
   if nargin < 3
      ksize = 3;
   else
      ksize = varargin{1};
      varargin = varargin(2:end);
   end
   
   switch dim
      case 2
         ker = fspecial2(ker, ksize, varargin{:});
      case 3
         ker = fspecial3(ker, ksize, varargin{:});
   end
   
   paddoff = 4;
   
else
   paddoff = 3;
end

if ndims(ker) ~= dim
   error('linearFilter: img and kernel dimensions do not agree')
end
if dim == 3 && size(img,3) == 3
   error('linearFilter: img must be gray scale img')
end

if nargin > paddoff && ischar(varargin{end})
   padding = varargin{end};
else
   padding = 'symmetric';
end


imgf = imfilter(img, ker, padding);

end