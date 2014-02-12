function out = linearFilter(image, ker, varargin)
%
% out = linearFilter(image, ker, varargin)
% out = linearFilter(image, type, ksize, varargin)
% out = linearFilter(image, ..., padding)
%
% description:
%    apply the filter kernel ker to image
%
% input:
%    image        image to be filtered
%    ker          filter kernel matrix or type
%    type         special type
%    ksize        h x w (x l) size of the fitler
%    padding      padding of array at borders as used in imfilter
%
% output:
%    out          filtered image
%
% See also: fspecial2, fspecial3, imfilter

dim = ndims(image);

if dim < 2 || dim > 4
   error('linearFilter: image must be a 2d or 3d gray scale image')
end

if ischar(ker)
   if nargin < 3
      ksize = 3;
   else
      ksize = varargin{1};
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
   error('linearFilter: image and kernel dimensions do not agree')
end
if dim == 3 && size(image,3) == 3
   error('linearFilter: image must gray scale image')
end

if nargin > paddoff && ischar(varargin{end})
   padding = varargin{end};
else
   padding = 'replicate';
end


out = imfilter(image, ker, padding);

end