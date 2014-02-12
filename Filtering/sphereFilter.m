function out = sphereFilter(image, varargin)
%
% out = sphereFilter(image, ksize, radius, padding)
%
% description:
%    apply a filter with spherical weights to image
%
% input:
%    image        image to be filtered
%    ksize        h x w ( x l) filter size
%    radius       radius of outer ring
%    padding      padding of array at borders
%
% output:
%    out          filtered image
%
% See also: fspecial, imfilter

out = linearFilter(image, 'sphere', varargin{:});

end