function out = laplacianFilter(image, varargin)
%
% out = laplacianFilter(image, alpha, padding)
%
% description:
%    apply a Laplacian filter to image
%
% input:
%    image        image to be filtered
%    alpha        space paramter ( 0 )
%    padding      padding of array at borders
%
% output:
%    out          filtered image
%
% See also: fspecial2, fspecial3, imfilter

out = linearFilter(image, 'laplacian', varargin{:});

end