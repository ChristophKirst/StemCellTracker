function out = logFilter(image, varargin)
%
% out = logFilter(img, radius, padding)
%
% description:
%    apply Laplacian of Gaussian filter to image 
%
% input:
%    image        image to be filtered
%    ksize        h x w (x l) of filer 
%    sigma        std of Gaussian = ( (radius-1)/3 )
%    padding      padding of array at borders
%
% output:
%    out          filtered image
%
% See also: fspecial2, fspecial3, imfilter

out = linearFilter(image, 'log', varargin{:});

end