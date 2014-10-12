function out = filterDoG(image, varargin)
%
% out = filterDoG(img, ksize, sigma_in, sigma_out, padding)
%
% description:
%    apply a difference of Gaussians filter with inner and outer radius radius_in and radius_out
%
% input:
%    image        image to be filtered
%    ksize        size of the fitler
%    sigma_in     std of inner Gaussian ([] = 1/1.5 * sigma_out)
%    sigma_out    std of outer Gaussian ([] = hsize / 2 / sqrt(2 log(2)) )
%    padding      padding of array at borders
%
% output:
%    out          filtered image
%
% See also: fspecial2, fspecial3, imfilter

out = filterLinear(image, 'dog', varargin{:});

end