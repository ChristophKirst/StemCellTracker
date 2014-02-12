function out = diskFilter(image, varargin)
%
% out = diskFilter(image, hsize, radius_in, radius_out, w_in, w_out, padding)
%
% description:
%    apply a double disk filter with weight w_in in radius_in and w_out
%    between radius_in and radius_out
%
% input:
%    image        image to be filtered
%    ksize        h x w (xl) bouding box and radius of outer ring
%    ring_width   width of outer ring (0)
%    w_disk       weight inner disk
%    w_ring       weight outer ring
%    padding      padding of array at borders
%
% output:
%    out          filtered image
%
% See also: fspecial2, fspecial3, imfilter

out = linearFilter(image, 'disk', varargin{:});

end