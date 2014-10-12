function out = filterDisk(image, varargin)
%
% out = filterDisk(image, ksize, ring_width, w_in, w_out, padding)
%
% description:
%    apply a double disk filter with weight w_in in inner disk and w_out
%    on out ring with width ring_width (inner radius is chooes to fit the 
%    filter kernel size (ring_width = 0 -> pure disk filter)
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

out = filterLinear(image, 'disk', varargin{:});

end