function [dim, s1, s2, minovl, maxovl, maxshift] = align2ImagesLeftRightParameter(img1, img2, varargin)
%
% [minovl, maxovl, maxshift] = align2ImagesLeftRightParameter(img1, img2, varargin)
%
% descritpion:
%   generate parameter for image alignment on gird routines align2ImagesLeftRightByXXX
%
% input:
%   imgs    images to align
%   param   parameter struct with entries:
%              .overlap.max   maximal overlap of images in primary direction ([] = max image width)
%              .overlap.min   minimal overlap if images in primary direction (1)
%              .shift.max     maximal shift away form border in secondary dimensions ([10, 10]) 
% output:
%   various parsed and checked parameter
%
% See also: align2ImagesLeftRightOrient

param = parseParameter(varargin{:});

% image sizes and dim
s1 = size(img1); s2 = size(img2);
dim = length(s1);
if dim ~= length(s2)
   warning('align2ImagesLeftRight: image dimension mismatch: %g ~= %g, padding with ones!', length(s1), length(s2));
   % fill smaller dim with ones
   dim = max(dim, length(s2));
   s1 = padright(s1, dim - length(s1), 1);
   s2 = padright(s2, dim - length(s2), 1);
end

% get overlaps and pad with zeros
minovl = getParameter(param, 'overlap.min', 1);
if minovl <= 0
   warning('align2ImagesLeftRight: overlap.min has to be at least 1!');
   minovl = 1;
end

maxovl = getParameter(param, 'overlap.max', []);
if isempty(maxovl) || isinf(maxovl)
   maxovl = min(s1(1), s2(1));
end
if maxovl > min(s1(1), s2(1))
   warning('align2ImagesLeftRight: overlap.max bigger than image sizes, setting it to %g !', min(s1(1), s2(1)));
   maxovl = min(s1(1), s2(1));
end

maxshift = getParameter(param, 'shift.max', 10);
maxshift = padright(maxshift, dim-1, maxshift(1));

end





