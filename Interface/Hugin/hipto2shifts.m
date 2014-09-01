function [shifts, varargout] = hipto2shifts(pto)
%
% shifts = hipto2shifts(pto)
% [shifts, isizes] = hipto2shifts(pto)
%
% description:
%   convert pto parameter to realtive shifts of images in pixel
%
% input:
%   pto     pto file or struct as obtained by hiparsepto
%
% output:
%   shifts  shifts of the images as cell array
%
% See also: hiparsepto

if ischar(pto)
   pto = hiparsepto(pto);
end

n = length(pto);
shifts = cell(1,n);

if n == 0
   return
end

% coordinates in panorama to pixel factor
hFoV = pto(1).v * 2*pi/360;
vFoV = 2 * atan(tan(hFoV/2)* pto(1).h / pto(1).w);

phfov = 2 * (1 - pto(1).TrZ) * sin(hFoV/2);
pano2pixel = pto(1).w / phfov;
%fprintf('p2p: %g\n', 1/pano2pixel)

% ref point of first image in pano coords
xy0 = - (1 - pto(1).TrZ) * sin([hFoV, vFoV]/2);

for i = 1:n
  
   % calculate vFoV from hFoV
   hFoV = pto(i).v  * 2*pi/360;
   vFoV = 2 * atan(tan(hFoV/2)* pto(i).h / pto(i).w);
   
   % calculate ref point of image
   xyi = - (1 - pto(i).TrZ) * sin([hFoV, vFoV]/2);
   xyi = xyi + [pto(i).TrX, pto(i).TrY];
   
   dxy = xyi - xy0;
   shifts{i} = circshift(round(pano2pixel .* dxy), [0, 1]);
   
end

if nargout > 1
   varargout{1} = hipto2sizes(pto);
end

end






