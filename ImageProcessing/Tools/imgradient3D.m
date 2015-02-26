function [imggrad, Gx, Gy, Gz] = imgradient3D(img, varargin)
%
% imggrad = imgradient3D(img)
%
% description:
%     calculates 3d gradient magnitudes of 3d image stack
%
% input:
%     img      3d grayscale image 
%     param    parameter struct with entries
%              .depth     aspect ratio assuming x,y is 1
%
% output:
%     imggrad  gradient magnitude

d = getParameter(parseParameter(varargin), 'depth', 1);

sx = [-1 0 1; -2 0 2; -1 0 1];
sy = [1 2 1; 0 0 0; -1 -2 -1];


szx = cat(3, sx, sx, sx);
szy = cat(3, sy, sy, sy);
szz = cat(3, ones(3,3), zeros(3,3), -ones(3,3));


Gx = convn(img, szx, 'same');
Gy = convn(img, szy, 'same');
Gz = 1/d * convn(img, szz, 'same');
imggrad = sqrt(Gx.^2 + Gy.^2 + Gz.^2);

end