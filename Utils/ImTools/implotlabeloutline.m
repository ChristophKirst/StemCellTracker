function implotlabeloutline(img, imglab, varargin)
%
% implotlabeloutline(img, imglab, param)
%
% description:
%     plots image as tiling together with pixel surface of the labeled image
%
% input:
%     img      grayscale image
%     imglab   labeled image
%     param    (optional) paramter struct for color map as in imlabelcolormap
%
% See also: impixelsurface, imsurfaceplot3d, imoverlay, imlabelcolormap

implottiling(imoverlaylabel(img, impixelsurface(imglab), false, varargin{:}))

end